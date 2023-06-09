"""construct (pseudo) MSA from processed LAST output"""


def setup_argparse(parser):
    parser.add_argument(
        "--coords",
        dest="coords",
        help="""process_alignments.py coordinate output""",
    )
    parser.add_argument(
        "--sequences",
        dest="sequences",
        help="""process_alignments.py sequence output""",
    )
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("--out", dest="out", help="""file with  read info""")
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000:-",
        help="""coordinates of switch region [chr14:106050000-106337000:-]""",
    )
    parser.add_argument(
        "--switch_annotation",
        dest="switch_annotation",
        help="""bed file with switch annotation""",
    )
    parser.add_argument("--nmax", dest="nmax", type=int, help="""use only nmax reads""")
    parser.add_argument(
        "--use_orientation",
        dest="use_orientation",
        action="store_true",
        default=False,
        help="""include alignment block orientation""",
    )


def run(args):

    import numpy as np
    import re
    import pandas as pd
    from Bio import SeqIO
    import scipy.sparse
    from logzero import logger
    from collections import defaultdict
    import gzip
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
        get_switch_coverage,
        decode_coords,
        decode_insert,
        shift_coord,
        ncodes,
    )

    (
        switch_chrom,
        switch_start,
        switch_end,
        switch_orientation,
    ) = parse_switch_coords(args.switch_coords)
    switch_anno = read_switch_anno(args.switch_annotation)
    cov_int, Ltot, eff_start, eff_end, anno_recs = get_switch_coverage(
        switch_anno, switch_chrom, switch_start, switch_end
    )

    logger.info("reading processed read coordinates from {0}".format(args.coords))

    read_mappings = {}
    read_isotypes = {}
    read_orientation = {}
    read_coverage = {}
    read_inserts = {}

    nreads = 0
    nignored = 0

    for line in open(args.coords):

        ls = line.strip("\n").split("\t")
        read = ls[0]
        isotype = ls[1]
        orientation = ls[2]
        cov = float(ls[3])

        mappings = []
        inserts = []
        for ll in ls[4:]:
            if "insert" not in ll:
                coords = decode_coords(ll)
                mappings.append(coords)
            elif "insert" in ll:
                inserts.append(decode_insert(ll))

        read_mappings[read] = mappings
        read_inserts[read] = inserts
        read_isotypes[read] = isotype
        read_orientation[read] = orientation
        read_coverage[read] = cov

        nreads += 1

    logger.info("{0} reads kept, {1} reads ignored".format(nreads, nignored))

    reads_all = list(read_mappings.keys())

    if args.nmax:
        reads = reads_all[: args.nmax]
    else:
        reads = reads_all

    nreads = len(reads)
    read_isotypes = pd.Series(read_isotypes)[reads]
    read_coverage = pd.Series(read_coverage)[reads]
    read_orientation = pd.Series(read_orientation)[reads]
    read_inserts = pd.Series(read_inserts)[reads]

    logger.info("reading processed read sequences from {0}".format(args.sequences))
    nt_ignored = defaultdict(int)
    nt_assigned = defaultdict(int)

    i, j, x = [], [], []
    for rec in SeqIO.parse(
        gzip.open("{0}".format(args.sequences), "rt")
        if args.sequences.endswith(".gz")
        else args.sequences,
        "fasta",
    ):
        read = rec.id.rsplit("@", 1)[0]
        if read not in reads:
            continue
        chrom, start, end, orientation = decode_coords(rec.id.rsplit("@", 1)[1])
        # find coordinates of this part of the read sequence
        ns = shift_coord(start, cov_int, ignore=True) - eff_start
        ne = shift_coord(end, cov_int, ignore=True) - eff_start
        if len(rec.seq) == ne - ns and ns >= 0 and ne <= Ltot:
            seq = str(rec.seq).upper()
            if args.use_orientation and orientation == "-":
                seq = seq.lower()
            nt_assigned[read] += len(seq)
        elif min(ne, Ltot) >= max(ns, 0):
            logger.warn("coordinates {0}:{1}-{2} partially outside range".format(chrom, start, end))
            ne = min(ne, Ltot)
            ns = max(ns, 0)
            if np.isfinite(shift_coord(start, cov_int)) or ne > Ltot:
                # left end outside bounds
                seq = str(rec.seq[(len(rec.seq) - ne + ns) :]).upper()
            else:
                # right end outside bounds
                seq = str(rec.seq[: (ne - ns)]).upper()
            if args.use_orientation and orientation == "-":
                seq = seq.lower()
            nt_assigned[read] += len(seq)
            nt_ignored[read] += len(rec.seq) - len(seq)
        else:
            logger.warn("coordinates {0}:{1}-{2} outside range".format(chrom, start, end))
            nt_ignored[read] += len(rec.seq)
            continue
        n = read_isotypes.index.get_loc(read)
        # find all the non-gap positions in this part of the read sequence
        for m in re.finditer("[ACGTacgt]", seq):
            i.append(n)
            j.append(ns + m.start())
            x.append(ncodes[m.group()])

    if np.unique(np.vstack([i, j]), axis=1).shape[1] < len(i):
        raise Exception("indices appear multiple times!\n")

    msa = scipy.sparse.csr_matrix((x, (i, j)), shape=(nreads, Ltot), dtype=np.int8)

    use = msa.sum(1).A1 > 0
    logger.info("removing {0} reads without coverage".format((~use).sum()))
    msa = msa[use].tocoo()

    logger.info("saving msa to {0}".format(args.msa))
    scipy.sparse.save_npz(args.msa, msa)

    logger.info("saving info to {0}".format(args.out))
    nt_assigned = pd.Series(nt_assigned)
    nt_ignored = pd.Series(nt_ignored)
    outside_range = (nt_ignored / (nt_assigned + nt_ignored)).fillna(0)
    read_info = pd.concat(
        [read_isotypes, read_orientation, read_coverage, outside_range],
        axis=1,
        keys=["isotype", "orientation", "coverage", "outside_range"],
    )
    read_info["insert"] = read_inserts.apply(
        lambda x: ",".join(
            map(
                lambda y: "{0}:{1}-{2}:{3}".format(
                    y.group("insert_chrom"),
                    y.group("insert_start"),
                    y.group("insert_end"),
                    y.group("orientation"),
                ),
                x,
            )
        )
    )

    read_info[use].to_csv(args.out)
