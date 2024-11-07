"""\
construct (pseudo) MSA from processed alignments: 
this script will take aligned segments and sequences from process_alignments and construct a MSA.
the MSA is stored as a sparse matrix of n_reads x n_positions, where the positions run over the concatenation of individual switch regions specified in a bed file.
matrix values m code for coverage and nucleotide identity:

- m % 10 gives nucleotide identity with A=1, C=2, G=3, T=4, gaps are zeros
- m // 10 gives coverage of genomic regions by one read (1x, 2x, ...)  and indicates tandem duplications; if ``--use_orientation`` is set, inversions are also considered and associated coverage values are negative
"""


def setup_argparse(parser):
    parser.add_argument(
        "--coords",
        dest="coords",
        help="""required: process_alignments coordinate output""",
    )
    parser.add_argument(
        "--sequences",
        dest="sequences",
        help="""required: process_alignments sequence output""",
    )
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""required: output file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("--out", dest="out", help="""required: output file with read info""")
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000:-",
        help="""coordinates of switch region [chr14:106050000-106337000:-]""",
    )
    parser.add_argument(
        "--switch_annotation",
        dest="switch_annotation",
        help="""required: bed file with coordinates of individual switch regions""",
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

    if args.nmax:
        process = pd.read_csv(args.coords, sep="\t", header=0, index_col=0, nrows=args.nmax)
    else:
        process = pd.read_csv(args.coords, sep="\t", header=0, index_col=0)
    nreads = process.shape[0]

    read_mappings = dict(
        (read, [decode_coords(m) for m in row.split(";")])
        for read, row in process["switch_mappings"].items()
    )

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
        if read not in read_mappings.keys():
            continue
        chrom, start, end, orientation = decode_coords(rec.id.rsplit("@", 1)[1])
        # find coordinates of this part of the read sequence
        ns = shift_coord(start, cov_int, ignore=True)
        ne = shift_coord(end, cov_int, ignore=True)
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
        n = process.index.get_loc(read)
        assert isinstance(n, int), "read index is not integer!"
        # find all the non-gap positions in this part of the read sequence
        for m in re.finditer("[ACGTacgt]", seq):
            i.append(n)
            j.append(ns + m.start())
            x.append(ncodes[m.group()])

    i = np.array(i)
    j = np.array(j)
    x = np.array(x)

    _, ind, inv = np.unique(np.vstack([i, j]), axis=1, return_index=True, return_inverse=True)

    cov = np.sign(x[ind])
    cons = np.abs(x[ind])
    # cov2 = np.sign(x[ind])
    # cons2 = np.abs(x[ind])

    if len(ind) < len(i):
        xx = scipy.sparse.csr_matrix((x, (np.arange(len(inv)), inv)))
        cts = (xx != 0).sum(0).A1
        p1 = np.where(cts > 1)[0]
        xx = xx[:, p1]
        p2 = np.where((xx != 0).sum(1).A1 > 0)[0]
        xx = xx[p2, :]
        logger.warn(
            "{0} positions in MSA are covered multiple times; getting consensus".format(len(p1))
        )
        # for k, p in enumerate(p1):
        #    cov2[p] = np.sum(np.sign(xx[:, k].data))
        #    cons2[p] = np.bincount(np.abs(xx[:, k].data)).argmax()
        cov[p1] = (np.sum(xx > 0, axis=0) - np.sum(xx < 0, axis=0)).A1
        cons[p1] = (
            np.vstack(
                [
                    np.bincount(xx.nonzero()[1][np.abs(xx.data) == k + 1], minlength=xx.shape[1])
                    for k in range(4)
                ]
            ).argmax(0)
            + 1
        )

    msa = scipy.sparse.csr_matrix(
        (10 * cov + cons, (i[ind], j[ind])), shape=(nreads, Ltot), dtype=np.int8
    )

    unique_entries = np.unique(msa.data % 10)
    assert np.max(unique_entries) < 5 and np.min(unique_entries > 0), "invalid entries in MSA"

    use = np.abs(msa).sum(1).A1 > 0
    logger.info("removing {0} reads without coverage".format((~use).sum()))
    msa = msa[use].tocoo()

    logger.info("saving msa to {0}".format(args.msa))
    scipy.sparse.save_npz(args.msa, msa)

    logger.info("saving info to {0}".format(args.out))
    nt_assigned = pd.Series(nt_assigned)
    nt_ignored = pd.Series(nt_ignored)
    process["frac_ignored"] = (nt_ignored / (nt_assigned + nt_ignored)).fillna(0)
    process["inserts"] = (
        process["inserts"]
        .dropna()
        .apply(
            lambda x: ",".join(
                map(
                    lambda y: "{0}:{1}-{2}:{3}".format(
                        y.group("insert_chrom"),
                        y.group("insert_start"),
                        y.group("insert_end"),
                        y.group("orientation"),
                    ),
                    [decode_insert(insert) for insert in x.split(";")],
                )
            )
        )
    )

    out_cols = [
        "isotype",
        "orientation",
        "frac_mapped",
        "frac_mapped_multi",
        "frac_ignored",
        "inserts",
    ]
    if "mate_breaks" in process.columns:
        process["mate_breaks"] = process["mate_breaks"].apply(
            lambda x: ";".join(str(shift_coord(int(y), cov_int)) for y in x.split(";"))
            if not pd.isnull(x)
            else ""
        )
        out_cols += ["mate_breaks"]
    process[use][out_cols].to_csv(args.out)
