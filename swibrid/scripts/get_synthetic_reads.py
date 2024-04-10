"""\
create synthetic reads from bed file
mutations, insertions and deletions will be added according to parameters estimated by LAST
clone sizes can be uniform or distributed as Poisson or Negative Binomial
additional homozygous / heterozygous / other variants can be added
"""


def setup_argparse(parser):
    parser.add_argument("-b", "--bed", dest="bed", help="""required: input bed file""")
    parser.add_argument(
        "-r", "--reference", dest="reference", help="""required: reference sequence"""
    )
    parser.add_argument("-p", "--par", dest="par", help="""required: file with LAST parameters""")
    parser.add_argument(
        "-o", "--out", dest="out", help="""required: output fasta / fastq (.gz) file"""
    )
    parser.add_argument("-i", "--info", dest="info", help="""output info file""")
    parser.add_argument(
        "-n",
        "--n",
        dest="n",
        type=int,
        help="""use n random entries of input bed file""",
    )
    parser.add_argument(
        "-k",
        "--k",
        dest="k",
        type=float,
        default=1,
        help="""number of mutated copies per sequence [1]""",
    )
    parser.add_argument(
        "-d",
        "--distribution",
        dest="distribution",
        default="poisson",
        help="""distribution type for k (delta | poisson | nbinom | P | NB | D) [nbinom]""",
    )
    parser.add_argument(
        "--nbinom_alpha",
        dest="nbinom_alpha",
        default=4,
        type=float,
        help="""dispersion alpha of negative binomial [4]""",
    )
    parser.add_argument(
        "--variants",
        dest="variants",
        nargs="?",
        help="""variant file (like output of find_variants) to add variants""",
    )
    parser.add_argument(
        "--rearrangements",
        dest="rearrangements",
        nargs="?",
        help="""rearrangement file (bed-like output of find_rearrangements) to add rearrangements""",
    )
    parser.add_argument("-s", "--seed", dest="seed", type=int, default=1)
    parser.add_argument(
        "--no-mutations",
        dest="no_mutations",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--no-deletions",
        dest="no_deletions",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--no-insertions",
        dest="no_insertions",
        action="store_true",
        default=False,
    )


def mutate_seq(seq, pars):
    import numpy as np

    def rand_string(size):
        return "".join("ACGT"[r] for r in np.random.randint(0, 4, size))

    # add deletions
    r = np.random.rand(len(seq))
    # positions where deletions occur
    j = np.where(r < pars["p_open_del"])[0]
    if len(j) > 0:
        # sizes of deletions
        s = np.random.poisson(pars["p_extend_del"], len(j)) + 1
        seq2 = (
            seq[: j[0]]
            + "".join(seq[j[i] + s[i] : j[i + 1]] for i in range(len(j) - 1))
            + seq[j[-1] + s[-1] :]
        )
    else:
        seq2 = seq

    # add insertions
    r = np.random.rand(len(seq2))
    # positions where insertions occur
    j = np.where(r < pars["p_open_ins"])[0]
    if len(j) > 0:
        # sizes of insertions
        s = np.random.poisson(pars["p_extend_ins"], len(j)) + 1
        seq3 = (
            seq2[: j[0]]
            + rand_string(s[0])
            + "".join(seq2[j[i] : j[i + 1]] + rand_string(s[i + 1]) for i in range(len(j) - 1))
            + seq2[j[-1] :]
        )
    else:
        seq3 = seq2

    # add mutations
    s = np.array(list(seq3))
    # get cumulative probs for mutated bases given bases in the sequence
    m = (pars["p_c"] / pars["p_c"].sum(1)[:, np.newaxis]).cumsum(1)
    # get random numbers for each position
    r = np.random.rand(len(seq3))
    for i, x in enumerate("ACGT"):
        # find all positions where seq = x
        j = np.where(s == x)[0]
        # choose new bases at those positions according to m[i]
        k = np.searchsorted(m[i], r[j])
        # replace bases at those positions
        s[j] = np.array(list("ACGT"))[k]

    return "".join(s)


def mutate_rec(rec, pars, add_id=""):
    import numpy as np
    from Bio import SeqRecord, Seq

    if pars is not None:
        seq = Seq.Seq(mutate_seq(str(rec.seq).upper(), pars))
    else:
        seq = Seq.Seq(str(rec.seq).upper())
    if np.random.rand() < 0.5:
        seq = seq.reverse_complement()
    return SeqRecord.SeqRecord(
        seq,
        id=rec.id + add_id,
        name=rec.name + add_id,
        description=rec.description + add_id,
    )


def run(args):
    import gzip
    import numpy as np
    import pandas as pd
    import pysam
    import scipy.stats
    from Bio import SeqIO, Seq, SeqRecord
    from logzero import logger
    from .utils import parse_LAST_pars, shift_coord, RC

    np.random.seed(args.seed)

    logger.info("reading LAST pars from " + args.par)
    pars = parse_LAST_pars(args.par)

    if args.no_deletions:
        pars["p_open_del"] = 0
        pars["p_extend_del"] = 0

    if args.no_insertions:
        pars["p_open_ins"] = 0
        pars["p_extend_ins"] = 0

    if args.no_mutations:
        pars = None

    logger.info("reading bed file from " + args.bed)
    bed = (
        pd.read_csv(args.bed, sep="\t", header=None, index_col=None)
        .sample(n=args.n)
        .reset_index(drop=True)
    )
    logger.info("loading reference from " + args.reference)
    reference = pysam.FastaFile(args.reference)

    if args.variants:
        logger.info("loading variants information from " + args.variants)
        variants = pd.read_csv(args.variants, sep="\t", header=0)
    else:
        variants = pd.DataFrame(
            [],
            columns=[
                "chrom",
                "position",
                "region",
                "rel_pos",
                "ref",
                "alt",
                "counts",
                "pval",
                "padj",
                "pval_clust",
                "padj_clust",
                "pval_strand",
                "type",
                "anno",
                "motif",
            ],
        )

    if args.rearrangements:
        logger.info("loading rearrangement information from " + args.rearrangements)
        rearrangements = pd.read_csv(
            args.rearrangements, sep="\t", index_col=None, header=None
        ).dropna()
    else:
        rearrangements = pd.DataFrame([], columns=["chrom", "start", "end", "type", "score"])

    if args.out.endswith("fastq.gz"):
        logger.info("writing fastq.gz output to " + args.out)
        outf = gzip.open(args.out, "wt")
    else:
        logger.info("writing fasta output to " + args.out)
        outf = open(args.out, "w")

    logger.info("copying reads" if args.no_mutations else "mutating reads")
    reads = []
    nreads = 0
    info = {}
    nvars_added = 0
    ninvs_added = 0
    ndups_added = 0

    for n, (
        chrom,
        start,
        end,
        name,
        score,
        strand,
        cstart,
        cend,
        color,
        nblocks,
        blocksizes,
        blockstarts,
    ) in bed.iterrows():
        blocksizes = list(map(int, blocksizes.strip(",").split(",")))
        blockstarts = list(map(int, blockstarts.strip(",").split(",")))

        cov_int = [
            (start + bstart, start + bstart + bsize)
            for bsize, bstart in zip(blocksizes, blockstarts)
        ]

        seq = list("".join(reference.fetch(chrom, s, e) for s, e in cov_int).upper())

        vars_added = []
        invs_added = []
        dups_added = []

        # loop over variants and add them to clone sequence
        for _, (
            _,
            pos,
            _,
            _,
            REF,
            ALT,
            counts,
            _,
            _,
            _,
            _,
            _,
            var_type,
            _,
            _,
        ) in variants.iterrows():
            rel_pos = shift_coord(pos - 1, cov_int)
            counts = np.array(list(map(int, counts.split(","))))
            # ignore variants outside the current region
            if pos < start or pos >= end or not np.isfinite(rel_pos) or seq[rel_pos] != REF:
                continue
            # ignore heterozygous variants for every other clone
            if (var_type == "het0" and n % 2 == 1) or (var_type == "het1" and n % 2 == 0):
                continue
            # ignore other variants depending on their frequency
            var_freq = counts["ACGT".index(ALT)] / np.sum(counts)
            if var_type == "n.d." and np.random.rand() > var_freq:
                continue

            nvars_added += 1
            vars_added.append("{0}:{1}{2}>{3}".format(chrom, pos, REF, ALT))
            seq[rel_pos] = ALT

        # loop over rearrangements ad add them to clone sequence
        for _, (_, rstart, rend, rtype, rfreq) in rearrangements.iterrows():
            rel_start = shift_coord(rstart, cov_int)
            rel_end = shift_coord(rend, cov_int)

            # ignore events outside the current regions
            if (
                rstart < start
                or rend >= end
                or not np.isfinite(rel_start)
                or not np.isfinite(rel_end)
            ):
                continue
            # ignore events depending on allele frequency
            if np.random.rand() > rfreq:
                continue

            if rtype == "inversion":
                ninvs_added += 1
                invs_added.append("{0}:{1}-{2}".format(chrom, rstart, rend))
                seq = (
                    seq[:rel_start] + [RC[s] for s in seq[rel_start:rel_end][::-1]] + seq[rel_end:]
                )
            if rtype == "duplication":
                ndups_added += 1
                dups_added.append("{0}:{1}-{2}".format(chrom, rstart, rend))
                seq = (
                    seq[:rel_start]
                    + seq[rel_start:rel_end]
                    + seq[rel_start:rel_end]
                    + seq[rel_end:]
                )

        rec = SeqRecord.SeqRecord(seq=Seq.Seq("".join(seq)), id=name, name=name, description=name)

        if args.distribution in ["poisson", "P"]:
            ncopies = np.random.poisson(args.k)
        elif args.distribution in ["nbinom", "NB"]:
            ncopies = scipy.stats.nbinom.rvs(
                1.0 / args.nbinom_alpha, 1.0 / (1.0 + args.nbinom_alpha * args.k), size=1
            )[0]
        else:
            ncopies = int(args.k)

        for k in range(ncopies):
            read = mutate_rec(
                rec,
                pars,
                add_id="::{0}:{1}-{2}_{3}".format(chrom, start, end, k + 1),
            )
            read.letter_annotations["phred_quality"] = [30] * len(read)
            reads.append(read)
            info[read.id] = {
                "length": len(read),
                "barcodes": "BC01@1-30:+",
                "primers": "primer_fw@0-50:+;primer_rv@{0}-{1}:+".format(len(read) - 50, len(read)),
                "internal": "",
                "true_cluster": "c{0}".format(n),
                "true_sequence": str(rec.seq),
                "true_haplotype": "het0" if n % 2 == 0 else "het1",
                "true_variants": ";".join(vars_added),
                "true_inversions": ";".join(invs_added),
                "true_duplications": ";".join(dups_added),
            }

        if args.out.endswith("fastq.gz"):
            SeqIO.write(reads, outf, "fastq")
        else:
            SeqIO.write(reads, outf, "fasta")

        nreads += len(reads)
        reads = []

    if args.info is not None:
        logger.info("writing read info to " + args.info)
        pd.DataFrame.from_dict(info, orient="index").to_csv(args.info)

    outf.close()

    logger.info(
        "{0} variants, {1} inversions, {2} duplications added in {3} reads".format(
            nvars_added, ninvs_added, ndups_added, nreads
        )
    )
