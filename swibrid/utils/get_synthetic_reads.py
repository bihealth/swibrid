"""create synthetic reads from bed file with coordinates"""


def setup_argparse(parser):

    parser.add_argument("-b", "--bed", dest="bed", help="""input bed file""")
    parser.add_argument(
        "-r", "--reference", dest="reference", help="""reference sequence"""
    )
    parser.add_argument(
        "-n",
        "--n",
        dest="n",
        type=int,
        help="""use n random entries of input bed file""",
    )
    parser.add_argument(
        "-p", "--par", dest="par", help="""file with LAST parameters"""
    )
    parser.add_argument(
        "-o", "--out", dest="out", help="""output fasta / fastq file"""
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
        "-i", "--info", dest="info", help="""output info file"""
    )
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
    parser.add_argument("--seed", dest="seed", type=int, default=1)
    parser.add_argument(
        "--poisson_n",
        dest="poisson_n",
        action="store_true",
        default=False,
        help="""treat n as mean of Poisson random variable""",
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
            + "".join(
                seq2[j[i] : j[i + 1]] + rand_string(s[i + 1])
                for i in range(len(j) - 1)
            )
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


def mutate_rec(rec, pars, add_id=None):

    import copy
    from Bio import SeqRecord, Seq

    seq = str(rec.seq).upper()
    if pars is None:
        res = copy.copy(rec)
    else:
        res = SeqRecord.SeqRecord(
            seq=Seq.Seq(mutate_seq(seq, pars)),
            id=rec.id,
            name=rec.name,
            description=rec.description,
        )
    if add_id:
        res.id += add_id
        res.name += add_id
        res.description += add_id
    return res


def run(args):

    import gzip
    import numpy as np
    import pandas as pd
    import pysam
    from Bio import SeqIO, Seq, SeqRecord
    from logzero import logger
    from .helpers import parse_LAST_pars

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
    bed = pd.read_csv(args.bed, sep="\t", header=None, index_col=None).sample(
        n=args.n
    )
    logger.info("loading reference from " + args.reference)
    reference = pysam.FastaFile(args.reference)

    logger.info("copying reads" if args.no_mutations else "mutating reads")
    reads = []
    info = {}
    for _, (
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

        seq = "".join(
            reference.fetch(chrom, start + bstart, start + bstart + bsize)
            for bsize, bstart in zip(blocksizes, blockstarts)
        )
        rec = SeqRecord.SeqRecord(
            seq=Seq.Seq(seq), id=name, name=name, description=name
        )

        if args.poisson_n:
            ncopies = np.random.poisson(args.k)
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
                "primers": "primer_fw@0-50:+;primer_rv@{0}-{1}:+".format(
                    len(read) - 50, len(read)
                ),
                "internal": "",
            }

    logger.info("writing sequences to " + args.out)
    if args.out.endswith("fastq.gz"):
        SeqIO.write(reads, gzip.open(args.out, "wt"), "fastq")
    else:
        SeqIO.write(reads, open(args.out, "w"), "fasta")
    if args.info is not None:
        logger.info("writing read info to " + args.info)
        pd.DataFrame.from_dict(info, orient="index").to_csv(args.info)
