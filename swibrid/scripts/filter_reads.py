"""filter reads"""


def setup_argparse(parser):
    """setup argparser subparser"""
    parser.add_argument("-i", "--input", dest="input", default="input", help="""input fastq""")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="output",
        help="""output fasta (.gz)""",
    )
    parser.add_argument("--info", dest="info", help="""read info""")
    parser.add_argument("--nmax", dest="nmax", type=int, help="""use only nmax reads""")
    parser.add_argument("--stats", dest="stats", help="""output stats""")
    parser.add_argument(
        "--min-length",
        dest="min_length",
        default=500,
        type=int,
        help="""minimum read length""",
    )
    parser.add_argument(
        "--only-complete",
        dest="only_complete",
        action="store_true",
        default=False,
        help="""use only complete reads with two primers""",
    )
    parser.add_argument(
        "--keep-internal",
        dest="keep_internal",
        action="store_true",
        default=False,
        help="""keep reads with internal primers""",
    )


def filter_reads(reads, info, args, stats):
    import pandas as pd
    import re

    for rec in reads:
        if args.nmax and stats["nreads_initial"] > args.nmax:
            break

        stats["nreads_initial"] += 1

        keep = True
        if len(rec) < args.min_length:
            stats["nreads_removed_short"] += 1
            keep = False

        if info is not None:
            if pd.isna(info.loc[rec.id, "primers"]):
                primers = set()
            else:
                primers = set(
                    map(
                        lambda x: re.sub("primer_|_[1-9]", "", x.split("@")[0]),
                        info.loc[rec.id, "primers"].split(";"),
                    )
                )

            if args.only_complete and not (
                any("rv" in p for p in primers) and any("fw" in p for p in primers)
            ):
                stats["nreads_removed_incomplete"] += 1
                keep = False

            internal = info.loc[rec.id, "internal"]
            if not args.keep_internal and not pd.isna(internal) and "primer" in internal:
                stats["nreads_removed_internal_primer"] += 1
                keep = False

        if keep:
            stats["nreads_to_map"] += 1
            yield rec


def run(args):
    from Bio import SeqIO
    import gzip
    import pandas as pd

    if args.info:
        read_info = pd.read_csv(args.info, header=0, index_col=0)

    stats = {"nreads_initial": 0, "nreads_removed_short": 0, "nreads_removed_incomplete": 0, "nreads_removed_internal_primer": 0, "nreads_to_map": 0}

    SeqIO.write(
        filter_reads(
            SeqIO.parse(gzip.open(args.input, "rt"), "fastq"),
            read_info,
            args,
            stats,
        ),
        gzip.open(args.output, "wt") if args.output.endswith(".gz") else args.output,
        "fastq",
    )

    if args.stats:
        pd.Series(stats).to_csv(args.stats, header=False)
