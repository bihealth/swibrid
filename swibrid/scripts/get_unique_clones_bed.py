"""\
produce a bed file with unique clones from the biggest clusters.
"""


def setup_argparse(parser):
    parser.add_argument("-b", "--bed", dest="bed", help="""required: input bed file""")
    parser.add_argument(
        "-c", "--clustering", dest="clustering", help="""required: clustering results"""
    )
    parser.add_argument("-o", "--output", dest="output", help="""output file [stdout]""")
    parser.add_argument(
        "--cut",
        dest="cut",
        type=float,
        default=0.8,
        help="""use only this fraction of the biggest clusters [.8]""",
    )


def run(args):
    import sys
    from logzero import logger
    import pandas as pd

    logger.info("reading clustering from " + args.clustering)
    clustering = pd.read_csv(args.clustering, header=0)
    clustering = clustering[~clustering["cluster"].isna()]
    clusters = clustering["filtered_cluster"]
    clones = clusters[clusters >= 0]

    reads = (
        clustering[clustering["cluster"].astype(int).isin(clones)]
        .groupby("cluster")
        .first()
        .iloc[:, 0]
        .values
    )

    if args.output:
        outf = open(args.output, "w")
    else:
        outf = sys.stdout

    logger.info("selecting entries from " + args.bed)
    for line in open(args.bed):
        ls = line.split()
        if ls[3] in reads and ls[0] == "chr14":
            outf.write(line)
