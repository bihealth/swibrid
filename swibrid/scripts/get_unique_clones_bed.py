""" get bed file with unique clones"""


def setup_argparse(parser):

    parser.add_argument("-b", "--bed", dest="bed", help="""input bed file""")
    parser.add_argument(
        "-c", "--clustering", dest="clustering", help="""clustering results"""
    )
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
    from .utils import get_eff_nclust

    logger.info("reading clustering from " + args.clustering)
    clustering = pd.read_csv(args.clustering, header=0)
    clustering = clustering[~clustering["cluster"].isna()]
    n_eff = get_eff_nclust(clustering["cluster"], cut=args.cut)
    clones = (
        clustering["cluster"].dropna().astype(int).value_counts().index[:n_eff]
    )
    reads = (
        clustering[clustering["cluster"].astype(int).isin(clones)]
        .groupby("cluster")
        .first()
        .iloc[:, 0]
        .values
    )

    logger.info("selecting entries from " + args.bed)
    for line in open(args.bed):
        ls = line.split()
        if ls[3] in reads and ls[0] == "chr14":
            sys.stdout.write(line)
