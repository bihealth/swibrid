"""collect stats from different samples"""


def setup_argparse(parser):
    parser.add_argument(
        "--samples", dest="samples", help="""comma-separated list of samples"""
    )
    parser.add_argument(
        "--sample_stats",
        dest="sample_stats",
        help="""output file with sample_stats""",
    )
    parser.add_argument(
        "--inserts", dest="inserts", help="""output file with inserts"""
    )
    parser.add_argument(
        "--cluster_stats",
        dest="cluster_stats",
        help="""output file with cluster stats""",
    )


def run(args):

    import os
    import pandas as pd
    import glob
    from logzero import logger
    from openpyxl import Workbook

    if args.samples is None:
        raise Exception("no list of samples given!")
    samples = args.samples.split(",")

    dfs = dict(
        (
            sample,
            pd.read_csv(
                "plots/" + sample + "_summary.csv", header=None, index_col=0
            )
            .squeeze()
            .dropna(),
        )
        for sample in samples
        if os.path.isfile("plots/" + sample + "_summary.csv")
    )
    dfs = dict((k, v[v.index.notnull()]) for k, v in dfs.items())
    df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1)
    logger.info("saving sample stats to " + args.sample_stats)
    df.T.to_csv(args.sample_stats)

    wb = Workbook()
    with pd.ExcelWriter(args.inserts, engine="openpyxl") as writer:
        writer.book = wb
        for rf in glob.glob(os.path.join("inserts", "*_inserts.tsv")):
            name = rf.split("/")[-1].split("_inserts.tsv")[0]
            try:
                tmp = pd.read_csv(rf, sep="\t", header=0)
                tmp.to_excel(writer, sheet_name=name, index=False)
                logger.info("adding inserts for {0}".format(name))
            except pd.errors.EmptyDataError:
                logger.warn("no inserts for {0}".format(name))
                pass

    wb = Workbook()
    with pd.ExcelWriter(args.cluster_stats, engine="openpyxl") as writer:
        writer.book = wb
        for rf in glob.glob(os.path.join("cluster", "*_analysis.csv")):
            name = rf.split("/")[-1].split("_analysis.csv")[0]
            try:
                tmp = pd.read_csv(rf, header=0, index_col=0).sort_values(
                    "size", ascending=False
                )
                tmp.to_excel(writer, sheet_name=name, index=True)
                logger.info("adding clustering stats for {0}".format(name))
            except pd.errors.EmptyDataError:
                logger.warn("no clustering stats for {0}".format(name))
                pass
