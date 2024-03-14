"""\
collect stats from different samples
this will collect all summary files for the listed samples and produce one big table
optionally, this will also collect all inserts and produce a multi-sheet excel file
clustering analysis results can also be collected into an excel file
"""


def setup_argparse(parser):
    parser.add_argument("--samples", dest="samples", help="""comma-separated list of samples""")
    parser.add_argument(
        "--sample_stats",
        dest="sample_stats",
        help="""required: output file with sample_stats""",
    )
    parser.add_argument("--summary_path", dest="summary_path", default="output/{sample}/{sample}_summary.csv",
                        help="""path pattern for summary files ["output/{sample}/{sample}_summary.csv"]""")
    parser.add_argument("--inserts", dest="inserts", help="""output file with inserts""")
    parser.add_argument("--inserts_path", dest="inserts_path", default="output/{sample}/{sample}_inserts.tsv",
                        help="""path pattern for insert tables ["output/{sample}/{sample}_inserts.tsv"]""")
    parser.add_argument(
        "--cluster_stats",
        dest="cluster_stats",
        help="""output file with cluster stats""",
    )
    parser.add_argument("--cluster_analysis_path", dest="cluster_analysis_path", default="output/{sample}/{sample}_cluster_analysis.csv",
                        help="""path pattern for cluster_analysis files ["output/{sample}/{sample}_cluster_analysis.csv"]""")


def run(args):
    import os
    import pandas as pd
    from logzero import logger
    import glob
    from openpyxl import Workbook

    if args.samples is None:
        raise Exception("no list of samples given!")
    samples = args.samples.split(",")

    dfs = dict(
        (
            sample,
            pd.read_csv(gl, header=None, index_col=0)
            .squeeze()
            .dropna(),
        )
        for sample in samples
        for gl in glob.glob(args.summary_path.format(sample=sample))
        )
    dfs = dict((k, v[v.index.notnull()]) for k, v in dfs.items())
    df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1)
    logger.info("saving sample stats to " + args.sample_stats)
    df.T.to_csv(args.sample_stats)

    if args.inserts is not None:
        wb = Workbook()
        with pd.ExcelWriter(args.inserts, engine="openpyxl") as writer:
            writer.workbook = wb
            for sample in samples:
                for gl in glob.glob(args.inserts_path.format(sample=sample)):
                    try:
                        tmp = pd.read_csv(gl,
                            sep="\t",
                            header=0,
                        )
                    except pd.errors.EmptyDataError:
                        pass
                    tmp.to_excel(writer, sheet_name=sample, index=False)
                    logger.info("adding inserts for {0}".format(sample))

    if args.cluster_stats is not None:
        wb = Workbook()
        with pd.ExcelWriter(args.cluster_stats, engine="openpyxl") as writer:
            writer.workbook = wb
            for sample in samples:
                for gl in glob.glob(args.cluster_analysis_path.format(sample=sample)):
                    tmp = pd.read_csv(
                        gl,
                        header=0,
                        index_col=0,
                    ).sort_values("size", ascending=False)
                    tmp.to_excel(writer, sheet_name=sample, index=True)
                    logger.info("adding clustering stats for {0}".format(sample))
