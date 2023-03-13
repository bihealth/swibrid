"""collect stats from different samples"""


def setup_argparse(parser):
    parser.add_argument(
        "--samples", dest="samples", help="""comma-separated list of samples"""
    )
    parser.add_argument(
        "-o", "--output", dest="outf", help="""output csv file"""
    )


def run(args):

    import os
    import pandas as pd

    if args.samples is None:
        raise Exception("no list of samples given!")
    samples = args.samples.split(",")

    dfs = dict(
        (
            sample,
            pd.read_csv("plots/" + sample + "_QC.csv", header=None, index_col=0)
            .squeeze()
            .dropna(),
        )
        for sample in samples
        if os.path.isfile("plots/" + sample + "_QC.csv")
    )
    dfs = dict((k, v[v.index.notnull()]) for k, v in dfs.items())
    df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1)
    df.T.to_csv(args.outf)
