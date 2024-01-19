"""downsample (aligned + processed) reads a sample"""


def setup_argparse(parser):
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="""input sample to downsample""",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output",
        help="""downsampled sample name""",
    )
    parser.add_argument(
        "-n", "--n", dest="n", type=int, default=5000, help="""use only n reads [5000]"""
    )
    parser.add_argument(
        "--aligner", dest="aligner", default="last", help="""alignment algorithm used [last]"""
    )


def run(args):
    import os
    import numpy as np
    import pandas as pd
    import scipy.sparse
    from logzero import logger

    assert args.input, "please provide an input sample"
    sample = args.input

    assert args.output, "please provide an output sample name"

    assert os.path.isfile("output/{0}/{0}_msa.npz".format(sample)), "please run pipeline on input sample at least up to construct_msa"

    logger.info("getting reads from MSA")
    logger.info("reading MSA for " + sample)
    msa = scipy.sparse.load_npz("output/{0}/{0}_msa.npz".format(sample))
    msa_df = pd.read_csv("output/{0}/{0}_msa.csv".format(sample), index_col=0, header=0)

    reads = msa_df.sample(n=min(args.n, msa_df.shape[0])).index
    inds = msa_df.index.get_indexer(reads)

    logger.info("using {0} of {1} reads".format(len(reads), msa_df.shape[0]))

    os.makedirs("output/{0}".format(args.output), exist_ok=True)

    logger.info("getting read info")
    info = pd.read_csv("input/{0}_info.csv".format(sample), header=0, index_col=0)
    info.loc[reads].to_csv("input/{0}_info.csv".format(args.output))

    logger.info("getting filter stats")
    filter_stats = pd.read_csv(
            "output/{0}/{0}_filter_stats.csv".format(sample), header=None, index_col=0
        ).squeeze()
    filter_stats.to_csv("output/{0}/{0}_filter_stats.csv".format(args.output), header=False)

    logger.info("getting alignment pars")
    pars=np.load("output/{0}/{0}_{1}_pars.npz".format(sample, args.aligner))
    np.savez("output/{0}/{0}_{1}_pars.npz".format(args.output, args.aligner), **pars)

    logger.info("getting output from process_alignments")
    process = pd.read_csv(
        "output/{0}/{0}_processed.out".format(sample), header=0, index_col=0, sep="\t"
    )
    process_stats = pd.read_csv(
        "output/{0}/{0}_process_stats.csv".format(sample), header=None, index_col=0
        ).squeeze()
    breakpoint_alignments = pd.read_csv(
            "output/{0}/{0}_breakpoint_alignments.csv".format(sample), header=0, index_col=0
            )

    process.loc[reads].to_csv("output/{0}/{0}_processed.out".format(args.output), sep="\t", index=True)
    process_stats.to_csv("output/{0}/{0}_process_stats.csv".format(args.output), header=False)
    breakpoint_alignments.loc[reads.intersection(breakpoint_alignments.index)].to_csv(
        "output/{0}/{0}_breakpoint_alignments.csv".format(args.output)
    )

    logger.info("getting inserts")
    try:
        inserts = pd.read_csv(
            "output/{0}/{0}_inserts.tsv".format(sample), index_col=0, header=0, sep="\t"
            )
    except (ValueError, pd.errors.EmptyDataError):
        inserts = pd.DataFrame(
            [],
            columns=[
                "isotype",
                "switch_coords",
                "insert_pos",
                "insert_coords",
                "insert_anno",
                "sequence",
                ],
            )
    bed = pd.read_csv(
        "output/{0}/{0}_inserts.bed".format(sample), header=None, sep="\t"
       )

    insert_reads = reads.intersection(inserts.index)
    inserts = inserts.loc[insert_reads].reset_index(names='read')
    inserts.to_csv("output/{0}/{0}_inserts.tsv".format(args.output), index=False, sep="\t")
    bed.set_index(3).loc[reads].reset_index()[[0, 1, 2, "read", 4, 5, 6, 7, 8, 9, 10, 11]].to_csv(
        "output/{0}/{0}_inserts.bed".format(args.output), header=False, index=False, sep="\t"
        )

    logger.info("writing MSA to output/{0}/{0}_msa.npz".format(args.output))
    msa_df.loc[reads].to_csv("output/{0}/{0}_msa.csv".format(args.output))
    scipy.sparse.save_npz(
        "output/{0}/{0}_msa.npz".format(args.output),
        msa.tocsr()[inds].tocoo(),
    )
