"""\
downsample (aligned + processed) reads for a sample
"""


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
        help="""downsampled sample name""",
    )
    parser.add_argument(
        "-n", "--n", dest="n", type=int, default=5000, help="""use only n reads [%(default)d]"""
    )
    parser.add_argument(
        "--aligner",
        dest="aligner",
        default="last",
        help="""alignment algorithm used [%(default)s]""",
    )
    parser.add_argument(
        "--msa_path",
        dest="msa_path",
        default="pipeline/{sample}/{sample}_msa.npz",
        help="""path pattern for msa files [%(default)s]""",
    )
    parser.add_argument(
        "--msa_csv_path",
        dest="msa_csv_path",
        default="pipeline/{sample}/{sample}_msa.csv",
        help="""path pattern for msa csv files [%(default)s]""",
    )
    parser.add_argument(
        "--info_path",
        dest="info_path",
        default="input/{sample}_info.csv",
        help="""path pattern for info csv files [%(default)s]""",
    )
    parser.add_argument(
        "--process_path",
        dest="process_path",
        default="pipeline/{sample}/{sample}_processed.out",
        help="""path pattern for processing output [%(default)s]""",
    )
    parser.add_argument(
        "--process_stats_path",
        dest="process_stats_path",
        default="pipeline/{sample}/{sample}_process_stats.csv",
        help="""path pattern for process stats csv [%(default)s]""",
    )
    parser.add_argument(
        "--alignment_pars_path",
        dest="alignment_pars_path",
        default="pipeline/{sample}/{sample}_{aligner}_pars.npz",
        help="""path pattern for alignment pars [%(default)s]""",
    )
    parser.add_argument(
        "--breakpoint_alignments_path",
        dest="breakpoint_alignments_path",
        default="pipeline/{sample}/{sample}_breakpoint_alignments.csv",
        help="""path pattern for breakpoint alignments [%(default)s]""",
    )
    parser.add_argument(
        "--inserts_tsv_path",
        dest="inserts_tsv_path",
        default="pipeline/{sample}/{sample}_inserts.tsv",
        help="""path pattern for inserts tsv [%(default)s]""",
    )
    parser.add_argument(
        "--bed_path",
        dest="bed_path",
        default="pipeline/{sample}/{sample}.bed",
        help="""path pattern for bed [%(default)s]""",
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

    assert os.path.isfile(
        args.msa_path.format(sample=sample)
    ), "please run pipeline on input sample at least up to construct_msa"

    logger.info("getting reads from MSA")
    logger.info("reading MSA for " + sample)
    msa = scipy.sparse.load_npz(args.msa_path.format(sample=sample))
    msa_df = pd.read_csv(args.msa_csv_path.format(sample=sample), index_col=0, header=0)

    reads = msa_df.sample(n=min(args.n, msa_df.shape[0])).index
    inds = msa_df.index.get_indexer(reads)

    logger.info("using {0} of {1} reads".format(len(reads), msa_df.shape[0]))

    logger.info("getting read info")
    info = pd.read_csv(args.info_path.format(sample=sample), header=0, index_col=0)
    os.makedirs(os.path.dirname(args.info_path.format(sample=args.output)), exist_ok=True)
    info.loc[reads].to_csv(args.info_path.format(sample=args.output))

    logger.info("getting alignment pars")
    pars = np.load(args.alignment_pars_path.format(sample=sample, aligner=args.aligner))
    os.makedirs(
        os.path.dirname(args.alignment_pars_path.format(sample=args.output, aligner=args.aligner)),
        exist_ok=True,
    )
    np.savez(args.alignment_pars_path.format(sample=args.output, aligner=args.aligner), **pars)

    logger.info("getting output from process_alignments")
    process = pd.read_csv(args.process_path.format(sample=sample), header=0, index_col=0, sep="\t")
    process_stats = pd.read_csv(
        args.process_stats_path.format(sample=sample), header=None, index_col=0
    ).squeeze()
    breakpoint_alignments = pd.read_csv(
        args.breakpoint_alignments_path.format(sample=sample), header=0, index_col=0
    )

    os.makedirs(os.path.dirname(args.process_path.format(sample=args.output)), exist_ok=True)
    process.loc[reads].to_csv(args.process_path.format(sample=args.output), sep="\t", index=True)
    os.makedirs(os.path.dirname(args.process_stats_path.format(sample=args.output)), exist_ok=True)
    process_stats.to_csv(args.process_stats_path.format(sample=args.output), header=False)

    os.makedirs(
        os.path.dirname(args.breakpoint_alignments_path.format(sample=args.output)), exist_ok=True
    )
    breakpoint_alignments.loc[reads.intersection(breakpoint_alignments.index)].to_csv(
        args.breakpoint_alignments_path.format(sample=args.output)
    )

    logger.info("getting inserts")
    try:
        inserts = pd.read_csv(
            args.inserts_tsv_path.format(sample=sample), index_col=0, header=0, sep="\t"
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
    bed = pd.read_csv(args.bed_path.format(sample=sample), header=None, sep="\t")

    insert_reads = reads.intersection(inserts.index)
    inserts = inserts.loc[insert_reads].reset_index(names="read")
    os.makedirs(os.path.dirname(args.inserts_tsv_path.format(sample=args.output)), exist_ok=True)
    inserts.to_csv(args.inserts_tsv_path.format(sample=args.output), index=False, sep="\t")
    os.makedirs(os.path.dirname(args.bed_path.format(sample=args.output)), exist_ok=True)
    bed.set_index(3).loc[reads].reset_index()[[0, 1, 2, "read", 4, 5, 6, 7, 8, 9, 10, 11]].to_csv(
        args.bed_path.format(sample=args.output), header=False, index=False, sep="\t"
    )

    logger.info("writing MSA to " + args.msa_path.format(sample=args.output))
    os.makedirs(os.path.dirname(args.msa_csv_path.format(sample=args.output)), exist_ok=True)
    msa_df.loc[reads].to_csv(args.msa_csv_path.format(sample=args.output))
    os.makedirs(os.path.dirname(args.msa_path.format(sample=args.output)), exist_ok=True)
    scipy.sparse.save_npz(
        args.msa_path.format(sample=args.output),
        msa.tocsr()[inds].tocoo(),
    )
