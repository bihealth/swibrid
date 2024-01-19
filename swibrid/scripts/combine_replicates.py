"""combine (aligned + processed) reads from replicates"""


def setup_argparse(parser):
    parser.add_argument(
        "-s",
        "--samples",
        dest="samples",
        help="""comma-separated list of samples to combine""",
    )
    parser.add_argument(
        "-c",
        "--combined",
        default="combined",
        dest="combined",
        help="""combined sample name""",
    )
    parser.add_argument(
        "--nmax", dest="nmax", type=int, default=50000, help="""use only nmax reads [50000]"""
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

    samples = args.samples.split(",")
    assert len(samples) > 1, "please provide a comma-separated list of more than one sample"

    assert all(
        os.path.isfile("output/{0}/{0}_msa.npz".format(s)) for s in samples
    ), "please run pipeline on invidividual samples at least up to construct_msa"

    logger.info("combining reads from MSA output")
    msa = {}
    msa_df = {}
    for sample in samples:
        logger.info("reading MSA for " + sample)
        msa[sample] = scipy.sparse.load_npz("output/{0}/{0}_msa.npz".format(sample))
        msa_df[sample] = pd.read_csv("output/{0}/{0}_msa.csv".format(sample), index_col=0, header=0)

    msa_df = pd.concat(msa_df.values(), axis=0, keys=msa_df.keys()).reset_index(
        level=0, names="sample"
    )

    reads = msa_df.sample(n=min(args.nmax, msa_df.shape[0])).index
    inds = msa_df.index.get_indexer(reads)

    logger.info("using {0} of {1} reads".format(len(reads), msa_df.shape[0]))

    os.makedirs("output/{0}".format(args.combined), exist_ok=True)

    logger.info("combining read info")
    info = {}
    for sample in samples:
        info[sample] = pd.read_csv("input/{0}_info.csv".format(sample), header=0, index_col=0)
    info = (
        pd.concat(info.values(), axis=0, keys=info.keys())
        .reset_index(level=0, names="sample")
        .loc[reads]
        .to_csv("input/{0}_info.csv".format(args.combined))
    )

    logger.info("combining filter stats")
    filter_stats = {}
    for sample in samples:
        filter_stats[sample] = pd.read_csv(
            "output/{0}/{0}_filter_stats.csv".format(sample), header=None, index_col=0
        ).squeeze()

    pd.concat(filter_stats.values(), axis=0, keys=filter_stats.keys()).groupby(level=0).agg(
        sum
    ).to_csv("output/{0}/{0}_filter_stats.csv".format(args.combined), header=False)

    logger.info("averaging alignment pars")
    pars = []
    for sample in samples:
        pars.append(np.load("output/{0}/{0}_{1}_pars.npz".format(sample, args.aligner)))
    combined_pars = {}
    for slot in set.intersection(*(set(p.files) for p in pars)):
        combined_pars[slot] = np.mean([p[slot] for p in pars], axis=0)
    np.savez("output/{0}/{0}_{1}_pars.npz".format(args.combined, args.aligner), **combined_pars)

    logger.info("combining output from process_alignments")
    process = {}
    process_stats = {}
    breakpoint_alignments = {}
    for sample in samples:
        process[sample] = pd.read_csv(
            "output/{0}/{0}_processed.out".format(sample), header=0, index_col=0, sep="\t"
        )
        process_stats[sample] = pd.read_csv(
            "output/{0}/{0}_process_stats.csv".format(sample), header=None, index_col=0
        ).squeeze()
        breakpoint_alignments[sample] = pd.read_csv(
            "output/{0}/{0}_breakpoint_alignments.csv".format(sample), header=0, index_col=0
        )

    pd.concat(process.values(), axis=0, keys=process.keys()).reset_index(
        level=0, names="sample"
    ).loc[reads].to_csv("output/{0}/{0}_processed.out".format(args.combined), sep="\t", index=True)
    pd.concat(process_stats.values(), axis=0, keys=process_stats.keys()).groupby(level=0).agg(
        sum
    ).to_csv("output/{0}/{0}_process_stats.csv".format(args.combined), header=False)
    breakpoint_alignments = pd.concat(
        breakpoint_alignments.values(), axis=0, keys=breakpoint_alignments.keys()
    ).reset_index(level=0, names="sample")
    breakpoint_alignments.loc[reads.intersection(breakpoint_alignments.index)].to_csv(
        "output/{0}/{0}_breakpoint_alignments.csv".format(args.combined)
    )

    logger.info("combining inserts")
    inserts = {}
    bed = {}
    for sample in samples:
        try:
            inserts[sample] = pd.read_csv(
                "output/{0}/{0}_inserts.tsv".format(sample), index_col=0, header=0, sep="\t"
            )
        except (ValueError, pd.errors.EmptyDataError):
            inserts[sample] = pd.DataFrame(
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
        bed[sample] = pd.read_csv(
            "output/{0}/{0}_inserts.bed".format(sample), header=None, sep="\t"
        )

    inserts = pd.concat(inserts.values(), axis=0)
    insert_reads = reads.intersection(inserts.index)
    inserts = inserts.loc[insert_reads].reset_index(names='read')
    inserts.to_csv("output/{0}/{0}_inserts.tsv".format(args.combined), index=False, sep="\t")
    bed = (
        pd.concat(bed.values(), axis=0)
        .set_index(3)
        .loc[reads]
        .reset_index()[[0, 1, 2, "read", 4, 5, 6, 7, 8, 9, 10, 11]]
    )
    bed.to_csv(
        "output/{0}/{0}_inserts.bed".format(args.combined), header=False, index=False, sep="\t"
    )

    logger.info("writing combined MSA to output/{0}/{0}_msa.npz".format(args.combined))
    msa_df.loc[reads].to_csv("output/{0}/{0}_msa.csv".format(args.combined))
    scipy.sparse.save_npz(
        "output/{0}/{0}_msa.npz".format(args.combined),
        scipy.sparse.vstack(msa.values()).tocsr()[inds].tocoo(),
    )
