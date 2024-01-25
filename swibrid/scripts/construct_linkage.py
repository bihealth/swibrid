"""construct hiearchical agglomerative clustering from MSA"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("--gaps", dest="gaps", help="""output of get_gaps.py""")
    parser.add_argument("--distances", dest="distances", help="""distance matrix""")
    parser.add_argument(
        "--max_gap",
        dest="max_gap",
        default=75,
        type=int,
        help="""max gap size to ignore [75]""",
    )
    parser.add_argument(
        "--metric",
        dest="metric",
        default="jaccard",
        help="""clustering metric [jaccard]""",
    )
    parser.add_argument(
        "--method",
        dest="method",
        default="average",
        help="""clustering method for hierarchical clustering [average]""",
    )
    parser.add_argument("--linkage", dest="linkage", help="""output file contains linkage""")
    parser.add_argument(
        "--use_sparse",
        dest="use_sparse",
        action="store_true",
        default=False,
        help="""EXPERIMENTAL: use sparse nearest-neighbor clustering""",
    )
    parser.add_argument(
        "--n_neighbors",
        dest="n_neighbors",
        default=100,
        type=int,
        help="""# of nearest neighbors for adjacency graph [100]""",
    )
    parser.add_argument(
        "--n_threads",
        dest="n_threads",
        default=1,
        type=int,
        help="""# of threads [1]""",
    )
    parser.add_argument(
        "--n_backup",
        dest="n_backup",
        type=int,
        default=0,
        help="""save full distances for these points in addition to the nearest neighbors [0]""",
    )
    parser.add_argument("--nmax", dest="nmax", type=int, help="""use only nmax reads""")
    parser.add_argument(
        "--ignore_unused_positions",
        dest="ignore_unused_positions",
        action="store_true",
        default=False,
        help="""ignore positions that have no coverage""",
    )
    parser.add_argument(
        "--downsample_nreads",
        dest="downsample_nreads",
        type=int,
        default=0,
        help="""construct additional linkage from downsampling to xx reads [0]""",
    )
    parser.add_argument(
        "--downsample_nreps",
        dest="downsample_nreps",
        type=int,
        default=0,
        help="""construct additional linkage from downsampling xx times [0]""",
    )
    parser.add_argument(
        "--downsampled_linkage",
        dest="downsampled_linkage",
        help="""output file for downsampled linkage(s)""",
    )


def run(args):
    import os
    import sys
    import numpy as np
    import scipy.sparse
    import scipy.spatial
    from logzero import logger
    from .utils import remove_gaps

    if not os.path.isfile(args.msa):
        logger.warn("no msa at {0}; run construct_msa.py first!".format(args.msa))
        sys.exit(1)

    logger.info("loading msa from {0}".format(args.msa))
    msa = scipy.sparse.load_npz(args.msa)

    if args.nmax is not None and msa.shape[0] > args.nmax:
        logger.info("restricting msa to {0} reads".format(args.nmax))
        msa = msa.tocsr()[: args.nmax].tocoo()

    logger.info("loading gaps from {0}".format(args.gaps))
    gaps = np.load(args.gaps)

    logger.info("removing gaps > {0} from msa".format(args.max_gap))
    msa_cleaned = remove_gaps(msa, gaps=gaps, max_gap=args.max_gap)
    if args.ignore_unused_positions:
        used = msa_cleaned.sum(0) > 0
        logger.info("ignoring {0}/{1} unused positions".format(np.sum(~used), len(used)))
        msa_cleaned = msa_cleaned[:, used]

    if args.use_sparse:
        import sparsecluster

        logger.info(
            "running sparsecluster hierarchical clustering of {3} reads with {0} metric and {1} method for {2} nearest neighbors".format(
                args.metric, args.method, args.n_neighbors, msa_cleaned.shape[0]
            )
        )
        if args.distances is not None and os.path.isfile(args.distances):
            logger.info("loading distance matrix from {0}".format(args.distances))
            dist = scipy.sparse.load_npz(args.distances)
            Z = sparsecluster.linkage(
                msa_cleaned,
                dist=dist,
                metric=args.metric,
                method=args.method,
                n_neighbors=args.n_neighbors,
                n_backup=args.n_backup,
                n_jobs=args.n_threads,
                verbose=True,
            )
        else:
            Z = sparsecluster.linkage(
                msa_cleaned,
                metric=args.metric,
                method=args.method,
                n_neighbors=args.n_neighbors,
                n_backup=args.n_backup,
                n_jobs=args.n_threads,
                verbose=True,
            )

    else:
        import fastcluster

        logger.info(
            "running fastcluster hierarchical clustering of {2} reads with {0} metric and {1} method".format(
                args.metric, args.method, msa_cleaned.shape[0]
            )
        )
        if args.distances is not None and os.path.isfile(args.distances):
            logger.info("loading distance matrix from {0}".format(args.distances))
            dist = np.load(args.distances)["dist"]
            Z = fastcluster.linkage(dist, metric=args.metric, method=args.method)
        else:
            Z = fastcluster.linkage(msa_cleaned.todense(), metric=args.metric, method=args.method)

    logger.info("saving linkage to {0}".format(args.linkage))
    np.savez(args.linkage, Z=Z)

    if (
        args.downsample_nreads > 0
        and args.downsample_nreps > 0
        and args.downsample_nreads < msa_cleaned.shape[0]
        and args.downsampled_linkage
    ):
        import fastcluster

        logger.info(
            "running {0} replicates of linkage construction for {1} downsampled reads".format(
                args.downsample_nreps, args.downsample_nreads
            )
        )

        Z = {}
        for n in range(args.downsample_nreps):
            ii = np.random.choice(msa_cleaned.shape[0], args.downsample_nreads, replace=False)
            msa_here = msa_cleaned[ii]
            logger.info(
                "rep {3}: running fastcluster hierarchical clustering of {2} reads with {0} metric and {1} method".format(
                    args.metric, args.method, msa_here.shape[0], n
                )
            )
            Z["Z" + str(n)] = fastcluster.linkage(
                msa_here.todense(), metric=args.metric, method=args.method
            )

        logger.info("saving downsampled linkages to {0}".format(args.downsampled_linkage))
        np.savez(args.downsampled_linkage, **Z)
