"""construct read-read graph from MSA"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("--gaps", dest="gaps", help="""output of get_gaps.py""")
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
        "--use_graph",
        dest="use_graph",
        action="store_true",
        default=False,
        help="""use nearest-neighbor graph""",
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
        "--graph",
        dest="graph",
        help="""output file contains adjacency graph or distance matrix""",
    )
    parser.add_argument(
        "--ignore_unused_positions",
        dest="ignore_unused_positions",
        action="store_true",
        default=False,
        help="""ignore positions that have no coverage""",
    )
    parser.add_argument(
        "--nmax_full",
        dest="nmax_full",
        type=int,
        default=50000,
        help="""switch to graph-based clustering if more than nmax_full reads [50000]""",
    )
    parser.add_argument(
        "--n_spike",
        dest="n_spike",
        type=int,
        default=0,
        help="""compute full distances for these points in addition to the nearest neighbors [0]"""
    )
    parser.add_argument(
        "--nmax", dest="nmax", type=int, help="""use only nmax reads"""
    )


def get_adjacency_graph(
    X,
    n_neighbors=10,
    random_state=0,
    metric="jaccard",
    n_jobs=-1
    ):

    # taken and modified from umap.umap_.nearest_neighbors

    from pynndescent import NNDescent
    import numpy as np

    # TODO: Hacked values for now
    n_trees = min(64, 5 + int(round((X.shape[0]) ** 0.5 / 20.0)))
    n_iters = max(5, int(round(np.log2(X.shape[0]))))

    knn_search_index = NNDescent(
        X,
        n_neighbors=n_neighbors,
        metric=metric,
        metric_kwds={},
        random_state=random_state,
        n_trees=n_trees,
        n_iters=n_iters,
        max_candidates=60,
        low_memory=True,
        n_jobs=n_jobs,
        verbose=False,
        compressed=False,
    )
    knn_indices, knn_dists = knn_search_index.neighbor_graph

    return knn_indices, knn_dists


def run(args):

    import os
    import sys
    import numpy as np
    import scipy.sparse
    import scipy.spatial
    from logzero import logger
    from .utils import remove_gaps

    if not os.path.isfile(args.msa):
        logger.warn(
            "no msa at {0}; run construct_msa.py first!".format(args.msa)
        )
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
        logger.info(
            "ignoring {0}/{1} unused positions".format(np.sum(~used), len(used))
        )
        msa_cleaned = msa_cleaned[:, used]

    args.use_graph = args.use_graph or msa_cleaned.shape[0] > args.nmax_full

    if args.use_graph:
        logger.info(
            "constructing nearest neighbor graph of {2} reads with {0} neighbors and {1} metric".format(
                args.n_neighbors, args.metric, msa_cleaned.shape[0]
            )
        )
        knn_inds, knn_dists = get_adjacency_graph(
            msa_cleaned,
            n_neighbors=args.n_neighbors,
            n_jobs=args.n_threads,
            metric=args.metric,
        )

        if args.n_spike > 0:

            logger.info(
                "computing full distance matrix for {0} additional reads".format(args.n_spike)
                )

            ii = np.random.choice(msa_cleaned.shape[0], args.n_spike)
            n = len(ii)
            dd = scipy.spatial.distance.pdist(msa_cleaned[ii], metric=args.metric)
            sinds = -1 * np.ones((msa_cleaned.shape[0],n)).astype(int)
            sinds[ii] = np.tile(ii, n).reshape(n, n)
            sdists = np.zeros((msa_cleaned.shape[0], n))
            sdists[ii] = scipy.spatial.distance.squareform(dd)
            
            knn_inds = np.hstack([knn_inds, sinds])
            knn_dists = np.hstack([knn_dists, sdists])

        logger.info("saving nearest neighbor graph to {0}".format(args.graph))
        np.savez(args.graph, knn_inds=knn_inds, knn_dists=knn_dists)
    else:
        logger.info(
            "constructing distance matrix of {0} reads with {1} metric".format(
                msa_cleaned.shape[0], args.metric
            )
        )
        dist = scipy.spatial.distance.pdist(msa_cleaned, metric=args.metric)
        logger.info("saving distance matrix to {0}".format(args.graph))
        np.savez(args.graph, dist=dist)
