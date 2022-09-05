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
        "--nmax", dest="nmax", type=int, help="""use only nmax reads"""
    )


def get_sparse_matrix_from_knn(knn_indices, knn_dists, min_dist=1.0e-10):

    import numpy as np
    import scipy.sparse

    nnodes, nneighbors = knn_indices.shape
    ii = np.repeat(np.arange(nnodes), nneighbors)
    jj = knn_indices.flatten("C")
    xx = knn_dists.flatten("C")
    ii, jj = np.minimum(ii, jj), np.maximum(ii, jj)
    # take unique non-diagonal non-negative elements
    take = np.intersect1d(
        np.unique(np.vstack([ii, jj]).T, axis=0, return_index=True)[1],
        np.where((ii < jj) & (ii >= 0))[0],
    )
    # make symmetric
    ii, jj = np.concatenate([ii[take], jj[take]]), np.concatenate(
        [jj[take], ii[take]]
    )
    xx = np.maximum(np.concatenate([xx[take], xx[take]]), min_dist)

    result = scipy.sparse.coo_matrix((xx, (ii, jj)), shape=(nnodes, nnodes))
    # result.eliminate_zeros()
    return result.tocsr()


def get_adjacency_graph(
    xx,
    n_neighbors=10,
    random_state=0,
    metric="jaccard",
    angular=False,
    n_jobs=-1,
    return_matrix=False,
):

    from umap.umap_ import nearest_neighbors

    knn_indices, knn_dists, forest = nearest_neighbors(
        xx,
        n_neighbors,
        random_state=random_state,
        metric=metric,
        metric_kwds={},
        angular=angular,
        n_jobs=n_jobs,
        verbose=False,
    )

    if return_matrix:
        return get_sparse_matrix_from_knn(knn_indices, knn_dists)
    else:
        return knn_indices, knn_dists


def run(args):

    import os
    import sys
    import numpy as np
    import scipy.sparse
    import scipy.spatial
    from logzero import logger
    from .helpers import remove_gaps

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
            return_matrix=False,
        )
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
