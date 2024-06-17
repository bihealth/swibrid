"""\
downsample the clustering to get more robust diversity measures
given the input MSA (and associated gaps), clustering is repeated `nreps` times on a subsample of `nreads` reads
using `fastcluster`  with (by defalt) cosine metric and average linkage and a cutoff derived from cluster_stats input
the following diversity measured are calculated:
- mean_cluster_size (as fraction of reads)
- std_cluster_size
- nclusters_final (number of clusters after filtering)
- nclusters_eff (from the entropy of the cluster size distribution pre-filtering)
- cluster_gini: gini coefficient (post-filtering)
- cluster_entropy: entropy (post-filtering)
- cluster_inverse_simpson: inverse simpson coefficient (post-filtering)
- top_clone_occupancy: relative fraction of reads in biggest cluster
- big_clones_occupancy: fraction of reads in clusters > 1% occupancy
"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""required: file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument(
        "--cluster_stats",
        dest="cluster_stats",
        help="""required: file with statistics from original clustering""",
    )
    parser.add_argument(
        "-s",
        "--stats",
        dest="stats",
        help="""required: output file with clustering statistics after downsampling""",
    )
    parser.add_argument(
        "--nreads", dest="nreads", type=int, default=1000, help="""number of reads used [1000]"""
    )
    parser.add_argument(
        "--nreps", dest="nreps", type=int, default=10, help="""number of replicates [10]"""
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
        default="cosine",
        help="""clustering metric [cosine]""",
    )
    parser.add_argument(
        "--method",
        dest="method",
        default="average",
        help="""clustering method for hierarchical clustering [average]""",
    )
    parser.add_argument(
        "--ignore_unused_positions",
        dest="ignore_unused_positions",
        action="store_true",
        default=False,
        help="""ignore positions that have no coverage""",
    )
    parser.add_argument(
        "--filtering_cutoff",
        dest="filtering_cutoff",
        default=0.95,
        type=float,
        help="""filter out smallest clusters, keeping at least this fraction of reads [.95]""",
    )
    parser.add_argument(
        "--big_clone_cutoff",
        dest="big_clone_cutoff",
        default=0.01,
        type=float,
        help="""cutoff to determine a "big" clone as fraction of clustered reads [.01]""",
    )


def run(args):
    import os
    import sys
    import numpy as np
    import pandas as pd
    import scipy.sparse
    import scipy.cluster.hierarchy
    import scipy.stats
    from logzero import logger
    from collections import defaultdict
    import fastcluster
    from .utils import remove_gaps, filter_clustering, calculate_gini

    if not os.path.isfile(args.msa):
        logger.warn("no msa at {0}; run construct_msa.py first!".format(args.msa))
        sys.exit(1)

    logger.info("loading msa from {0}".format(args.msa))
    msa = scipy.sparse.load_npz(args.msa)

    logger.info("loading gaps from {0}".format(args.gaps))
    gaps = np.load(args.gaps)

    logger.info("removing gaps > {0} from msa".format(args.max_gap))
    msa_cleaned = remove_gaps(msa, gaps=gaps, max_gap=args.max_gap)
    if args.ignore_unused_positions:
        used = msa_cleaned.sum(0) > 0
        logger.info("ignoring {0}/{1} unused positions".format(np.sum(~used), len(used)))
        msa_cleaned = msa_cleaned[:, used]

    logger.info("reading original cluster stats file {0}".format(args.cluster_stats))
    cluster_stats = pd.read_csv(args.cluster_stats, index_col=0, header=None).squeeze()
    cutoff = cluster_stats["c_opt"]

    if args.nreads > msa_cleaned.shape[0]:
        nreps = 1
        nreads = msa_cleaned.shape[0]
    else:
        nreps = args.nreps
        nreads = args.nreads

    logger.info(
        "running {0} replicates of linkage ({2} metric, {3} method), clustering (cutoff: {4:.2g}) and diversity analysis for {1}/{5} downsampled reads".format(
            nreps, nreads, args.metric, args.method, cutoff, msa_cleaned.shape[0]
        )
    )

    stats = defaultdict(list)
    for n in range(nreps):
        logger.info("rep {0}".format(n + 1))
        ii = np.random.choice(msa_cleaned.shape[0], nreads, replace=False)
        msa_here = msa_cleaned[ii]
        Z = fastcluster.linkage(msa_here.todense(), metric=args.metric, method=args.method)
        cc = scipy.cluster.hierarchy.cut_tree(Z, height=[cutoff])
        clustering = filter_clustering(Z, cc[:, 0].flatten(), p=args.filtering_cutoff)
        clusters, cinv, csize = np.unique(clustering, return_inverse=True, return_counts=True)
        use = clusters >= 0
        rel_size = csize[use] / csize[use].sum()
        if (rel_size.sum() > 0.999) | (rel_size.sum() < 1.001):
            continue

        stats["mean_cluster_size_downsampled"].append(rel_size.mean())
        stats["std_cluster_size_downsampled"].append(rel_size.std())
        stats["nclusters_final_downsampled"].append(len(rel_size))
        stats["nclusters_eff_downsampled"] = np.exp(scipy.stats.entropy(csize))
        stats["cluster_gini_downsampled"].append(calculate_gini(rel_size))
        stats["cluster_entropy_downsampled"].append(
            scipy.stats.entropy(rel_size) / np.log(len(rel_size))
        )
        stats["cluster_inverse_simpson_downsampled"].append(1.0 / (rel_size**2).sum())
        stats["top_clone_occupancy_downsampled"].append(rel_size.max())
        stats["big_clones_occupancy_downsampled"].append(
            rel_size[rel_size > args.big_clone_cutoff].sum()
        )

    if args.stats:
        logger.info("saving results to {0}".format(args.stats))
        pd.DataFrame.from_dict(stats).to_csv(args.stats)
