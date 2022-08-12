"""construct linkage from graph"""


def setup_argparse(parser):
    parser.add_argument(
        "--method",
        dest="method",
        default="average",
        help="""clustering method for hierarchical clustering [average]""",
    )
    parser.add_argument(
        "--linkage", dest="linkage", help="""output file contains linkage"""
    )
    parser.add_argument(
        "--graph",
        dest="graph",
        help="""input file containing kNN graph or distance matrix""",
    )
    parser.add_argument(
        "--ignore_unused_positions",
        dest="ignore_unused_positions",
        action="store_true",
        default=False,
        help="""ignore positions that have no coverage""",
    )
    parser.add_argument(
        "--gbbs_bin",
        dest="gbbs_bin",
        help="""path to directory with gbbs binaries, e.g., bazel-bin""",
    )


def write_gbbs_graph(knn_indices, knn_dists, outfile, min_dist=1.0e-10):

    import numpy as np

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
    # sort by node indices
    o = np.lexsort([jj, ii])
    edges = jj[o]
    weights = xx[o]
    # determine node offsets
    offsets = np.concatenate([[0], np.cumsum(np.bincount(ii))[:-1]])
    nedges = len(xx)

    with open(outfile, "w") as outf:
        outf.write("WeightedAdjacencyGraph\n{0}\n{1}\n".format(nnodes, nedges))
        outf.write("\n".join("{0}".format(o) for o in offsets) + "\n")
        outf.write("\n".join("{0}".format(e) for e in edges) + "\n")
        outf.write("\n".join("{0:.10g}".format(w) for w in weights) + "\n")


def convert_gbbs_linkage(linkage_file, max_dist=1):

    import numpy as np
    import pandas as pd
    import scipy.cluster.hierarchy

    Ls = pd.read_csv(
        linkage_file,
        index_col=0,
        header=None,
        names=["parent", "wgh"],
        delim_whitespace=True,
    )
    # find nodes with only one child
    num_children = Ls["parent"].value_counts()
    nonbinary = num_children.index[num_children == 1]
    if len(nonbinary) > 1:
        raise Exception(
            "convert_gbbs_linkage: not sure what to do if there is more than one node with only one child!"
        )
    elif len(nonbinary) == 1:
        for nb in nonbinary:
            child = Ls.index[Ls["parent"] == nb]
            if nb in Ls.index:
                # change parent node of child
                Ls.loc[child, "parent"] = Ls.loc[nb, "parent"]
                # remove this node
                Ls = Ls.drop(nb, axis=0)
            else:
                # if nb has no parent, drop child
                Ls = Ls.drop(child, axis=0)
        # adjust node indices
        ii = np.array(Ls.index)
        pp = np.array(Ls["parent"])
        ii[ii >= nb] = ii[ii >= nb] - len(nonbinary)
        pp[pp >= nb] = pp[pp >= nb] - len(nonbinary)
        Ls.index = ii
        Ls["parent"] = pp
    # sort by parent node index and by distance
    o = np.lexsort([Ls["wgh"].values, Ls["parent"].values])
    zz = []
    i = 0
    while i < len(o) - 1:
        if Ls["parent"][o[i]] == Ls["parent"][o[i + 1]]:
            zz.append([o[i] + 1, o[i + 1] + 1, min(Ls["wgh"][o[i]], max_dist)])
            i += 2
        else:
            # this should not happen
            i += 1
            raise Exception(
                "convert_gbbs_linkage: nodes with single children left in input array!"
            )
    return scipy.cluster.hierarchy.from_mlab_linkage(np.array(zz))


def run(args):

    import os
    import sys
    import numpy as np
    import tempfile
    from logzero import logger

    logger.info("reading graph from {0}".format(args.graph))
    tmp = np.load(args.graph)
    if "dist" in tmp:
        import fastcluster

        logger.info(
            "running hierarchical clustering of with {0} method".format(
                args.method
            )
        )
        Z = fastcluster.linkage(tmp["dist"], method=args.method)
    elif "knn_inds" in tmp and "knn_dists" in tmp:
        knn_inds = tmp["knn_inds"]
        knn_dists = tmp["knn_dists"]
        with tempfile.TemporaryDirectory() as tmpdir:
            graph_file = os.path.join(tmpdir, "graph.txt")
            linkage_file = os.path.join(tmpdir, "linkage.txt")
            write_gbbs_graph(knn_inds, knn_dists, graph_file)
            logger.info(
                "running gbbs HierarchicalAgglomerativeClustering with {0} linkage".format(
                    args.method
                )
            )
            os.system(
                args.gbbs_bin
                + "/benchmarks/Clustering/SeqHAC/HACDissimilarity -s -of {0} -linkage {1} {2} > /dev/null".format(
                    linkage_file, args.method, graph_file
                )
            )
            if not os.path.isfile(linkage_file):
                raise Exception(
                    "gbbs clustering failed - linkage file missing!\n"
                )
            Z = convert_gbbs_linkage(linkage_file)
    else:
        logger.warn("invalid input {0}".format(args.graph))
        sys.exit(1)

    logger.info("saving linkage to {0}\n".format(args.linkage))
    np.savez(args.linkage, Z=Z)
