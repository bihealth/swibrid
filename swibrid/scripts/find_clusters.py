"""find clusters by cutting linkage dendrogram"""


def setup_argparse(parser):
    parser.add_argument(
        "-l",
        "--linkage",
        dest="linkage",
        help="""output of construct_linkage.py""",
    )
    parser.add_argument(
        "-f",
        "--fix_cutoff",
        dest="cutoff",
        type=float,
        help="""use fixed cutoff instead of data-derived""",
    )
    parser.add_argument(
        "-c",
        dest="cmin",
        default=0.001,
        type=float,
        help="""minimum cutoff value [.001]""",
    )
    parser.add_argument(
        "-C",
        dest="cmax",
        default=0.1,
        type=float,
        help="""maximum cutoff value [.1]""",
    )
    parser.add_argument(
        "-n",
        dest="nc",
        default=100,
        type=int,
        help="""number of cutoff values between cmin and cmax [100]""",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="""read info (output of construct_msa.py)""",
    )
    parser.add_argument("-o", "--output", dest="output", help="""output file with clustering""")
    parser.add_argument(
        "--scanning",
        dest="scanning",
        help="""file with scanning data""",
    )
    parser.add_argument("-s", "--stats", dest="stats", help="""file with statistics""")
    parser.add_argument(
        "--fit-method",
        dest="fit_method",
        default="distance",
        help="""cutoff determination method: "trend" or "distance" [distance]""",
    )


def fit_cutoff(cc, nn, method="fit"):
    import scipy
    import pandas as pd
    import numpy as np

    cc = np.asarray(cc)
    nn = np.asarray(nn)
    cmin = np.min(cc)
    cmax = np.max(cc)
    nmin = np.min(nn)
    nmax = np.max(nn)
    if method == "trend":
        from .utils import res2

        # p0 = [nmax, 0, nmin, np.log(nmax / nmin) / (cmax - cmin)]
        p0 = [0.75 * (nmax - nmin), 10, 0.25 * (nmax - nmin), 200]
        opt2 = scipy.optimize.least_squares(
            res2,
            1.0e-8 + np.array(p0),
            args=(cc, nn),
            method="dogbox",
            bounds=(0, np.inf),
        )
        if opt2.x[1] > opt2.x[3]:
            res = [opt2.x[2], opt2.x[3], opt2.x[0], opt2.x[1]]
        else:
            res = list(opt2.x)
        diff = pd.Series(np.abs(nn - res[0]), index=cc)
        cutoff = diff.sort_values().index[0]
        return res, cutoff
    elif method == "distance":
        diff = pd.Series(
            np.abs((cmax - cmin) * (nmax - nn) - (cmin - cc) * (nmin - nmax))
            / np.sqrt((cmax - cmin) ** 2 + (nmax - nmin) ** 2),
            index=cc,
        )
        cutoff = diff.sort_values().index[-1]
        return [], cutoff


def run(args):

    import numpy as np
    import pandas as pd
    import scipy.sparse
    import scipy.stats
    import scipy.cluster.hierarchy
    import scipy.optimize
    from logzero import logger
    from .utils import filter_clustering

    cutoffs = np.linspace(args.cmin, args.cmax, args.nc)
    if args.cutoff is not None:
        cutoffs = np.sort(np.unique(np.concatenate([cutoffs, [args.cutoff]])))
    ncutoffs = len(cutoffs)

    logger.info("loading linkage from {0}".format(args.linkage))
    Z = np.load(args.linkage)["Z"]

    logger.info("scanning {0} cutoffs on linkage".format(ncutoffs))
    cc = scipy.cluster.hierarchy.cut_tree(Z, height=cutoffs)
    nclust = np.max(cc, 0) + 1
    entropy = np.array(
        [scipy.stats.entropy(np.bincount(cc[:, i])) / np.log(nclust[i]) for i in range(ncutoffs)]
    )

    if args.scanning:
        logger.info("saving scanning data to {0}".format(args.scanning))
        pd.DataFrame({"nclusters": nclust, "entropy": entropy}, index=cutoffs).to_csv(
            args.scanning, header=True, index=True
        )

    if args.cutoff is None:

        logger.info("finding optimal cutoff")
        res, c_opt = fit_cutoff(cutoffs, nclust, method=args.fit_method)

    else:

        logger.info("fixing cutoff at {0}".format(args.cutoff))
        res, c_opt = [], args.cutoff

    clustering = cc[:, cutoffs == c_opt].flatten()

    logger.info("filtering clusters")
    filtered_clustering = filter_clustering(Z, clustering)

    if args.stats:
        logger.info("saving stats to {0}".format(args.stats))
        df = pd.Series(
            res
            + [
                c_opt,
                len(np.unique(clustering)),
                sum(np.unique(filtered_clustering) >= 0),
                np.sum(np.bincount(clustering) > 1),
                np.mean(np.bincount(clustering)[np.unique(clustering)] == 1),
            ],
            index=(["p0", "p1", "p2", "p3"] if len(res) == 4 else [])
            + [
                "c_opt",
                "nclusters",
                "eff_nclusters",
                "nclusters_multi",
                "frac_singletons",
            ],
        )
        df.to_csv(args.stats)

    if args.input:
        logger.info("reading read info from {0}".format(args.input))
        df = pd.read_csv(args.input, header=0, index_col=0)
        if len(clustering) == df.shape[0]:
            df["cluster"] = clustering
            df["filtered_cluster"] = filtered_clustering
        else:
            df["cluster"] = np.concatenate([clustering, [np.nan] * (df.shape[0] - len(clustering))])
            df["filtered_cluster"] = np.concatenate(
                [filtered_clustering, [np.nan] * (df.shape[0] - len(filtered_clustering))]
            )
    else:
        df = pd.DataFrame(dict(cluster=clustering, filtered_cluster=filtered_clustering))

    if args.output:
        logger.info("saving read info to {0}".format(args.output))
        df.to_csv(args.output)
