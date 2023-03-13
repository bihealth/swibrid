"""get sample summary and plot"""


def setup_argparse(parser):

    parser.add_argument("--sample", dest="sample", help="""sample""")
    parser.add_argument("--figure", dest="figure", help="""output figure""")
    parser.add_argument("--stats", dest="stats", help="""output stats""")
    parser.add_argument(
        "--filter", dest="filter", help="""stats file from filter_reads.py"""
    )
    parser.add_argument(
        "--process",
        dest="process",
        help="""stats file from process_last_output.py""",
    )
    parser.add_argument("--info", dest="info", help="""read info""")
    parser.add_argument("--gaps", dest="gaps", help="""gap distribution""")
    parser.add_argument(
        "--gap_stats", dest="gap_stats", help="""gap statistics"""
    )
    parser.add_argument(
        "--max_gap",
        dest="max_gap",
        default=75,
        type=int,
        help="""max gap size to ignore [75]""",
    )
    parser.add_argument(
        "--clustering",
        dest="clustering",
        help="""find_clusters output file""",
    )
    parser.add_argument(
        "--extrapolation",
        dest="extrapolation",
        help="""find_clusters.py cluster extrapolation""",
    )
    parser.add_argument(
        "--cluster_stats",
        dest="cluster_stats",
        help="""find_clusters.py stats""",
    )
    parser.add_argument(
        "--clustering_analysis",
        dest="clustering_analysis",
        help="""analyze_clusters.py output""",
    )
    parser.add_argument(
        "--mutations",
        dest="mutations",
        help="""file with mutation info (from find_mutations.py)""",
    )
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000:-",
        help="""coordinates of switch region [chr14:106050000-106337000:-]""",
    )
    parser.add_argument(
        "--switch_annotation",
        dest="switch_annotation",
        help="""bed file with switch annotation""",
    )
    parser.add_argument(
        "--use_clones",
        dest="use_clones",
        help="""comma-separated list of clones to use""",
    )
    parser.add_argument(
        "--use_weights",
        dest="use_weights",
        default='inverse',
        help="""use different cluster weights ("inverse" | "uniform" | "adjusted") [inverse]"""
    )


def run(args):

    import os
    import re
    import numpy as np
    import scipy.stats
    import pandas as pd
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    import seaborn as sns

    from logzero import logger

    from .utils import (
        parse_switch_coords,
        get_switch_coverage,
        decode_coords,
        read_switch_anno,
        merge_intervals,
        shift_coord,
        calculate_gini,
        weighted_avg_and_std,
        ncodes,
        f2,
    )

    matplotlib.rcParams.update({"font.size": 8})
    matplotlib.rcParams.update({"axes.linewidth": 0.5})
    matplotlib.rcParams.update({"xtick.major.width": 0.5})
    matplotlib.rcParams.update({"ytick.major.width": 0.5})

    (
        switch_chrom,
        switch_start,
        switch_end,
        switch_orientation,
    ) = parse_switch_coords(args.switch_coords)
    switch_anno = read_switch_anno(args.switch_annotation)
    cov_int, Ltot, eff_start, eff_end, anno_recs = get_switch_coverage(
        switch_anno, switch_chrom, switch_start, switch_end
    )

    stats = {}
    logger.info("reading filter stats")
    stats["filter"] = pd.read_csv(args.filter, index_col=0, header=0).squeeze()
    logger.info("reading processing stats")
    stats["process"] = (
        pd.read_csv(args.process, index_col=0, header=0)
        .squeeze()
        .rename({"nreads": "mapped"})
    )
    logger.info("reading clustering, clustering stats, clustering analysis")
    clustering = pd.read_csv(args.clustering, index_col=0, header=0)
    clustering_analysis = pd.read_csv(
        args.clustering_analysis, index_col=0, header=0
    )
    stats["cluster"] = pd.read_csv(
        args.cluster_stats, index_col=0, header=0
    ).squeeze()
    stats = pd.concat(stats.values(), axis=0)

    neff = stats["eff_nclusters"].astype(int)
    if args.use_clones:
        clones = list(map(int, args.use_clones.split(",")))
    else:
        clones = clustering_analysis.sort_values("size", ascending=False).index[
            :neff
        ]
    logger.info("using {0} clones".format(len(clones)))

    if args.use_weights == 'inverse':
        logger.info('using inverse weights')
        w = np.ones(len(clones))
    elif args.use_weights == 'uniform':
        logger.info('using uniform weights')
        w = clustering_analysis.loc[clones, 'size']
    elif args.use_weights == 'adjusted':
        logger.info('using adjusted weights')
        w = clustering_analysis.loc[clones, 'adj_size']

    stats["nclusters_used"] = len(clones)
    stats["clustered"] = sum(~clustering["cluster"].isna())
    stats["mean_length"], stats["std_length"] = weighted_avg_and_std (clustering_analysis.loc[clones, "length"], w)
    stats["mean_GC"], stats["std_GC"] = weighted_avg_and_std (clustering_analysis.loc[clones, "GC"], w)
    stats["mean_cluster_size"] = clustering_analysis.loc[clones, "size"].mean()
    stats["std_cluster_size"] = clustering_analysis.loc[clones, "size"].std()
    stats["mean_adj_cluster_size"] = clustering_analysis.loc[clones, "adj_size"].mean()
    stats["std_adj_cluster_size"] = clustering_analysis.loc[clones, "adj_size"].std()

    stats["cluster_gini"] = calculate_gini(
        clustering["cluster"][clustering["cluster"].isin(clones)],
    )
    stats["cluster_entropy"] = scipy.stats.entropy(
        clustering_analysis.loc[clones, "size"]
    ) / np.log(len(clones))
    stats["cluster_inverse_simpson"] = 1. / ((clustering_analysis.loc[clones, "size"] / clustering_analysis.loc[clones, "size"].sum())**2).sum()

    stats["PCR_length_bias"] = scipy.stats.linregress(
        clustering_analysis.loc[clones, "length"],
        np.log(clustering_analysis.loc[clones, "size"]),
    )[0]

    stats["PCR_GC_bias"] = scipy.stats.linregress(
        clustering_analysis.loc[clones, "GC"],
        np.log(clustering_analysis.loc[clones, "size"]),
    )[0]

    logger.info("reading read info")
    read_info = pd.read_csv(args.info, header=0, index_col=0)

    barcode_locs = (
        [
            (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2]))
            / (2 * l)
            for bc, l in read_info[["barcodes", "length"]].values
            for x in bc.split(";")
        ]
        if "barcodes" in read_info.columns
        else []
    )

    primer_locs = (
        [
            (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2]))
            / (2 * l)
            for pr, l in read_info[["primers", "length"]].dropna().values
            for x in pr.split(";")
        ]
        if "primers" in read_info.columns
        else []
    )

    internal_primer_locs = (
        [
            (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2]))
            / (2 * l)
            for pr, l in read_info[["internal", "length"]].dropna().values
            if pr.count("primer") > 0
            for x in pr.split(";")
            if "primer" in x
        ]
        if "internal" in read_info.columns
        else []
    )

    internal_barcode_locs = (
        [
            (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2]))
            / (2 * l)
            for bc, l in read_info[["internal", "length"]].dropna().values
            for x in bc.split(";")
            if "BC" in x
        ]
        if "internal" in read_info.columns
        else []
    )

    read_isotype_count = clustering["isotype"].value_counts()
    cluster_isotype_count = clustering_analysis.loc[
        clones, "isotype"
    ].value_counts()
    insert_isotype_count = (
        clustering_analysis["insert_pos_isotype"].dropna().value_counts()
    )

    nreads = read_isotype_count.sum()
    nclusters = cluster_isotype_count.sum()
    isotype_fracs = pd.DataFrame(
        {
            "reads\n(n={0})".format(nreads): read_isotype_count / nreads,
            "clusters\n(n={0})".format(nclusters): cluster_isotype_count
            / nclusters,
        }
    )

    read_isotype_count.index = "nreads_" + read_isotype_count.index
    cluster_isotype_count.index = "nclusters_" + cluster_isotype_count.index
    insert_isotype_count.index = (
        "ninserts_" + insert_isotype_count.index.astype(str)
    )
    take = ~clustering["cluster"].isna()
    inserts = [
        decode_coords(m)
        for insert in clustering[take]["insert"].dropna()
        for m in insert.split(",")
    ]
    cluster_stats = {
        "n_inserts": len(inserts),
        "n_unique_inserts": len(merge_intervals(inserts)),
        "n_clusters_inserts": len(
            clustering_analysis["insert_pos_isotype"].dropna()
        ),
    }
    insert_stats = clustering_analysis.agg(
        {
            "insert_frequency": "mean",
            "insert_overlap": "mean",
            "insert_pos_overlap": "mean",
            "insert_length": "mean",
            "insert_gap_length": "mean",
        },
        axis=0,
    )

    stats = pd.concat(
        [
            stats,
            read_isotype_count,
            cluster_isotype_count,
            pd.Series(cluster_stats),
            insert_stats,
            insert_isotype_count,
        ],
        axis=0,
    )

    extrapolation = pd.read_csv(
        args.extrapolation, header=None, index_col=0
    ).squeeze()

    if args.gap_stats and os.path.isfile(args.gap_stats):
        logger.info("loading gap stats from " + args.gap_stats)
        gap_stats = pd.read_csv(
            args.gap_stats, header=None, index_col=0
        ).squeeze()
        stats = pd.concat([stats, gap_stats], axis=0)

    if args.mutations and os.path.isfile(args.mutations):

        logger.info("loading mutations from " + args.mutations)
        mutations = pd.read_csv(args.mutations, sep="\t", header=0)

        ref = []
        alt = []
        padj = []
        take = (mutations["strand_bias"] - 0.5).abs() < 0.25
        for _, row in mutations[take].iterrows():
            r = ncodes[row["ref"]] - 1
            ref.append(r)
            c = list(map(int, row["alt_counts"].split(",")))
            c[r] = 0
            alt.append(np.argmax(c))
            padj.append(row["padj_clust"])

        mm = pd.crosstab(ref,alt).reindex(index=range(4),columns=range(4)).fillna(0).astype(int).values

        mut_stats = pd.Series(
            {
                "num_mutated_pos": len(ref),
                "fraction_cluster_specific_mutations": np.mean(
                    np.array(padj) < 0.05
                ),
                "fraction_transitions": mm[(0, 2, 1, 3), (2, 0, 3, 1)].sum()
                / len(ref),
                "fraction_C>A": mm[(1,2), (0,3)].sum() / len(ref),
                "fraction_C>G": mm[(1,2), (2,1)].sum() / len(ref),
                "fraction_C>T": mm[(1,2), (3,0)].sum() / len(ref),
                "fraction_T>A": mm[(3,0), (0,3)].sum() / len(ref),
                "fraction_T>C": mm[(3,0), (1,2)].sum() / len(ref),
                "fraction_T>G": mm[(3,0), (2,1)].sum() / len(ref),
            }
        )

        stats = pd.concat([stats, mut_stats], axis=0)

    if args.stats is not None:
        logger.info("saving summary to " + args.stats)
        stats.to_csv(args.stats)

    if args.figure is not None:

        fig, axs = plt.subplots(4, 2, figsize=(6, 8))
        fig.subplots_adjust(hspace=0.5, wspace=0.4, bottom=0.1, top=0.95)

        fig.text(0.5, 0.98, args.sample, size="medium", ha="center")

        ax = axs[0, 0]
        ax.hist(
            read_info["length"],
            bins=np.geomspace(
                read_info["length"].min(), read_info["length"].max(), 100
            ),
            histtype="step",
        )
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylabel("# reads")
        ax.set_xlabel("read length")
        sns.despine(ax=ax)

        ax = axs[0, 1]
        ax.hist(
            [barcode_locs, primer_locs],
            label=["barcode", "primer"],
            bins=np.linspace(0, 1, 101),
            histtype="step",
            color=("#1f77b4", "#ff7f0e"),
        )
        ax.hist(
            [internal_barcode_locs, internal_primer_locs],
            bins=np.linspace(0, 1, 101),
            histtype="step",
            color=("#1f77b4", "#ff7f0e"),
            linestyle=("dashed"),
        )
        ax.legend(loc=9, handlelength=1, frameon=False, prop={"size": "small"})
        ax.set_yscale("log")
        ax.set_ylabel("# reads")
        ax.set_xlabel("relative position")
        sns.despine(ax=ax)

        ax = axs[1, 0]
        isotype_fracs.T.plot(kind="barh", stacked=True, ax=ax)
        ax.set_xlabel("fraction")
        ax.legend(
            ncol=3,
            handlelength=1,
            columnspacing=0.5,
            loc=2,
            frameon=False,
            prop={"size": "small"},
        )
        ax.set_ylim([-0.5, 2])
        sns.despine(ax=ax)

        ax = axs[1, 1]
        cc = extrapolation.index
        nc = extrapolation.values
        ax.plot(cc, nc, ".", color="#1f77b4", markersize=1)

        if "p0" in stats.index:
            pars = stats[["p0", "p1", "p2", "p3"]].values
            ax.plot(cc, f2(pars, cc), "-", color="k", lw=1)
            ax.plot(
                cc, f2((pars[0], pars[1], 0, 0), cc), "--", color="k", lw=0.5
            )
        ax.vlines(
            stats["c_opt"],
            0,
            extrapolation.loc[stats["c_opt"]],
            color="r",
            linestyle="dashed",
            lw=0.5,
        )
        ax.hlines(
            extrapolation.loc[stats["c_opt"]],
            0,
            stats["c_opt"],
            color="r",
            linestyle="dashed",
            lw=0.5,
        )
        ax.set_yscale("log")
        ax.set_xlabel("cutoff")
        ax.set_ylabel("# clusters")
        ax.set_ylim([max(1, ax.get_ylim()[0]), None])
        sns.despine(ax=ax)

        major_ticks = []
        minor_ticks = []
        minor_labels = []
        for rec in anno_recs:
            start = shift_coord(int(rec[3][1]), cov_int) - eff_start
            end = shift_coord(int(rec[3][2]), cov_int) - eff_start
            major_ticks += [start, end]
            minor_ticks.append((start + end) / 2)
            minor_labels.append(rec[3][3])

        ax = axs[2, 0]
        ax.hist(
            clustering["cluster"].dropna().value_counts(),
            bins=np.geomspace(1, stats["clustered"], 50),
            histtype="step",
            label="all",
        )
        ax.hist(
            clustering["cluster"].dropna().value_counts()[:neff],
            bins=np.geomspace(1, stats["clustered"], 50),
            histtype="stepfilled",
            label="top 95%",
        )
        ax.hist(
            clustering_analysis["adj_size"].dropna(), 
            bins=np.geomspace(1, stats["clustered"], 50),
            histtype="stepfilled",
            label="adjusted",
        )
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylabel("# clusters")
        ax.set_xlabel("cluster size")
        ax.legend(loc=1, handlelength=1, prop={"size": "small"}, frameon=False)
        sns.despine(ax=ax)

        ax = axs[2, 1]
        ax.plot(
            clustering_analysis["size"],
            clustering_analysis["break_spread"],
            ".",
        )
        ax.set_xscale("log")
        ax.set_ylabel("break spread")
        ax.set_xlabel("cluster size")
        sns.despine(ax=ax)

        logger.info("loading gaps from " + args.gaps)
        gaps = np.load(args.gaps)

        ax = axs[3, 0]
        ax.hist(
            gaps["gap_size"], bins=np.geomspace(1, 3.0e4, 50), histtype="step"
        )
        ax.axvline(args.max_gap, lw=0.5, color="r", ls="dashed")
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylabel("# events")
        ax.set_xlabel("gap size")
        sns.despine(ax=ax)

        ax = axs[3, 1]
        use = gaps["gap_size"] > args.max_gap
        ax.hist(
            [gaps["pos_right"][use], gaps["pos_left"][use]],
            label=["left", "right"],
            # bins=np.arange(Ltot),
            bins=np.linspace(0, Ltot, 100),
            histtype="step",
        )
        ax.set_xlim([Ltot, 0])
        ax.set_yscale("log")
        ax.set_ylabel("# events")
        ax.set_xlabel("breakpoint position")
        ax.set_xticks(np.unique(major_ticks))
        ax.set_xticklabels([])
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_xticklabels(minor_labels, minor=True, size="small")
        ax.tick_params(which="minor", length=0)
        ax.legend(
            loc=9, handlelength=1, frameon=False, prop={"size": "small"}, ncol=2
        )
        sns.despine(ax=ax)

        fig.savefig(args.figure)
