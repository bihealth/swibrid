"""\
produce a summary of features derived from a sample and a plot.
collects statistics produced by process_alignments, find_clusters, and downsample_clustering.
produces cluster distribution statistics similar to downsample_clustering (the latter are averaged):

- mean_cluster_size (as fraction of reads)
- std_cluster_size
- nclusters_final (number of clusters after filtering)
- nclusters_eff (from the entropy of the cluster size distribution pre-filtering)
- cluster_gini: gini coefficient (post-filtering)
- cluster_entropy: entropy (post-filtering)
- cluster_inverse_simpson: inverse simpson coefficient (post-filtering)
- top_clone_occupancy: relative fraction of reads in biggest cluster
- big_clones_occupancy: fraction of reads in clusters > 1% occupancy
- size_length_bias: regression coefficient of log(cluster_size) ~ length
- size_GC_bias: regression coefficient of log(cluster_size) ~ GC

averages cluster-specific features from analyze_clustering and get_breakpoint_stats over clusters.
collects number of reads / clusters per isotype.
gets statistics on variants (germline vs. somatic, transitions vs. transversions, etc.)
"""


def setup_argparse(parser):
    parser.add_argument("--sample", dest="sample", help="""sample name""")
    parser.add_argument("--figure", dest="figure", help="""output figure""")
    parser.add_argument("--stats", dest="stats", help="""output stats""")
    parser.add_argument(
        "--process",
        dest="process",
        help="""stats file from process_alignments""",
    )
    parser.add_argument("--info", dest="info", help="""read info""")
    parser.add_argument("--gaps", dest="gaps", help="""gap distribution""")
    parser.add_argument(
        "--breakpoint_stats", dest="breakpoint_stats", help="""breakpoint statistics"""
    )
    parser.add_argument(
        "--max_gap",
        dest="max_gap",
        default=75,
        type=int,
        help="""max gap size to ignore [%(default)d]""",
    )
    parser.add_argument(
        "--clustering",
        dest="clustering",
        help="""find_clusters output file""",
    )
    parser.add_argument(
        "--scanning",
        dest="scanning",
        help="""find_clusters dendrogram scanning""",
    )
    parser.add_argument(
        "--cluster_stats",
        dest="cluster_stats",
        help="""find_clusters stats""",
    )
    parser.add_argument(
        "--cluster_analysis",
        dest="cluster_analysis",
        help="""analyze_clusters output""",
    )
    parser.add_argument(
        "--cluster_downsampling",
        dest="cluster_downsampling",
        help="""downsample_clusters output""",
    )
    parser.add_argument(
        "--variants",
        dest="variants",
        nargs='?',
        help="""file with variant table (from find_variants)""",
    )
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000:-",
        help="""coordinates of switch region [%(default)s]""",
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
        "--weights",
        dest="weights",
        default="cluster",
        help="""specify weights ("cluster" | "reads" | "adjusted") [%(default)s]""",
    )
    parser.add_argument(
        "--max_n_blunt",
        dest="max_n_blunt",
        default=0.5,
        type=float,
        help="""max avg. number of untemplated / homologous nucleotides of reads in cluster to be considered a blunt end [%(default).1f]""",
    )
    parser.add_argument(
        "--big_clone_cutoff",
        dest="big_clone_cutoff",
        default=0.01,
        type=float,
        help="""cutoff to determine a "big" clone as fraction of clustered reads [%(default).2f]""",
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
        isotype_colors,
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
    logger.info("reading processing stats")
    stats["process"] = pd.read_csv(args.process, index_col=0, header=None).squeeze()
    stats["process"]["nreads_unmapped"] = (
        stats["process"]["nreads_initial"] - stats["process"]["nreads_mapped"]
    )

    logger.info("reading clustering, clustering stats, clustering analysis")
    clustering = pd.read_csv(args.clustering, index_col=0, header=0)
    cluster_analysis = pd.read_csv(args.cluster_analysis, index_col=0, header=0)
    stats["cluster"] = (
        pd.read_csv(args.cluster_stats, index_col=0, header=None).squeeze().drop("eff_nclusters")
    )
    stats = pd.concat(stats.values(), axis=0)

    stats["frac_reads_unused"] = (
        stats[["nreads_removed_incomplete", "nreads_unmapped", "nreads_removed_no_switch"]].sum()
        / stats["nreads_initial"]
    )

    if args.use_clones:
        clones = list(map(int, args.use_clones.split(",")))
    else:
        clusters = clustering["filtered_cluster"].dropna()
        clones = clusters[clusters >= 0].astype(int).unique()

    logger.info("using {0} clones".format(len(clones)))

    if args.weights == "cluster":
        logger.info("using uniform weights per cluster")
        w = np.ones(len(clones))
    elif args.weights == "reads":
        logger.info("using uniform weights per read")
        w = cluster_analysis.loc[clones, "size"]
    elif args.weights == "adjusted":
        logger.info("using adjusted weights per cluster")
        w = cluster_analysis.loc[clones, "adj_size"]
    else:
        raise ValueError("invalid value {0} for args.weights!".format(args.weights))

    stats["nclusters_final"] = len(clones)
    stats["nclusters_eff"] = np.exp(scipy.stats.entropy(cluster_analysis["size"]))
    stats["nreads_final"] = sum(~clustering["cluster"].isna())

    if len(clones) > 0:
        rel_size = cluster_analysis.loc[clones, "size"] / cluster_analysis.loc[clones, "size"].sum()

        assert (rel_size.sum() > 0.999) & (
            rel_size.sum() < 1.001
        ), "relative sizes don't sum up to 1"

        stats["mean_length"], stats["std_length"] = weighted_avg_and_std(
            cluster_analysis.loc[clones, "length"], w
        )

        stats["mean_GC"], stats["std_GC"] = weighted_avg_and_std(
            cluster_analysis.loc[clones, "GC"], w
        )
        stats["mean_cluster_size"] = rel_size.mean()
        stats["std_cluster_size"] = rel_size.std()

        stats["cluster_gini"] = calculate_gini(rel_size)
        stats["cluster_entropy"] = scipy.stats.entropy(rel_size) / np.log(len(rel_size))
        stats["cluster_inverse_simpson"] = 1.0 / (rel_size**2).sum()

        stats["top_clone_occupancy"] = rel_size.max()
        stats["big_clones_occupancy"] = rel_size[rel_size > args.big_clone_cutoff].sum()

        stats["mean_frac_mapped"] = cluster_analysis.loc[clones, "frac_mapped"].mean()
        stats["mean_frac_mapped_multi"] = cluster_analysis.loc[clones, "frac_mapped_multi"].mean()
        stats["mean_frac_ignored"] = cluster_analysis.loc[clones, "frac_ignored"].mean()

        stats["median_inversion_size"] = cluster_analysis.loc[clones, "inversion_size"][
            cluster_analysis.loc[clones]["inversion_size"] > 0
        ].median()
        stats["median_duplication_size"] = cluster_analysis.loc[clones, "duplication_size"][
            cluster_analysis.loc[clones]["duplication_size"] > 0
        ].median()

        try:
            stats["size_length_bias"] = scipy.stats.linregress(
                cluster_analysis.loc[clones, "length"],
                np.log(cluster_analysis.loc[clones, "size"]),
            )[0]
        except ValueError:
            stats["size_length_bias"] = 0

        try:
            stats["size_GC_bias"] = scipy.stats.linregress(
                cluster_analysis.loc[clones, "GC"],
                np.log(cluster_analysis.loc[clones, "size"]),
            )[0]
        except ValueError:
            stats["size_GC_bias"] = 0

    if args.cluster_downsampling:
        logger.info(
            "reading cluster downsampling results from {0}".format(args.cluster_downsampling)
        )
        downsampling = pd.read_csv(
            args.cluster_downsampling,
            header=0,
            index_col=0,
        ).mean(0)
        stats = pd.concat([stats, downsampling], axis=0)

    isotype_read_count = clustering["isotype"].dropna().value_counts()
    isotype_cluster_count = cluster_analysis.loc[clones, "isotype"].dropna().value_counts()
    isotype_insert_count = cluster_analysis["insert_pos_isotype"].dropna().value_counts()

    isotype_read_fraction = isotype_read_count / isotype_read_count.sum()
    isotype_cluster_fraction = isotype_cluster_count / isotype_cluster_count.sum()

    isotype_read_fraction.index = "frac_reads_" + isotype_read_fraction.index
    isotype_cluster_fraction.index = "frac_clusters_" + isotype_cluster_fraction.index
    isotype_insert_count.index = "ninserts_" + isotype_insert_count.index.astype(str)
    take = ~clustering["cluster"].isna()
    inserts = [
        decode_coords(m)
        for insert in clustering[take]["inserts"].dropna()
        for m in insert.split(",")
    ]
    cluster_stats = {
        "n_inserts": len(inserts),
        "n_unique_inserts": len(merge_intervals(inserts)),
        "n_clusters_inserts": len(cluster_analysis["insert_pos_isotype"].dropna()),
        "insert_frequency": len(merge_intervals(inserts)) / len(clones) if len(clones) > 0 else 0,
    }
    insert_stats = cluster_analysis.agg(
        {
            "insert_frequency": "mean",
            "insert_overlap": "mean",
            "insert_pos_overlap": "mean",
            "insert_length": "mean",
            "insert_gap_length": "mean",
        },
        axis=0,
    )
    insert_stats.index = "mean_" + insert_stats.index.astype(str)

    cluster_analysis["isotype_simple"] = cluster_analysis["isotype"].str[:2].str.upper()
    alpha_clusters = cluster_analysis.loc[clones, "isotype_simple"] == "SA"
    alpha_ratio = pd.Series(
        {
            "alpha_ratio_reads": cluster_analysis.loc[clones][alpha_clusters]["size"].sum()
            / cluster_analysis.loc[clones]["size"].sum(),
            "alpha_ratio_clusters": alpha_clusters.mean(),
        }
    )

    if len(clones) > 0:
        cluster_analysis["weights"] = np.nan
        cluster_analysis.loc[clones, "weights"] = w
        isotype_length = (
            cluster_analysis.loc[clones]
            .groupby("isotype_simple")[["length", "weights"]]
            .apply(
                lambda df: pd.DataFrame(
                    weighted_avg_and_std(df["length"], df["weights"]), index=["mean", "std"]
                )
            )
            .squeeze()
        )
        isotype_length.index = [y + "_length_" + x for x, y in isotype_length.index.tolist()]
    else:
        isotype_length = None

    realignment_stats = (
        cluster_analysis.loc[clones]
        .groupby("isotype_simple")[["n_untemplated_switch", "n_homology_switch"]]
        .mean()
        .unstack()
    )
    realignment_stats.index = ["_".join(i) for i in realignment_stats.index.tolist()]
    blunt_ends = (
        cluster_analysis.loc[
            clones, ["isotype_simple", "n_untemplated_switch", "n_homology_switch"]
        ]
        .groupby("isotype_simple")
        .apply(lambda x: (x < args.max_n_blunt).all(axis=1).mean())
    )
    blunt_ends.index = "frac_blunt_" + blunt_ends.index
    if not isinstance(blunt_ends, pd.Series):
        blunt_ends = pd.Series(pd.NA, index=[])

    blunt_ends["frac_blunt"] = (
        (
            cluster_analysis.loc[clones][["n_untemplated_switch", "n_homology_switch"]]
            < args.max_n_blunt
        )
        .all(axis=1)
        .mean()
    )

    realignment_stats = pd.concat(
        [
            cluster_analysis.loc[clones][["n_untemplated_switch", "n_homology_switch"]].mean(),
            realignment_stats,
            blunt_ends,
        ],
        axis=0,
    )

    stats = pd.concat(
        [
            stats,
            isotype_length,
            isotype_read_fraction,
            isotype_cluster_fraction,
            alpha_ratio,
            pd.Series(cluster_stats),
            insert_stats,
            isotype_insert_count,
            realignment_stats,
        ],
        axis=0,
    )

    cutoff_scanning = pd.read_csv(args.scanning, header=0, index_col=0)

    if args.breakpoint_stats and os.path.isfile(args.breakpoint_stats):
        logger.info("loading breakpoint stats from " + args.breakpoint_stats)
        breakpoint_stats = pd.read_csv(args.breakpoint_stats, header=None, index_col=0).squeeze()
        stats = pd.concat([stats, breakpoint_stats], axis=0)

    if args.variants and os.path.isfile(args.variants):
        logger.info("loading variants from " + args.variants)
        variants = pd.read_csv(args.variants, sep="\t", header=0)

        ref = variants["ref"].apply(lambda x: ncodes[x] - 1)
        alt = variants["alt"].apply(lambda x: ncodes[x] - 1)
        germline = (variants["type"] != "n.d.") | ~variants["anno"].isnull()

        mm = (
            pd.crosstab(ref, alt)
            .reindex(index=range(4), columns=range(4))
            .fillna(0)
            .astype(int)
            .values
        )

        var_stats = pd.Series(
            {
                "num_variants": len(ref),
                "frac_variants_germline": np.mean(germline),
                "frac_variants_transitions": mm[(0, 2, 1, 3), (2, 0, 3, 1)].sum() / len(ref),
                "frac_variants_C>A": mm[(1, 2), (0, 3)].sum() / len(ref),
                "frac_variants_C>G": mm[(1, 2), (2, 1)].sum() / len(ref),
                "frac_variants_C>T": mm[(1, 2), (3, 0)].sum() / len(ref),
                "frac_variants_T>A": mm[(3, 0), (0, 3)].sum() / len(ref),
                "frac_variants_T>C": mm[(3, 0), (1, 2)].sum() / len(ref),
                "frac_variants_T>G": mm[(3, 0), (2, 1)].sum() / len(ref),
            }
        )

        if (~variants["motif"].isnull()).sum() > 0:
            for mot in set(",".join(variants["motif"].dropna()).split(",")):
                var_stats["num_variants_" + mot] = (variants["motif"].str.contains(mot)).sum()

        stats = pd.concat([stats, var_stats], axis=0)

    if args.stats is not None:
        logger.info("saving summary to " + args.stats)
        stats.to_csv(args.stats)

    if args.figure is not None:
        logger.info("reading read info")
        read_info = pd.read_csv(args.info, header=0, index_col=0)
        assert read_info.index.is_unique, "index of info file is not unique!"

        barcode_locs = (
            [
                (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2])) / (2 * l)
                for bc, l in read_info[["barcodes", "length"]].dropna().values
                for x in bc.split(";")
                if "BC" in bc
            ]
            if "barcodes" in read_info.columns
            else []
        )

        primer_locs = (
            [
                (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2])) / (2 * l)
                for pr, l in read_info[["primers", "length"]].dropna().values
                for x in pr.split(";")
                if "primer" in pr
            ]
            if "primers" in read_info.columns
            else []
        )

        internal_primer_locs = (
            [
                (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2])) / (2 * l)
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
                (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2])) / (2 * l)
                for bc, l in read_info[["internal", "length"]].dropna().values
                for x in bc.split(";")
                if "BC" in x
            ]
            if "internal" in read_info.columns
            else []
        )

        fig, axs = plt.subplots(4, 2, figsize=(6, 8))
        fig.subplots_adjust(hspace=0.5, wspace=0.4, bottom=0.1, top=0.95)

        fig.text(0.5, 0.98, args.sample, size="medium", ha="center")

        ax = axs[0, 0]
        ax.hist(
            read_info["length"],
            bins=np.geomspace(
                0.95 * read_info["length"].min(), 1.05 * read_info["length"].max(), 100
            ),
            histtype="step",
        )
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylabel("# reads")
        ax.set_xlabel("read length")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

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
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax = axs[1, 0]
        isotype_read_fraction.index = isotype_read_fraction.index.str.split("_").str[-1]
        isotype_cluster_fraction.index = isotype_cluster_fraction.index.str.split("_").str[-1]
        isotype_fracs = pd.DataFrame(
            {
                "reads\n(n={0})".format(isotype_read_count.sum()): isotype_read_fraction,
                "clusters\n(n={0})".format(isotype_cluster_count.sum()): isotype_cluster_fraction,
            }
        )
        isotype_fracs.T.plot(
            kind="barh",
            stacked=True,
            ax=ax,
            color=[
                isotype_colors[it] if it in isotype_colors else "gray" for it in isotype_fracs.index
            ],
        )
        ax.set_xlabel("fraction")
        ax.legend(
            ncol=4,
            handlelength=1,
            columnspacing=0.5,
            loc=2,
            frameon=False,
            prop={"size": "small"},
        )
        ax.set_ylim([-0.5, 2])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax = axs[1, 1]
        cc = cutoff_scanning.index
        nc = cutoff_scanning["nclusters"].values
        ax.plot(cc, nc, "-", color="#e41a1c")

        if "p0" in stats.index:
            pars = stats[["p0", "p1", "p2", "p3"]].values
            ax.plot(cc, f2(pars, cc), "-", color="k", lw=1)
            ax.plot(cc, f2((pars[0], pars[1], 0, 0), cc), "--", color="k", lw=0.5)
        ax.axvline(
            stats["c_opt"],
            color="k",
            linestyle="dashed",
            lw=0.5,
        )
        ax.set_yscale("log")
        ax.set_xlabel("cutoff")
        ax.set_ylabel("# clusters", color="#e41a1c")
        ax.set_ylim([max(1, ax.get_ylim()[0]), None])
        ax.spines["top"].set_visible(False)

        ax2 = ax.twinx()
        entropy = cutoff_scanning["entropy"].values
        ax2.plot(cc, entropy, "-", color="#377eb8")
        ax2.set_ylabel("entropy", color="#377eb8")
        ax2.spines["top"].set_visible(False)

        major_ticks = []
        minor_ticks = []
        minor_labels = []
        for rec in anno_recs:
            start = shift_coord(int(rec[3][1]), cov_int)
            end = shift_coord(int(rec[3][2]), cov_int)
            major_ticks += [start, end]
            minor_ticks.append((start + end) / 2)
            minor_labels.append(rec[3][3])

        ax = axs[2, 0]
        ax.hist(
            clustering["cluster"].dropna().value_counts(),
            bins=np.geomspace(1, stats["nreads_final"], 50),
            histtype="step",
            label="all",
        )
        ax.hist(
            clustering["filtered_cluster"][clustering["filtered_cluster"] >= 0].value_counts(),
            bins=np.geomspace(1, stats["nreads_final"], 50),
            histtype="stepfilled",
            label="filtered",
            alpha=0.5,
        )
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylabel("# clusters")
        ax.set_xlabel("cluster size")
        ax.legend(loc=1, handlelength=1, prop={"size": "small"}, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax = axs[2, 1]
        ax.plot(
            cluster_analysis["size"],
            cluster_analysis["length"],
            ".",
        )
        ax.set_xscale("log")
        ax.set_ylabel("length")
        ax.set_xlabel("cluster size")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        logger.info("loading gaps from " + args.gaps)
        gaps = np.load(args.gaps)

        ax = axs[3, 0]
        ax.hist(gaps["gap_size"], bins=np.geomspace(1, 3.0e4, 50), histtype="step")
        ax.axvline(args.max_gap, lw=0.5, color="r", ls="dashed")
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylabel("# events")
        ax.set_xlabel("gap size")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax = axs[3, 1]
        use = gaps["gap_size"] > args.max_gap
        ax.hist(
            [gaps["gap_right"][use], gaps["gap_left"][use]],
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
        ax.legend(loc=9, handlelength=1, frameon=False, prop={"size": "small"}, ncol=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.savefig(args.figure)
