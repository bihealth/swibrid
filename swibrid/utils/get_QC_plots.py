"""plot QC stats"""


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
        help="""find_clusters.py output file""",
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


def run(args):

    import os
    import re
    import numpy as np
    import pandas as pd
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    import seaborn as sns

    from logzero import logger

    from .helpers import (
        parse_switch_coords,
        get_switch_coverage,
        decode_coords,
        read_switch_anno,
        merge_intervals,
        shift_coord,
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
    stats["filter"] = pd.read_csv(
        args.filter, index_col=0, header=None
    ).squeeze()
    logger.info("reading processing stats")
    stats["process"] = (
        pd.read_csv(args.process, index_col=0, header=0)
        .squeeze()
        .rename({"nreads": "mapped"})
    )
    logger.info("reading clustering stats")
    stats["cluster"] = pd.read_csv(
        args.cluster_stats, index_col=0, header=0
    ).squeeze()
    stats = pd.concat(stats.values(), axis=0)

    logger.info("reading clustering")
    clustering = pd.read_csv(args.clustering, index_col=0, header=0)
    clustering_analysis = pd.read_csv(
        args.clustering_analysis, index_col=0, header=0
    )
    stats["clustered"] = sum(~clustering["cluster"].isna())

    logger.info("reading read info")
    read_info = pd.read_csv(args.info, header=0, index_col=0)

    barcode_locs = [
        (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2])) / (2 * l)
        for bc, l in read_info[["barcodes", "length"]].values
        for x in bc.split(";")
    ]

    primer_locs = [
        (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2])) / (2 * l)
        for pr, l in read_info[["primers", "length"]].dropna().values
        for x in pr.split(";")
    ]

    internal_primer_locs = [
        (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2])) / (2 * l)
        for pr, l in read_info[["internal", "length"]].dropna().values
        if pr.count("primer") > 0
        for x in pr.split(";")
        if "primer" in x
    ]

    internal_barcode_locs = [
        (int(re.split("[@:-]", x)[1]) + int(re.split("[@:-]", x)[2])) / (2 * l)
        for bc, l in read_info[["internal", "length"]].dropna().values
        for x in bc.split(";")
        if "BC" in x
    ]

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

    isotype_read_count = clustering["isotype"].value_counts()
    isotype_cluster_count = clustering_analysis["isotype"].value_counts()

    nreads = isotype_read_count.sum()
    nclusters = isotype_cluster_count.sum()
    df = pd.DataFrame(
        {
            "reads\n(n={0})".format(nreads): isotype_read_count / nreads,
            "clusters\n(n={0})".format(nclusters): isotype_cluster_count
            / nclusters,
        }
    )
    ax = axs[1, 0]
    df.T.plot(kind="barh", stacked=True, ax=ax)
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

    isotype_read_count.index = "nreads_" + isotype_read_count.index
    isotype_cluster_count.index = "nclusters_" + isotype_cluster_count.index
    stats = pd.concat(
        [stats, isotype_read_count, isotype_cluster_count], axis=0
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
    }
    stats = pd.concat([stats, pd.Series(cluster_stats)], axis=0)

    extrapolation = pd.read_csv(
        args.extrapolation, header=None, index_col=0
    ).squeeze()

    ax = axs[1, 1]
    cc = extrapolation.index
    nc = extrapolation.values
    ax.plot(cc, nc, ".", color="#1f77b4", markersize=1)

    if "p0" in stats.index:
        pars = stats[["p0", "p1", "p2", "p3"]].values
        ax.plot(cc, f2(pars, cc), "-", color="k", lw=1)
        ax.plot(cc, f2((pars[0], pars[1], 0, 0), cc), "--", color="k", lw=0.5)
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
    neff = stats["eff_nclusters"].astype(int)
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
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("# clusters")
    ax.set_xlabel("cluster size")
    ax.legend(loc=1, handlelength=1, prop={"size": "small"}, frameon=False)
    sns.despine(ax=ax)

    logger.info("loading gaps from " + args.gaps)
    gaps = np.load(args.gaps)

    ax = axs[2, 1]
    ax.plot(
        clustering_analysis["size"], clustering_analysis["break_spread"], "."
    )
    ax.set_xscale("log")
    ax.set_ylabel("break spread")
    ax.set_xlabel("cluster size")
    sns.despine(ax=ax)

    ax = axs[3, 0]
    ax.hist(gaps["gap_size"], bins=np.geomspace(1, 3.0e4, 50), histtype="step")
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

    if args.gap_stats and os.path.isfile(args.gap_stats):
        gap_stats = pd.read_csv(
            args.gap_stats, header=None, index_col=0
        ).squeeze()
        stats = pd.concat([stats, gap_stats], axis=0)

    if args.stats is not None:
        stats.to_csv(args.stats)
