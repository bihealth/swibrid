"""\
given a MSA and clustering results, this will display individual reads mapping over the switch regions
reads will be ordered as dictated by the linkage, or simply by isotype and cluster value
reads can be colored by different variables (isotype, cluster, haplotype, coverage, strand, orientation or other columns present in the info file)
an additional sidebar can display additional variables
variant positions or breakpoint realignment statistics can also be indicated
"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""required: file(s) with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("--figure", dest="figure", help="""required: output figure""")
    parser.add_argument(
        "--clustering_results",
        dest="clustering_results",
        help="""required: file contains clustering results""",
    )
    parser.add_argument(
        "--switch_annotation",
        dest="switch_annotation",
        help="""required: bed file with switch annotation""",
    )
    parser.add_argument("--linkage", dest="linkage", help="""file containing flinkage""")
    parser.add_argument("--sample", dest="sample", help="""sample name (for figure title)""")
    parser.add_argument("--info", dest="info", help="""file with read info""")
    parser.add_argument(
        "--clustering_stats",
        dest="clustering_stats",
        help="""file contains clustering stats with cutoff value""",
    )
    parser.add_argument(
        "--cutoff",
        dest="cutoff",
        type=float,
        help="""show cutoff line""",
    )
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000:-",
        help="""coordinates of switch region [chr14:106050000-106337000:-]""",
    )
    parser.add_argument(
        "--annotation",
        dest="annotation",
        nargs='?',
        help="""bed file with gene annotation""",
    )
    parser.add_argument(
        "--color_by",
        dest="color_by",
        default="isotype",
        help="""color reads by isotype, cluster, haplotype or by other columns present in the "info" file [isotype]""",
    )
    parser.add_argument(
        "--sidebar_color_by",
        dest="sidebar_color_by",
        help="""color sidebar reads by (comma-separated list of) isotype, cluster, haplotype or by other columns present in the "info" file [none]""",
    )
    parser.add_argument(
        "--show_inserts",
        dest="show_inserts",
        action="store_true",
        default=False,
        help="""show insert locations""",
    )
    parser.add_argument(
        "--coords",
        dest="coords",
        help="""file with processed read coordinates (required if inserts are shown)""",
    )
    parser.add_argument(
        "--no_x_ticks",
        dest="no_x_ticks",
        action="store_true",
        default=False,
        help="""omit the x axis ticks""",
    )
    parser.add_argument(
        "--omit_scale_bar",
        dest="omit_scale_bar",
        action="store_true",
        default=False,
        help="""omit scale bar"""
    )
    parser.add_argument(
        "--fig_width",
        dest="fig_width",
        default=6,
        type=float,
        help="""width of figure in inches [6]""",
    )
    parser.add_argument(
        "--fig_height",
        dest="fig_height",
        type=float,
        default=8,
        help="""width of figure in inches [8]""",
    )
    parser.add_argument(
        "--linkage_border",
        dest="linkage_border",
        type=float,
        default=0.1,
        help="""fraction of figure used for dendrogram [.1]""",
    )
    parser.add_argument(
        "--cutoff_color",
        dest="cutoff_color",
        default="r",
        help="""color of dendrogram cutoff line [r]""",
    )
    parser.add_argument(
        "--variants_table",
        dest="variants_table",
        help="""show variant positions from variant table (from find_variants)""",
    )
    parser.add_argument(
        "--variants_matrix",
        dest="variants_matrix",
        help="""indicate variant occurrences from variant matrix (from find_variants)""",
    )
    parser.add_argument(
        "--haplotypes",
        dest="haplotypes",
        help="""file with haplotype clustering (from find_variants)""",
    )
    parser.add_argument(
        "--realignments", dest="realignments", help="""indicate realignment scores with breakpoint realignments"""
    )
    parser.add_argument(
        "--dpi",
        dest="dpi",
        default=300,
        type=int,
        help="""dpi of output figure [300]""",
    )
    parser.add_argument(
        "--chunksize",
        dest="chunksize",
        default=5000,
        type=int,
        help="""plot image in chunks of reads of this size [5000]""",
    )
    parser.add_argument(
        "--cmax",
        dest="cmax",
        type=float,
        help="""maximum height in dendrogram plot""",
    )


def rand_cmap(
    nlabels,
    type="bright",
    first_color_black=False,
    last_color_black=False,
    verbose=False,
    bad="white",
    under="black",
    over="black",
    seed=0,
):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    (from https://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib/14720445)
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    np.random.seed(seed)

    # Generate color map for bright colors, based on hsv
    if type == "bright":
        randHSVcolors = [
            (
                np.random.uniform(low=0.0, high=1),
                np.random.uniform(low=0.2, high=1),
                np.random.uniform(low=0.9, high=1),
            )
            for i in range(nlabels)
        ]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

    # Generate soft pastel colors, by limiting the RGB spectrum
    elif type == "soft":
        low = 0.6
        high = 0.95
        randRGBcolors = [
            (
                np.random.uniform(low=low, high=high),
                np.random.uniform(low=low, high=high),
                np.random.uniform(low=low, high=high),
            )
            for i in range(nlabels)
        ]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

    else:
        print('Please choose "bright" or "soft" for type')
        return

    random_colormap = LinearSegmentedColormap.from_list("new_map", randRGBcolors, N=nlabels)

    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        colorbar.ColorbarBase(
            ax,
            cmap=random_colormap,
            norm=norm,
            spacing="proportional",
            ticks=None,
            boundaries=bounds,
            format="%1i",
            orientation="horizontal",
        )

    random_colormap.set_bad(bad)
    random_colormap.set_under(under)
    random_colormap.set_over(over)

    return random_colormap


def run(args):
    import os
    import sys
    import numpy as np
    import scipy.sparse
    import scipy.cluster.hierarchy
    import re
    import pandas as pd
    from collections import defaultdict
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.collections import LineCollection
    from logzero import logger
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
        intersect_intervals,
        interval_length,
        get_switch_coverage,
        shift_coord,
        decode_insert,
        merge_intervals,
    )

    matplotlib.rcParams.update({"font.size": 8})
    matplotlib.rcParams.update({"axes.linewidth": 0.5})
    matplotlib.rcParams.update({"xtick.major.width": 0.5})
    matplotlib.rcParams.update({"ytick.major.width": 0.5})

    # increase recursion limit for plotting large dendrograms
    sys.setrecursionlimit(100000)

    scale_bar_x_length = 5000
    scale_bar_x_legend = "5kb"

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

    if args.annotation:
        logger.info("using annotation from " + args.annotation)
        annotation = defaultdict(list)
        for line in open(args.annotation):
            ls = line.strip("\n").split("\t")
            chrom = ls[0]
            start = int(ls[1])
            end = int(ls[2])
            annotation[chrom].append((chrom, start, end) + tuple(ls[3:]))

    logger.info("loading msa from {0}".format(args.msa))
    msa = scipy.sparse.load_npz(args.msa).tocsr()
    nreads = msa.shape[0]

    logger.info("loading clustering from {0}".format(args.clustering_results))
    clustering = pd.read_csv(args.clustering_results, index_col=0, header=0)
    if clustering.shape[0] != nreads:
        logger.warn("# of reads in clustering different from # of reads in msa!")
        exit(1)
    reads = clustering["cluster"].dropna().index
    clustering = clustering.loc[reads]
    nreads = len(reads)
    if len(reads) < msa.shape[0]:
        logger.info("restricting msa to {0} reads in clustering\n".format(nreads))
        msa = msa[:nreads]

    if args.info:
        logger.info("reading read info from {0}".format(args.info))
        read_info = pd.read_csv(args.info, header=0, index_col=0)

    if args.variants_table:
        logger.info("reading variant table from {0}".format(args.variants_table))
        variants = pd.read_csv(args.variants_table, sep="\t", header=0)
    if args.variants_matrix:
        logger.info("reading variant matrix from {0}".format(args.variants_matrix))
        vmat = scipy.sparse.load_npz(args.variants_matrix)

    if args.haplotypes:
        logger.info("reading haplotypes from {0}".format(args.haplotypes))
        haplotypes = pd.read_csv(args.haplotypes, index_col=0, header=0)

    logger.info("coloring msa by {0}".format(args.color_by))
    if args.color_by == "cluster":
        csize = clustering["cluster"].value_counts()
        nclust = len(csize)
        crank = pd.Series(range(nclust), csize.index)
        cmap = rand_cmap(
            max(2, nclust),
            type="bright",
            first_color_black=False,
            last_color_black=False,
            verbose=False,
            seed=10,
        )
        values = crank[clustering["cluster"].astype(int)].values % nclust
    elif args.color_by == "isotype":
        cmap = plt.cm.Set1
        values = clustering["isotype"].astype("category").cat.codes.values % 9
    elif args.info and args.color_by in read_info.columns:
        info_col = read_info.loc[reads, args.color_by]
        if pd.api.types.is_numeric_dtype(read_info[args.color_by]):
            values = (info_col - info_col.min()) / (info_col.max() - info_col.min())
            cmap = plt.cm.cool
        else:
            if args.color_by in ["primers", "barcodes"]:
                info_col = info_col.apply(
                    lambda x: "+".join(
                        set(
                            map(
                                lambda y: re.sub("primer_", "", y.split("@")[0]),
                                x.split(";"),
                            )
                        )
                    )
                )
            ncats = len(info_col.astype("category").cat.categories)
            if ncats < 9:
                values = info_col.astype("category").cat.codes % 9
                cmap = plt.cm.Set1
            else:
                values = info_col.astype("category").cat.codes % 20
                cmap = plt.cm.tab20
    elif args.color_by == "orientation":
        cmap = plt.cm.PiYG
    elif args.color_by == "coverage":
        cmap = plt.cm.rainbow_r
    elif args.color_by == "strand" and args.info:
        values = (clustering["orientation"] == "+").astype(int) + 1
        cmap = plt.cm.PiYG
    elif args.color_by == "haplotype" and args.haplotypes is not None:
        values = haplotypes.loc[clustering["cluster"].astype(int).values, "haplotype"].values
        cmap = plt.cm.coolwarm
    else:
        logger.warn("no info on {0}!".format(args.color_by))
        values = [0.3] * nreads
        cmap = plt.cm.Greys

    if args.cutoff is not None:
        logger.info("setting clustering cutoff at {0:.2f}".format(args.cutoff))
        cutoff = args.cutoff
    elif args.clustering_stats is not None:
        logger.info("reading clustering stats from {0}".format(args.clustering_stats))
        clustering_stats = pd.read_csv(args.clustering_stats, header=None, index_col=0).squeeze()
        cutoff = clustering_stats["c_opt"]
    else:
        cutoff = None

    logger.info("creating figure")

    figsize = (args.fig_width, args.fig_height + 1)
    fig = plt.figure(figsize=figsize)
    fig.clf()

    bottom = 0.2 / figsize[1]
    height = args.fig_height / figsize[1]
    linkage_border = args.linkage_border if args.linkage else 0.01
    insert_border = 0.2 if args.show_inserts else 0.01
    if args.sidebar_color_by is not None:
        sidebar_colors = args.sidebar_color_by.split(',')
        valid_cols = set(['strand','isotype','cluster','haplotype'])
        if args.info:
            valid_cols |= set(read_info.columns)
        sidebar_colors = [s for s in sidebar_colors if s in valid_cols]
        num_sidebars = len(sidebar_colors)
    else:
        num_sidebars = 0
    sidebar_width = .01
    left = .01 + linkage_border + num_sidebars * sidebar_width
    width = .99 - linkage_border - insert_border - num_sidebars * sidebar_width

    if args.linkage:
        logger.info("reading linkage from {0} and creating dendrogram".format(args.linkage))
        Z = np.load(args.linkage)["Z"]
        ax = fig.add_axes([0.01, bottom, left - num_sidebars * sidebar_width - .015, height])
        lw = fig.bbox_inches.height * ax.get_position().height * 72 / nreads
        with plt.rc_context({"lines.linewidth": min(lw, 0.5)}):
            L = scipy.cluster.hierarchy.dendrogram(
                Z,
                orientation="left",
                no_labels=True,
                link_color_func=lambda x: "k",
                ax=ax,
                color_threshold=0,
            )

            if cutoff:
                ax.axvline(cutoff, color=args.cutoff_color, lw=0.5)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            if args.cmax is None:
                ax.set_xlim([Z[:, 2].max(), 0])
                ax.set_xscale("function", functions=(lambda x: x ** (1 / 2), lambda x: x**2))
            else:
                ax.set_xlim([args.cmax, 0])
            ax.set_ylim(ax.get_ylim()[::-1])
            ax.set_xticks([])
            ax.set_yticks([])
            order = L["leaves"]
            if len(order) != nreads:
                raise Exception("linkage does not fit to clustering and/or msa!")
    else:
        order = np.lexsort((clustering["cluster"], clustering["isotype"]))[::-1]

    if num_sidebars > 0:
        logger.info("adding {0} sidebars".format(len(sidebar_colors)))
        for nsc, sidebar_color in enumerate(sidebar_colors):
            ax = fig.add_axes([left - (num_sidebars - nsc) * sidebar_width, 
                               bottom, .9 * sidebar_width, height])
            if sidebar_color == "isotype":
                sidebar_values = clustering["isotype"].astype("category").cat.codes.values % 9
                sidebar_cmap = plt.cm.Set1
            elif sidebar_color == "cluster":
                csize = clustering["cluster"].value_counts()
                nclust = len(csize)
                crank = pd.Series(range(nclust), csize.index)
                sidebar_cmap = rand_cmap(
                    max(2, nclust),
                    type="bright",
                    first_color_black=False,
                    last_color_black=False,
                    verbose=False,
                    seed=10,
                )
                sidebar_values = crank[clustering["cluster"].astype(int)].values % nclust
            elif sidebar_color == "haplotype":
                sidebar_values = haplotypes.loc[
                    clustering["cluster"].astype(int).values, "haplotype"
                ].values
                sidebar_cmap = plt.cm.coolwarm
            elif args.info and sidebar_color in read_info.columns:
                info_col = read_info.loc[reads, sidebar_color]
                if pd.api.types.is_numeric_dtype(read_info[sidebar_color]):
                    sidebar_values = (
                        (info_col - info_col.min()) / (info_col.max() - info_col.min())
                    ).values
                    sidebar_cmap = plt.cm.cool
                else:
                    if sidebar_color in ["primers", "barcodes"]:
                        info_col = info_col.apply(
                            lambda x: "+".join(
                                set(
                                    map(
                                        lambda y: re.sub("primer_", "", y.split("@")[0]),
                                        x.split(";"),
                                    )
                                )
                            )
                        )
                    ncats = len(info_col.astype("category").cat.categories)
                    if ncats < 9:
                        sidebar_values = info_col.astype("category").cat.codes.values % 9
                        sidebar_cmap = plt.cm.Set1
                    else:
                        sidebar_values = info_col.astype("category").cat.codes.values % 20
                        sidebar_cmap = plt.cm.tab20
            elif sidebar_color == "strand":
                sidebar_values = (clustering["orientation"].values == "+").astype(int) + 1
                sidebar_cmap = plt.cm.PiYG
            else:
                logger.warn("no info on {0}!".format(sidebar_color))
                sidebar_values = np.array([0.3] * nreads)
                sidebar_cmap = plt.cm.Greys

            ax.imshow(
                sidebar_values[order, np.newaxis],
                aspect="auto",
                interpolation="nearest",
                cmap=sidebar_cmap,
                vmin=0,
                vmax=max(1, sidebar_values.max()),
            )
            ax.set_axis_off()
            ax.set_title(sidebar_color, rotation=45, size=4, ha='left')

    logger.info("plotting MSA")
    ax = fig.add_axes([left, bottom, width, height])
    ax.set_facecolor('w')

    # plot the image in chunks
    nchunks = np.ceil(nreads / args.chunksize).astype(int)
    for n in range(nchunks):
        order_chunk = order[n * args.chunksize : min((n + 1) * args.chunksize, nreads)]
        extent = [
            0,
            Ltot,
            nreads - min((n + 1) * args.chunksize, nreads),
            nreads - n * args.chunksize,
        ]
        msa_chunk = msa[order_chunk]
        im = np.empty(msa_chunk.shape)
        im[:] = np.nan
        if args.color_by == "orientation":
            pos = msa_chunk.data > 0
            neg = msa_chunk.data < 0
            i, j = np.nonzero(msa_chunk)
            im[(i[pos], j[pos])] = 0.8
            im[(i[neg], j[neg])] = 0.1
            values = np.array([0, 1])
        elif args.color_by == "coverage":
            im[np.nonzero(msa_chunk)] = (msa_chunk.data // 10 + 2) / 4
            values = np.array([0, 1])
        else:
            tmp = np.broadcast_to(values[order_chunk], msa_chunk.T.shape).T
            im[np.nonzero(msa_chunk)] = tmp[np.nonzero(msa_chunk)]

        if args.variants_matrix and vmat.nnz > 0:
            ax.imshow(
                im,
                aspect="auto",
                interpolation="nearest",
                cmap=cmap,
                alpha=0.5,
                extent=extent,
                vmin=0,
                vmax=max(1, values.max()),
            )
            x, y = vmat[order_chunk].nonzero()
            use = np.isin(y, variants["rel_pos"])
            ax.scatter(
                y[use],
                nreads - (x[use] + n * args.chunksize) - .5,
                marker=".",
                lw=0,
                s=0.1,
                zorder=2,
                color="k",
                edgecolor=None,
            )
        else:
            ax.imshow(
                im,
                aspect="auto",
                interpolation="nearest",
                cmap=cmap,
                extent=extent,
                vmin=0,
                vmax=max(1, values.max()),
            )

    if args.variants_table:
        logger.info("adding variant positions")
        for typ, c in zip(("n.d.", "het0", "het1", "hom"), ("gray", "b", "r", "k")):
            take = variants["type"] == typ
            ax.scatter(
                variants[take]["rel_pos"].values,
                nreads * np.ones(take.sum()),
                s=6,
                c=c,
                marker=7,
                lw=0,
                zorder=2,
                edgecolor=None,
                clip_on=False,
            )

    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([])

    major_ticks = []
    minor_ticks = []
    minor_labels = []
    for rec in anno_recs:
        start = shift_coord(int(rec[3][1]), cov_int) 
        end = shift_coord(int(rec[3][2]), cov_int) 
        major_ticks += [start, end]
        minor_ticks.append((start + end) / 2)
        minor_labels.append(rec[3][3])
    ax.set_xticks(np.unique(major_ticks))
    ax.set_xticklabels([])
    ax.set_xticks(minor_ticks, minor=True)
    if not args.no_x_ticks:
        ax.set_xticklabels(minor_labels, minor=True, size="small")
    ax.tick_params(which="minor", length=0)

    if not args.omit_scale_bar:
        # add genomic scale bar
        ax.hlines(
            nreads + 0.01 * nreads,
            0,
            scale_bar_x_length,
            color="k",
            lw=2,
            clip_on=False,
        )
        ax.text(
            0.5 * scale_bar_x_length,
            nreads + 0.015 * nreads,
            scale_bar_x_legend,
            color="k",
            clip_on=False,
            size="xx-small",
            ha="center",
            va="bottom",
        )

    xlim = [0, Ltot]
    ylim = [0, nreads]
    ax.set_xlim(xlim if switch_orientation == "+" else xlim[::-1])
    ax.set_ylim(ylim)

    if args.show_inserts:
        assert args.coords is not None, "need processed read coordinates to show inserts!"

        logger.info("reading processed read coordinates from {0}".format(args.coords))
        process = pd.read_csv(args.coords, sep="\t", header=0, index_col=0)
        read_inserts = dict(
            (read, [decode_insert(insert) for insert in inserts.split(";")])
            for read, inserts in process["inserts"].dropna().items()
        )

        stat_string = (
            "{0}: {1} reads "
            "({2} inserts; {3} unique inserts; "
            "{4} clusters; {5} eff. clusters)\n"
        )

        inserts = [
            (
                m.group("insert_chrom"),
                int(m.group("insert_start")),
                int(m.group("insert_end")),
            )
            for read, inserts in read_inserts.items()
            for m in inserts
        ]
        ninserts = len(inserts)
        unique_inserts = merge_intervals(inserts)
        ninserts_unique = len(unique_inserts)
        nclusts = len(np.unique(clustering["cluster"].astype(int)))
        nclusts_eff = len(
            np.unique(clustering["filtered_cluster"][clustering["filtered_cluster"] >= 0])
        )
        stats = stat_string.format(
            args.sample, nreads, ninserts, ninserts_unique, nclusts, nclusts_eff
        )

        fig.text(0.01, 0.99, stats, size="x-small", ha="left", va="top")

        dx = 0.01 * Ltot

        insert_pos = {}
        dn = nreads / float(ninserts_unique) if ninserts_unique > 0 else 1
        n = nreads + 0.5 * dn

        arrows = []
        for read in clustering.iloc[order].dropna(subset=["inserts"]).index:
            p = nreads - reads[order].get_loc(read) - 1
            for m in read_inserts[read]:
                insert_chrom = m.group("insert_chrom")
                insert_start = int(m.group("insert_start"))
                insert_end = int(m.group("insert_end"))
                switch_left = int(m.group("switch_left"))
                switch_right = int(m.group("switch_right"))
                if (switch_orientation == "+" and switch_right < switch_left) or (
                    switch_orientation == "-" and switch_right > switch_left
                ):
                    switch_right, switch_left = switch_left, switch_right
                insert = (insert_chrom, insert_start, insert_end)
                ax.plot(
                    shift_coord(switch_left, cov_int),
                    p + 0.5,
                    ">",
                    color="k",
                    markersize=0.5,
                )
                ax.plot(
                    shift_coord(switch_right, cov_int),
                    p + 0.5,
                    "<",
                    color="k",
                    markersize=0.5,
                )
                arrows.append([(0, p + 0.5), (-dx, p + 0.5)])
                uinsert = [
                    ins
                    for ins in unique_inserts
                    if interval_length(intersect_intervals([insert], [ins])) > 0
                ][0]
                if uinsert in insert_pos:
                    arrows.append([(-dx, p + 0.5), (-3 * dx, insert_pos[uinsert])])
                    continue
                n -= dn
                arrows.append([(-dx, p + 0.5), (-3 * dx, n + 0.5)])
                insert_pos[uinsert] = n + 0.5
                insert_anno = set()
                if args.annotation:
                    for rec in intersect_intervals([uinsert], annotation[uinsert[0]], loj=True):
                        insert_anno.add(re.sub(".exon[0-9]*", "", rec[3][3]))
                if len(insert_anno) == 0 or insert_anno == set("."):
                    insert_anno = "{0}:{1}-{2}".format(*uinsert)
                else:
                    insert_anno = "|".join(insert_anno)
                # insert_anno='{0}:{1}-{2}'.format(*uinsert)
                arrows.append([(-3 * dx, n + 0.5), (-4 * dx, n + 0.5)])
                ax.text(
                    -4.2 * dx,
                    n + 0.5,
                    insert_anno,
                    size="xx-small",
                    color="k",
                    clip_on=False,
                    ha="left",
                    va="center",
                )

        ax.add_collection(LineCollection(arrows, linewidths=0.25, colors="k", clip_on=False))
    else:
        stat_string = "{0}: {1} reads " "({2} clusters; {3} eff. clusters)\n"
        nclusts = len(np.unique(clustering["cluster"].astype(int)))
        nclusts_eff = len(
            np.unique(clustering["filtered_cluster"][clustering["filtered_cluster"] >= 0])
        )
        stats = stat_string.format(args.sample, nreads, nclusts, nclusts_eff)

        fig.text(0.01, 0.99, stats, size="x-small", ha="left", va="top")

    if args.realignments is not None and os.path.isfile(args.realignments):
        logger.info("adding breakpoint realignment results")
        realignments = pd.read_csv(args.realignments, header=0, index_col=0)
        realignments = realignments.loc[realignments.index.intersection(clustering.index)]

        if realignments.shape[0] > 0:
            realignments = realignments[realignments["type"] == "switch"]
            realignments["cluster"] = clustering.loc[realignments.index, "cluster"]

            pleft = (
                realignments["pos_left"]
                .apply(lambda x: shift_coord(int(x.split(":")[1]), cov_int))
                .values
            )
            pright = (
                realignments["pos_right"]
                .apply(lambda x: shift_coord(int(x.split(":")[1]), cov_int))
                .values
            )
            nh = realignments["n_homology"].values
            nu = realignments["n_untemplated"].values

            pos = clustering.iloc[order].index.get_indexer(realignments.index)

            ax.scatter(
                np.maximum(pleft, pright),
                nreads - pos,
                s=0.1 * (1 + nu),
                c=nh,
                lw=0,
                zorder=2,
                edgecolor=None,
                clip_on=False,
                vmin=0,
                vmax=10,
                cmap=plt.cm.winter,
            )

    if args.figure is not None:
        logger.info("saving figure to {0}".format(args.figure))
        fig.savefig(args.figure, dpi=args.dpi)
