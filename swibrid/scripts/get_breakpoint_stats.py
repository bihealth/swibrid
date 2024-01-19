"""analyze breakpoint statistics"""


def setup_argparse(parser):
    parser.add_argument(
        "-g",
        "--gaps",
        dest="gaps",
        help="""file with gap positions (output of get_gaps.py)""",
    )
    parser.add_argument(
        "-r",
        "--rearrangements",
        dest="rearrangements",
        help="""file with inversion/duplication positions (output of find_variants.py)""",
    )
    parser.add_argument(
        "-c",
        "--clustering",
        dest="clustering",
        help="""file with clustering results""",
    )
    parser.add_argument(
        "-a",
        "--clustering_analysis",
        dest="clustering_analysis",
        help="""file with clustering analysis""",
    )
    parser.add_argument(
        "-b",
        "--binsize",
        dest="binsize",
        default=50,
        type=int,
        help="""binsize [50]""",
    )
    parser.add_argument(
        "--scale_factor",
        dest="scale_factor",
        default=10,
        type=int,
        help="""factor to increase binsize for 2D histogram [10]""",
    )
    parser.add_argument(
        "--max_gap",
        dest="max_gap",
        default=75,
        type=int,
        help="""max gap size to ignore [75]""",
    )
    parser.add_argument(
        "--homology",
        dest="homology",
        help="""file with homology values (output of get_switch_homology.py)""",
    )
    parser.add_argument(
        "--motifs",
        dest="motifs",
        help="""file with motif counts (output of get_switch_motifs.py)""",
    )
    parser.add_argument("-o", "--out", help="""output file""")
    parser.add_argument("-p", "--plot", help="""plot file""")
    parser.add_argument("--sample", dest="sample", help="""sample name (for figure title)""")
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
        "--range",
        dest="range",
        default="5",
        help="""range of kmer sizes, e.g., 3,5-7 [5]""",
    )
    parser.add_argument(
        "-n",
        "--ntop",
        dest="ntop",
        type=int,
        default=10,
        help="""number of top bins to select""",
    )
    parser.add_argument("--reference", dest="reference", help="""genome fasta file""")
    parser.add_argument(
        "--top_donor",
        dest="top_donor",
        help="""output file with top donor sequences""",
    )
    parser.add_argument(
        "--top_receiver",
        dest="top_receiver",
        help="""output file with top receiver sequences""",
    )
    parser.add_argument(
        "--use_clones",
        dest="use_clones",
        help="""comma-separated list of clones to use or 'all' (default: filtered clusters, excluding singletons)""",
    )
    parser.add_argument(
        "--weights",
        dest="weights",
        default="cluster",
        help="""use different weights ("cluster" | "reads" | "adjusted") [cluster]""",
    )


def run(args):
    import numpy as np
    import pandas as pd
    import scipy.sparse
    import scipy.stats
    import pysam
    from logzero import logger
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
        shift_coord,
        get_switch_coverage,
        parse_range,
        get_switch_iis,
    )

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

    binsize = args.binsize
    switch_iis = get_switch_iis(anno_recs, cov_int, eff_start, binsize)

    logger.info("reading gaps from " + args.gaps)
    gaps = np.load(args.gaps)
    nreads = gaps["read_idx"].max() + 1

    if args.clustering:
        logger.info("reading clustering from " + args.clustering)
        clustering = pd.read_csv(args.clustering, header=0, index_col=0)
        clustering = clustering[~clustering["cluster"].isna()]
        nreads = clustering.shape[0]
        logger.info("reading clustering analysis from " + args.clustering_analysis)
        analysis = pd.read_csv(args.clustering_analysis, header=0, index_col=0)
        if args.use_clones:
            if args.use_clones == "all":
                clones = clustering["cluster"].astype(int).unique()
            else:
                clones = list(map(int, args.use_clones.split(",")))
        else:
            clusters = clustering["filtered_cluster"].dropna()
            clones = clusters[clusters >= 0].astype(int).unique()
        logger.info("using {0} clones".format(len(clones)))
        nc = clustering["cluster"].dropna().astype(int).value_counts()
        singletons = nc.index[nc == 1]
        if args.weights == "cluster":
            logger.info("using uniform weights per cluster")
            w = 1.0 / nc.loc[clustering["cluster"].values]
        elif args.weights == "reads":
            logger.info("using uniform weights per read")
            w = pd.Series(1, index=np.arange(nreads))
        elif args.weights == "adjusted":
            logger.info("using adjusted weights per cluster")
            w = (
                analysis.loc[clustering["cluster"].values, "adj_size"]
                / analysis.loc[clustering["cluster"].values, "size"]
            )
        else:
            raise ValueError("invalid value {0} for args.weights!".format(args.weights))
        w[~w.index.isin(clones) | w.index.isin(singletons)] = 0
        weights = pd.Series(w.values / w.sum(), index=np.arange(nreads))
    else:
        weights = pd.Series(np.ones(nreads) / nreads, index=np.arange(nreads))

    gap_read = gaps["read_idx"]
    gap_left = gaps["gap_left"]
    gap_right = gaps["gap_right"]
    gap_size = gaps["gap_size"]

    # breaks in consistent orientation
    Leff = Ltot // binsize
    take = (gap_size >= args.max_gap) & (gap_left // binsize < Leff) & (gap_right // binsize < Leff)
    bp_hist = scipy.sparse.csr_matrix(
        (
            weights[gap_read[take]],
            (
                np.minimum(gap_left[take], gap_right[take]) // binsize,
                np.maximum(gap_left[take], gap_right[take]) // binsize,
            ),
        ),
        shape=(Leff, Leff),
    ).todense()
    np.nan_to_num(bp_hist, copy=False)

    if args.rearrangements:
        logger.info("reading rearrangements from " + args.rearrangements)
        rearrangements = np.load(args.rearrangements)
        inv_read = rearrangements["inv_read"]
        inv_left = rearrangements["inv_left"]
        inv_right = rearrangements["inv_right"]
        inv_size = rearrangements["inv_size"]

        # inversions
        take = (
            (inv_size >= args.max_gap)
            & (inv_left // binsize < Leff)
            & (inv_right // binsize < Leff)
        )
        bp_hist_inv = scipy.sparse.csr_matrix(
            (
                weights[inv_read[take]],
                (
                    np.maximum(inv_left[take], inv_right[take]) // binsize,
                    np.minimum(inv_left[take], inv_right[take]) // binsize,
                ),
            ),
            shape=(Leff, Leff),
        ).todense()
        np.nan_to_num(bp_hist_inv, copy=False)

        dup_read = rearrangements["dup_read"]
        dup_left = rearrangements["dup_left"]
        dup_right = rearrangements["dup_right"]
        dup_size = rearrangements["dup_size"]

        # duplications
        take = (
            (dup_size >= args.max_gap)
            & (dup_left // binsize < Leff)
            & (dup_right // binsize < Leff)
        )
        bp_hist_dup = scipy.sparse.csr_matrix(
            (
                weights[dup_read[take]],
                (
                    np.maximum(dup_left[take], dup_right[take]) // binsize,
                    np.minimum(dup_left[take], dup_right[take]) // binsize,
                ),
            ),
            shape=(Leff, Leff),
        ).todense()
        np.nan_to_num(bp_hist_dup, copy=False)

    else:
        bp_hist_inv = np.matrix(np.zeros((Leff, Leff)))
        bp_hist_dup = np.matrix(np.zeros((Leff, Leff)))

    stats = {}

    xx, yy = np.meshgrid(np.arange(Leff), np.arange(Leff), indexing="ij")

    nbreaks = bp_hist.sum()
    ninversions = bp_hist_inv.sum()
    nduplications = bp_hist_dup.sum()
    stats["breaks_normalized"] = nbreaks
    stats["frac_breaks_inversions"] = ninversions / (nbreaks + ninversions + nduplications)
    stats["frac_breaks_duplications"] = nduplications / (nbreaks + ninversions + nduplications)

    single_event = ((switch_iis[xx] == "SM") | (switch_iis[yy] == "SM")) & (
        switch_iis[xx] != switch_iis[yy]
    )
    multiple_event = (
        (switch_iis[xx] != "SM") & (switch_iis[yy] != "SM") & (switch_iis[xx] != switch_iis[yy])
    )
    within_event = switch_iis[xx] == switch_iis[yy]

    stats["frac_breaks_single"] = bp_hist[single_event].sum() / nbreaks
    stats["frac_breaks_multiple"] = bp_hist[multiple_event].sum() / nbreaks
    stats["frac_breaks_within"] = bp_hist[within_event].sum() / nbreaks
    stats["frac_breaks_inversions_within"] = bp_hist_inv[within_event].sum() / ninversions
    stats["frac_breaks_duplications_within"] = bp_hist_dup[within_event].sum() / nduplications

    # collapse different gamma and alpha isotypes for frac_breaks, spread and homology scores
    switch_iis = np.array([si[:2] for si in switch_iis])
    regions = np.unique(switch_iis)
    if switch_orientation == "-":
        regions = regions[::-1]

    for i in range(len(regions)):
        r1 = regions[i]
        for j in range(i + 1):
            r2 = regions[j]
            take = (switch_iis[xx] == r1) & (switch_iis[yy] == r2) | (switch_iis[yy] == r1) & (
                switch_iis[xx] == r2
            )
            stats["frac_breaks_{1}_{0}".format(r1, r2)] = bp_hist[take].sum() / (nbreaks)

    bps = (bp_hist + bp_hist.T).sum(1).A1
    for sr in np.unique(switch_iis):
        bps_here = bps[switch_iis == sr]
        pos = binsize * np.arange(np.sum(switch_iis == sr))
        m = np.sum(pos * bps_here) / np.sum(bps_here)
        m2 = np.sum(pos**2 * bps_here) / np.sum(bps_here)
        stats["spread_" + sr] = np.sqrt(m2 - m**2)

    if args.homology:
        homology = np.load(args.homology)
        for n in parse_range(args.range):
            stats["homology_fw"] = (
                np.asarray(bp_hist) * homology["fw_{0}".format(n)]
            ).sum() / np.sum(bp_hist)
            stats["homology_rv"] = (
                np.asarray(bp_hist) * homology["rv_{0}".format(n)]
            ).sum() / np.sum(bp_hist)
        r1 = regions[0]
        for r2 in regions[1:]:
            take = (switch_iis[xx] == r1) & (switch_iis[yy] == r2) | (switch_iis[yy] == r1) & (
                switch_iis[xx] == r2
            )
            for n in parse_range(args.range):
                stats["homology_fw_{0}_{1}".format(r1, r2)] = np.sum(
                    bp_hist[take].A1 * homology["fw_{0}".format(n)][take]
                ) / np.sum(bp_hist[take])
                stats["homology_rv_{0}_{1}".format(r1, r2)] = np.sum(
                    bp_hist[take].A1 * homology["rv_{0}".format(n)][take]
                ) / np.sum(bp_hist[take])

    if args.motifs:
        motif_counts = np.load(args.motifs)
        for motif in motif_counts.files:
            counts = motif_counts[motif]
            stats["donor_score_{0}".format(motif)] = np.sum(bp_hist.sum(0).A1 * counts) / (
                bp_hist.sum(0).A1.sum() * counts.sum()
            )
            stats["receiver_score_{0}".format(motif)] = np.sum(bp_hist.sum(1).A1 * counts) / (
                bp_hist.sum(1).A1.sum() * counts.sum()
            )
        r1 = regions[0]
        take1 = switch_iis == r1
        for r2 in regions[1:]:
            take2 = switch_iis == r2
            for motif in motif_counts.files:
                counts = motif_counts[motif]
                stats["donor_score_{0}_{1}_{2}".format(motif, r1, r2)] = np.sum(
                    bp_hist.sum(0).A1[take1] * counts[take1]
                ) / (bp_hist.sum(0).A1[take1].sum() * counts[take1].sum())
                stats["donor_score_{0}_{1}_{2}".format(motif, r1, r2)] = np.sum(
                    bp_hist.sum(1).A1[take2] * counts[take2]
                ) / (bp_hist.sum(1).A1[take2].sum() * counts[take2].sum())

    logger.info("saving results to {0}\n".format(args.out))
    pd.Series(stats).to_csv(args.out, header=False)

    if args.plot:
        logger.info("creating figure and saving to {0}\n".format(args.plot))
        from matplotlib import pyplot as plt
        from matplotlib.markers import MarkerStyle

        major_ticks = []
        minor_ticks = []
        minor_labels = []
        for rec in anno_recs:
            start = shift_coord(int(rec[3][1]), cov_int) - eff_start
            end = shift_coord(int(rec[3][2]), cov_int) - eff_start
            major_ticks += [start, end]
            minor_ticks.append((start + end) / 2)
            minor_labels.append(rec[3][3])

        fig = plt.figure(figsize=(4, 5.5))

        fig.text(0.535, 0.99, args.sample, size="large", ha="center", va="top")

        scale_factor = args.scale_factor
        assert Leff % scale_factor == 0, "Leff is not a multiple of scale_factor"
        bph_p = scipy.sparse.csr_matrix(
            np.asarray(bp_hist.T)
            .reshape((Leff, Leff // scale_factor, scale_factor))
            .sum(-1)
            .reshape((Leff // scale_factor, scale_factor, Leff // scale_factor))
            .sum(1)
        )
        bph_p.eliminate_zeros()
        bph_i = scipy.sparse.csr_matrix(
            np.asarray(bp_hist_inv.T)
            .reshape((Leff, Leff // scale_factor, scale_factor))
            .sum(-1)
            .reshape((Leff // scale_factor, scale_factor, Leff // scale_factor))
            .sum(1)
        )
        bph_i.eliminate_zeros()
        bph_d = scipy.sparse.csr_matrix(
            np.asarray(bp_hist_dup.T)
            .reshape((Leff, Leff // scale_factor, scale_factor))
            .sum(-1)
            .reshape((Leff // scale_factor, scale_factor, Leff // scale_factor))
            .sum(1)
        )
        bph_d.eliminate_zeros()

        ax = fig.add_axes([0.12, 0.32, 0.83, 0.6])
        ax.scatter(
            bph_p.nonzero()[1],
            bph_p.nonzero()[0],
            c=np.log(bph_p.data),
            cmap=plt.cm.Greens,
            marker="s",
            s=12,
            linewidths=0.2,
            edgecolors="k",
        )
        ax.scatter(
            bph_i.nonzero()[1],
            bph_i.nonzero()[0],
            c=np.log(bph_i.data),
            cmap=plt.cm.Reds,
            marker=MarkerStyle("o", fillstyle="right"),
            s=12,
            linewidths=0.2,
            edgecolors="k",
        )
        ax.scatter(
            bph_d.nonzero()[1],
            bph_d.nonzero()[0],
            c=np.log(bph_d.data),
            cmap=plt.cm.Blues,
            marker=MarkerStyle("o", fillstyle="left"),
            s=12,
            linewidths=0.2,
            edgecolors="k",
        )
        ax.set_xlim([Leff // scale_factor, 0])
        ax.set_ylim([Leff // scale_factor, 0])
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="gray", lw=0.5)
        ax.set_yticks(np.array(major_ticks) // (binsize * scale_factor))
        ax.set_xticks(np.array(major_ticks) // (binsize * scale_factor))
        ax.set_yticks(np.array(minor_ticks) // (binsize * scale_factor), minor=True)
        ax.set_xticks(np.array(minor_ticks) // (binsize * scale_factor), minor=True)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_yticklabels(minor_labels, minor=True)
        ax.tick_params(which="minor", length=0)
        ax.grid(alpha=0.5, color="lightgrey", which="major")
        ax.set_title("2D breakpoint histogram", size="medium")

        ax = fig.add_axes([0.12, 0.07, 0.83, 0.15])
        ax.plot(np.arange(Leff), (bp_hist + bp_hist.T).mean(0).A1, "g-", lw=0.5)
        ax.plot(np.arange(Leff), -(bp_hist_inv + bp_hist_inv.T).mean(0).A1, "r-", lw=0.5)
        ax.plot(np.arange(Leff), (bp_hist_dup + bp_hist_dup.T).mean(0).A1, "b-", lw=0.5)
        ax.set_xlim([Leff, 0])
        ax.set_xticks(np.array(major_ticks) // binsize)
        ax.set_xticks(np.array(minor_ticks) // binsize, minor=True)
        ax.set_xticklabels([])
        ax.set_xticklabels(minor_labels, rotation=90, minor=True)
        ax.set_yticks([])
        ax.set_ylabel("frequency")
        ax.tick_params(which="minor", length=0)
        ax.grid(alpha=0.5, color="lightgrey", which="major")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_title("1D breakpoint histogram", size="medium")

        fig.savefig(args.plot, dpi=300)

    if args.reference and args.top_donor and args.top_receiver:
        genome = pysam.FastaFile(args.reference)

        switch_seqs = [genome.fetch(switch_chrom, ci[0], ci[1]).upper() for ci in cov_int]
        stot = "".join(switch_seqs)

        top_donor_bins = np.argsort(bp_hist.sum(0).A1)[-args.ntop :]
        top_donor_bins.sort()
        top_donor_seqs = [stot[k * binsize : (k + 1) * binsize] for k in top_donor_bins]
        logger.info("saving top donor sequences to " + args.top_donor)
        with open(args.top_donor, "w") as outf:
            outf.write(
                "\n".join(">donor_{0}\n{1}".format(k + 1, s) for k, s in enumerate(top_donor_seqs))
                + "\n"
            )

        top_receiver_bins = np.argsort(bp_hist.sum(1).A1)[-args.ntop :]
        top_receiver_bins.sort()
        top_receiver_seqs = [stot[k * binsize : (k + 1) * binsize] for k in top_receiver_bins]
        logger.info("saving top receiver sequences to " + args.top_receiver)
        with open(args.top_receiver, "w") as outf:
            outf.write(
                "\n".join(
                    ">receiver_{0}\n{1}".format(k + 1, s) for k, s in enumerate(top_receiver_seqs)
                )
                + "\n"
            )
