"""analyze gap statistics"""


def setup_argparse(parser):

    parser.add_argument(
        "-g",
        "--gaps",
        dest="gaps",
        help="""file with gap positions (output of get_gaps.py)""",
    )
    parser.add_argument(
        "-c",
        "--clustering",
        dest="clustering",
        help="""file with clustering results""",
    )
    parser.add_argument(
        "-s",
        "--clustering_stats",
        dest="clustering_stats",
        help="""file with clustering stats""",
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
        default="4-7",
        help="""range of kmer sizes, e.g., 3,5-7 [4-7]""",
    )
    parser.add_argument(
        "-n",
        "--ntop",
        dest="ntop",
        type=int,
        default=10,
        help="""number of top bins to select""",
    )
    parser.add_argument(
        "--reference", dest="reference", help="""genome fasta file"""
    )
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
        help="""comma-separated list of clones to use or 'all' (default: top 95% of clusters excluding singletons)""",
    )
    parser.add_argument(
        "--use_weights",
        dest="use_weights",
        default='inverse',
        help="""use different cluster weights ("inverse" | "uniform" | "adjusted") [inverse]"""
    )


def run(args):

    import numpy as np
    import pandas as pd
    import scipy.sparse
    import scipy.stats
    import pysam
    import re
    from logzero import logger
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
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

    gaps = np.load(args.gaps)
    nreads = gaps["read_idx"].max() + 1

    if args.clustering:
        logger.info("reading clustering from " + args.clustering)
        clustering = pd.read_csv(args.clustering, header=0, index_col=0)
        clustering = clustering[~clustering["cluster"].isna()]
        nreads = clustering.shape[0]
        logger.info("reading clustering stats from " + args.clustering_stats)
        stats = pd.read_csv(
            args.clustering_stats, header=0, index_col=0
        ).squeeze()
        logger.info("reading clustering analysis from " + args.clustering_analysis)
        analysis = pd.read_csv(
            args.clustering_analysis, header=0, index_col=0
        )
        neff = stats["eff_nclusters"].astype(int)
        if args.use_clones:
            if args.use_clones == "all":
                clones = clustering['cluster'].astype(int).unique()
            else:
                clones = list(map(int, args.use_clones.split(",")))
        else:
            clones = (
                clustering["cluster"]
                .dropna()
                .astype(int)
                .value_counts()
                .index[:neff]
            )
        logger.info("using {0} clones".format(len(clones)))
        nc = clustering["cluster"].dropna().astype(int).value_counts()
        singletons = nc.index[nc == 1]
        if args.use_weights == 'inverse':
            logger.info('using inverse weights')
            w = 1.0 / nc.loc[clustering["cluster"].values]
        elif args.use_weights == 'uniform':
            logger.info('using uniform weights')
            w = pd.Series(1, index=np.arange(nreads))
        elif args.use_weights == 'adjusted':
            logger.info('using adjusted weights')
            w = analysis.loc[clustering['cluster'].values,'adj_size'] / analysis.loc[clustering['cluster'].values,'size']
        w[~w.index.isin(clones) | w.index.isin(singletons)] = 0
        weights = pd.Series(w.values / w.sum(), index=np.arange(nreads))
    else:
        weights = pd.Series(np.ones(nreads) / nreads, index=np.arange(nreads))

    gap_reads = gaps["read_idx"]
    gap_left = gaps["pos_left"]
    gap_right = gaps["pos_right"]
    gap_size = gaps["gap_size"]
    Leff = Ltot // binsize
    take = (
        (gap_size >= args.max_gap)
        & (gap_left // binsize < Leff)
        & (gap_right // binsize < Leff)
    )
    bp_hist = scipy.sparse.csr_matrix(
        (
            weights[gap_reads[take]],
            (gap_left[take] // binsize, gap_right[take] // binsize),
        ),
        shape=(Leff, Leff),
    ).todense()

    stats = {}

    xx, yy = np.meshgrid(np.arange(Leff), np.arange(Leff), indexing="ij")
    iu = np.triu_indices(Leff)

    nbreaks = bp_hist[iu].A1.sum()
    stats["breaks_normalized"] = nbreaks

    single_event = ((switch_iis[xx] == "SM") | (switch_iis[yy] == "SM")) & (
        switch_iis[xx] != switch_iis[yy]
    )
    multiple_event = (
        (switch_iis[xx] != "SM")
        & (switch_iis[yy] != "SM")
        & (switch_iis[xx] != switch_iis[yy])
    )
    within_event = switch_iis[xx] == switch_iis[yy]

    stats["frac_single_break"] = bp_hist[single_event].sum() / nbreaks
    stats["frac_multiple_breaks"] = bp_hist[multiple_event].sum() / nbreaks
    stats["frac_within_break"] = bp_hist[within_event].sum() / nbreaks

    # check if certain regions in SM break to different isotypes with different frequencies
    df = {}
    for isotype in np.unique(switch_iis):
        for k, b in enumerate(np.where(switch_iis == "SM")[0]):
            df[(isotype, k)] = bp_hist[switch_iis == isotype][:, b].sum()
    m = pd.Series(df).unstack(level=0)
    stats["SM_downstream_bias"] = (
        np.nansum(
            [
                m.loc[i].sum()
                * scipy.stats.entropy(m.loc[i], m.sum(0))
                / np.log(m.shape[1])
                for i in m.index
            ]
        )
        / m.sum().sum()
    )

    if args.homology:
        homology = np.load(args.homology)
        for n in parse_range(args.range):
            stats["homology_fw_{0}".format(n)] = (
                np.sum(bp_hist[iu].A1 * homology["fw_{0}".format(n)][iu])
                / nbreaks
            )
            stats["homology_rv_{0}".format(n)] = (
                np.sum(bp_hist[iu].A1 * homology["rv_{0}".format(n)][iu])
                / nbreaks
            )

    bps = (bp_hist + bp_hist.T).sum(1).A1
    for sr in np.unique(switch_iis):
        bps_here = bps[switch_iis == sr]
        pos = binsize * np.arange(np.sum(switch_iis == sr))
        m = np.sum(pos * bps_here) / np.sum(bps_here)
        m2 = np.sum(pos**2 * bps_here) / np.sum(bps_here)
        stats["spread_" + sr] = np.sqrt(m2 - m**2)

    if args.motifs:

        motif_counts = np.load(args.motifs)

        for motif in motif_counts.files:
            counts = motif_counts[motif]
            stats["donor_score_" + motif] = np.sum(
                bp_hist.sum(0).A1 * counts
            ) / (bp_hist.sum(0).A1.sum() * counts.sum())
            stats["receiver_score_" + motif] = np.sum(
                bp_hist.sum(1).A1 * counts
            ) / (bp_hist.sum(1).A1.sum() * counts.sum())

    logger.info("saving results to {0}\n".format(args.out))
    pd.Series(stats).to_csv(args.out, header=False)

    if args.reference and args.top_donor and args.top_receiver:

        genome = pysam.FastaFile(args.reference)

        switch_seqs = [
            genome.fetch(switch_chrom, ci[0], ci[1]).upper() for ci in cov_int
        ]
        stot = "".join(switch_seqs)

        top_donor_bins = np.argsort(bp_hist.sum(0).A1)[-args.ntop :]
        top_donor_bins.sort()
        top_donor_seqs = [
            stot[k * binsize : (k + 1) * binsize] for k in top_donor_bins
        ]
        logger.info("saving top donor sequences to " + args.top_donor)
        with open(args.top_donor, "w") as outf:
            outf.write(
                "\n".join(
                    ">donor_{0}\n{1}".format(k + 1, s)
                    for k, s in enumerate(top_donor_seqs)
                )
                + "\n"
            )

        top_receiver_bins = np.argsort(bp_hist.sum(1).A1)[-args.ntop :]
        top_receiver_bins.sort()
        top_receiver_seqs = [
            stot[k * binsize : (k + 1) * binsize] for k in top_receiver_bins
        ]
        logger.info("saving top receiver sequences to " + args.top_receiver)
        with open(args.top_receiver, "w") as outf:
            outf.write(
                "\n".join(
                    ">receiver_{0}\n{1}".format(k + 1, s)
                    for k, s in enumerate(top_receiver_seqs)
                )
                + "\n"
            )
