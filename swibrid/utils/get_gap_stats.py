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
    parser.add_argument('-n',
                        '--ntop',
                        dest='ntop',
                        type=int,
                        default=10,
                        help="""number of top bins to select""")
    parser.add_argument(
        '--reference',
        dest='reference',
        help="""genome fasta file"""
    )
    parser.add_argument(
        '--top_donor',
        dest='top_donor',
        help="""output file with top donor sequences"""
    )
    parser.add_argument(
        '--top_receiver',
        dest='top_receiver',
        help="""output file with top receiver sequences"""
    )


def run(args):

    import numpy as np
    import pandas as pd
    import scipy.sparse
    import pysam
    import re
    from logzero import logger
    from .helpers import (
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
        clustering = pd.read_csv(args.clustering, header=0, index_col=0)
        clustering = clustering[~clustering["cluster"].isna()]
        stats = pd.read_csv(
            args.clustering_stats, header=None, index_col=0
        ).squeeze()
        n_eff = stats["eff_nclusters"].astype(int)
        clones = (
            clustering["cluster"]
            .dropna()
            .astype(int)
            .value_counts()
            .index[:n_eff]
        )
        nc = clustering["cluster"].dropna().astype(int).value_counts()
        singletons = nc.index[nc == 1]
        w = 1.0 / nc.loc[clustering["cluster"].values]
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

    if args.reference:

        reference = pysam.FastaFile(args.reference)
        switch_seq = ''.join([
            reference.fetch(switch_chrom, ci[0], ci[1]).upper() for ci in cov_int
        ])

        simple_repeats = ['G_C','A_T','TG','CA','TCCA','CAGCC','CAGCT','CACAT','CAGAGA']
    
        for repeat in simple_repeats:
            counts = np.array([sum(1 for _ in re.finditer(repeat.replace("_",'|'), switch_seq[k*binsize:(k+1)*binsize]))
                               for k in range(len(switch_seq)//binsize)])
            stats["donor_score_" + repeat] = np.sum(bp_hist.sum(0).A1*counts) / (bp_hist.sum(0).A1.sum() * counts.sum())
            stats["receiver_score_" + repeat] = np.sum(bp_hist.sum(1).A1*counts) / (bp_hist.sum(1).A1.sum() * counts.sum())

    logger.info("saving results to {0}\n".format(args.out))
    pd.Series(stats).to_csv(args.out, header=False)

    if args.top_donor and args.top_receiver:

        top_donor_bins = np.argsort(bp_hist.sum(0).A1)[-args.ntop:]
        top_donor_bins.sort()
        top_donor_seqs=[switch_seq[k*binsize:(k+1)*binsize] for k in top_donor_bins]
        logger.info('saving top donor sequences to ' + args.top_donor)
        with open(args.top_donor,'w') as outf:
            outf.write('\n'.join('>donor_{0}\n{1}'.format(k+1,s) for k,s in enumerate(top_donor_seqs))+'\n')

        top_receiver_bins = np.argsort(bp_hist.sum(1).A1)[-args.ntop:]
        top_receiver_bins.sort()
        top_receiver_seqs=[switch_seq[k*binsize:(k+1)*binsize] for k in top_receiver_bins]
        logger.info('saving top receiver sequences to ' + args.top_receiver)
        with open(args.top_receiver,'w') as outf:
            outf.write('\n'.join('>receiver_{0}\n{1}'.format(k+1,s) for k,s in enumerate(top_receiver_seqs))+'\n')
    
