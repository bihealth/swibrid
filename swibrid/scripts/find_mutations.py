"""find mutations in MSA"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""file with  (pseudo) multiple alignment of read sequences for clustering""",
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
        "--clustering",
        dest="clustering",
        help="""find_clusters.py output file""",
    )
    parser.add_argument(
        "--reference",
        dest="reference",
        help="""reference sequence (switch chromosome)""",
    )
    parser.add_argument(
        "--last",
        dest="last",
        help="""LAST parameters to estimate mutation probability""",
    )
    parser.add_argument(
        "--fdr", dest="fdr", default=0.05, help="""FDR for mutation detection"""
    )
    parser.add_argument(
        "--min_cov",
        dest="min_cov",
        default=50,
        type=int,
        help="""minimum read coverage at potentially mutated positions [50]""",
    )
    parser.add_argument(
        "--gap_window_size",
        dest="gap_window_size",
        default=10,
        type=int,
        help="""size of window to compute local gap frequency [10]""",
    )
    parser.add_argument(
        "--max_local_gap_freq",
        dest="max_local_gap_freq",
        default=0.7,
        type=float,
        help="""max. local gap frequency in window [.7]""",
    )
    parser.add_argument("-o", "--out", dest="out", help="""output file""")


def run(args):
    import numpy as np
    import pandas as pd
    import scipy.cluster.hierarchy
    import scipy.optimize
    import scipy.stats
    import pysam
    from logzero import logger
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
        get_switch_coverage,
        intersect_intervals,
        parse_LAST_pars,
        p_adjust_bh,
        ncodes,
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

    cov_map = np.zeros(Ltot, dtype=int)
    n = 0
    for start, end in cov_int:
        d = end - start
        cov_map[np.arange(d) + n] = np.arange(start, end).astype(int)
        n += d

    clustering = pd.read_csv(args.clustering, index_col=0, header=0)

    reads = clustering["cluster"].dropna().index
    nreads = len(reads)
    clustering = clustering.loc[reads]

    logger.info("loading msa from {0}".format(args.msa))
    msa = scipy.sparse.load_npz(args.msa).tocsr()
    if nreads < msa.shape[0]:
        logger.info("restricting msa to {0} reads in clustering".format(nreads))
        msa = msa[:nreads]

    # get LAST pars
    logger.info("reading LAST pars from {0}".format(args.last))
    pars = parse_LAST_pars(args.last)
    m0 = pars["p_c"] / pars["p_c"].sum(1)[:, np.newaxis]

    logger.info("looking for mutations")
    reference = pysam.FastaFile(args.reference)
    ref_seq = "".join(
        reference.fetch(switch_chrom, s, e) for s, e in cov_int
    ).upper()
    # get reference in integer coding array format
    rr = np.array([ncodes[c] for c in ref_seq])

    # get average gap frequency in a window around each nuclotide
    gap_freq = (
        np.apply_along_axis(
            np.convolve,
            1,
            np.asarray(~(msa != 0).todense()),
            np.ones(args.gap_window_size),
            mode="same",
        )
        / args.gap_window_size
    )

    # find all positions where reads are different from reference and not gap, and also not near gaps (<70% local gap frequency)
    all_diff = (
        (msa != 0).todense()
        & (msa != rr[np.newaxis, :])
        & (gap_freq < args.max_local_gap_freq)
    )

    # compare distribution of nucleotides at each position against reference and mutation prob under a chisquare model
    nuc_dist = np.apply_along_axis(
        lambda x: np.bincount(np.abs(x), minlength=5),
        axis=0,
        arr=np.asarray(msa.todense()),
    ).T[:, 1:]
    # number of non-gaps at each position
    nr = nuc_dist.sum(1)
    pmut = scipy.stats.chisquare(
        nuc_dist, nr[:, np.newaxis] * m0[rr - 1], axis=1
    )[1]
    # exclude positions with less than min_cov non-gaps
    pmut[nr < args.min_cov] = 1

    # now let's check how many mutations we see and how many we expect under the null model
    nmut = all_diff.sum(0).A1
    nmut_exp = nr * (1 - m0[rr - 1].max(1))

    # set pmut to 1 for all positions where we see fewer mutations than expected
    pmut[nmut <= nmut_exp] = 1

    # multiple testing correction using Benjamini-Hochberg
    padj = p_adjust_bh(pmut)

    SNP_pos = np.where(padj < args.fdr)[0]

    # check strand bias
    npos = (
        all_diff[clustering["orientation"].values == "+", :][:, SNP_pos]
        .sum(0)
        .A1
    )
    nneg = (
        all_diff[clustering["orientation"].values == "-", :][:, SNP_pos]
        .sum(0)
        .A1
    )
    strand_bias = npos / (npos + nneg)

    # check distribution across clusters and find clusters where most reads are mutated
    pclust = np.zeros(len(SNP_pos))
    mut_clust = {}
    for n, pos in enumerate(SNP_pos):
        take = (msa[:, pos] != 0).todense().A1
        cdist = np.bincount(clustering["cluster"][take])
        p = cdist > 0
        cdist_mut = np.bincount(
            clustering["cluster"][all_diff[:, pos].A1], minlength=len(cdist)
        )
        cdist_ref = cdist * cdist_mut.sum() / nr[pos]
        pclust[n] = scipy.stats.chisquare(cdist_mut[p], cdist_ref[p])[1]
        mut_clust[pos] = np.where(
            (cdist_mut > cdist_ref)
            & p
            & (
                p_adjust_bh(1 - scipy.stats.poisson.cdf(cdist_mut, cdist_ref))
                < 0.05
            )
        )[0]

    pclust_adj = p_adjust_bh(pclust)

    # test distribution of inversions across clusters
    inversions = (msa < 0).todense().any(1).A1
    cdist = np.bincount(clustering["cluster"])
    cdist_inv = np.bincount(
        clustering["cluster"][inversions], minlength=len(cdist)
    )
    cdist_ref = cdist * cdist_inv.sum() / cdist.sum()
    pinv = scipy.stats.chisquare(cdist_inv, cdist_ref)[1]

    logger.info("saving results to {0}".format(args.out))
    with open(args.out, "w") as outf:
        msa = msa.todense()
        outf.write(
            "chrom\tposition\tregion\trel_pos\tref\talt_counts\tpval\tpadj\tpval_clust\tpadj_clust\tstrand_bias\tclusts\talt_A\talt_C\talt_G\talt_T\n"
        )
        lines = []
        line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.3g}\t{7:.3g}\t{8:.3g}\t{9:.3g}\t{10:.3g}\t{11}\t{12}\t{13}\t{14}\t{15}\n"
        for n, pos in enumerate(SNP_pos):
            ref = ref_seq[pos]
            real_pos = cov_map[pos]
            region = []
            for _, _, _, hit in intersect_intervals(
                [("chr14", real_pos, real_pos + 1)], switch_anno, loj=True
            ):
                region.append(hit[3])
            region = ",".join(region)
            vals = [
                switch_chrom,
                real_pos,
                region,
                pos,
                ref,
                ",".join(map(str, nuc_dist[pos])),
                pmut[pos],
                padj[pos],
                pclust[n],
                pclust_adj[n],
                strand_bias[n],
                ",".join(map(str, mut_clust[pos])),
            ] + [
                ",".join(
                    map(
                        str,
                        np.where((all_diff[:, pos] & (msa[:, pos] == k)))[0],
                    )
                )
                for k in range(1, 5)
            ]
            lines.append(line.format(*vals))
        outf.writelines(lines)
