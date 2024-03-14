"""\
EXPERIMENTAL: find single-nucleotide variants in MSA and determine haplotypes

single-nucleotide variants: variable positions in the input MSA (excluding positions at or near gaps) where
nucleotide distributions are different than what's expected given mutation frequencies (estimated at the initial alignment step)
potential variants are removed if
- there's less than min_cov non-gaps
- fewer variants than expected
- high strand bias
- no cluster with at least min_cluster_cov reads and allele frequency > min_freq
variants are aggregated over clusters, and the distribution across clusters is tested for evenness
variants can be annotated by dbSNP id
motif occurrences around variants are also scored

haplotypes are determined by performing a weighted NMF on a matrix of allele frequencies per cluster and variant

output is a vcf-style file with genomic (1-based) and relative (0-based) coordinates, and a matrix (sparse integer array, same shape as MSA) indicating which read contains a variant at which position 
"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""required: file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument(
        "--clustering",
        dest="clustering",
        help="""required: clustering (find_clusters output)""",
    )
    parser.add_argument(
        "--stats",
        dest="stats",
        help="""required: clustering stats (find_clusters output)""",
    )
    parser.add_argument(
        "--reference",
        dest="reference",
        help="""required: reference sequence (switch chromosome)""",
    )
    parser.add_argument(
        "--pars",
        dest="pars",
        help="""required: alignment parameters to estimate mutation probability""",
    )
    parser.add_argument("-o", "--out", dest="out", help="""required: output file (text)""")
    parser.add_argument("-m", "--mat", dest="mat", help="""required: output file (matrix)""")
    parser.add_argument("--out_complete", dest="out_complete", help="""required: output file with all variants (text)""")
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
    parser.add_argument("--fdr", dest="fdr", default=0.05, help="""FDR for variant calling""")
    parser.add_argument(
        "--min_cov",
        dest="min_cov",
        default=50,
        type=int,
        help="""minimum read coverage at potentially variable positions [50]""",
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
    parser.add_argument(
        "--min_cluster_cov",
        dest="min_cluster_cov",
        default=10,
        type=int,
        help="""minimum read coverage at potentially variable positions per cluster [10]""",
    )
    parser.add_argument(
        "--min_freq",
        dest="min_freq",
        default=0.4,
        type=float,
        help="""minimum allele frequency at potentially variable positions (in total or per cluster) [.4]""",
    )
    parser.add_argument(
        "--variant_annotation",
        dest="variant_annotation",
        nargs='?',
        help="""variant annotation file (vcf.gz; e.g., from 1000Genomes project)""",
    )
    parser.add_argument(
        "--haplotypes",
        dest="haplotypes",
        help="""cluster haplotypes""",
    )
    parser.add_argument(
        "--motifs",
        dest="motifs",
        default="Cg,wrCy,Tw",
        help="""comma-separated list of sequence motifs to look up at variant positions [Cg,wrCy,Tw]""",
    )


def run(args):
    import os
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
        p_adjust_bh,
        vrange,
        nmf_consistency,
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

    # remove coverage and orientation info from msa
    msa.data = msa.data % 10

    # get alignment pars
    logger.info("reading alignment pars from {0}".format(args.pars))
    pars = np.load(args.pars)
    m0 = pars["p_c"] / pars["p_c"].sum(1)[:, np.newaxis]

    logger.info("looking for variants")
    reference = pysam.FastaFile(args.reference)
    ref_seq = "".join(reference.fetch(switch_chrom, s, e) for s, e in cov_int).upper()
    # get reference in integer coding array format
    ref = np.array([ncodes[c] for c in ref_seq])

    # find all positions where reads are different from reference and not gap
    # and also not near gaps (<70% local gap frequency)

    all_diff = scipy.sparse.csr_matrix(
        ((msa.data != ref[msa.nonzero()[1]]) & (msa.data != 0), msa.nonzero()),
        shape=msa.shape,
    )
    all_diff.eliminate_zeros()
    i, j = all_diff.nonzero()
    ii = np.repeat(i, args.gap_window_size)
    j[j < args.gap_window_size // 2] = args.gap_window_size // 2
    j[j >= msa.shape[1] - args.gap_window_size // 2] = msa.shape[1] - args.gap_window_size // 2
    jj = vrange(j - args.gap_window_size // 2, j + args.gap_window_size // 2)

    gap_freq = 1 - (msa != 0)[(ii, jj)].reshape(len(i), args.gap_window_size).mean(1).A1
    all_diff.data[gap_freq >= args.max_local_gap_freq] = False
    all_diff.eliminate_zeros()

    logger.info(
        "{0} positions ignored due to high gap frequency".format(
            sum(gap_freq >= args.max_local_gap_freq)
        )
    )

    nuc_dist = np.vstack([np.sum(msa == k, 0).A1 for k in range(1, 5)]).T
    # number of non-gaps at each position
    nr = nuc_dist.sum(1)
    # find most frequent alternative allele at each position
    nd = nuc_dist.copy()
    nd[(np.arange(nuc_dist.shape[0]), ref - 1)] = 0
    alt = nd.argmax(1) + 1

    #  compare distribution against reference-dependent mutation probability under a chisquare model
    pvar = scipy.stats.chisquare(nuc_dist, nr[:, np.newaxis] * m0[ref - 1], axis=1)[1]

    # check how many variants we see and how many we expect under the null model
    nvar = all_diff.sum(0).A1
    nvar_exp = nr * m0[ref - 1, alt - 1]

    # calculate another p-value for having more variants than expected
    pp = scipy.stats.poisson.sf(nvar, nvar_exp)
    pp_adj = p_adjust_bh(pp)

    # multiple testing correction using Benjamini-Hochberg
    padj = p_adjust_bh(pvar)

    # check strand bias
    npos = all_diff[clustering["orientation"].values == "+", :].sum(0).A1
    nneg = all_diff[clustering["orientation"].values == "-", :].sum(0).A1
    pstrand = scipy.stats.chisquare(np.array([npos, nneg]), axis=0)[1]

    logger.info("aggregating counts over clusters")
    i, j = all_diff.multiply(msa == alt).nonzero()
    # create matrix with only positions that have these frequent non-reference nucleotides
    mat = scipy.sparse.csc_matrix((np.ones_like(i), (i, j)), shape=msa.shape, dtype=np.int8)

    clusters, cind, cinv, csize = np.unique(
        clustering["cluster"].dropna(), return_index=True, return_inverse=True, return_counts=True
    )
    clones = np.unique(clustering["filtered_cluster"][clustering["filtered_cluster"] >= 0])
    data = np.ones_like(cinv)
    data[~np.isin(cinv, clones) | (csize[cinv] < args.min_cluster_cov)] = 0

    mm = scipy.sparse.csr_matrix(
        (data, (cinv, np.arange(len(cinv)))),
        shape=(len(clusters), len(cinv)),
    )
    mm.eliminate_zeros()

    # A is counts of alternative allele per cluster, D is counts of reference allele per cluster
    A = mm.dot(mat)
    D = mm.dot(msa != 0)

    # set padj to 1 for all positions that have:
    # - less than min_cov non-gaps
    # - fewer variants than expected
    # - high strand bias
    # - no cluster with at least min_cluster_cov reads and allele frequency > args.min_freq
    logger.info("filtering variants")
    low_cov = nr < args.min_cov
    few_var = pp_adj > args.fdr
    stranded = pstrand < 0.05
    if D.nnz > 0:
        dnz = D.nonzero()
        no_clust_data = (A[dnz] > args.min_freq * D[dnz]).A1 & (D[dnz] >= args.min_cluster_cov).A1
        no_clust = ~(scipy.sparse.csr_matrix((no_clust_data, dnz), shape=D.shape).sum(0) > 0).A1
    else:
        no_clust = np.ones_like(padj).astype(bool)
    freq = mat.sum(0).A1 / (msa != 0).sum(0).A1
    no_clust = no_clust & (freq <= args.min_freq)

    padj_complete = padj.copy()
    padj_complete[few_var] = 1
    SNP_pos_complete = np.where(padj_complete < args.fdr)[0]

    padj[low_cov | few_var | stranded | no_clust] = 1
    SNP_pos = np.where(padj < args.fdr)[0]
    logger.info("{0} variants removed due to low coverage".format(sum(low_cov)))
    logger.info("{0} variants removed due to low # of events".format(sum(few_var)))
    logger.info("{0} variants removed due to strand bias".format(sum(stranded)))
    logger.info("{0} variants removed due to low allele frequency".format(sum(no_clust)))
    logger.info("{0} variants kept".format(len(SNP_pos)))

    logger.info("testing distribution across clusters")
    # check distribution across clusters with at least args.min_cluster_cov reads
    pclust = np.zeros(len(SNP_pos))

    for n, pos in enumerate(SNP_pos):
        take = (msa[:, pos] != 0).todense().A1
        cdist = np.bincount(clustering["cluster"][take])
        p = cdist > args.min_cluster_cov
        cdist_var = np.bincount(
            clustering["cluster"][all_diff[:, pos].todense().A1], minlength=len(cdist)
        )
        cdist_ref = cdist * cdist_var.sum() / nr[pos]
        cdist_ref = cdist[p] * cdist_var[p].sum() / cdist[p].sum()
        pclust[n] = scipy.stats.chisquare(cdist_var[p], cdist_ref)[1]

    pclust_adj = p_adjust_bh(pclust)

    # test distribution of inversions across clusters
    # inversions = (msa < 0).sum(1).A1 > 0
    # cdist = np.bincount(clustering["cluster"])
    # cdist_inv = np.bincount(clustering["cluster"][inversions], minlength=len(cdist))
    # cdist_ref = cdist * cdist_inv.sum() / cdist.sum()
    # pinv = scipy.stats.chisquare(cdist_inv, cdist_ref)[1]

    anno = [""] * len(SNP_pos)
    if args.variant_annotation and os.path.isfile(args.variant_annotation):
        logger.info("annotating variants with " + args.variant_annotation)

        var_anno = pd.read_csv(
            args.variant_annotation,
            sep="\t",
            comment="#",
            header=None,
            index_col=1,
            names=["CHROM", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
        )
        # keep only PASS variants and remove polyclonal variants
        var_anno = var_anno[
            (var_anno["FILTER"] == "PASS")
            & (var_anno["REF"].str.len() == 1)
            & (var_anno["ALT"].str.len() == 1)
        ]

        for n, pos in enumerate(SNP_pos):
            real_pos = cov_map[pos] + 1
            if (
                real_pos in var_anno.index
                and var_anno.loc[real_pos, "REF"] == "ACGT"[ref[pos] - 1]
                and var_anno.loc[real_pos, "ALT"] == "ACGT"[alt[pos] - 1]
            ):
                anno[n] = var_anno.loc[real_pos, "ID"]

    mtype = np.array(["n.d."] * len(SNP_pos))

    if args.haplotypes is not None:
        use_clusts = np.isin(clusters, clones) & (csize >= args.min_cluster_cov)

        if use_clusts.sum() > 2 and len(SNP_pos) > 1:
            logger.info(
                "clustering haplotypes for {0} clusters via wNMF using {1} SNPs".format(
                    sum(use_clusts), len(SNP_pos)
                )
            )

            from wNMF import wNMF

            model = wNMF(
                n_components=2,
                beta_loss="kullback-leibler",
                max_iter=1000,
                verbose=False,
                track_error=True,
            )
            fit = model.fit(
                X=np.asarray(A[use_clusts][:, SNP_pos].todense().astype(float)),
                W=np.asarray(D[use_clusts][:, SNP_pos].todense().astype(float)),
                n_run=100,
            )

            haplotypes = 0.5 * np.ones_like(clusters)
            haplotypes[use_clusts] = np.argmax(fit.U, axis=1)
            consistency = np.zeros_like(clusters)
            consistency[use_clusts] = nmf_consistency(fit.U, fit.U_all)

            af0 = A[haplotypes == 0].sum(0).A1[SNP_pos] / D[haplotypes == 0].sum(0).A1[SNP_pos]
            af1 = A[haplotypes == 1].sum(0).A1[SNP_pos] / D[haplotypes == 1].sum(0).A1[SNP_pos]

            mtype[(af0 > 0.6) & (af1 < 0.6)] = "het0"
            mtype[(af0 < 0.6) & (af1 > 0.6)] = "het1"
            mtype[(af0 > 0.6) & (af1 > 0.6)] = "hom"

        else:
            haplotypes = 0.5 * np.ones_like(clusters)
            consistency = np.zeros_like(clusters)

        logger.info("saving haplotypes to " + args.haplotypes)
        pd.DataFrame({"haplotype": haplotypes, "consistency": consistency}, index=clusters).to_csv(
            args.haplotypes
        )

    if args.motifs is not None:
        logger.info("checking motifs")

        import re
        from collections import defaultdict
        from .utils import IUPAC, RC

        motif_anno = defaultdict(set)

        for mot in args.motifs.split(","):
            rx = r"".join(IUPAC[m.upper()] for m in mot)
            i = re.search(r"[A-Z]", mot).start()
            for p in SNP_pos:
                seq = ref_seq[(p - i) : (p - i + len(mot))]
                seq_rc = "".join(RC[s] for s in ref_seq[(p - len(mot) + i + 1) : (p + i + 1)][::-1])
                if re.match(rx, seq) or re.match(rx, seq_rc):
                    motif_anno[p].add(mot)

        motif_anno = np.array([",".join(motif_anno[p]) for p in SNP_pos])

    else:
        motif_anno = np.array([""] * len(SNP_pos))

    logger.info("saving variant table to {0}".format(args.out))
    with open(args.out, "w") as outf:
        outf.write(
            "chrom\tposition\tregion\trel_pos\tref\talt\tcounts\tpval\tpadj\tpval_clust\tpadj_clust\tpval_strand\ttype\tanno\tmotif\n"
        )
        lines = []
        line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:.3g}\t{8:.3g}\t{9:.3g}\t{10:.3g}\t{11:.3g}\t{12}\t{13}\t{14}\n"
        for n, pos in enumerate(SNP_pos):
            ref = ref_seq[pos]
            real_pos = cov_map[pos] + 1
            region = []
            for _, _, _, hit in intersect_intervals(
                [("chr14", real_pos - 1, real_pos)], switch_anno, loj=True
            ):
                region.append(hit[3])
            region = ",".join(region)
            vals = [
                switch_chrom,
                real_pos,
                region,
                pos,
                ref,
                "ACGT"[alt[pos] - 1],
                ",".join(map(str, nuc_dist[pos])),
                pvar[pos],
                padj[pos],
                pclust[n],
                pclust_adj[n],
                pstrand[pos],
                mtype[n],
                anno[n],
                motif_anno[n],
            ]
            lines.append(line.format(*vals))
        outf.writelines(lines)

    if args.mat is not None:
        logger.info("saving variant matrix to {0}".format(args.mat))
        scipy.sparse.save_npz(args.mat, mat)

    if args.out_complete:
        logger.info("saving complete variant table to {0}".format(args.out_complete))
        with open(args.out_complete, "w") as outf:
            outf.write(
                "chrom\tposition\tregion\trel_pos\tref\talt\tcounts\tpval\tpadj\tcov\tpp_adj\tpstrand\tfreq\tused\n"
            )
            lines = []
            line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:.3g}\t{8:.3g}\t{9}\t{10:.3g}\t{11:.3g}\t{12:.3g}\t{13}\n"
            for n, pos in enumerate(SNP_pos_complete):
                ref = ref_seq[pos]
                real_pos = cov_map[pos] + 1
                region = []
                for _, _, _, hit in intersect_intervals(
                    [("chr14", real_pos - 1, real_pos)], switch_anno, loj=True
                ):
                    region.append(hit[3])
                region = ",".join(region)
                vals = [
                    switch_chrom,
                    real_pos,
                    region,
                    pos,
                    ref,
                    "ACGT"[alt[pos] - 1],
                    ",".join(map(str, nuc_dist[pos])),
                    pvar[pos],
                    padj_complete[pos],
                    nr[pos],
                    pp_adj[pos],
                    pstrand[pos],
                    freq[pos],
                    pos in SNP_pos,
                ]
                lines.append(line.format(*vals))
            outf.writelines(lines)

