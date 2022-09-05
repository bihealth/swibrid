"""analyze clustering"""


def setup_argparse(parser):

    parser.add_argument("-o", "--output", dest="output", help="""output""")
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""file(s) with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("--gaps", dest="gaps", help="""gap distribution""")
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


def run(args):

    import numpy as np
    import scipy.sparse
    import scipy.cluster.hierarchy
    import pandas as pd

    from logzero import logger

    from .helpers import (
        parse_switch_coords,
        get_switch_coverage,
        read_switch_anno,
        construct_mut_matrix,
        remove_gaps,
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

    logger.info("reading clustering")
    clustering = pd.read_csv(args.clustering, index_col=0, header=0)
    clusters, cinv, csize = np.unique(
        clustering["cluster"].dropna(), return_inverse=True, return_counts=True
    )
    avg_isotype = (
        clustering.groupby("cluster")["isotype"]
        .agg(pd.Series.mode)
        .loc[clusters]
    )

    logger.info("loading gaps from " + args.gaps)
    gaps = np.load(args.gaps)
    logger.info("loading MSA from " + args.msa)
    msa = scipy.sparse.load_npz(args.msa)
    logger.info("removing gaps > {0} from MSA".format(args.max_gap))
    msa_cleaned = remove_gaps(msa, gaps=gaps, max_gap=args.max_gap)
    logger.info("adding mutations from {0}".format(args.mutations))
    mutations = pd.read_csv(args.mutations, sep="\t", header=0)
    mut = construct_mut_matrix(mutations, msa.shape[0], msa.shape[1])

    mm = scipy.sparse.csr_matrix(
        (1.0 / csize[cinv], (cinv, np.arange(len(cinv)))),
        shape=(len(clusters), len(cinv)),
    )

    logger.info("averaging MSA")
    avg_msa = mm.dot(np.abs(msa_cleaned))
    cluster_length = avg_msa.sum(1)
    break_spread = np.sum((avg_msa > 0) & (avg_msa < 0.95), 1) / cluster_length

    logger.info("averaging mutations")
    avg_mut = mm.dot(mut > 0)
    mut_positions = (avg_mut > 0).sum(1).A1
    tmp = avg_mut.copy()
    tmp.data[(tmp.data <= 0.05) | (tmp.data >= 0.95)] = 0
    tmp.eliminate_zeros()
    mut_spread = np.sum(tmp > 0, 1).A1 / mut_positions

    df = pd.DataFrame(
        {
            "size": csize,
            "isotype": avg_isotype,
            "length": cluster_length,
            "break_spread": break_spread,
            "mut_spread": mut_spread,
        },
        index=clusters,
    )

    switch_iis = get_switch_iis(anno_recs, cov_int, eff_start, 1)
    for isotype in np.unique(switch_iis):
        df["length_" + isotype] = avg_msa[:, switch_iis == isotype].sum(1)

    logger.info("saving results to " + args.output)
    df.to_csv(args.output)
