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
    parser.add_argument(
        "--inserts", dest="inserts", help="""results.tsv file with inserts"""
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
        decode_coords,
        interval_length,
        merge_intervals,
        intersect_intervals,
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
        .agg(lambda x: pd.Series.mode(x)[0])
        .loc[clusters]
    )

    logger.info("loading gaps from " + args.gaps)
    gaps = np.load(args.gaps)
    logger.info("loading MSA from " + args.msa)
    msa = scipy.sparse.load_npz(args.msa)
    logger.info("removing gaps > {0} from MSA".format(args.max_gap))
    msa_cleaned = remove_gaps(msa, gaps=gaps, max_gap=args.max_gap)

    mm = scipy.sparse.csr_matrix(
        (1.0 / csize[cinv], (cinv, np.arange(len(cinv)))),
        shape=(len(clusters), len(cinv)),
    )

    logger.info("averaging MSA")
    avg_msa = mm.dot(np.abs(msa_cleaned))
    cluster_length = avg_msa.sum(1)
    break_spread = np.sum((avg_msa > 0) & (avg_msa < 0.95), 1) / cluster_length

    if args.mutations is not None:
        logger.info("adding mutations from {0}".format(args.mutations))
        mutations = pd.read_csv(args.mutations, sep="\t", header=0)
        mut = construct_mut_matrix(mutations, msa.shape[0], msa.shape[1])
        logger.info("averaging mutations")
        avg_mut = mm.dot(mut > 0)
        mut_positions = (avg_mut > 0).sum(1).A1
        tmp = avg_mut.copy()
        tmp.data[(tmp.data <= 0.05) | (tmp.data >= 0.95)] = 0
        tmp.eliminate_zeros()
        mut_spread = np.sum(tmp > 0, 1).A1 / mut_positions

    else:
        mut_spread = np.zeros(len(clusters))

    logger.info("analyzing GC content")
    GC = ((np.abs(msa) == 2).sum(1).A1 + (np.abs(msa) == 3).sum(1).A1) / (
        msa != 0
    ).sum(1).A1
    cluster_GC = mm.dot(GC)

    try:
        inserts = pd.read_csv(args.inserts, index_col=0, header=0, sep="\t")
    except pd.errors.EmptyDataError:
        inserts = None

    insert_stats = pd.DataFrame(
        [],
        columns=[
            "insert_overlap",
            "insert_pos_overlap",
            "insert_length",
            "insert_gap_length",
            "insert_pos_isotype",
        ],
    )
    insert_frequency = 0

    if inserts is not None:

        logger.info("checking inserts")

        def get_insert_isotype(x):
            left_isotype = ",".join(
                rec[3][3]
                for rec in intersect_intervals(
                    [
                        (
                            switch_chrom,
                            x["insert_pos_left"],
                            x["insert_pos_left"] + 1,
                        )
                    ],
                    switch_anno,
                    loj=True,
                )
            )
            right_isotype = ",".join(
                rec[3][3]
                for rec in intersect_intervals(
                    [
                        (
                            switch_chrom,
                            x["insert_pos_right"],
                            x["insert_pos_right"] + 1,
                        )
                    ],
                    switch_anno,
                    loj=True,
                )
            )
            return left_isotype + "_" + right_isotype

        def aggregate_inserts(x):
            import functools

            insert_coords = x[["insert_left", "insert_right"]].drop_duplicates()
            insert_union = interval_length(
                merge_intervals(
                    [
                        ("", y["insert_left"], y["insert_right"])
                        for _, y in insert_coords.iterrows()
                    ]
                )
            )
            insert_intersection = interval_length(
                functools.reduce(
                    intersect_intervals,
                    (
                        [("", y["insert_left"], y["insert_right"])]
                        for _, y in insert_coords.iterrows()
                    ),
                )
            )
            insert_overlap = insert_intersection / insert_union

            insert_pos = x[
                ["insert_pos_left", "insert_pos_right"]
            ].drop_duplicates()
            pos_union = (
                interval_length(
                    merge_intervals(
                        [
                            ("", y["insert_pos_left"], y["insert_pos_right"])
                            for _, y in insert_pos.iterrows()
                        ]
                    )
                )
                + 1
            )
            pos_intersection = (
                interval_length(
                    functools.reduce(
                        intersect_intervals,
                        (
                            [("", y["insert_pos_left"], y["insert_pos_right"])]
                            for _, y in insert_pos.iterrows()
                        ),
                    )
                )
                + 1
            )
            pos_overlap = pos_intersection / pos_union

            insert_len = (x["insert_right"] - x["insert_left"]).mean()
            gap_len = (x["insert_pos_right"] - x["insert_pos_left"]).mean()
            main_isotype = x["insert_isotype"].mode()[0]
            return pd.Series(
                {
                    "insert_overlap": insert_overlap,
                    "insert_pos_overlap": pos_overlap,
                    "insert_length": insert_len,
                    "insert_gap_length": gap_len,
                    "insert_pos_isotype": main_isotype,
                }
            )

        insert_coords = pd.DataFrame(
            inserts["insert_coords"].apply(decode_coords).tolist(),
            columns=["chrom", "insert_left", "insert_right", "."],
            index=inserts.index,
        ).drop(["chrom", "."], axis=1)
        insert_pos = pd.DataFrame(
            inserts["insert_pos"].apply(decode_coords).tolist(),
            columns=["chrom", "insert_pos_left", "insert_pos_right", "."],
            index=inserts.index,
        ).drop(["chrom", "."], axis=1)
        (
            insert_pos["insert_pos_left"],
            insert_pos["insert_pos_right"],
        ) = insert_pos.min(axis=1) - 1, insert_pos.max(axis=1)

        inserts = pd.concat([inserts, insert_coords, insert_pos], axis=1)
        inserts["insert_isotype"] = inserts.apply(get_insert_isotype, axis=1)

        tmp = pd.concat([clustering, inserts], axis=1)
        insert_frequency = (
            tmp.groupby("cluster")["insert"]
            .agg(lambda x: np.mean(~pd.isnull(x)))
            .loc[clusters]
        )
        if tmp.dropna().shape[0] > 0:
            insert_stats = (
                tmp.dropna().groupby("cluster").apply(aggregate_inserts)
            )

    df = pd.DataFrame(
        {
            "size": csize,
            "isotype": avg_isotype,
            "length": cluster_length,
            "GC": cluster_GC,
            "break_spread": break_spread,
            "mut_spread": mut_spread,
            "insert_frequency": insert_frequency,
        },
        index=clusters,
    )
    df = pd.concat([df, insert_stats], axis=1)

    switch_iis = get_switch_iis(anno_recs, cov_int, eff_start, 1)
    for isotype in np.unique(switch_iis):
        df["length_" + isotype] = avg_msa[:, switch_iis == isotype].sum(1)

    logger.info("saving results to " + args.output)
    df.to_csv(args.output)
