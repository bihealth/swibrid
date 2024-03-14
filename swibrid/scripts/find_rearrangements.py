"""\
find rearrangements (inversions/duplications) in MSA: 
regions with coverage larger or smaller than 1 and without gaps bigger than 25 nt
output files
- a .npz file with read indices, left/right positions and sizes for inversions and duplications per read
- a bed file with consensus rearrangement coordinates (by merging individual events), type and allele frequency
"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""required: file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("-o", "--out", dest="out", help="""required: output file (bed)""")
    parser.add_argument("-m", "--mat", dest="mat", help="""required: output file (matrix)""")
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
        "--max_rearrangement_gap",
        dest="max_rearrangement_gap",
        default=25,
        type=int,
        help="""max gap within rearrangements to ignore [25]""",
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
        merge_intervals,
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

    logger.info("loading msa from {0}".format(args.msa))
    msa = scipy.sparse.load_npz(args.msa).tocsr()
    nreads = msa.shape[0]

    logger.info("finding rearrangements (inversions / duplications)")

    # get indices of nonzero elements
    ii, jj = msa.nonzero()
    # get coverage of nonzero entries
    cov = msa.data // 10

    # find positions where coverage is larger or smaller than 1
    # and within that, get intervals of consecutive positions containing only gaps smaller than max_rearrangement_gap
    def get_intervals(pp):
        if len(pp) == 0:
            return [], [], []
        rr = ii[pp]
        pos = jj[pp] + 1
        inds = np.concatenate(
            [
                [0],
                np.where((np.diff(pos) > args.max_rearrangement_gap) | (np.diff(rr) != 0))[0]
                + 1,
                [len(pos)],
            ]
        )
        return (
            np.asarray(rr[inds[:-1]]),
            np.asarray(pos[inds[:-1]]),
            np.asarray(pos[(inds[1:] - 1)]),
        )

    inv_read, inv_left, inv_right = get_intervals(np.where(cov < 1)[0])
    dup_read, dup_left, dup_right = get_intervals(np.where(cov > 1)[0])
    inv_size = inv_right - inv_left if len(inv_read) > 0 else np.array([])
    dup_size = dup_right - dup_left if len(dup_read) > 0 else np.array([])

    assert np.all(inv_size >= 0) & np.all(dup_size >= 0), "negative inv/dup sizes!"

    if len(inv_read) >= 0:
        msa_inv = msa.tocsr()[inv_read]
        inv_cov = np.array(
            [
                (msa_inv[k, pl:pr].data // 10).mean()
                for k, (pl, pr) in enumerate(zip(inv_left, inv_right))
            ]
        )

        if np.any(inv_cov > -0.5):
            logger.warn(
                "{0} of {1} inversions have mean coverage > -0.5".format(
                    np.sum(inv_cov > -0.5), len(inv_read)
                )
            )

    if len(dup_read) >= 0:
        msa_dup = msa.tocsr()[dup_read]
        dup_cov = np.array(
            [
                (msa_dup[k, pl:pr].data // 10).mean()
                for k, (pl, pr) in enumerate(zip(dup_left, dup_right))
            ]
        )

        if np.any(dup_cov < 1.5):
            logger.warn(
                "{0} of {1} duplications have mean coverage < 1.5".format(
                    np.sum(dup_cov < 1.5), len(dup_read)
                )
            )

    logger.info("saving rearrangement positions to {0}".format(args.mat))
    np.savez(
        args.mat,
        inv_read=inv_read,
        inv_left=inv_left,
        inv_right=inv_right,
        inv_size=inv_size,
        dup_read=dup_read,
        dup_left=dup_left,
        dup_right=dup_right,
        dup_size=dup_size,
    )

    logger.info("getting consensus events and saving to {0}".format(args.out))

    merged_inversions = merge_intervals((('inversion',l,r) for l,r in zip(inv_left,inv_right)), add_count=True)
    merged_duplications = merge_intervals((('duplication',l,r) for l,r in zip(dup_left,dup_right)), add_count=True)
        
    consensus = pd.DataFrame(merged_inversions + merged_duplications,
                             columns=['name', 'left', 'right', 'nreads'])
    consensus['chrom'] = switch_chrom
    consensus['left'] = np.maximum(0, np.minimum(consensus['left'], Ltot - 1))
    consensus['right'] = np.maximum(0, np.minimum(consensus['right'], Ltot - 1))
    consensus_nreads = [(msa[:,l:r]!=0).sum(0).mean() for _,(l,r) in consensus[['left','right']].iterrows()]
    consensus['score'] = consensus['nreads'] / consensus_nreads
    consensus['start'] = cov_map[consensus['left']]
    consensus['end'] = cov_map[consensus['right']]
    consensus[['chrom','start','end','name','score']].sort_values(['chrom','start','end']).to_csv(args.out, header=False, index=False, sep='\t')
