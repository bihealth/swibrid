"""get gaps in (pseudo)MSA"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("-o", "--out", help="""output file (.npz)""")


def run(args):
    import os
    import sys
    import numpy as np
    import scipy.sparse
    from logzero import logger
    from .utils import get_gap_positions

    assert os.path.isfile(args.msa), "no msa at {0}; run construct_msa.py first!".format(args.msa)

    logger.info("loading msa from {0}".format(args.msa))
    msa = scipy.sparse.load_npz(args.msa)

    read_idx, pos_left, pos_right, gap_size, inversion = get_gap_positions(msa)

    logger.info("saving results to {0}".format(args.out))
    np.savez(
        args.out,
        read_idx=read_idx,
        pos_left=pos_left,
        pos_right=pos_right,
        gap_size=gap_size,
        inversion=inversion,
    )
