"""\
find gaps in (pseudo)MSA
takes the MSA matrix from construct_MSA and finds all gap positions
output is a .npz file containing several integer arrays
- read_idx: index of read
- gap_left: left end (within the MSA; 0-based)
- gap_right: right end 
- gap_size: size of gap
"""


def setup_argparse(parser):
    parser.add_argument(
        "--msa",
        dest="msa",
        help="""required: .npz file with  (pseudo) multiple alignment of read sequences for clustering""",
    )
    parser.add_argument("-o", "--out", help="""required: output file (.npz)""")


def run(args):
    import os
    import numpy as np
    import scipy.sparse
    from logzero import logger
    from .utils import get_gap_positions

    assert os.path.isfile(args.msa), "no msa at {0}; run construct_msa.py first!".format(args.msa)

    logger.info("loading msa from {0}".format(args.msa))
    msa = scipy.sparse.load_npz(args.msa)

    read_idx, gap_left, gap_right, gap_size = get_gap_positions(msa)

    logger.info("saving results to {0}".format(args.out))
    np.savez(
        args.out,
        read_idx=read_idx,
        gap_left=gap_left,
        gap_right=gap_right,
        gap_size=gap_size,
    )
