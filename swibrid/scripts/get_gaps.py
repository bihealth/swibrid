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
    parser.add_argument(
        "--paired_end_mode",
        dest="paired_end_mode",
        action="store_true",
        default=False,
        help="""use paired-end mode (--raw_reads needs to be a comma-separated list of mates)""",
    )
    parser.add_argument("--msa_csv", dest="msa_csv", help="""msa read info (for paire_end_mode)""")


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

    if args.paired_end_mode:
        import pandas as pd

        logger.info("paired_end_mode: removing mate breaks using {0}".format(args.msa_csv))
        msa_csv = pd.read_csv(args.msa_csv, index_col=None, header=0)
        mate_breaks = dict(
            (i, (int(mb.split(";")[0]), int(mb.split(";")[1])))
            for i, mb in msa_csv["mate_breaks"].dropna().items()
            if not "nan" in mb
        )
        use = [
            i not in mate_breaks or mate_breaks[i] != (l, r)
            for i, l, r in zip(read_idx, gap_left, gap_right)
        ]
        read_idx = read_idx[use]
        gap_left = gap_left[use]
        gap_right = gap_right[use]
        gap_size = gap_size[use]

    logger.info("saving results to {0}".format(args.out))
    np.savez(
        args.out,
        read_idx=read_idx,
        gap_left=gap_left,
        gap_right=gap_right,
        gap_size=gap_size,
    )
