"""\
get LAST-like parameters from a SAM or MAF file
this will go through the alignments and tally mutation and indel frequencies
"""


def setup_argparse(parser):
    parser.add_argument(
        "-i", "--inf", dest="inf", help="""input .bam or .maf.gz file, or .par from LAST"""
    )
    parser.add_argument("-r", "--ref", dest="ref", help="""reference genome""")
    parser.add_argument("-o", "--out", dest="out", help="""output file (.npz)""")


def run(args):
    from collections import defaultdict
    import numpy as np
    import pandas as pd
    import pysam
    from Bio import AlignIO
    import gzip
    from logzero import logger
    from .utils import RC

    stats = defaultdict(int)
    ins_sizes = []
    del_sizes = []

    if args.inf.endswith(".par"):
        from .utils import parse_LAST_pars

        logger.info("reading parameter file " + args.inf)
        pars = parse_LAST_pars(args.inf)

        res = {
            "p_c": pars["p_c"],
            "p_open_del": pars["p_open_del"],
            "p_open_ins": pars["p_open_ins"],
            "p_extend_del": pars["p_extend_del"],
            "p_extend_ins": pars["p_extend_ins"],
        }

        logger.info("saving output to " + args.out)
        np.savez(args.out, **res)

        return

    elif args.inf.endswith(".bam"):
        logger.info("loading reference from " + args.ref)
        genome = pysam.FastaFile(args.ref)

        logger.info("reading bam file " + args.inf)

        for rec in pysam.Samfile(args.inf):
            if rec.is_unmapped or rec.seq is None or rec.is_secondary:
                continue

            pairs = rec.get_aligned_pairs()
            ref_seq = genome.fetch(
                rec.reference_name, rec.reference_start, rec.reference_end + 1
            ).upper()

            ref_gap_size = 0
            read_gap_size = 0
            seq_started = False

            for read_pos, ref_pos in pairs:
                if read_pos is None and seq_started:
                    read_gap_size += 1

                elif ref_pos is None and seq_started:
                    ref_gap_size += 1

                elif read_pos is not None and ref_pos is not None:
                    seq_started = True

                    if read_gap_size > 0:
                        del_sizes.append(read_gap_size)
                        read_gap_size = 0

                    if ref_gap_size > 0:
                        ins_sizes.append(ref_gap_size)
                        ref_gap_size = 0

                    ref_nuc = ref_seq[ref_pos - rec.reference_start]
                    read_nuc = str(rec.seq[read_pos]).upper()

                    if rec.is_reverse:
                        stats[(RC[read_nuc], RC[ref_nuc])] += 1
                    else:
                        stats[(read_nuc, ref_nuc)] += 1

    elif args.inf.endswith(".maf") or args.inf.endswith(".maf.gz"):
        logger.info("loading reference from " + args.ref)
        genome = pysam.FastaFile(args.ref)

        logger.info("reading maf file " + args.inf)

        for ref, read in AlignIO.parse(
            gzip.open(args.inf, "rt") if args.inf.endswith(".gz") else open(args.inf),
            "maf",
            seq_count=2,
        ):
            ref_gap_size = 0
            read_gap_size = 0
            seq_started = False

            for read_nuc, ref_nuc in zip(str(read.seq), str(ref.seq)):
                if read_nuc == "-" and seq_started:
                    read_gap_size += 1

                elif ref_nuc == "-" and seq_started:
                    ref_gap_size += 1

                elif read_nuc != "-" and ref_nuc != "-":
                    seq_started = True

                    if read_gap_size > 0:
                        del_sizes.append(read_gap_size)
                        read_gap_size = 0

                    if ref_gap_size > 0:
                        ins_sizes.append(ref_gap_size)
                        ref_gap_size = 0

                    stats[(read_nuc.upper(), ref_nuc.upper())] += 1

    stats = pd.Series(stats).unstack(level=1).values

    # remove trailing gaps
    if read_gap_size > 0:
        del_sizes.pop()
    if ref_gap_size > 0:
        ins_sizes.pop()

    ntot = stats.sum().sum()

    res = {
        "p_c": stats,
        "p_open_del": len(del_sizes) / ntot,
        "p_open_ins": len(ins_sizes) / ntot,
        "p_extend_del": 1.0 - 1.0 / np.mean(del_sizes),
        "p_extend_ins": 1.0 - 1.0 / np.mean(ins_sizes),
    }

    logger.info("saving output to " + args.out)
    np.savez(args.out, **res)
