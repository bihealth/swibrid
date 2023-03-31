""" analyze sequence content within switch region bins by a number of motifs """


def setup_argparse(parser):

    parser.add_argument("-g", "--genome", dest="genome", help="""genome fasta file""")
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
        "-b",
        "--binsize",
        dest="binsize",
        type=int,
        default=100,
        help="""binsize [100]""",
    )
    parser.add_argument("-o", "--output", dest="output", help="""output file (npz)""")
    parser.add_argument(
        "--motifs",
        dest="motifs",
        default="TG,CA,TCCA,CAGCC,CAGCT",
        help="""comma-separated list of sequence motifs [TG,CA,TCCA,CAGCC,CAGCT]""",
    )
    parser.add_argument("--figure", dest="figure", help="""output figure""")


def run(args):

    import pysam
    import numpy as np
    import re
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
        get_switch_coverage,
        shift_coord,
    )

    genome = pysam.FastaFile(args.genome)

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

    switch_seqs = [genome.fetch(switch_chrom, ci[0], ci[1]).upper() for ci in cov_int]
    stot = "".join(switch_seqs)

    binsize = args.binsize

    motif_counts = {}
    motif_counts["G_C"] = (
        np.array(
            [
                sum(1 for _ in re.finditer("G|C", stot[k * binsize : (k + 1) * binsize]))
                for k in range(len(stot) // binsize)
            ]
        )
        / binsize
    )
    motif_counts["A_T"] = (
        np.array(
            [
                sum(1 for _ in re.finditer("A|T", stot[k * binsize : (k + 1) * binsize]))
                for k in range(len(stot) // binsize)
            ]
        )
        / binsize
    )
    for motif in args.motifs.split(","):
        motif_counts[motif] = np.array(
            [
                sum(1 for _ in re.finditer(motif, stot[k * binsize : (k + 1) * binsize]))
                for k in range(len(stot) // binsize)
            ]
        ) / (binsize / len(motif))

    if args.output:
        np.savez(args.output, **motif_counts)

    if args.figure:

        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt

        major_ticks = []
        minor_ticks = []
        minor_labels = []
        for rec in anno_recs:
            start = shift_coord(int(rec[3][1]), cov_int) - eff_start
            end = shift_coord(int(rec[3][2]), cov_int) - eff_start
            major_ticks += [start, end]
            minor_ticks.append((start + end) / 2)
            minor_labels.append(rec[3][3])
        minor_ticks = np.array(minor_ticks)
        major_ticks = np.array(np.unique(major_ticks))

        fig, axs = plt.subplots(len(motif_counts), 1, figsize=(4, 2 * len(motif_counts)))
        fig.subplots_adjust(top=0.97, bottom=0.05)

        for k, motif in enumerate(motif_counts.keys()):
            ax = axs[k]
            ax.plot(np.arange(len(stot) // binsize), motif_counts[motif], "-")
            ax.set_xlim([len(stot) // binsize, 0])
            ax.set_xticks(np.array(major_ticks) // binsize)
            ax.set_xticks(np.array(minor_ticks) // binsize, minor=True)
            ax.set_xticklabels([])
            if k == len(motif_counts.keys()) - 1:
                ax.set_xticklabels(minor_labels, rotation=90, minor=True)
            ax.tick_params(which="minor", length=0)
            ax.grid(alpha=0.5, color="lightgrey", which="major", axis="x")
            ax.set_title(motif, size="medium")

        fig.savefig(args.figure)
