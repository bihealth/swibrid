""" analyze sequence homology within switch region bins by jaccard similarity of k-mers  """


def setup_argparse(parser):

    parser.add_argument(
        "-g", "--genome", dest="genome", help="""genome fasta file"""
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
        "-b",
        "--binsize",
        dest="binsize",
        type=int,
        default=100,
        help="""binsize [100]""",
    )
    parser.add_argument(
        "-o", "--output", dest="output", help="""output file (npz)"""
    )
    parser.add_argument(
        "--range",
        dest="range",
        default="4-7",
        help="""range of kmer sizes, e.g., 3,5-7 [4-7]""",
    )
    parser.add_argument("--figure", dest="figure", help="""output figure""")


def get_rc(seq):
    rc = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(rc[c] for c in seq[::-1])


def jaccard_distance(s1, s2, n=5, rc=False):
    k1 = set(s1[i : i + n] for i in range(len(s1) - n))
    if rc:
        s2r = get_rc(s2)
        k2 = set(s2r[i : i + n] for i in range(len(s2r) - n))
    else:
        k2 = set(s2[i : i + n] for i in range(len(s2) - n))
    return len(k1 & k2) / len(k1 | k2)


def run(args):

    import pysam
    import numpy as np
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
        get_switch_coverage,
        shift_coord,
        parse_range,
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

    switch_seqs = [
        genome.fetch(switch_chrom, ci[0], ci[1]).upper() for ci in cov_int
    ]
    stot = "".join(switch_seqs)

    jd = dict()
    binsize = args.binsize
    for n in parse_range(args.range):
        jd["fw_{0}".format(n)] = np.array(
            [
                [
                    jaccard_distance(
                        stot[
                            max(binsize * k - n, 0) : min(
                                binsize * (k + 1) + n, len(stot)
                            )
                        ],
                        stot[
                            max(binsize * j - n, 0) : min(
                                binsize * (j + 1) + n, len(stot)
                            )
                        ],
                        n=n,
                    )
                    for j in range(len(stot) // binsize)
                ]
                for k in range(len(stot) // binsize)
            ]
        )
        jd["rv_{0}".format(n)] = np.array(
            [
                [
                    jaccard_distance(
                        stot[
                            max(binsize * k - n, 0) : min(
                                binsize * (k + 1) + n, len(stot)
                            )
                        ],
                        stot[
                            max(binsize * j - n, 0) : min(
                                binsize * (j + 1) + n, len(stot)
                            )
                        ],
                        n=n,
                        rc=True,
                    )
                    for j in range(len(stot) // binsize)
                ]
                for k in range(len(stot) // binsize)
            ]
        )

    if args.output:
        np.savez(args.output, **jd)

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

        fig, axs = plt.subplots(1, len(jd) // 2, figsize=(2 * len(jd), 4))

        for k, n in enumerate(parse_range(args.range)):
            ax = axs[k]
            jj = np.zeros_like(jd["fw_{0}".format(n)])
            iu = np.triu_indices(len(stot) // binsize)
            il = np.tril_indices(len(stot) // binsize)
            jj[il] = jd["fw_{0}".format(n)][il]
            jj[iu] = jd["rv_{0}".format(n)][iu]
            ax.imshow(
                jj, cmap=plt.cm.Reds, origin="lower", interpolation="none"
            )
            ax.set_xlim([len(stot) // binsize, 0])
            ax.set_ylim([len(stot) // binsize, 0])
            ax.plot(
                [0, 1], [0, 1], transform=ax.transAxes, color="gray", lw=0.5
            )
            ax.set_yticks(np.array(major_ticks) // binsize)
            ax.set_xticks(np.array(major_ticks) // binsize)
            ax.set_yticks(np.array(minor_ticks) // binsize, minor=True)
            ax.set_xticks(np.array(minor_ticks) // binsize, minor=True)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticklabels(minor_labels, rotation=90, minor=True)
            ax.set_yticklabels(minor_labels if k == 0 else [], minor=True)
            ax.tick_params(which="minor", length=0)
            ax.grid(alpha=0.5, color="lightgrey", which="major")
            ax.set_title("k={0}".format(n), size="medium")

        fig.savefig(args.figure)
