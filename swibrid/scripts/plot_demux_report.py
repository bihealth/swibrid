"""plot demultiplexing report"""


def setup_argparse(parser):
    parser.add_argument("-o", "--outdir", dest="outdir", help="""output directory""")
    parser.add_argument("-f", "--figure", dest="figure", help="""summary figure""")
    parser.add_argument(
        "-s",
        "--sample-sheet",
        dest="sample_sheet",
        help="""sample sheet (barcode <tab> sample_name, no header)""",
    )


def run(args):
    import os
    import numpy as np
    from collections import Counter
    import pandas as pd
    import re
    from logzero import logger
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    logger.info("reading sample sheet from " + args.sample_sheet)
    sample_sheet = pd.read_csv(args.sample_sheet, sep="\t", index_col=0, header=None).squeeze()
    whitelist = sample_sheet.index

    logger.info("reading read info from " + args.outdir)
    read_info = {}
    for comb in list(whitelist):
        fname = os.path.join(args.outdir, sample_sheet[comb] + "_info.csv")
        if os.path.isfile(fname):
            read_info[comb] = pd.read_csv(fname)

    fname = os.path.join(args.outdir, "undetermined_info.csv")
    if os.path.isfile(fname):
        read_info["undetermined"] = pd.read_csv(fname)

    nreads = pd.Series(dict((k, v.shape[0]) for k, v in read_info.items())).sort_values(
        ascending=False
    )

    logger.info("creating report figure in " + args.figure)

    fig = plt.figure(figsize=(8, 8))
    fig.clf()

    ax = fig.add_axes([0.25, 0.55, 0.5, 0.4])
    labels = []
    for n, comb in enumerate(nreads.index):
        ax.barh([n], [nreads[comb]], lw=0.5, edgecolor="k")
        ax.text(
            nreads[comb],
            n,
            str(nreads[comb]),
            size="x-small",
            ha="left",
            va="center",
        )
        if args.sample_sheet is not None and comb in sample_sheet.index:
            labels.append(sample_sheet[comb])
        else:
            labels.append(comb)
    ax.set_xticks([])
    ax.set_yticks(range(len(nreads.index)))
    ax.set_yticklabels(labels)
    ax.set_title("{0} reads".format(nreads.sum()))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    ax = fig.add_axes([0.8, 0.7, 0.15, 0.2])

    def classify_barcodes(bc, whitelist):
        if pd.isnull(bc):
            return "none"
        bcs = sorted(set(map(lambda x: x.split("@")[0], bc.split(";"))))
        if len(bcs) > 1:
            return "multiple"
        else:
            return ";".join(bcs)

    ignored = pd.Series(
        Counter(classify_barcodes(bc, whitelist) for bc in read_info["undetermined"]["ignored"])
    ).sort_values(ascending=False)
    labels = [i if ignored[i] >= 0.01 * ignored.sum() else "" for i in ignored.index]
    colors = [
        plt.cm.Set2(k) if ignored[i] >= 0.01 * ignored.sum() else (0.5, 0.5, 0.5, 1.0)
        for k, i in enumerate(ignored.index)
    ]
    ax.pie(ignored.values, labels=labels, colors=colors, textprops={"size": "x-small"})
    ax.set_title("undetermined reads", size="small")

    ax = fig.add_axes([0.13, 0.1, 0.35, 0.4])
    for comb in nreads.index:
        bc_pos = [
            np.mean(list(map(int, re.split("[-:]", p.split("@")[1])[:2]))) / r["length"]
            for _, r in read_info[comb][["length", "barcodes"]].dropna().iterrows()
            if len(r["barcodes"]) > 0
            for p in r["barcodes"].split(";")
            if comb in p
        ]
        ax.hist(bc_pos, bins=np.linspace(0, 1, 50), histtype="step")

    ax.set_xticks([0, 0.5, 1])
    ax.set_ylabel("# of reads")
    ax.set_title("barcode position")
    ax.set_xlabel("rel pos in read")

    ax = fig.add_axes([0.58, 0.1, 0.35, 0.4])
    for comb in nreads.index:
        ax.hist(
            np.log10(read_info[comb]["length"]),
            bins=np.linspace(2, 4, 50),
            histtype="step",
        )

    ax.set_xticks([2, 3, 4])
    ax.set_xticklabels([100, 1000, 10000])
    ax.set_title("length distribution")
    ax.set_xlabel("read length")

    fig.savefig(args.figure)
