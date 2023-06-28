"""demultiplex reads using output of BLAST against primers and barcodes"""


def setup_argparse(parser):
    parser.add_argument("-i", "--input", dest="input", help="""input fastq(.gz)""")
    parser.add_argument("-b", "--blast", dest="blast", help="""BLAST output file (or stdin)""")
    parser.add_argument(
        "-c",
        "--cutoff",
        dest="cutoff",
        default=70,
        type=float,
        help="""%% identify cutoff [70]""",
    )
    parser.add_argument(
        "-d",
        "--dist",
        dest="dist",
        default=100,
        type=float,
        help="""max distance of match to end of sequence""",
    )
    parser.add_argument(
        "-m",
        "--min_reads",
        dest="min_reads",
        default=1000,
        type=float,
        help="""min # of reads required to create separate file [1000]""",
    )
    parser.add_argument("-o", "--outdir", dest="outdir", help="""output directory""")
    parser.add_argument("-f", "--figure", dest="figure", help="""summary figure""")
    parser.add_argument("-r", "--report", dest="report", help="""summary report""")
    parser.add_argument(
        "-s",
        "--sample-sheet",
        dest="sample_sheet",
        help="""sample sheet (barcode <tab> sample_name, no header)""",
    )
    parser.add_argument(
        "--collapse",
        dest="collapse",
        default=False,
        action="store_true",
        help="""count reads regardless of whether one or both barcodes are present (for identical barcodes)? [False]""",
    )
    parser.add_argument(
        "--split-reads",
        dest="split_reads",
        default=False,
        action="store_true",
        help="""split reads with adjacent internal primer binding sites""",
    )
    parser.add_argument(
        "--max-split-dist",
        dest="max_split_dist",
        default=200,
        help="""maximum distance between internal primer matches to split read""",
    )
    parser.add_argument("--nmax", dest="nmax", type=int, help="""maximum number of reads""")


def run(args):
    import os
    import sys
    import gzip
    import time
    import numpy as np
    from collections import defaultdict, Counter
    import pandas as pd
    import re
    from Bio import SeqIO
    from logzero import logger

    if args.sample_sheet is not None:
        logger.info("reading sample sheet from " + args.sample_sheet)
        sample_sheet = pd.read_csv(args.sample_sheet, sep="\t", index_col=0, header=None).squeeze()
        whitelist = sample_sheet.index

    if args.blast is not None:
        logger.info("reading and filtering BLAST output from " + args.blast)
        inf = open(args.blast)
    else:
        logger.info("reading and filtering BLAST output from stdin")
        inf = sys.stdin

    all_hits = defaultdict(list)
    n = 0
    m = 0
    for line in inf:
        n += 1
        if n % 1000 == 0:
            logger.info("{0}k hits read, {1} kept".format(n // 1000, m))
        ls = line.strip("\n").split("\t")
        read = ls[1]
        hit = dict(
            query=ls[0],
            query_len=int(ls[2]),
            seq_len=int(ls[3]),
            pid=float(ls[4]),
            matchlen=int(ls[5]),
            orientation="+" if int(ls[7]) > int(ls[6]) else "-",
            sstart=min(int(ls[6]), int(ls[7])),
            send=max(int(ls[6]), int(ls[7])),
            evalue=float(ls[8]),
            pideff=int(ls[5]) * float(ls[4]) / float(ls[2]),
            relpos=0.5 * (int(ls[6]) + int(ls[7])) / float(ls[3]),
        )
        if hit["pideff"] > args.cutoff:
            m += 1
            all_hits[read].append(hit)

    if len(all_hits) == 0:
        logger.warn("no hits in BLAST output file! exiting ...")
        sys.exit()
    else:
        logger.info("{0} hits read, {1} kept".format(n, m))

    logger.info("aggregating reads")
    nreads = Counter()
    for read, hits in all_hits.items():
        barcodes = []
        nbarcodes = 0
        for hit in hits:
            if not hit["query"].startswith("primer_"):
                nbarcodes += 1
                barcodes.append(hit["query"])
        if args.collapse:
            barcodes = set(barcodes)
        if args.sample_sheet is not None:
            barcodes = [bc for bc in barcodes if bc in whitelist]
        if len(barcodes) > 0:
            barcodes = "+".join(sorted(barcodes))
            nreads[barcodes] += 1

    read_stats = pd.Series(nreads).sort_values(ascending=False)
    combs = read_stats.index[read_stats > args.min_reads]
    if args.collapse:
        combs = pd.Index(list(set([bc for c in combs for bc in c.split("+")])))
    if args.sample_sheet is not None:
        combs = whitelist

    if args.outdir is not None:
        logger.info("using these barcode combinations for demultiplexing: " + ", ".join(combs))

        outf = {}
        for comb in combs:
            if args.sample_sheet is not None and comb in sample_sheet.index:
                outf[comb] = gzip.open(
                    os.path.join(args.outdir, sample_sheet[comb] + ".fastq.gz"),
                    "wt",
                )
            else:
                outf[comb] = gzip.open(os.path.join(args.outdir, comb + ".fastq.gz"), "wt")

        outf["undetermined"] = gzip.open(os.path.join(args.outdir, "undetermined.fastq.gz"), "wt")

        read_info = defaultdict(dict)

        logger.info("demultiplexing reads from " + args.input)
        if args.input.endswith(".gz"):
            inf = gzip.open(args.input, "rt")
        else:
            inf = open(args.input)

        nsplit = 0
        for n, rec in enumerate(SeqIO.parse(inf, "fastq")):
            if n % 1000 == 0:
                logger.info("{0:4d}k reads processed".format(n // 1000))
            if args.nmax and n >= args.nmax:
                break

            # get positions of all primer and barcode occurrences in the read
            read_hits = sorted(
                set(
                    (
                        re.sub("_[0-9]*$", "", h["query"]),
                        h["sstart"],
                        h["send"],
                        h["orientation"],
                    )
                    for h in all_hits[rec.id]
                ),
                key=lambda x: x[1],
            )

            # check if there are clusters of multiple primers or barcodes within args.max_split_dist
            clusters = []
            ii = []
            for i, rh in enumerate(read_hits):
                rh = read_hits[i]
                if not rh[0].startswith("primer"):
                    continue
                if ii:
                    rh2 = read_hits[ii[-1]]
                    dist = rh[2] - rh2[1]
                    overlap = min(rh[2], rh2[2]) >= max(rh[1], rh2[1])
                    if dist < args.max_split_dist and not overlap:
                        ii.append(i)
                    else:
                        ii = [i]
                else:
                    ii = [i]
                if len(ii) > 1:
                    clusters.append(ii)

            # don't split reads or no clusters
            if not args.split_reads or len(clusters) == 0:
                barcodes = [
                    rh
                    for rh in read_hits
                    if "primer" not in rh[0]
                    and (rh[1] <= args.dist or rh[2] >= len(rec) - args.dist)
                ]
                primers = [
                    rh
                    for rh in read_hits
                    if "primer" in rh[0] and (rh[1] <= args.dist or rh[2] >= len(rec) - args.dist)
                ]
                internal = [
                    rh for rh in read_hits if (rh[1] > args.dist and rh[2] < len(rec) - args.dist)
                ]
                if args.collapse:
                    comb = "+".join(set(bc[0] for bc in barcodes))
                else:
                    comb = "+".join(bc[0] for bc in barcodes)
                barcodes = ";".join("{0}@{1}-{2}:{3}".format(*bc) for bc in barcodes)
                primers = ";".join("{0}@{1}-{2}:{3}".format(*pr) for pr in primers)
                internal = ";".join("{0}@{1}-{2}:{3}".format(*h) for h in internal)
                if comb in combs:
                    comb = comb
                else:
                    comb = "undetermined"
                SeqIO.write(rec, outf[comb], "fastq")
                read_info[comb][rec.id] = dict(
                    map(lambda x: x.split("="), rec.description.split()[2:])
                )
                read_info[comb][rec.id].update(
                    dict(
                        meanQ=np.mean(rec.letter_annotations["phred_quality"]),
                        length=len(rec),
                        barcodes=barcodes,
                        primers=primers,
                        internal=internal,
                    )
                )
            else:
                nsplit += 1
                npart = 0
                sp = [0]
                for cl in clusters:
                    rhs = [read_hits[i] for i in range(cl[0], cl[-1] + 1)]
                    is_bc = [not lb.startswith("primer") for lb, _, _, _ in rhs]
                    if sum(is_bc) == 2:
                        # if there are two barcodes, split between them
                        si = min(i for i, t in enumerate(is_bc) if t)
                    else:
                        # otherwise split in the largest gap
                        gaps = [rhs[i][1] - rhs[i - 1][2] for i in range(1, len(rhs))]
                        si = np.argmax(gaps)
                    p = (rhs[si + 1][1] + rhs[si][2]) // 2
                    if p > sp[-1]:
                        sp.append(p)
                sp.append(len(rec))
                if len(sp) > len(set(sp)):
                    raise Exception("stop")
                for k in range(len(sp) - 1):
                    npart += 1
                    rsub = rec[sp[k] : sp[k + 1]]
                    rsub.id += "_part{0}".format(npart)
                    rsub.description += " part={0}".format(k + 1)
                    barcodes = [
                        (lb, st - sp[k], en - sp[k], o)
                        for lb, st, en, o in read_hits
                        if not lb.startswith("primer")
                        and (st >= sp[k] and en <= sp[k + 1])
                        and (st <= sp[k] + args.dist or en >= sp[k + 1] - args.dist)
                    ]
                    primers = [
                        (lb, st - sp[k], en - sp[k], o)
                        for lb, st, en, o in read_hits
                        if lb.startswith("primer")
                        and (st >= sp[k] and en <= sp[k + 1])
                        and (st <= sp[k] + args.dist or en >= sp[k + 1] - args.dist)
                    ]
                    internal = [
                        (lb, st - sp[k], en - sp[k], o)
                        for lb, st, en, o in read_hits
                        if (st >= sp[k] and en <= sp[k + 1])
                        and (st > sp[k] + args.dist and en < sp[k + 1] - args.dist)
                    ]
                    if args.collapse:
                        comb = "+".join(set(bc[0] for bc in barcodes))
                    else:
                        comb = "+".join(bc[0] for bc in barcodes)
                    barcodes = ";".join("{0}@{1}-{2}:{3}".format(*bc) for bc in barcodes)
                    primers = ";".join("{0}@{1}-{2}:{3}".format(*pr) for pr in primers)
                    internal = ";".join("{0}@{1}-{2}:{3}".format(*h) for h in internal)
                    if comb in combs:
                        comb = comb
                    else:
                        comb = "undetermined"
                    SeqIO.write(rsub, outf[comb], "fastq")
                    read_info[comb][rsub.id] = dict(
                        map(lambda x: x.split("="), rsub.description.split()[2:])
                    )
                    read_info[comb][rsub.id].update(
                        dict(
                            meanQ=np.mean(rsub.letter_annotations["phred_quality"]),
                            length=len(rsub),
                            barcodes=barcodes,
                            primers=primers,
                            internal=internal,
                        )
                    )

        logger.info("{0} reads split ... collecting info".format(nsplit))
        for comb in read_info.keys():
            read_info[comb] = pd.DataFrame.from_dict(read_info[comb], orient="index")
            if "start_time" in read_info[comb].columns:
                read_info[comb]["start_time"] = pd.to_datetime(read_info[comb]["start_time"]).apply(
                    lambda x: time.mktime(x.timetuple())
                )
                read_info[comb]["start_time"] = (
                    read_info[comb]["start_time"] - read_info[comb]["start_time"].min()
                )

            if args.sample_sheet is not None and comb in sample_sheet.index:
                read_info[comb].to_csv(os.path.join(args.outdir, sample_sheet[comb] + "_info.csv"))
            else:
                read_info[comb].to_csv(os.path.join(args.outdir, comb + "_info.csv"))

        logger.info("done")

    nreads = pd.Series(dict((k, v.shape[0]) for k, v in read_info.items())).sort_values(
        ascending=False
    )
    if args.report is not None:
        nreads.to_csv(args.report)

    if args.figure is not None:
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt

        logger.info("creating report figure in " + args.figure)

        fig = plt.figure(figsize=(8, 8))
        fig.clf()

        ax = fig.add_axes([0.25, 0.55, 0.6, 0.4])
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
