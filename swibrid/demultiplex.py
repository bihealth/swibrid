"""\
demultiplex reads using output of BLAST against primers and barcodes

BLAST output is created with ``blastn -db {barcodes_primers} -query {query} -task blastn-short -max_target_seqs 50 -outfmt "6 saccver qaccver slen qlen pident length qstart qend evalue" -gapopen 5 -gapextend 2 -reward 3 -penalty -4 -evalue 1 -num_threads 1 -perc_identity 50`` where {barcodes_primers} is a database containing barcode and primer sequences and {query} is the input fasta file

optional (but recommended) input is a sample sheet (tab delimited, no header or row names, 1st column barcode, 2nd column sample name)
"""


def setup_argparse(parser):
    parser.add_argument(
        "-i", "--input", dest="input", help="""required: input fastq(.gz)""", required=True
    )
    parser.add_argument(
        "-b",
        "--blast",
        dest="blast",
        required=True,
        help="""required: BLAST output file (or stdin)""",
    )
    parser.add_argument(
        "-o", "--outdir", dest="outdir", help="""required: output directory""", required=True
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        dest="cutoff",
        default=70,
        type=float,
        help="""%% identify cutoff [%(default).1f]""",
    )
    parser.add_argument(
        "-d",
        "--dist",
        dest="dist",
        default=100,
        type=float,
        help="""max distance of match to end of sequence [%(default)d]""",
    )
    parser.add_argument(
        "-m",
        "--min_reads",
        dest="min_reads",
        default=1000,
        type=float,
        help="""min # of reads required to create separate file if no sample sheet is given [%(default).0f]""",
    )
    parser.add_argument("-f", "--figure", dest="figure", help="""summary figure""")
    parser.add_argument("-r", "--report", dest="report", help="""summary report""")
    parser.add_argument(
        "-s",
        "--sample_sheet",
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
        "--split_reads",
        dest="split_reads",
        default=False,
        action="store_true",
        help="""split reads with adjacent internal primer binding sites""",
    )
    parser.add_argument(
        "--max_split_dist",
        dest="max_split_dist",
        default=200,
        help="""maximum distance between internal primer matches to split read [%(default)d]""",
    )
    parser.add_argument("--nmax", dest="nmax", type=int, help="""maximum number of reads""")
    parser.add_argument(
        "--add_umi",
        dest="add_umi",
        action="store_true",
        default=False,
        help="""add UMI to read name""",
    )
    parser.add_argument(
        "--umi_regex",
        dest="umi_regex",
        default=r"[ACGTN]{12}",
        help="""regex to find UMI in read description [%(default)s]""",
    )
    parser.add_argument(
        "--write_info_chunks",
        dest="write_info_chunks",
        action="store_true",
        default=False,
        help="""write info files in chunks""",
    )


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
        sample_sheet = pd.read_csv(args.sample_sheet, sep="\t", index_col=0, header=None).squeeze(
            axis=1
        )
        whitelist = sample_sheet.index
        assert len(whitelist) == len(
            whitelist.unique()
        ), "sample sheet contains repeated barcodes. exiting ..."

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
        if n % 1000000 == 0:
            logger.info("{0}M hits read, {1} kept".format(n // 1000000, m))
        ls = line.strip("\n").split("\t")
        read = ls[1]
        hit = dict(
            query=ls[0],
            # seq_len=int(ls[2]),
            # query_len=int(ls[3]),
            # pid=float(ls[4]),
            # matchlen=int(ls[5]),
            orientation="+" if int(ls[7]) > int(ls[6]) else "-",
            sstart=min(int(ls[6]), int(ls[7])),
            send=max(int(ls[6]), int(ls[7])),
            # evalue=float(ls[8]),
            pideff=int(ls[5]) * float(ls[4]) / float(ls[2]),
        )
        if hit["pideff"] > args.cutoff:
            m += 1
            all_hits[read].append(hit)

        if args.nmax and len(all_hits) > args.nmax:
            break

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
    else:
        whitelist = combs

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
        read_info_created = set()
        nreads = defaultdict(int)

        logger.info("demultiplexing reads from " + args.input)
        if args.input.endswith(".gz"):
            inf = gzip.open(args.input, "rt")
        else:
            inf = open(args.input)

        nsplit = 0
        for n, rec in enumerate(SeqIO.parse(inf, "fastq")):
            if args.nmax and n >= args.nmax:
                break

            # get positions of all primer and barcode occurrences in the read
            read_hits = sorted(
                set(
                    (
                        # collapse different primer variants
                        re.sub(r"^(primer.*)_[0-9]*$", r"\1", h["query"]),
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
                    and rh[0] in whitelist
                ]
                primers = [
                    rh
                    for rh in read_hits
                    if "primer" in rh[0] and (rh[1] <= args.dist or rh[2] >= len(rec) - args.dist)
                ]
                internal = [
                    rh for rh in read_hits if (rh[1] > args.dist and rh[2] < len(rec) - args.dist)
                ]
                ignored = [
                    rh
                    for rh in read_hits
                    if "primer" not in rh[0]
                    and (rh[1] <= args.dist or rh[2] >= len(rec) - args.dist)
                    and rh[0] not in whitelist
                ]
                if args.collapse:
                    comb = "+".join(set(bc[0] for bc in barcodes))
                else:
                    comb = "+".join(bc[0] for bc in barcodes)
                barcodes = ";".join("{0}@{1}-{2}:{3}".format(*bc) for bc in barcodes)
                primers = ";".join("{0}@{1}-{2}:{3}".format(*pr) for pr in primers)
                internal = ";".join("{0}@{1}-{2}:{3}".format(*h) for h in internal)
                ignored = ";".join("{0}@{1}-{2}:{3}".format(*h) for h in ignored)
                if comb in combs:
                    comb = comb
                else:
                    comb = "undetermined"

                if args.add_umi:
                    umi_match = re.search(args.umi_regex, rec.description)
                    if umi_match:
                        rec.id += "_" + umi_match.group()

                SeqIO.write(rec, outf[comb], "fastq")
                nreads[comb] += 1
                read_info[comb][rec.id] = dict(
                    map(lambda x: x.split("="), rec.description.split()[2:])
                )
                read_info[comb][rec.id]["part"] = pd.NA
                read_info[comb][rec.id].update(
                    dict(
                        meanQ=np.mean(rec.letter_annotations["phred_quality"]),
                        length=len(rec),
                        barcodes=barcodes,
                        primers=primers,
                        internal=internal,
                        ignored=ignored,
                    )
                )
            else:
                nsplit += 1
                npart = 0
                sp = [0]
                for cl in clusters:
                    rhs = [read_hits[i] for i in range(cl[0], cl[-1] + 1)]
                    is_bc = [not lb.startswith("primer") for lb, _, _, _ in rhs if lb in whitelist]
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
                        and lb in whitelist
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
                    ignored = [
                        (lb, st - sp[k], en - sp[k], o)
                        for lb, st, en, o in read_hits
                        if not lb.startswith("primer")
                        and (st >= sp[k] and en <= sp[k + 1])
                        and (st <= sp[k] + args.dist or en >= sp[k + 1] - args.dist)
                        and lb not in whitelist
                    ]
                    if args.collapse:
                        comb = "+".join(set(bc[0] for bc in barcodes))
                    else:
                        comb = "+".join(bc[0] for bc in barcodes)
                    barcodes = ";".join("{0}@{1}-{2}:{3}".format(*bc) for bc in barcodes)
                    primers = ";".join("{0}@{1}-{2}:{3}".format(*pr) for pr in primers)
                    internal = ";".join("{0}@{1}-{2}:{3}".format(*h) for h in internal)
                    ignored = ";".join("{0}@{1}-{2}:{3}".format(*h) for h in ignored)
                    if comb in combs:
                        comb = comb
                    else:
                        comb = "undetermined"
                    SeqIO.write(rsub, outf[comb], "fastq")
                    nreads[comb] += 1
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
                            ignored=ignored,
                        )
                    )

            if n % 1000 == 0 and n > 0:
                logger.info("{0:4d}k reads processed".format(n // 1000))

                if args.write_info_chunks:
                    for comb in read_info.keys():
                        read_info[comb] = pd.DataFrame.from_dict(read_info[comb], orient="index")
                        if "start_time" in read_info[comb].columns:
                            read_info[comb]["start_time"] = pd.to_datetime(
                                read_info[comb]["start_time"]
                            ).apply(lambda x: time.mktime(x.timetuple()))
                            read_info[comb]["start_time"] = (
                                read_info[comb]["start_time"] - read_info[comb]["start_time"].min()
                            )

                        if args.sample_sheet is not None and comb in sample_sheet.index:
                            info_file = os.path.join(args.outdir, sample_sheet[comb] + "_info.csv")
                        else:
                            info_file = os.path.join(args.outdir, comb + "_info.csv")

                        if read_info[comb].shape[0] == 0:
                            continue

                        if n == 1000 or comb not in read_info_created:
                            read_info[comb].to_csv(info_file, mode="w", header=True)
                            read_info_created.add(comb)
                        else:
                            read_info[comb].to_csv(info_file, mode="a", header=False)

                    read_info = defaultdict(dict)

        logger.info("{0:4d} reads processed".format(n))
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
                info_file = os.path.join(args.outdir, sample_sheet[comb] + "_info.csv")
            else:
                info_file = os.path.join(args.outdir, comb + "_info.csv")

            if read_info[comb].shape[0] == 0:
                continue

            if comb not in read_info_created:
                read_info[comb].to_csv(info_file, mode="w", header=True)
                read_info_created.add(comb)
            else:
                read_info[comb].to_csv(info_file, mode="a", header=False)

        for f in outf.values():
            f.close()

    if args.report is not None:
        pd.Series(nreads).sort_values(ascending=False).to_csv(args.report)

    if args.figure is not None:
        from .plot_demux_report import run as plot_demux_report

        plot_demux_report(args)
