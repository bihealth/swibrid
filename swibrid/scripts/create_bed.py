"""\
create bed files for all reads and summary table for insert-containing reads, including read sequence with insert part highlighted in uppercase
"""


def setup_argparse(parser):
    parser.add_argument(
        "--processed_alignments",
        dest="processed_alignments",
        help="""required: output of process_alignments""",
    )
    parser.add_argument(
        "--bed",
        dest="bed",
        help="""required: bed output with insert + switch coordinates for selected reads""",
    )
    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="""required: output table with insert information""",
    )
    parser.add_argument("--raw_reads", dest="raw_reads", help="""fasta file with minION reads""")
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000",
        help="""coordinates of switch region [chr14:106050000-106337000]""",
    )
    parser.add_argument(
        "--annotation",
        dest="annotation",
        nargs='?',
        help="""bed file with gene annotation""",
    )
    parser.add_argument(
        "--switch_annotation",
        dest="switch_annotation",
        help="""bed file with switch annotation""",
    )


def run(args):
    import re
    import pandas as pd
    from collections import defaultdict
    from Bio import SeqIO
    from logzero import logger
    import gzip
    from .utils import (
        parse_switch_coords,
        intersect_intervals,
        decode_coords,
        decode_insert,
    )

    (
        switch_chrom,
        switch_start,
        switch_end,
        switch_orientation,
    ) = parse_switch_coords(args.switch_coords)

    if args.annotation is not None:
        logger.info("adding annotation from " + args.annotation)
        annotation = defaultdict(list)
        for line in open(args.annotation):
            ls = line.strip("\n").split("\t")
            annotation[ls[0]].append((ls[0], int(ls[1]), int(ls[2])) + tuple(ls[3:]))

    logger.info("reading processed alignments from " + args.processed_alignments)
    process = pd.read_csv(args.processed_alignments, sep="\t", header=0, index_col=0)

    logger.info("writing read alignment coordinates to " + args.bed)
    logger.info("writing summary table to " + args.outfile)

    bed = open(args.bed, "w")
    inserts = defaultdict(list)
    for read, row in process.iterrows():
        mappings = [decode_coords(m) for m in row["switch_mappings"].split(";")]

        tstart = mappings[0][1]
        tend = mappings[-1][2]

        bed_entries = [
            switch_chrom,
            tstart,
            tend,
            read,
            0,
            ".",
            tstart,
            tstart,
            "0,0,0",
            len(mappings),
            ",".join(map(str, (sm[2] - sm[1] for sm in mappings))) + ",",
            ",".join(map(str, (sm[1] - tstart for sm in mappings))) + ",",
        ]

        bed.write("\t".join(map(str, bed_entries)) + "\n")

        if pd.isnull(row["inserts"]):
            continue

        for insert in (decode_insert(i) for i in row["inserts"].split(";")):
            insert_chrom = insert.group("insert_chrom")
            insert_start = int(insert.group("insert_start"))
            insert_end = int(insert.group("insert_end"))
            insert_len = insert_end - insert_start

            insert_anno = set()
            if args.annotation is not None:
                for rec in intersect_intervals(
                    [(insert_chrom, insert_start, insert_end)],
                    annotation[insert_chrom],
                    loj=True,
                ):
                    insert_anno.add(re.sub(".exon[0-9]*", "", rec[3][3]))
            if len(insert_anno) == 0 or insert_anno == set("."):
                insert_anno = "{0}:{1}-{2}".format(insert_chrom, insert_start, insert_end)
            else:
                insert_anno = "|".join(insert_anno)

            bed_str = "{0}\t{1}\t{2}\t{3}\t0\t.\t{1}\t{2}\t0,0,0\t1\t{4},\t0,\n"
            bed.write(bed_str.format(insert_chrom, insert_start, insert_end, read, insert_len))

            inserts[read].append((insert, insert_anno, tstart, tend, row["isotype"]))

    bed.close()

    insert_table = {}
    logger.info("parsing raw reads from " + args.raw_reads)
    for rec in SeqIO.parse(
        gzip.open(args.raw_reads, "rt") if args.raw_reads.endswith(".gz") else open(args.raw_reads),
        "fastq",
    ):
        if rec.id not in inserts:
            continue

        for insert, insert_anno, tstart, tend, isotype in inserts[rec.id]:
            insert_chrom = insert.group("insert_chrom")
            insert_start = int(insert.group("insert_start"))
            insert_end = int(insert.group("insert_end"))
            insert_len = insert_end - insert_start
            cstart = int(insert.group("switch_left")) - 1
            cend = int(insert.group("switch_right")) + 1
            if cend < cstart:
                cstart, cend = cend, cstart
            istart = int(insert.group("istart"))
            iend = int(insert.group("iend"))

            seq = rec.seq[:istart].lower() + rec.seq[istart:iend].upper() + rec.seq[iend:].lower()

            insert_table[rec.id] = dict(
                read=rec.id,
                isotype=isotype,
                switch_coords="{0}:{1}-{2}".format(switch_chrom, tstart + 1, tend),
                insert_pos="{0}:{1}-{2}".format(switch_chrom, cstart + 1, cend),
                insert_coords="{0}:{1}-{2}".format(insert_chrom, insert_start + 1, insert_end),
                insert_anno=insert_anno,
                sequence=seq,
            )

    pd.DataFrame.from_dict(insert_table, orient="index").to_csv(
        args.outfile, header=True, index=False, sep="\t"
    )
