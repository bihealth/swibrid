"""create bed files and summary table for insert-containing reads"""


def setup_argparse(parser):

    parser.add_argument("--raw_reads", dest="raw_reads", help="""fasta file with minION reads""")
    parser.add_argument(
        "--processed_reads",
        dest="processed_reads",
        help="""output of process_alignments""",
    )
    parser.add_argument(
        "--bed",
        dest="bed",
        help="""bed output with insert + switch coordinates for selected reads""",
    )
    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="""output table with all information""",
    )
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000",
        help="""coordinates of switch region [chr14:106050000-106337000]""",
    )
    parser.add_argument(
        "--annotation",
        dest="annotation",
        help="""bed file with gene annotation""",
    )
    parser.add_argument(
        "--switch_annotation",
        dest="switch_annotation",
        help="""bed file with switch annotation""",
    )
    parser.add_argument(
        "--interrupt_for_read",
        dest="interrupt_for_read",
        help="""interrupt for specified read(s)""",
    )


def run(args):

    import re
    import pandas as pd
    from collections import defaultdict
    from Bio import SeqIO
    from logzero import logger
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

    logger.info("parsing raw reads file " + args.raw_reads)
    raw_reads = SeqIO.to_dict(SeqIO.parse(args.raw_reads, "fasta"))

    logger.info("adding annotation from " + args.annotation)
    annotation = defaultdict(list)
    for line in open(args.annotation):
        ls = line.strip("\n").split("\t")
        annotation[ls[0]].append((ls[0], int(ls[1]), int(ls[2])) + tuple(ls[3:]))

    logger.info("reading filtered read matches from " + args.processed_reads)
    logger.info("writing coordinates for selected reads to " + args.bed)
    logger.info("writing summary table to " + args.outfile)

    bed = open(args.bed, "w")
    out_table = dict()
    for line in open(args.processed_reads):

        ls = line.strip("\n").split("\t")
        read = ls[0]
        isotype = ls[1]
        mappings = []
        inserts = []
        read_orientation = ls[2]
        for ll in ls[4:]:
            if "insert" not in ll:
                coords = decode_coords(ll)
                mappings.append(coords)
            elif "insert" in ll:
                inserts.append(decode_insert(ll))

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

        for insert in inserts:

            insert_chrom = insert.group("insert_chrom")
            insert_start = int(insert.group("insert_start"))
            insert_end = int(insert.group("insert_end"))
            insert_len = insert_end - insert_start

            insert_anno = set()
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

            cstart = int(insert.group("switch_left")) - 1
            cend = int(insert.group("switch_right")) + 1
            if cend < cstart:
                cstart, cend = cend, cstart

            seq = str(raw_reads[read].seq)
            istart = int(insert.group("istart"))
            iend = int(insert.group("iend"))
            seq = seq[:istart].lower() + seq[istart:iend].upper() + seq[iend:].lower()
            out_table[read] = dict(
                read=read,
                isotype=isotype,
                switch_coords="{0}:{1}-{2}".format(switch_chrom, tstart + 1, tend),
                insert_pos="{0}:{1}-{2}".format(switch_chrom, cstart + 1, cend),
                insert_coords="{0}:{1}-{2}".format(insert_chrom, insert_start + 1, insert_end),
                insert_anno=insert_anno,
                sequence=seq,
            )

            if args.interrupt_for_read is not None and read in args.interrupt_for_read.split(","):
                return raw_reads[read].seq, ls

    bed.close()

    pd.DataFrame.from_dict(out_table, orient="index").to_csv(
        args.outfile, header=True, index=False, sep="\t"
    )
