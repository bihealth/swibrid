"""process LAST output"""


def setup_argparse(parser):
    parser.add_argument(
        "--last",
        dest="last",
        help="""MAF file with aligned reads (LAST output)""",
    )
    parser.add_argument(
        "--telo",
        dest="telo",
        help="""output of blasting reads against telomer repeats""",
    )
    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="""output table with processed reads""",
    )
    parser.add_argument("--stats", dest="stats", help="""output stats""")
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000",
        help="""coordinates of switch region [chr14:106050000-106337000]""",
    )
    parser.add_argument(
        "--switch_annotation",
        dest="switch_annotation",
        help="""bed file with switch annotation""",
    )
    parser.add_argument(
        "--min_switch_match",
        dest="min_switch_match",
        default=100,
        type=int,
        help="""min match length to switch region on both sides [100]""",
    )
    parser.add_argument(
        "--max_switch_overlap",
        dest="max_switch_overlap",
        default=0,
        type=int,
        help="""max overlap between switch parts of different orientation [0]""",
    )
    parser.add_argument(
        "--max_gap",
        dest="max_gap",
        default=100,
        type=int,
        help="""max gap between insert and switch parts (on the read) [100]""",
    )
    parser.add_argument(
        "--max_overlap",
        dest="max_overlap",
        default=50,
        type=int,
        help="""max overlap between switch and insert parts [50]""",
    )
    parser.add_argument(
        "--min_insert_match",
        dest="min_insert_match",
        default=50,
        type=int,
        help="""min length of insert [50]""",
    )
    parser.add_argument(
        "--min_num_switch_matches",
        dest="min_num_switch_matches",
        default=2,
        type=int,
        help="""min number of switch matches in read [2]""",
    )
    parser.add_argument(
        "--min_cov",
        dest="min_cov",
        default=0.9,
        type=float,
        help="""coverage cutoff on reads [0.95 but 0.4 for plasmid]""",
    )
    parser.add_argument(
        "--min_gap",
        dest="min_gap",
        default=-1,
        type=float,
        help="""ignore reads containing small gaps [-1]""",
    )
    parser.add_argument(
        "--isotype_extra",
        dest="isotype_extra",
        default=200,
        type=int,
        help="""extra space added to define isotypes""",
    )
    parser.add_argument(
        "--telo_cutoff",
        dest="telo_cutoff",
        default=90,
        type=float,
        help="""%% identify cutoff for telomer match [90]""",
    )
    parser.add_argument(
        "--min_telo_matchlen",
        dest="min_telo_matchlen",
        type=int,
        default=25,
        help="""minimum matchlen for telomeric repeat matches [25]""",
    )
    parser.add_argument(
        "--sequences",
        dest="sequences",
        help="""save sequence alignment to reference for pseudo-multiple alignment""",
    )
    parser.add_argument(
        "--interrupt_for_read",
        dest="interrupt_for_read",
        help="""interrupt for specified read(s)""",
    )
    parser.add_argument(
        "--blacklist_regions",
        dest="blacklist_regions",
        help="""ignore inserts from blacklist regions defined in this bed file""",
    )


def parse_maf(maf_input, min_gap=50):

    from collections import defaultdict
    import re
    from logzero import logger
    from Bio import AlignIO

    read_matches = []
    read_id = ''
    for ref, read_part in AlignIO.parse(maf_input, "maf", seq_count=2):
        read_start = read_part.annotations["start"]
        ref_start = ref.annotations["start"]
        pos = 0
        num_ref_gaps = 0
        num_read_gaps = 0
        chunks = []
        # find ref gaps of at least min_gap
        for m in re.finditer(r"\-{{{0},}}".format(min_gap), str(ref.seq)):

            new_ref_gaps = ref.seq[pos : m.start()].count("-")
            ref_start = ref.annotations["start"] + pos - num_ref_gaps
            ref_end = ref.annotations["start"] + m.start() - num_ref_gaps - new_ref_gaps

            new_read_gaps = read_part.seq[pos : m.start()].count("-")
            read_start = read_part.annotations["start"] + pos - num_read_gaps
            read_end = read_part.annotations["start"] + m.start() - num_read_gaps - new_read_gaps

            # get sequence of read at non-gap positions in ref
            ref_seq = ref.seq[pos : m.start()]
            aligned_seq = "".join(
                x for k, x in enumerate(read_part.seq[pos : m.start()]) if ref_seq[k] != "-"
            )

            chnk = (
                read_start,
                read_end - read_start,
                read_part.annotations["srcSize"],
                ref.id,
                ref_start,
                ref_end - ref_start,
                ref.annotations["size"],
                read_part.annotations["strand"],
                str(aligned_seq),
            )

            chunks.append(chnk)

            num_read_gaps += new_read_gaps
            num_ref_gaps += m.end() - m.start() + new_ref_gaps
            pos = m.end()

        new_ref_gaps = ref.seq[pos:].count("-")
        ref_start = ref.annotations["start"] + pos - num_ref_gaps
        ref_end = ref.annotations["start"] + len(ref) - num_ref_gaps - new_ref_gaps

        new_read_gaps = read_part.seq[pos:].count("-")
        read_start = read_part.annotations["start"] + pos - num_read_gaps
        read_end = read_part.annotations["start"] + len(ref) - num_read_gaps - new_read_gaps

        # get sequence of read at non-gap positions in ref
        ref_seq = ref.seq[pos:]
        aligned_seq = "".join(x for k, x in enumerate(read_part.seq[pos:]) if ref_seq[k] != "-")

        chnk = (
            read_start,
            read_end - read_start,
            read_part.annotations["srcSize"],
            ref.id,
            ref_start,
            ref_end - ref_start,
            ref.annotations["size"],
            read_part.annotations["strand"],
            str(aligned_seq),
        )

        chunks.append(chnk)

        if read_part.id != read_id:
            if len(read_matches) > 0:
                yield read_id, read_matches
            read_id = read_part.id
            read_matches = []

        read_matches += chunks

    yield read_id, read_matches


def run(args):

    import operator
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from Bio import SeqIO, Seq, SeqRecord
    import gzip
    from logzero import logger
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
        intersect_intervals,
        interval_length,
        merge_intervals,
    )

    if args.telo:
        telo = pd.read_csv(args.telo, index_col=0, sep="\t")

    if args.blacklist_regions is not None:
        blacklist_regions = sorted(
            [
                (line.split()[0],) + tuple(map(int, line.split()[1:3]))
                for line in open(args.blacklist_regions)
            ],
            key=lambda x: x[1],
        )
    else:
        blacklist_regions = []

    (
        switch_chrom,
        switch_start,
        switch_end,
        switch_orientation,
    ) = parse_switch_coords(args.switch_coords)
    switch_anno = read_switch_anno(args.switch_annotation)

    stats = {
        "nreads": 0,
        "no_switch": 0,
        "length_mismatch": 0,
        "orientation_mismatch": 0,
        "inversions": 0,
        "low_cov": 0,
        "switch_order": 0,
        "small_gap": 0,
    }

    logger.info("processing reads from MAF file " + args.last)

    outf = open(args.outfile, "w")
    seq_out = gzip.open(args.sequences, "wt") if args.sequences.endswith(".gz") else open(args.sequences, "w")

    for read, matches in parse_maf(args.last):

        use = True
        stats["nreads"] += 1

        # collect parts matching to switch region
        switch_matches = []
        # collect parts matching elsewhere
        inserts = []
        # get total read coverage

        # get positions of matches to switch region and to elsewhere ("inserts")
        for match in matches:
            (
                read_start,
                read_len,
                tot_read_len,
                ref_chrom,
                ref_start,
                ref_len,
                tot_ref_len,
                orientation,
                aligned_seq,
            ) = match
            if orientation == -1:
                read_start = tot_read_len - read_start - read_len
            read_end = read_start + read_len
            ref_end = ref_start + ref_len

            # find read parts that map to the switch region
            if (
                ref_chrom == switch_chrom
                and ref_start >= switch_start
                and ref_end <= switch_end
                and ref_len >= args.min_switch_match
            ):
                switch_matches.append(
                    (
                        read_start,
                        read_end,
                        ref_chrom,
                        ref_start,
                        ref_end,
                        orientation,
                        aligned_seq,
                    )
                )
            # find read parts that map elsewhere (including parts on the switch chromosome itself)
            elif (
                ref_chrom != switch_chrom or ref_end < switch_start or ref_start > switch_end
            ) and ref_len > args.min_insert_match:
                inserts.append(
                    (
                        read_start,
                        read_end,
                        ref_chrom,
                        ref_start,
                        ref_end,
                        orientation,
                    )
                )

        # check for telomer matches from demultiplexing
        if args.telo and read in telo.index:
            telo_matches = []
            for _, df in telo.loc[[read]].iterrows():
                if df["pident"] > args.telo_cutoff and df["length"] > args.min_telo_matchlen:
                    telo_matches.append(("telo", int(df["qstart"]) - 1, int(df["qend"])))
            for _, read_start, read_end in merge_intervals(telo_matches):
                tot_len = read_end - read_start
                inserts.append((read_start, read_end, "telomer", 0, tot_len, 1))

        # check that matches have consistent orientation (choose the main orientation from read coverage)
        orientations = {"1": 0, "-1": 0}
        for read_start, read_end, _, _, _, orientation, _ in switch_matches:
            orientations[str(orientation)] += read_end - read_start
        main_orientation = max([1, -1], key=lambda x: orientations[str(x)])
        read_overlap = interval_length(
            intersect_intervals(
                sorted(
                    [("", m[0], m[1]) for m in switch_matches if m[5] == 1],
                    key=lambda x: x[1],
                ),
                sorted(
                    [("", m[0], m[1]) for m in switch_matches if m[5] == -1],
                    key=lambda x: x[1],
                ),
            )
        )
        ref_overlap = interval_length(
            intersect_intervals(
                sorted(
                    [("", m[3], m[4]) for m in switch_matches if m[5] == 1],
                    key=lambda x: x[1],
                ),
                sorted(
                    [("", m[3], m[4]) for m in switch_matches if m[5] == -1],
                    key=lambda x: x[1],
                ),
            )
        )

        if read_overlap > args.max_switch_overlap or ref_overlap > args.max_switch_overlap:
            stats["orientation_mismatch"] += 1
            use = False

        # sort switch matches by their order along the read
        switch_matches = sorted(switch_matches, key=operator.itemgetter(0))

        filtered_inserts = []
        # focus now on reads with structure switch - insert - switch
        if len(switch_matches) >= args.min_num_switch_matches:

            for insert in inserts:

                # ignore inserts at beginning or ends of reads
                if insert[0] < min(sw[0] for sw in switch_matches) or insert[1] > max(
                    sw[1] for sw in switch_matches
                ):
                    continue

                # ignore inserts in blacklist regions
                if interval_length(intersect_intervals([insert[2:5]], blacklist_regions)) > 0:
                    continue

                # get adjoining switch matches (along the read)
                # left: last match starting before insert
                switch_left = [sm for sm in switch_matches if sm[0] <= insert[0]][-1]
                # right: first match ending after insert
                switch_right = [sm for sm in switch_matches if sm[1] >= insert[1]][0]

                # separate match if an insert falls completely within a switch match
                if (
                    switch_left == switch_right
                    and switch_left[1] > insert[1]
                    and switch_right[0] < insert[0]
                ):
                    # split length mismatch between read and genomic coordinates evenly
                    extra = (
                        switch_left[4]
                        - switch_left[3]
                        - (insert[0] - switch_left[0])
                        - (switch_left[1] - insert[1])
                    )
                    switch_left = (
                        switch_left[0],
                        insert[0],
                        switch_left[2],
                        switch_left[3],
                        switch_left[3] + insert[0] - switch_left[0] + extra // 2,
                        switch_left[5],
                    )
                    switch_right = (
                        insert[1],
                        switch_left[1],
                        switch_left[2],
                        switch_left[4] - (switch_left[1] - insert[1]) - extra // 2,
                        switch_left[4],
                        switch_left[5],
                    )
                elif switch_left[1] > insert[1] or switch_right[0] < insert[0]:
                    raise Exception("interval error for insert-adjoining switch matches!")

                # number of unmapped bases in read between switch parts and inserts
                gaps = (switch_right[0] - insert[1], insert[0] - switch_left[1])

                # amount of overlapping alignment length for this insert
                overlap = sum(
                    e - s
                    for _, s, e in intersect_intervals(
                        [("", m[0], m[1]) for m in switch_matches],
                        [("", insert[0], insert[1])],
                    )
                )

                # arrange left and right switch matches by genomic orientation
                switch_left, switch_right = sorted(
                    [switch_left, switch_right], key=operator.itemgetter(3)
                )

                if (
                    switch_left
                    and switch_right
                    and all(gap < args.max_gap for gap in gaps)
                    and overlap <= args.max_overlap
                ):
                    filtered_inserts.append(
                        (
                            insert[2],
                            insert[3],
                            insert[4],
                            switch_left[4],
                            gaps[0],
                            switch_right[3],
                            gaps[1],
                            insert[0],
                            insert[1],
                            insert[5],
                        )
                    )

        # take only reads matching to switch region (only once enough?)
        if len(switch_matches) == 0:
            stats["no_switch"] += 1
            use = False

        # merge the switch mappings and sort them by genomic coordinate
        read_mappings = [
            (rm[0][:-1], rm[1], rm[2], rm[0][-1])
            for rm in sorted(
                merge_intervals(
                    (
                        x[2] + ("+" if x[5] == main_orientation else "-"),
                        x[3],
                        x[4],
                    )
                    for x in switch_matches
                ),
                key=lambda x: x[1],
            )
        ]

        # check that mapped regions on the genome are not much shorter than mapped parts of the read
        if sum(sm[4] - sm[3] for sm in switch_matches) < 0.7 * sum(
            sm[1] - sm[0] for sm in switch_matches
        ):
            stats["length_mismatch"] += 1
            use = False

        # check that there are no small gaps < args.min_gap
        switch_maps = defaultdict(list)
        for rec in intersect_intervals(read_mappings, switch_anno, loj=False):
            if rec[3][3].startswith("S"):
                switch_maps[rec[3][3]].append((rec[1], rec[2]))

        gaps = [
            val[i][1] - val[i - 1][0] for val in switch_maps.values() for i in range(1, len(val))
        ]
        if len(gaps) > 0 and min(gaps) < args.min_gap:
            stats["small_gap"] += 1
            use = False

        # check that more than args.min_cov of the read maps to the genome
        if any(fi[8] < fi[7] for fi in filtered_inserts):
            raise Exception("stop")
        tot_cov = sum(sm[1] - sm[0] for sm in switch_matches) + sum(
            fi[8] - fi[7] for fi in filtered_inserts
        )
        if tot_cov / tot_read_len < args.min_cov:
            stats["low_cov"] += 1
            use = False

        # check that switch matches are in consistent order along the read and the genome (but ignore inverted segments)
        read_order = np.argsort([sm[0] for sm in switch_matches if sm[5] == main_orientation])
        genomic_order = np.argsort([sm[3] for sm in switch_matches if sm[5] == main_orientation])
        if main_orientation == -1:
            genomic_order = genomic_order[::-1]
        if not (read_order == genomic_order).all():
            stats["switch_order"] += 1
            use = False

        if args.interrupt_for_read is not None and read in args.interrupt_for_read.split(","):
            raise Exception("breaking for read {0}".format(read))

        if not use:
            continue

        if orientations["1"] > 0 and orientations["-1"] > 0:
            stats["inversions"] += 1

        # isotype
        last_map = read_mappings[-1] if switch_orientation == "+" else read_mappings[0]
        tmp = intersect_intervals(
            [
                (
                    last_map[0],
                    last_map[1] - args.isotype_extra,
                    last_map[2] + args.isotype_extra,
                )
            ],
            switch_anno,
            loj=True,
        )
        isotype = ",".join(rec[3][3] for rec in tmp)
        if len(isotype) == 0:
            isotype = "none"

        outf.write(
            "{0}\t{1}\t{2}\t{3:.4f}".format(
                read,
                isotype,
                "+" if main_orientation == 1 else "-",
                tot_cov / tot_read_len,
            )
        )
        for chrom, start, end, orientation in read_mappings:
            outf.write("\t{0}:{1}-{2}:{3}".format(chrom, start, end, orientation))
        for (
            chrom,
            start,
            end,
            left,
            gapl,
            right,
            gapr,
            istart,
            iend,
            orientation,
        ) in filtered_inserts:
            outf.write(
                ("\tinsert_{0}:{1}_{2}_{3}-{4}_{5}_" "{6}:{7}-{8}_{9}_{10}:{11}").format(
                    switch_chrom,
                    left,
                    gapl,
                    istart,
                    iend,
                    "+" if orientation == 1 else "-",
                    chrom,
                    start,
                    end,
                    gapr,
                    switch_chrom,
                    right,
                )
            )

        outf.write("\n")

        if args.sequences is not None:
            seq = ["-" * (rm[2] - rm[1]) for rm in read_mappings]
            for (
                read_start,
                read_end,
                ref_chrom,
                ref_start,
                ref_end,
                orientation,
                aligned_seq,
            ) in switch_matches:
                for k, rm in enumerate(read_mappings):
                    if ref_start >= rm[1] and ref_end <= rm[2]:
                        seq[k] = (
                            seq[k][: (ref_start - rm[1])]
                            + aligned_seq
                            + seq[k][(ref_end - rm[1]) :]
                        )
            sequences = []
            for k, rm in enumerate(read_mappings):
                sequences.append(
                    SeqRecord.SeqRecord(
                        Seq.Seq(seq[k]),
                        id="{0}@{1}:{2}-{3}:{4}".format(read, *rm),
                    )
                )
            SeqIO.write(sequences, seq_out, "fasta")

    outf.close()
    seq_out.close()

    logger.info('done. saving stats to ' + args.stats)
    if args.stats:
        pd.Series(stats).to_csv(args.stats)

        
