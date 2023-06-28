"""process alignments"""


def setup_argparse(parser):
    parser.add_argument(
        "--alignments",
        dest="alignments",
        help="""MAF file (LAST output) or SAM file (minimap2 output) with aligned reads (LAST output)""",
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
        default=75,
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
    parser.add_argument(
        "--realign_breakpoints",
        dest="realign_breakpoints",
        help="""re-align reads around breakpoints and save results to this file (requires --raw_reads and --genome)""",
    )
    parser.add_argument("--raw_reads", dest="raw_reads", help="""fasta file with unaligned reads""")
    parser.add_argument("--genome", dest="genome", help="""reference genome file""")


def parse_sam(sam_input, min_gap=75):
    import pysam

    read_matches = []
    read_id = ""

    for rec in pysam.Samfile(sam_input):
        if rec.is_unmapped or rec.seq is None or rec.is_secondary:
            continue

        clip_start = rec.cigartuples[0][1] if rec.cigartuples[0][0] == 5 else 0
        clip_end = rec.cigartuples[-1][1] if rec.cigartuples[-1][0] == 5 else 0
        read_len = len(rec.seq) + clip_start + clip_end
        aligned_seq = ""
        al_pos = 0
        read_pos = 0
        ref_pos = 0
        read_start_chunk = clip_start
        ref_start_chunk = 0
        insertions = 0
        chunks = []
        for op, n in rec.cigartuples:
            if op in [0, 7, 8]:  # alignment match
                aligned_seq += str(rec.seq[al_pos : al_pos + n])
                al_pos += n
                read_pos += n
                ref_pos += n

            elif op == 1:  # insertion
                if n > min_gap:  # split insertions of more than min_gap nucleotides
                    chnk = (
                        read_start_chunk,
                        clip_start + read_pos - read_start_chunk,
                        read_len,
                        rec.reference_name,
                        rec.reference_start + ref_start_chunk,
                        ref_pos - ref_start_chunk,
                        -1 if rec.is_reverse else 1,
                        str(aligned_seq),
                    )
                    chunks.append(chnk)
                    aligned_seq = ""
                    read_start_chunk = clip_start + read_pos + n
                    ref_start_chunk = ref_pos
                    insertions += n

                al_pos += n
                read_pos += n

            elif op == 2:  # deletion
                aligned_seq += "-" * n
                ref_pos += n

            elif op == 4:  # soft clip
                if read_pos == 0:  # add to clip_start if at beginning
                    read_start_chunk += n
                    clip_start += n
                else:
                    clip_end += n
                al_pos += n
        
            elif op != 5:
                raise ValueError("CIGAR op {0} not implemented".format(op))

        chnk = (
            read_start_chunk,
            clip_start + read_pos - read_start_chunk,
            read_len,
            rec.reference_name,
            rec.reference_start + ref_start_chunk,
            ref_pos - ref_start_chunk,
            -1 if rec.is_reverse else 1,
            str(aligned_seq),
        )
        chunks.append(chnk)

        if rec.query_name != read_id:
            if len(read_matches) > 0:
                yield read_id, read_matches
            read_id = rec.query_name
            read_matches = []

        assert sum(c[5] for c in chunks) == rec.reference_length, 'ref lengths dont match!'
        assert sum(c[1] for c in chunks) + clip_start + clip_end + insertions == read_len, 'read lengths dont match!'

        read_matches += chunks

    yield read_id, read_matches


def parse_maf(alignments, min_gap=75):
    import re
    from Bio import AlignIO

    read_matches = []
    read_id = ""

    for ref, read in AlignIO.parse(alignments, "maf", seq_count=2):
        read_name = read.id
        read_start = read.annotations["start"]
        read_len = read.annotations["srcSize"]
        read_seq = str(read.seq)
        ref_name = ref.id
        ref_start = ref.annotations["start"]
        ref_seq = str(ref.seq)
        orientation = read.annotations["strand"]

        pos = 0
        num_ref_gaps = 0
        num_read_gaps = 0
        chunks = []
        # find ref gaps of at least min_gap
        for m in re.finditer(r"\-{{{0},}}".format(min_gap), ref_seq):
            new_ref_gaps = ref_seq[pos : m.start()].count("-")
            ref_start_chunk = ref_start + pos - num_ref_gaps
            ref_end_chunk = ref_start + m.start() - num_ref_gaps - new_ref_gaps

            new_read_gaps = read_seq[pos : m.start()].count("-")
            read_start_chunk = read_start + pos - num_read_gaps
            read_end_chunk = read_start + m.start() - num_read_gaps - new_read_gaps

            # get sequence of read at non-gap positions in ref
            ref_seq_chunk = ref_seq[pos : m.start()]
            aligned_seq_chunk = "".join(
                x for k, x in enumerate(read_seq[pos : m.start()]) if ref_seq_chunk[k] != "-"
            )

            chnk = (
                read_start_chunk,
                read_end_chunk - read_start_chunk,
                read_len,
                ref_name,
                ref_start_chunk,
                ref_end_chunk - ref_start_chunk,
                orientation,
                str(aligned_seq_chunk),
            )

            chunks.append(chnk)

            num_read_gaps += new_read_gaps
            num_ref_gaps += m.end() - m.start() + new_ref_gaps
            pos = m.end()

        new_ref_gaps = ref_seq[pos:].count("-")
        ref_start_chunk = ref_start + pos - num_ref_gaps
        ref_end_chunk = ref_start + len(ref_seq) - num_ref_gaps - new_ref_gaps

        new_read_gaps = read_seq[pos:].count("-")
        read_start_chunk = read_start + pos - num_read_gaps
        read_end_chunk = read_start + len(ref_seq) - num_read_gaps - new_read_gaps

        # get sequence of read at non-gap positions in ref
        ref_seq_chunk = ref_seq[pos:]
        aligned_seq_chunk = "".join(
            x for k, x in enumerate(read_seq[pos:]) if ref_seq_chunk[k] != "-"
        )

        chnk = (
            read_start_chunk,
            read_end_chunk - read_start_chunk,
            read_len,
            ref_name,
            ref_start_chunk,
            ref_end_chunk - ref_start_chunk,
            orientation,
            str(aligned_seq_chunk),
        )

        chunks.append(chnk)

        if read_name != read_id:
            if len(read_matches) > 0:
                yield read_id, read_matches
            read_id = read_name
            read_matches = []

        read_matches += chunks

    yield read_id, read_matches


def combine_alignments(al1, al2, pad):

    s1 = ""
    s0 = ""
    s2 = ""
    m1 = ""
    m2 = ""
    i = 0
    j = 0
    s1A = al1[0]
    s1B = al1[1]
    s2A = al2[0]
    s2B = al2[1]
    while i < len(s1A) and j < len(s2A):
        if s1A[i] == "-":
            s1 += s1B[i].lower()
            m1 += " "
            s0 += "-"
            s2 += "-"
            m2 += " "
            i += 1
        elif s2A[j] == "-":
            s1 += "-"
            m1 += " "
            s0 += "-"
            m2 += " "
            s2 += s2B[j].lower()
            j += 1
        elif s1A[i] == s2A[j]:
            s1 += s1B[i]
            m1 += "|" if s1A[i] == s1B[i] else (" " if s1B[i] == "-" else ".")
            s0 += s1A[i]
            m2 += "|" if s2A[j] == s2B[j] else (" " if s2B[j] == "-" else ".")
            s2 += s2B[j]
            i += 1
            j += 1
        else:
            raise Exception("stop")
        if len(s1A[:i].replace("-", "")) != len(s2A[:j].replace("-", "")):
            raise Exception("stop")
        if (
            len(s1A[:i].replace("-", "")) == pad or len(s1A[i:].replace("-", "")) == pad
        ) and not (s1A[i] == "-" or s2A[j] == "-"):
            s1 += "/"
            s0 += "/"
            s2 += "/"
            m1 += "/"
            m2 += "/"
    return s1, m1, s0, m2, s2


def realign_breakpoints(matches, genome, read_seq, pad=20):
    from Bio import Align

    def revcomp(seq):
        from swibrid.scripts.utils import RC

        return "".join(RC[s] for s in seq[::-1])

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -4
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -2

    res = []
    for i in range(len(matches) - 1):
        rseq = read_seq[
            min(matches[i][1], matches[i + 1][0])
            - pad : max(matches[i][1], matches[i + 1][0])
            + pad
        ]
        orientation = (matches[i][5], matches[i + 1][5])
        if orientation[0] == -1:
            pos_left = matches[i][2] + ":" + str(matches[i][3])
            gseq1 = revcomp(
                genome.fetch(matches[i][2], matches[i][3] - pad, matches[i][3] + pad).upper()
            )
        else:
            pos_left = matches[i][2] + ":" + str(matches[i][4])
            gseq1 = genome.fetch(matches[i][2], matches[i][4] - pad, matches[i][4] + pad).upper()
        if orientation[1] == -1:
            pos_right = matches[i + 1][2] + ":" + str(matches[i + 1][4])
            gseq2 = revcomp(
                genome.fetch(
                    matches[i + 1][2], matches[i + 1][4] - pad, matches[i + 1][4] + pad
                ).upper()
            )
        else:
            pos_right = matches[i + 1][2] + ":" + str(matches[i + 1][3])
            gseq2 = genome.fetch(
                matches[i + 1][2], matches[i + 1][3] - pad, matches[i + 1][3] + pad
            ).upper()

        #al1 = pairwise2.align.globalms(rseq, gseq1, 2, -4, -4, -2)
        #al2 = pairwise2.align.globalms(rseq, gseq2, 2, -4, -4, -2)

        al1 = aligner.align(rseq, gseq1)
        al2 = aligner.align(rseq, gseq2)

        s1, m1, s0, m2, s2 = combine_alignments(al1[0], al2[0], pad)

        m1s = m1.split("/")
        m2s = m2.split("/")
        s0s = s0.split("/")

        n_homology = 0
        n_untemplated = 0

        # check for homology on the left side of the breakpoint
        k = m2s[0].rfind(" ") + 1
        n_homology += sum(x == "|" and y == "|" for x, y in zip(m1s[0][k:], m2s[0][k:]))
        # check for homology on the right side of the breakpoint
        k = m1s[-1].find(" ")
        n_homology += sum(x == "|" and y == "|" for x, y in zip(m1s[-1][:k], m2s[-1][:k]))

        # check for untemplated or homologous nucleotides between the breakpoints on the read
        if len(m1s) > 2 and len(m2s) > 2:
            n_homology += sum(
                x == "|" and y == "|" and z != "-" for x, y, z in zip(m1s[1], m2s[1], s0s[1])
            )
            n_untemplated = sum(
                x == " " and y == " " and z != "-" for x, y, z in zip(m1s[1], m2s[1], s0s[1])
            )

        r = {
            "pos_left": pos_left,
            "pos_right": pos_right,
            "n_homology": n_homology,
            "n_untemplated": n_untemplated,
            "orientation": orientation,
            "type": "insert" if "insert" in [matches[i][-1], matches[i + 1][-1]] else "switch",
            "left_seq": s1,
            "match_left": m1,
            "read_seq": s0,
            "match_right": m2,
            "right_seq": s2,
        }

        res.append(r)

    return res


def run(args):
    import os
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

    if args.realign_breakpoints is not None:
        import pysam

        logger.info("getting raw reads from " + args.raw_reads)
        raw_reads = SeqIO.index(args.raw_reads, "fasta")
        logger.info("loading genomes from " + args.genome)
        genome = pysam.FastaFile(args.genome)

    if args.telo:
        telo = pd.read_csv(args.telo, index_col=0, sep="\t")

    if args.blacklist_regions is not None and os.path.isfile(args.blacklist_regions):
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
        "overlap_mismatch": 0,
        "inversions": 0,
        "low_cov": 0,
        "switch_order": 0,
        "small_gap": 0,
    }

    logger.info("processing reads from " + args.alignments)
    if args.alignments.endswith(".sam"):
        alignments = parse_sam(args.alignments, min_gap=args.min_gap)
    else:
        alignments = parse_maf(args.alignments, min_gap=args.min_gap)

    outf = open(args.outfile, "w")
    seq_out = (
        gzip.open(args.sequences, "wt")
        if args.sequences.endswith(".gz")
        else open(args.sequences, "w")
    )

    realignments = {}

    for read, matches in alignments:
        use = True
        stats["nreads"] += 1

        # collect parts matching to switch region
        switch_matches = []
        # collect parts matching elsewhere
        inserts = []
        # get total read coverage

        # get positions of matches to switch region and to elsewhere ("inserts")
        for (
            read_start,
            read_len,
            tot_read_len,
            ref_chrom,
            ref_start,
            ref_len,
            orientation,
            aligned_seq,
        ) in matches:

            assert read_len > 0 and ref_len > 0 and tot_read_len > 0, 'negative lengths in parsed alignments!'

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

        # check overlap of aligned segments
        read_overlap = sum(
            [
                switch_matches[i][1] - switch_matches[i + 1][0]
                for i in range(len(switch_matches) - 1)
            ]
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
            stats["overlap_mismatch"] += 1
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
                    raise ValueError("interval error for insert-adjoining switch matches!")

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
        #switch_maps = defaultdict(list)
        #for rec in intersect_intervals(read_mappings, switch_anno, loj=False):
        #    if rec[3][3].startswith("S"):
        #        switch_maps[rec[3][3]].append((rec[1], rec[2]))

        #gaps = [
        #    val[i][1] - val[i - 1][0] for val in switch_maps.values() for i in range(1, len(val))
        #]
        #if len(gaps) > 0 and min(gaps) < args.min_gap:
        #    stats["small_gap"] += 1
        #    use = False

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

        # re-align breakpoints
        if args.realign_breakpoints is not None:
            processed_matches = sorted(
                [sm[:6] + ("switch",) for sm in switch_matches]
                + [(f[7], f[8], f[0], f[1], f[2], f[9], "insert") for f in filtered_inserts if f[0] != 'telomer'],
                key=operator.itemgetter(0),
            )
            realignments[read] = realign_breakpoints(
                processed_matches, genome, str(raw_reads[read].seq)
            )

    outf.close()
    seq_out.close()

    logger.info("done. saving stats to " + args.stats)
    if args.stats:
        pd.Series(stats).to_csv(args.stats)

    if args.realign_breakpoints is not None:
        logger.info("saving breakpoint re-alignments to " + args.realign_breakpoints)
        df = pd.DataFrame(
            [pd.Series(al, name=read) for read, als in realignments.items() for al in als]
        )
        df.to_csv(args.realign_breakpoints)
