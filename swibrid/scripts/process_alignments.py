"""\
process alignments:
for an input file with aligned reads (MAF if LAST output, SAM if minimap2 output),
get alignments to switch regions and elsewhere (potential inserts)
output table contains
- isotype
- read orientation
- read coverage
- fraction of read sequence mapping to the same genomic regions
- mapped switch region segments
- inserts
aligned sequences can be saved separately (necessary to then construct a pseudo MSA)
reads are removed if
- they are too short
- they don't contain a forward and reverse primer
- they contain internal primers
- they contain no switch region
- if mapped regions on the genome are much shorter than mapped parts of the read
- there's too much overlap between different alignments on the read, or on the genome
- too little of the read maps
- alignments to the switch region are in the wrong order
- the isotype cannot be determined
if `--realign_breakpoints` is set, 20nt on each side of a breakpoint are re-aligned, and statistics like number of homologous or untemplated bases are extracted
"""


def setup_argparse(parser):
    parser.add_argument(
        "--alignments",
        dest="alignments",
        help="""required: MAF file (LAST output) or SAM file (minimap2 output) with aligned reads (LAST output)""",
    )
    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="""required: output table with processed reads""",
    )
    parser.add_argument("--info", dest="info", help="""required: csv file with read info""")
    parser.add_argument(
        "--sequences",
        dest="sequences",
        help="""save sequence alignment to reference for pseudo-multiple alignment""",
    )
    parser.add_argument("--stats", dest="stats", help="""output stats on removed reads""")
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
        "--min_length",
        dest="min_length",
        default=500,
        type=int,
        help="""minimum read length""",
    )
    parser.add_argument(
        "--only_complete",
        dest="only_complete",
        action="store_true",
        default=False,
        help="""use only complete reads with two primers""",
    )
    parser.add_argument(
        "--keep_internal",
        dest="keep_internal",
        action="store_true",
        default=False,
        help="""keep reads with internal primers""",
    )
    parser.add_argument(
        "--remove_duplicates",
        dest="remove_duplicates",
        action="store_true",
        default=False,
        help="""remove duplicate reads (based on UMIs in read IDs)""",
    )
    parser.add_argument(
        "--ignore_order",
        dest="ignore_order",
        action="store_true",
        default=False,
        help="""don't check order of matches along read""",
    )
    parser.add_argument(
        "--telo",
        dest="telo",
        default="",
        help="""output of blasting reads against telomer repeats""",
    )
    parser.add_argument(
        "--min_switch_match",
        dest="min_switch_match",
        default=100,
        type=int,
        help="""min match length to switch region""",
    )
    parser.add_argument(
        "--max_switch_overlap",
        dest="max_switch_overlap",
        default=0,
        type=int,
        help="""max overlap between switch parts of different orientation [0]""",
    )
    parser.add_argument(
        "--max_insert_gap",
        dest="max_insert_gap",
        default=100,
        type=int,
        help="""max gap between insert and switch parts (on the read) [100]""",
    )
    parser.add_argument(
        "--max_insert_overlap",
        dest="max_insert_overlap",
        default=50,
        type=int,
        help="""max overlap between switch and insert parts (on the read) [50]""",
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
        default=1,
        type=int,
        help="""min number of switch matches in read [1]""",
    )
    parser.add_argument(
        "--min_cov",
        dest="min_cov",
        default=0.9,
        type=float,
        help="""coverage cutoff on reads [0.95 but 0.4 for plasmid]""",
    )
    parser.add_argument(
        "--max_gap",
        dest="max_gap",
        default=75,
        type=float,
        help="""split alignments with gaps beyond that value [75]""",
    )
    parser.add_argument(
        "--isotype_extra",
        dest="isotype_extra",
        default=200,
        type=int,
        help="""extra space added to define isotypes [200]""",
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
        "--blacklist_regions",
        dest="blacklist_regions",
        nargs="?",
        help="""ignore inserts from blacklist regions defined in this bed file""",
    )
    parser.add_argument(
        "--realign_breakpoints",
        dest="realign_breakpoints",
        help="""re-align reads around breakpoints and save results to this file (requires --raw_reads and --genome)""",
    )
    parser.add_argument(
        "--realignment_penalties",
        dest="realignment_penalties",
        default="ont",
        help="""re-alignment penalties ("ont", "hifi" or "sr") [ont]""",
    )
    parser.add_argument(
        "--raw_reads",
        dest="raw_reads",
        help="""fasta file with unaligned reads (comma-separated list of mates for paired-end mode)""",
    )
    parser.add_argument("--genome", dest="genome", help="""reference genome file""")
    parser.add_argument(
        "--keep_secondary_alignments",
        dest="keep_secondary_alignments",
        action="store_true",
        default=False,
        help="""keep secondary alignments in SAM file (but conflicting alignments are not resolved!)""",
    )
    parser.add_argument(
        "--paired_end_mode",
        dest="paired_end_mode",
        action="store_true",
        default=False,
        help="""use paired-end mode (--raw_reads needs to be a comma-separated list of mates)""",
    )
    parser.add_argument(
        "--interrupt_for_read",
        dest="interrupt_for_read",
        help="""(for debugging) interrupt for specified read(s)""",
    )


def parse_sam(sam_input, max_gap=75, keep_secondary=False):
    import pysam

    read_matches = []
    read_id = ""

    for rec in pysam.Samfile(sam_input):
        if rec.is_unmapped or rec.seq is None:
            continue
        if not keep_secondary and rec.is_secondary:
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
        deletions = 0
        chunks = []
        for op, n in rec.cigartuples:
            if op in [0, 7, 8]:  # alignment match
                aligned_seq += str(rec.seq[al_pos : al_pos + n])
                al_pos += n
                read_pos += n
                ref_pos += n

            elif op == 1:  # insertion
                if n >= max_gap:  # split insertions beyond max_gap nucleotides
                    chnk = (
                        read_start_chunk,
                        clip_start + read_pos - read_start_chunk,
                        read_len,
                        rec.reference_name,
                        rec.reference_start + ref_start_chunk,
                        ref_pos - ref_start_chunk,
                        -1 if rec.is_reverse else 1,
                        str(aligned_seq),
                        "R1" if rec.is_read1 else "R2",
                    )
                    chunks.append(chnk)
                    aligned_seq = ""
                    read_start_chunk = clip_start + read_pos + n
                    ref_start_chunk = ref_pos
                    insertions += n

                al_pos += n
                read_pos += n

            elif op in [2, 3]:  # deletion / skip
                if n >= max_gap:  # split deletions beyond max_gap nucleotides
                    chnk = (
                        read_start_chunk,
                        clip_start + read_pos - read_start_chunk,
                        read_len,
                        rec.reference_name,
                        rec.reference_start + ref_start_chunk,
                        ref_pos - ref_start_chunk,
                        -1 if rec.is_reverse else 1,
                        str(aligned_seq),
                        "R1" if rec.is_read1 else "R2",
                    )
                    chunks.append(chnk)
                    aligned_seq = ""
                    read_start_chunk = clip_start + read_pos
                    ref_start_chunk = ref_pos + n
                    deletions += n

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
            "R1" if rec.is_read1 else "R2",
        )
        chunks.append(chnk)

        if rec.query_name != read_id:
            if len(read_matches) > 0:
                yield read_id, read_matches
            read_id = rec.query_name
            read_matches = []

        assert (
            sum(c[5] for c in chunks) + deletions == rec.reference_length
        ), "ref lengths dont match!"
        assert (
            sum(c[1] for c in chunks) + clip_start + clip_end + insertions == read_len
        ), "read lengths dont match!"

        read_matches += chunks

    yield read_id, read_matches


def parse_maf(alignments, max_gap=75):
    import re
    import gzip
    from Bio import AlignIO

    read_matches = []
    read_id = ""

    for ref, read in AlignIO.parse(
        gzip.open(alignments, "rt") if alignments.endswith(".gz") else open(alignments),
        "maf",
        seq_count=2,
    ):
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
        # find ref gaps beyond max_gap
        for m in re.finditer(r"\-{{{0},}}".format(max_gap), ref_seq):
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
                "R1",
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
            "R1",
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
            if s1B[i] != "N":
                s1 += s1B[i].lower()
                m1 += " "
                s0 += "-"
                s2 += "-"
                m2 += " "
            i += 1
        elif s2A[j] == "-":
            if s2B[j] != "N":
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
            raise ValueError("no match in realignment!")

        assert len(s1A[:i].replace("-", "")) == len(
            s2A[:j].replace("-", "")
        ), "mismatched sequence lengths in combine_alignments"

        if (len(s1A[:i].replace("-", "")) == pad or len(s1A[i:].replace("-", "")) == pad) and not (
            s1A[i] == "-" or s2A[j] == "-"
        ):
            s1 += "/"
            s0 += "/"
            s2 += "/"
            m1 += "/"
            m2 += "/"

    return s1, m1, s0, m2, s2


def realign_breakpoints(matches, genome_dict, read_seq, pad=20, penalties="ont"):
    import numpy as np
    from Bio import Align

    def revcomp(seq):
        from swibrid.scripts.utils import RC

        return "".join(RC[s] for s in seq[::-1])

    def fetch(match, genome_dict, ind, pad):
        if "insert" in match[6]:
            seq = (
                genome_dict["genome"]
                .fetch(
                    match[2],
                    max(0, match[ind] - pad),
                    min(genome_dict["genome"].get_reference_length(match[2]), match[ind] + pad),
                )
                .upper()
            )
        else:
            seq = genome_dict["switch_genome"][
                max(0, match[ind] - pad - genome_dict["offset"]) : min(
                    len(genome_dict["switch_genome"]), match[ind] + pad - genome_dict["offset"]
                )
            ].upper()
        assert len(seq) <= 2 * pad, "trying to align a sequence of length {0}".format(len(seq))

        return seq

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    if penalties == "ont":  # use minimap2 presets
        aligner.match_score = 2  # -A
        aligner.mismatch_score = -4  # -B
        aligner.open_gap_score = -4  # -O
        aligner.extend_gap_score = -2  # -E
    elif penalties == "hifi":
        aligner.match_score = 1  # -A
        aligner.mismatch_score = -4  # -B
        aligner.open_gap_score = -6  # -O
        aligner.extend_gap_score = -2  # -E
    elif penalties == "sr":
        aligner.match_score = 2  # -A
        aligner.mismatch_score = -8  # -B
        aligner.open_gap_score = -12  # -O
        aligner.extend_gap_score = -2  # -E
    else:
        raise ValueError("unknown gap penalties {0}".format(penalties))

    # genome = genome_dict['genome']
    res = []
    for i in range(len(matches) - 1):
        # ignore breaks between different read mates
        if matches[i][7] != matches[i + 1][7]:
            continue

        if type(read_seq) is tuple:
            read_seq_here = read_seq[0 if matches[i][7] else 1]
        else:
            read_seq_here = read_seq

        rseq = read_seq_here[
            min(matches[i][1], matches[i + 1][0])
            - pad : max(matches[i][1], matches[i + 1][0])
            + pad
        ]
        orientation = (matches[i][5], matches[i + 1][5])
        if orientation[0] == -1:
            pos_left = matches[i][2] + ":" + str(matches[i][3])
            gseq1 = revcomp(fetch(matches[i], genome_dict, 3, pad))
        else:
            pos_left = matches[i][2] + ":" + str(matches[i][4])
            gseq1 = fetch(matches[i], genome_dict, 4, pad)

        if orientation[1] == -1:
            pos_right = matches[i + 1][2] + ":" + str(matches[i + 1][4])
            gseq2 = revcomp(fetch(matches[i + 1], genome_dict, 4, pad))
        else:
            pos_right = matches[i + 1][2] + ":" + str(matches[i + 1][3])
            gseq2 = fetch(matches[i + 1], genome_dict, 3, pad)

        if gseq1 == "" or gseq2 == "" or rseq == "":
            continue

        # mask identical bases for small gaps
        lchrom = pos_left.split(":")[0]
        lpos = int(pos_left.split(":")[1])
        rchrom = pos_right.split(":")[0]
        rpos = int(pos_right.split(":")[1])
        posdiff = np.abs(rpos - lpos)
        if lchrom == rchrom and posdiff < pad:
            gseq1 = gseq1[: -(pad - posdiff)] + "N" * (pad - posdiff)
            gseq2 = "N" * (pad - posdiff) + gseq2[(pad - posdiff) :]

        al1 = aligner.align(rseq, gseq1)
        al2 = aligner.align(rseq, gseq2)

        s1, m1, s0, m2, s2 = combine_alignments(al1[0], al2[0], pad)

        m1s = m1.split("/")
        m2s = m2.split("/")
        s0s = s0.split("/")

        n_homology = 0
        n_untemplated = 0

        # check for homology on the left side of the breakpoint
        # k = m2s[0].rfind(" ") + 1
        # n_homology += sum(x == "|" and y == "|" for x, y in zip(m1s[0][k:], m2s[0][k:]))
        for k in range(len(m2s[0]) - 1, -1, -1):
            if m1s[0][k] == "|" and m2s[0][k] == "|":
                n_homology += 1
            if m2s[0][k] == " " and s0s[0][k] != "-":
                break
        # check for homology on the right side of the breakpoint
        # k = m1s[-1].find(" ")
        # n_homology += sum(x == "|" and y == "|" for x, y in zip(m1s[-1][:k], m2s[-1][:k]))
        for k in range(len(m1s[-1])):
            if m1s[-1][k] == "|" and m2s[-1][k] == "|":
                n_homology += 1
            if m1s[-1][k] == " " and s0s[-1][k] != "-":
                break

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
            "type": "insert" if "insert" in [matches[i][6], matches[i + 1][6]] else "switch",
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
    from Bio import SeqIO, Seq, SeqRecord
    import gzip
    import re
    from collections import defaultdict
    from logzero import logger
    from .utils import (
        parse_switch_coords,
        read_switch_anno,
        intersect_intervals,
        interval_length,
        merge_intervals,
    )

    if os.path.isfile(args.telo):
        telo = pd.read_csv(args.telo, index_col=0, sep="\t")
    else:
        telo = None

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
        "nreads_initial": 0,
        "nreads_mapped": 0,
        "nreads_removed_short": 0,
        "nreads_removed_incomplete": 0,
        "nreads_removed_no_info": 0,
        "nreads_removed_internal_primer": 0,
        "nreads_removed_no_switch": 0,
        "nreads_removed_length_mismatch": 0,
        "nreads_removed_overlap_mismatch": 0,
        "nreads_inversions": 0,
        "nreads_removed_low_cov": 0,
        "nreads_removed_switch_order": 0,
        "nreads_removed_no_isotype": 0,
        "nreads_removed_duplicate": 0,
    }

    if args.realign_breakpoints is not None:
        processed_matches = {}

    logger.info("reading read info from " + args.info)
    read_info = pd.read_csv(args.info, header=0, index_col=0)
    assert read_info.index.is_unique, "index of info file is not unique!"

    stats["nreads_initial"] = read_info.shape[0]

    logger.info("processing reads from " + args.alignments)
    if args.alignments.endswith(".bam") or args.alignments.endswith(".bam"):
        alignments = parse_sam(
            args.alignments, max_gap=args.max_gap, keep_secondary=args.keep_secondary_alignments
        )
    else:
        alignments = parse_maf(args.alignments, max_gap=args.max_gap)

    outf = open(args.outfile, "w")
    header = "read\tisotype\torientation\tfrac_mapped\tfrac_mapped_multi\tswitch_mappings\tinserts"
    if args.paired_end_mode:
        header += "\tmate_breaks"
    outf.write(header + "\n")

    known_reads = set()

    if args.sequences is not None:
        logger.info("saving aligned sequences to {0}".format(args.sequences))
        seq_out = (
            gzip.open(args.sequences, "wt")
            if args.sequences.endswith(".gz")
            else open(args.sequences, "w")
        )

    if args.remove_duplicates:
        processed_umis = defaultdict(dict)

    for read, matches in alignments:
        if read in known_reads:
            logger.warn("read {0} already encountered; skipping".format(read))
            continue
        known_reads.add(read)

        use = True
        stats["nreads_mapped"] += 1

        if read not in read_info.index:
            stats["nreads_removed_no_info"] += 1
            continue

        # filter reads based on length
        if read_info.loc[read, "length"] < args.min_length:
            stats["nreads_removed_short"] += 1
            use = False

        # check for primers
        if pd.isna(read_info.loc[read, "primers"]):
            primers = set()
        else:
            primers = set(
                map(
                    lambda x: re.sub("primer_|_[1-9]", "", x.split("@")[0]),
                    read_info.loc[read, "primers"].split(";"),
                )
            )

        if args.only_complete and not (
            any("rv" in p for p in primers) and any("fw" in p for p in primers)
        ):
            stats["nreads_removed_incomplete"] += 1
            use = False

        internal = read_info.loc[read, "internal"]
        if not args.keep_internal and not pd.isna(internal) and "primer" in internal:
            stats["nreads_removed_internal_primer"] += 1
            use = False

        if not use:
            continue

        # collect parts matching to switch region
        switch_matches = []
        # collect parts matching elsewhere
        inserts = []

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
            mate,
        ) in matches:
            if read_len <= 0 or ref_len <= 0:
                continue

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
                        mate,
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
                        mate,
                    )
                )

        # check for telomer matches from demultiplexing
        if telo is not None and read in telo.index:
            telo_matches = []
            for _, df in telo.loc[[read]].iterrows():
                if df["pident"] > args.telo_cutoff and df["length"] > args.min_telo_matchlen:
                    telo_matches.append(("telo", int(df["qstart"]) - 1, int(df["qend"])))
            for _, read_start, read_end in merge_intervals(telo_matches):
                tot_len = read_end - read_start
                inserts.append((read_start, read_end, "telomer", 0, tot_len, 1))

        # check that matches have consistent orientation (choose the main orientation from read coverage)
        orientations = {"1": 0, "-1": 0}
        for read_start, read_end, _, _, _, orientation, _, mate in switch_matches:
            orientations[str(orientation if mate == "R1" else -orientation)] += (
                read_end - read_start
            )
        main_orientation = max([1, -1], key=lambda x: orientations[str(x)])

        # sort switch matches by their order along the read (mates)
        switch_matches = sorted(switch_matches, key=lambda x: (x[0], x[-1]))

        # check overlap of aligned segments on the read, or of different read mates in paired-end mode
        read_overlap = sum(
            [
                switch_matches[i][1] - switch_matches[i + 1][0]
                for i in range(len(switch_matches) - 1)
                if switch_matches[i][-1] == switch_matches[i + 1][-1]
            ]
        )
        ref_overlap = interval_length(
            intersect_intervals(
                sorted(
                    [
                        ("", m[3], m[4])
                        for m in switch_matches
                        if (m[5] if m[-1] == "R1" else -m[5]) == 1
                    ],
                    key=lambda x: x[1],
                ),
                sorted(
                    [
                        ("", m[3], m[4])
                        for m in switch_matches
                        if (m[5] if m[-1] == "R1" else -m[5]) == -1
                    ],
                    key=lambda x: x[1],
                ),
            )
        )
        mate_overlap = interval_length(
            intersect_intervals(
                sorted(
                    [("", m[4], m[4] + m[5]) for m in matches if m[-1] == "R1"],
                    key=lambda x: x[1],
                ),
                sorted(
                    [("", m[4], m[4] + m[5]) for m in matches if m[-1] == "R2"],
                    key=lambda x: x[1],
                ),
            )
        )

        if (
            read_overlap > args.max_switch_overlap
            or ref_overlap > args.max_switch_overlap
            or (args.paired_end_mode and mate_overlap > args.max_switch_overlap)
        ):
            stats["nreads_removed_overlap_mismatch"] += 1
            use = False

        filtered_inserts = []
        # focus now on reads with structure switch - insert - switch
        if len(switch_matches) >= 2:
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
                    logger.warn(
                        "interval error for insert-adjoining switch matches; discarding read "
                        + read
                    )
                    use = False
                    continue

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
                    and all(gap < args.max_insert_gap for gap in gaps)
                    and overlap <= args.max_insert_overlap
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
                            insert[-1],
                        )
                    )

        # take only reads matching to switch region (only once enough?)
        if len(switch_matches) < args.min_num_switch_matches:
            stats["nreads_removed_no_switch"] += 1
            use = False

        # merge the switch mappings and sort them by genomic coordinate
        read_mappings = [
            (rm[0][:-1], rm[1], rm[2], rm[0][-1])
            for rm in sorted(
                merge_intervals(
                    (
                        x[2]
                        + ("+" if (x[5] if x[-1] == "R1" else -x[5]) == main_orientation else "-"),
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
            stats["nreads_removed_length_mismatch"] += 1
            use = False

        # check that more than args.min_cov of the read maps to the genome
        assert all(
            fi[8] >= fi[7] for fi in filtered_inserts
        ), "coordinate mismatch for filtered inserts"
        tot_cov = sum(sm[1] - sm[0] for sm in switch_matches) + sum(
            fi[8] - fi[7] for fi in filtered_inserts
        )
        # adjust total read length for paired-end reads
        if args.paired_end_mode:
            tot_read_len = sum(x[0] for x in set([(m[2], m[-1]) for m in matches]))

        if tot_cov / tot_read_len < args.min_cov:
            stats["nreads_removed_low_cov"] += 1
            use = False

        # check how much genomic sequence is covered multiple times
        unmerged_switch_coverage = sum(sm[4] - sm[3] for sm in switch_matches) + sum(
            fi[8] + fi[7] for fi in filtered_inserts
        )
        merged_switch_coverage = sum(rm[2] - rm[1] for rm in read_mappings) + sum(
            fi[8] + fi[7] for fi in filtered_inserts
        )

        # check that switch matches are in consistent order along the read and the genome (but ignore inverted segments)
        read_order = np.argsort(
            [
                sm[0]
                for sm in switch_matches
                if (sm[5] if sm[-1] == "R1" else -sm[5]) == main_orientation
            ]
        )
        genomic_order = np.argsort(
            [
                sm[3]
                for sm in switch_matches
                if (sm[5] if sm[-1] == "R1" else -sm[5]) == main_orientation
            ]
        )
        if main_orientation == -1 or args.paired_end_mode:
            genomic_order = genomic_order[::-1]
        if not (read_order == genomic_order).all() and not args.ignore_order:
            stats["nreads_removed_switch_order"] += 1
            use = False

        # isotype
        isotype = ""
        i = 0
        while len(isotype) == 0 and i < len(read_mappings):
            last_map = read_mappings[-(i + 1)] if switch_orientation == "+" else read_mappings[i]
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
            i += 1

        if isotype == "":
            use = False
            stats["nreads_removed_no_isotype"] += 1

        if orientations["1"] > 0 and orientations["-1"] > 0:
            stats["nreads_inversions"] += 1

        if args.interrupt_for_read and read in args.interrupt_for_read:
            print(read, [sm[:7] for sm in matches])
            raise Exception("stop")

        if not use:
            continue

        out_string = "{0}\t{1}\t{2:.4f}\t{3:.4f}\t".format(
            isotype,
            "+" if main_orientation == 1 else "-",
            tot_cov / tot_read_len,
            1 - merged_switch_coverage / unmerged_switch_coverage,
        )

        out_string += ";".join(
            "{0}:{1}-{2}:{3}".format(chrom, start, end, orientation)
            for chrom, start, end, orientation in read_mappings
        )
        out_string += "\t"
        out_string += ";".join(
            "insert_{0}:{1}_{2}_{3}-{4}_{5}{12}_{6}:{7}-{8}_{9}_{10}:{11}".format(
                switch_chrom,
                left,
                gapl,
                istart,
                iend,
                "+" if (orientation if mate == "R1" else -orientation) == 1 else "-",
                chrom,
                start,
                end,
                gapr,
                switch_chrom,
                right,
                mate if args.paired_end_mode else "",
            )
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
                mate,
            ) in filtered_inserts
        )

        # indicate which breakpoints occur between read mates (and should not be counted)
        if args.paired_end_mode and len(set(sm[-1] for sm in switch_matches)) > 1:
            mb = [
                k
                for k in range(len(switch_matches) - 1)
                if switch_matches[k][-1] != switch_matches[k + 1][-1]
            ][0]
            coords = sorted(switch_matches[mb][3:5] + switch_matches[mb + 1][3:5])
            out_string += "\t{0};{1}".format(coords[1], coords[2])

        if args.remove_duplicates:
            umi = read.split("_")[-1]
            if out_string in processed_umis[umi]:
                stats["nreads_removed_duplicate"] += 1
                processed_umis[umi][out_string] += 1
                continue
            else:
                processed_umis[umi][out_string] = 1

        outf.write(read + "\t" + out_string + "\n")

        # write genomic alignments to fasta file for MSA construction
        if args.sequences is not None:
            sequences = (
                SeqRecord.SeqRecord(
                    Seq.Seq(aligned_seq),
                    id="{0}@{1}:{2}-{3}:{4}".format(
                        read,
                        ref_chrom,
                        ref_start,
                        ref_end,
                        "+"
                        if ((orientation if mate == "R1" else -orientation) == main_orientation)
                        else "-",
                    ),
                )
                for (
                    read_start,
                    read_end,
                    ref_chrom,
                    ref_start,
                    ref_end,
                    orientation,
                    aligned_seq,
                    mate,
                ) in switch_matches
            )
            SeqIO.write(sequences, seq_out, "fasta")

        # store coordinates for breakpoint realignment
        if args.realign_breakpoints is not None:
            processed_matches[read] = sorted(
                [sm[:6] + ("switch", sm[-1]) for sm in switch_matches]
                + [
                    (f[7], f[8], f[0], f[1], f[2], f[9], "insert", f[-1])
                    for f in filtered_inserts
                    if f[0] != "telomer"
                ],
                key=lambda x: (x[-1], x[0]),
            )

    outf.close()
    if args.sequences is not None:
        seq_out.close()

    logger.info("done. saving stats to " + args.stats)
    if args.stats:
        pd.Series(stats).to_csv(args.stats, header=False)

    if args.realign_breakpoints is not None:
        import pysam

        logger.info("loading genome from " + args.genome)
        genome_dict = {"genome": pysam.FastaFile(args.genome)}

        # pre-load switch region genome to maybe speed up look-ups
        genome_dict["switch_genome"] = genome_dict["genome"].fetch(
            switch_chrom, switch_start, switch_end
        )
        genome_dict["offset"] = switch_start

        realignments = {}
        if args.paired_end_mode:
            logger.info(
                "re-aligning breakpoints using raw reads from "
                + " and ".join(args.raw_reads.split(","))
            )
            for rec1, rec2 in zip(
                SeqIO.parse(
                    gzip.open(args.raw_reads.split(",")[0], "rt")
                    if args.raw_reads.split(",")[0].endswith(".gz")
                    else open(args.raw_reads.split(",")[0]),
                    "fastq",
                ),
                SeqIO.parse(
                    gzip.open(args.raw_reads.split(",")[1], "rt")
                    if args.raw_reads.split(",")[1].endswith(".gz")
                    else open(args.raw_reads.split(",")[1]),
                    "fastq",
                ),
            ):
                if (
                    rec1.id != rec2.id
                    or rec1.id not in processed_matches
                    or rec2.id not in processed_matches
                ):
                    continue

                realignments[rec1.id] = realign_breakpoints(
                    processed_matches[rec1.id],
                    genome_dict,
                    (str(rec1.seq), str(rec2.seq)),
                    penalties=args.realignment_penalties,
                )

        else:
            logger.info("re-aligning breakpoints using raw reads from " + args.raw_reads)
            for rec in SeqIO.parse(
                gzip.open(args.raw_reads, "rt")
                if args.raw_reads.endswith(".gz")
                else open(args.raw_reads),
                "fastq",
            ):
                if rec.id not in processed_matches:
                    continue

                realignments[rec.id] = realign_breakpoints(
                    processed_matches[rec.id],
                    genome_dict,
                    str(rec.seq),
                    penalties=args.realignment_penalties,
                )

        logger.info("saving breakpoint re-alignments to " + args.realign_breakpoints)
        df = pd.DataFrame(
            [pd.Series(al, name=read) for read, als in realignments.items() for al in als]
        )
        df.to_csv(args.realign_breakpoints)
