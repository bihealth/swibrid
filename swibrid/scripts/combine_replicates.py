"""combine (aligned + processed) reads from replicates"""


def setup_argparse(parser):
    parser.add_argument(
        "-s",
        "--samples",
        dest="samples",
        help="""comma-separated list of samples to combine""",
    )
    parser.add_argument(
        "-c",
        "--combined",
        default="combined",
        dest="combined",
        help="""combined sample name""",
    )
    parser.add_argument(
        "--nmax", dest="nmax", type=int, default=50000, help="""use only nmax reads [50000]"""
    )
    parser.add_argument(
        "--aligner", dest="aligner", default="last", help="""alignment algorithm used [last]"""
    )
    parser.add_argument(
        "--msa_path",
        dest="msa_path",
        default="pipeline/{sample}/{sample}_msa.npz",
        help="""path pattern for msa files ["pipeline/{sample}/{sample}_msa.npz"]""",
    )
    parser.add_argument(
        "--msa_csv_path",
        dest="msa_csv_path",
        default="pipeline/{sample}/{sample}_msa.csv",
        help="""path pattern for msa csv files ["pipeline/{sample}/{sample}_msa.csv"]""",
    )
    parser.add_argument(
        "--info_path",
        dest="info_path",
        default="input/{sample}_info.csv",
        help="""path pattern for info csv files ["input/{sample}_info.csv"]""",
    )
    parser.add_argument(
        "--process_path",
        dest="process_path",
        default="pipeline/{sample}/{sample}_processed.out",
        help="""path pattern for processing output ["pipeline/{sample}/{sample}_processed.out"]""",
    )
    parser.add_argument(
        "--process_stats_path",
        dest="process_stats_path",
        default="pipeline/{sample}/{sample}_process_stats.csv",
        help="""path pattern for process stats csv ["pipeline/{sample}/{sample}_process_stats.csv"]""",
    )
    parser.add_argument(
        "--alignment_pars_path",
        dest="alignment_pars_path",
        default="pipeline/{sample}/{sample}_{aligner}_pars.npz",
        help="""path pattern for alignment pars ["pipeline/{sample}/{sample}_{aligner}_pars.npz"]""",
    )
    parser.add_argument(
        "--aligned_seqs_path",
        dest="aligned_seqs_path",
        default="pipeline/{sample}/{sample}_aligned.fasta.gz",
        help="""path pattern for aligned sequences ["pipeline/{sample}/{sample}_aligned.fasta.gz""",
    )
    parser.add_argument(
        "--breakpoint_alignments_path",
        dest="breakpoint_alignments_path",
        default="pipeline/{sample}/{sample}_breakpoint_alignments.csv",
        help="""path pattern for breakpoint alignments ["pipeline/{sample}/{sample}_breakpoint_alignments.csv"]""",
    )
    parser.add_argument(
        "--inserts_tsv_path",
        dest="inserts_tsv_path",
        default="pipeline/{sample}/{sample}_inserts.tsv",
        help="""path pattern for inserts tsv ["pipeline/{sample}/{sample}_inserts.tsv"]""",
    )
    parser.add_argument(
        "--bed_path",
        dest="bed_path",
        default="pipeline/{sample}/{sample}.bed",
        help="""path pattern for bed ["pipeline/{sample}/{sample}.bed"]""",
    )


def run(args):
    import os
    import numpy as np
    import pandas as pd
    import scipy.sparse
    from Bio import SeqIO
    import gzip
    from logzero import logger

    samples = args.samples.split(",")
    assert len(samples) > 1, "please provide a comma-separated list of more than one sample"

    assert all(
        os.path.isfile(args.msa_path.format(sample=sample)) for sample in samples
    ), "please run pipeline on invidividual samples at least up to construct_msa"

    logger.info("combining reads from MSA output")
    msa = {}
    msa_df = {}
    for sample in samples:
        logger.info("reading MSA for " + sample)
        msa[sample] = scipy.sparse.load_npz(args.msa_path.format(sample=sample))
        msa_df[sample] = pd.read_csv(args.msa_csv_path.format(sample=sample), index_col=0, header=0)

    msa_df = pd.concat(msa_df.values(), axis=0, keys=msa_df.keys()).reset_index(
        level=0, names="sample"
    )

    reads = msa_df.sample(n=min(args.nmax, msa_df.shape[0])).index
    inds = msa_df.index.get_indexer(reads)

    logger.info("using {0} of {1} reads".format(len(reads), msa_df.shape[0]))

    os.makedirs(os.path.dirname(args.msa_path.format(sample=args.combined)), exist_ok=True)

    logger.info("combining reads and read info")
    info = {}
    for sample in samples:
        info[sample] = pd.read_csv("input/{0}_info.csv".format(sample), header=0, index_col=0)
    info = (
        pd.concat(info.values(), axis=0, keys=info.keys())
        .reset_index(level=0, names="sample")
        .loc[reads]
        .to_csv("input/{0}_info.csv".format(args.combined))
    )
    with gzip.open("input/{0}.fastq.gz".format(args.combined), "wt") as outf:
        for sample in samples:
            for rec in SeqIO.parse(gzip.open("input/{0}.fastq.gz".format(sample), "rt"), "fastq"):
                if rec.id in reads:
                    SeqIO.write(rec, outf, "fastq")

    logger.info("averaging alignment pars")
    pars = []
    for sample in samples:
        pars.append(np.load(args.alignment_pars_path.format(sample=sample, aligner=args.aligner)))
    combined_pars = {}
    for slot in set.intersection(*(set(p.files) for p in pars)):
        combined_pars[slot] = np.mean([p[slot] for p in pars], axis=0)
    np.savez(
        args.alignment_pars_path.format(sample=args.combined, aligner=args.aligner), **combined_pars
    )

    logger.info("combining output from process_alignments")
    process = {}
    process_stats = {}
    breakpoint_alignments = {}
    for sample in samples:
        process[sample] = pd.read_csv(
            args.process_path.format(sample=sample), header=0, index_col=0, sep="\t"
        )
        process_stats[sample] = pd.read_csv(
            args.process_stats_path.format(sample=sample), header=None, index_col=0
        ).squeeze()
        breakpoint_alignments[sample] = pd.read_csv(
            args.breakpoint_alignments_path.format(sample=sample), header=0, index_col=0
        )

    pd.concat(process.values(), axis=0, keys=process.keys()).reset_index(
        level=0, names="sample"
    ).loc[reads].to_csv(args.process_path.format(sample=args.combined), sep="\t", index=True)
    pd.concat(process_stats.values(), axis=0, keys=process_stats.keys()).groupby(level=0).agg(
        sum
    ).to_csv(args.process_stats_path.format(sample=args.combined), header=False)
    breakpoint_alignments = pd.concat(
        breakpoint_alignments.values(), axis=0, keys=breakpoint_alignments.keys()
    ).reset_index(level=0, names="sample")
    breakpoint_alignments.loc[reads.intersection(breakpoint_alignments.index)].to_csv(
        args.breakpoint_alignments_path.format(sample=args.combined)
    )

    # update time stamp on (probably empty) file with aligned sequences
    aligned_seqs = args.aligned_seqs_path.format(sample=args.combined)
    with open(aligned_seqs, "a"):
        os.utime(aligned_seqs, None)

    logger.info("combining inserts")
    inserts = {}
    bed = {}
    for sample in samples:
        try:
            inserts[sample] = pd.read_csv(
                args.inserts_tsv_path.format(sample=sample), index_col=0, header=0, sep="\t"
            )
        except (ValueError, pd.errors.EmptyDataError):
            inserts[sample] = pd.DataFrame(
                [],
                columns=[
                    "isotype",
                    "switch_coords",
                    "insert_pos",
                    "insert_coords",
                    "insert_anno",
                    "sequence",
                ],
            )
        bed[sample] = pd.read_csv(args.bed_path.format(sample=sample), header=None, sep="\t")

    inserts = pd.concat(inserts.values(), axis=0)
    insert_reads = reads.intersection(inserts.index)
    inserts = inserts.loc[insert_reads].reset_index(names="read")
    inserts.to_csv(args.inserts_tsv_path.format(sample=args.combined), index=False, sep="\t")
    bed = (
        pd.concat(bed.values(), axis=0)
        .set_index(3)
        .loc[reads]
        .reset_index()[[0, 1, 2, "read", 4, 5, 6, 7, 8, 9, 10, 11]]
    )
    bed.to_csv(args.bed_path.format(sample=args.combined), header=False, index=False, sep="\t")

    logger.info("writing combined MSA to " + args.msa_path.format(sample=args.combined))
    msa_df.loc[reads].to_csv(args.msa_csv_path.format(sample=args.combined))
    scipy.sparse.save_npz(
        args.msa_path.format(sample=args.combined),
        scipy.sparse.vstack(msa.values()).tocsr()[inds].tocoo(),
    )
