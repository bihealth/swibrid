"""Command line interface for swibrid.
This is the main program entry point for the ``swibrid`` executable and its sub commands.
"""

import argparse
import logging
from textwrap import dedent
from logzero import logger, loglevel
from pathlib import Path
from ruamel.yaml import YAML

from .scripts import (
    demultiplex,
    plot_demux_report,
    process_alignments,
    construct_msa,
    get_gaps,
    construct_linkage,
    find_clusters,
    plot_clustering,
    find_variants,
    find_rearrangements,
    get_summary,
    get_breakpoint_stats,
    get_switch_homology,
    get_switch_motifs,
    create_bed,
    get_synthetic_reads,
    get_alignment_pars,
    get_unique_clones_bed,
    combine_replicates,
    downsample,
    analyze_clustering,
    downsample_clustering,
    collect_results,
)

from swibrid import __version__


def run_nocmd(parser):
    """No command given, print help and ``exit(1)``."""
    parser.print_help()
    parser.exit(1)


def run_pipeline(args, snake_options):
    """run the pipeline"""

    import os

    if not args.snakefile:
        snakefile = Path(__file__).parent.joinpath("pipeline.snake")
    else:
        snakefile = args.snakefile

    Path("logs").mkdir(exist_ok=True)
    config = YAML(typ="safe").load(open(args.config))
    s_command = """snakemake --snakefile {sfile} --configfile {cfile} {snakemake_options}""".format(
        sfile=snakefile, cfile=args.config, snakemake_options=config["SNAKEOPTS"]
    )
    s_command += " " + " ".join(snake_options)

    if args.slurm and args.sge:
        logger.error("cannot submit to slurm and SGE at the same time!")

    if args.slurm:
        s_command = """#!/bin/bash\n""" + s_command + " --profile=cubi-dev\n"
        logger.info("run script content:\n" + s_command)
        run_script = Path("run_pipeline.sh")
        run_script.write_text(s_command)
        command = "sbatch -t 48:00:00 --mem=4G -n 1 -p {queue} run_pipeline.sh".format(
            queue=args.queue
        )
    elif args.sge:
        # drmaa_string = "\" --mem={resources.mem_mb} --time={resources.time} -n {threads}\""
        s_command = """#!/bin/bash\n""" + s_command  # + " --drmaa=" + drmaa_string + "\n"
        logger.info("run script content:\n" + s_command)
        run_script = Path("run_pipeline.sh")
        run_script.write_text(s_command)
        command = "qsub -cwd -V -pe make 16 -l m_mem_free=32G,h_rt=48:00:00 -j y -o swibrid.log run_pipeline.sh"
    else:
        command = s_command
    # run
    logger.info("executing\n" + command)
    os.system(command)


def test_pipeline(args, snake_options):
    """test the pipeline"""

    import os

    if not os.path.isfile("pipeline.snake") or args.overwrite:
        logger.info("copying pipeline.snake")
        snakefile = Path(__file__).parent.joinpath("pipeline.snake")
        os.system("cp {0} .".format(snakefile))

    if not os.path.isfile("test_config.yaml") or args.overwrite:
        logger.info("copying test_config.yaml")
        configfile = Path(__file__).parent.parent.joinpath("test_data/test_config.yaml")
        os.system("cp {0} .".format(configfile))
    
    Path("index").mkdir(exist_ok=True)
    if not os.path.isfile("index/hs1_chr14_99-110MB.fa") or args.overwrite:
        logger.info("setting up index for switch region of hs1")
        reference = Path(__file__).parent.parent.joinpath("test_data/hs1_chr14_99-110MB.fa")
        os.system("cp {0} index".format(reference))
        os.system("samtools faidx index/hs1_chr14_99-110MB.fa")
        os.system("lastdb index/hs1_chr14_99-110MBdb index/hs1_chr14_99-110MB.fa")

    switch_regions = Path(__file__).parent.parent.joinpath("test_data/hs1_chr14_99-110MB_switch_regions.bed")
    if not os.path.isfile("index/hs1_chr14_99-110MB_switch_regions.bed") or args.overwrite:
        os.system("cp {0} index".format(switch_regions))

    Path("input").mkdir(exist_ok=True)
    if not os.path.isfile("input/input_clones.bed") or args.overwrite:
        logger.info("copying input clones")
        clones = Path(__file__).parent.parent.joinpath("test_data/input_clones.bed")
        os.system("cp {0} input".format(clones))
    if not os.path.isfile("input/input_variants.txt") or args.overwrite:
        logger.info("copying input variants")
        variants = Path(__file__).parent.parent.joinpath("test_data/input_variants.txt")
        os.system("cp {0} input".format(variants))
    if not os.path.isfile("input/input_pars.par") or args.overwrite:
        logger.info("copying input parameters")
        pars = Path(__file__).parent.parent.joinpath("test_data/input_pars.par")
        os.system("cp {0} input".format(pars))

    run_pipeline(args, snake_options)


def run_setup(args):
    """setup new directory"""

    import os

    snakefile = Path(__file__).parent.joinpath("pipeline.snake")
    configfile = Path(__file__).parent.joinpath("config.yaml")

    if not os.path.isfile("pipeline.snake") or args.overwrite:
        logger.info("copying pipeline.snake")
        os.system("cp {0} .".format(snakefile))
    else:
        logger.warn("refusing to overwrite pipeline.snake")
    if not os.path.isfile("config.yaml") or args.overwrite:
        logger.info("copying config.yaml")
        os.system("cp {0} .".format(configfile))
    else:
        logger.warn("refusing to overwrite config.yaml")


def main(argv=None):
    """Main entry point before parsing command line arguments."""

    description = """
    main commands:

    swibrid setup                  copy config and snakefile to current directory
    swibrid run                    run the entire pipeline
    swibrid test                   perform a test run on synthetic reads

    subcommands to run pipeline steps individually:

    swibrid demultiplex            demultiplex minION run
    swibrid process_alignments     process and filter alignment output
    swibrid get_alignment_pars     get alignment parameters (mismatch + gap rates) from align output
    swibrid create_bed             create bed files and summary tables for insert-containing reads
    swibrid construct_msa          construct a (pseudo) MSA
    swibrid get_gaps               find gaps / breakpoints in MSA
    swibrid construct_linkage      construct hierarchical clustering linkage
    swibrid find_clusters          get read clustering from linkage
    swibrid find_variants          variant calling
    swibrid find_rearrangements    detect rearrangements (structural variants)
    swibrid plot_clustering        plot the read clustering
    swibrid get_breakpoint_stats   get breakpoint histograms and stats
    swibrid get_switch_homology    get switch sequence homology
    swibrid get_switch_motifs      get switch sequence motifs
    swibrid analyze_clustering     analyze clustering results
    swibrid downsample_clustering  downsample and analyze clustering
    swibrid get_summary            get sample summary stats and plot
    swibrid collect_results        collect results for multiple samples

    additional subcommands:

    swibrid get_unique_clones_bed  get bed file with unique clones
    swibrid get_synthetic_reads    create synthetic reads from bed file
    swibrid plot_demux_report      make a graphical summary of demultiplexing output
    swibrid combine_replicate      combine (aligned & processed) reads from replicates
    swibrid downsample             downsample (aligned & processed) reads from a sample
    """

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=dedent(__doc__))
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Increase verbosity.",
    )
    parser.add_argument("--version", action="version", version="%%(prog)s %s" % __version__)

    subparsers = parser.add_subparsers(
        dest="cmd",
        title="commands",
        description=description,
        help=argparse.SUPPRESS,
    )

    setup_parser = subparsers.add_parser(
        "setup", description="main command: copy config and snakefile to current directory"
    )
    setup_parser.add_argument(
        "-f",
        "--overwrite",
        dest="overwrite",
        action="store_true",
        default=False,
        help="""overwrite files if present [no]""",
    )

    pipeline_parser = subparsers.add_parser("run", formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
        main command: run the entire pipeline
        run `swibrid setup` and edit config.yaml first!
        
        this will call snakemake; additional options are passed to snakemake
        e.g., for a dry run use `swibrid run -n`, or `swibrid run --unlock` after a failed run
        """))
    pipeline_parser.add_argument(
        "-s",
        "--snakefile",
        dest="snakefile",
        default="pipeline.snake",
        help="snakefile [pipeline.snake]",
    )
    pipeline_parser.add_argument(
        "-c",
        "--config",
        dest="config",
        default="config.yaml",
        help="config file [config.yaml]",
    )
    pipeline_parser.add_argument(
        "--slurm",
        dest="slurm",
        action="store_true",
        help="submit to slurm",
    )
    pipeline_parser.add_argument(
        "--sge",
        dest="sge",
        action="store_true",
        help="submit to SGE",
    )
    pipeline_parser.add_argument(
        "--queue",
        dest="queue",
        default="medium",
        help="slurm queue [medium]",
    )
    pipeline_parser.add_argument(
        "snake_options",
        nargs=argparse.REMAINDER,
        help="pass options to snakemake (...)",
    )

    test_parser = subparsers.add_parser("test", formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
        test the pipline
        this will create synthetic reads in `input` and run the pipeline on this data,
        using a reduced hs1 genome in `index` with only the switch region (chr14:99000000-110000000)
        this will probably take about 30-60 minutes and will call snakemake, passing on additional options
        e.g., for a dry run use `swibrid test -n`, and do `swibrid test --unlock` after a failed run
        """))
    test_parser.add_argument(
        "-f",
        "--overwrite",
        dest="overwrite",
        action="store_true",
        default=False,
        help="""overwrite files if present [no]""",
    )
    test_parser.add_argument(
        "-s",
        "--snakefile",
        dest="snakefile",
        default="pipeline.snake",
        help="snakefile [pipeline.snake]",
    )
    test_parser.add_argument(
        "-c",
        "--config",
        dest="config",
        default="test_config.yaml",
        help="config file [test_config.yaml]",
    )
    test_parser.add_argument(
        "--slurm",
        dest="slurm",
        action="store_true",
        help="submit to slurm",
    )
    test_parser.add_argument(
        "--sge",
        dest="sge",
        action="store_true",
        help="submit to SGE",
    )
    test_parser.add_argument(
        "--queue",
        dest="queue",
        default="medium",
        help="slurm queue [medium]",
    )
    test_parser.add_argument(
        "snake_options",
        nargs=argparse.REMAINDER,
        help="pass options to snakemake (...)",
    )

    demultiplex.setup_argparse(subparsers.add_parser("demultiplex", formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent(demultiplex.__doc__)))
    plot_demux_report.setup_argparse(
        subparsers.add_parser("plot_demux_report", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(plot_demux_report.__doc__))
    )
    process_alignments.setup_argparse(
        subparsers.add_parser("process_alignments", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(process_alignments.__doc__))
    )
    get_alignment_pars.setup_argparse(
        subparsers.add_parser("get_alignment_pars", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(get_alignment_pars.__doc__))
    )
    create_bed.setup_argparse(subparsers.add_parser("create_bed", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(create_bed.__doc__)))
    construct_msa.setup_argparse(
        subparsers.add_parser("construct_msa", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(construct_msa.__doc__))
    )
    get_gaps.setup_argparse(subparsers.add_parser("get_gaps", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(get_gaps.__doc__)))
    construct_linkage.setup_argparse(
        subparsers.add_parser("construct_linkage", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(construct_linkage.__doc__))
    )
    find_clusters.setup_argparse(
        subparsers.add_parser("find_clusters", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(find_clusters.__doc__))
    )
    find_variants.setup_argparse(
        subparsers.add_parser("find_variants", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(find_variants.__doc__))
    )
    find_rearrangements.setup_argparse(
        subparsers.add_parser("find_rearrangements", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(find_rearrangements.__doc__))
    )
    plot_clustering.setup_argparse(
        subparsers.add_parser("plot_clustering", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(plot_clustering.__doc__))
    )
    get_summary.setup_argparse(
        subparsers.add_parser("get_summary", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(get_summary.__doc__))
    )
    get_breakpoint_stats.setup_argparse(
        subparsers.add_parser("get_breakpoint_stats", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(get_breakpoint_stats.__doc__))
    )
    get_switch_homology.setup_argparse(
        subparsers.add_parser(
            "get_switch_homology", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(get_switch_homology.__doc__)
        )
    )
    get_switch_motifs.setup_argparse(
        subparsers.add_parser("get_switch_motifs", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(get_switch_motifs.__doc__))
    )
    analyze_clustering.setup_argparse(
        subparsers.add_parser("analyze_clustering", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(analyze_clustering.__doc__))
    )
    downsample_clustering.setup_argparse(
        subparsers.add_parser("downsample_clustering", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(downsample_clustering.__doc__))
    )
    get_synthetic_reads.setup_argparse(
        subparsers.add_parser("get_synthetic_reads", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(get_synthetic_reads.__doc__))
    )
    get_unique_clones_bed.setup_argparse(
        subparsers.add_parser("get_unique_clones_bed", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(get_unique_clones_bed.__doc__))
    )
    combine_replicates.setup_argparse(
        subparsers.add_parser("combine_replicates", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(combine_replicates.__doc__))
    )
    downsample.setup_argparse(
        subparsers.add_parser("downsample", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(downsample.__doc__))
    )
    collect_results.setup_argparse(
        subparsers.add_parser("collect_results", formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent(collect_results.__doc__))
    )

    args, extra = parser.parse_known_args(argv)

    # Setup logging verbosity.
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    loglevel(level=level)

    steps = {
        "setup": run_setup,
        "demultiplex": demultiplex.run,
        "plot_demux_report": plot_demux_report.run,
        "process_alignments": process_alignments.run,
        "get_alignment_pars": get_alignment_pars.run,
        "create_bed": create_bed.run,
        "construct_msa": construct_msa.run,
        "get_gaps": get_gaps.run,
        "construct_linkage": construct_linkage.run,
        "find_clusters": find_clusters.run,
        "find_variants": find_variants.run,
        "find_rearrangements": find_rearrangements.run,
        "plot_clustering": plot_clustering.run,
        "get_summary": get_summary.run,
        "get_breakpoint_stats": get_breakpoint_stats.run,
        "get_switch_homology": get_switch_homology.run,
        "get_switch_motifs": get_switch_motifs.run,
        "analyze_clustering": analyze_clustering.run,
        "downsample_clustering": downsample_clustering.run,
        "get_synthetic_reads": get_synthetic_reads.run,
        "get_unique_clones_bed": get_unique_clones_bed.run,
        "combine_replicates": combine_replicates.run,
        "downsample": downsample.run,
        "collect_results": collect_results.run,
    }

    if args.cmd == "run":
        return run_pipeline(args, extra)
    elif args.cmd == "test":
        return test_pipeline(args, extra)
    elif args.cmd in steps.keys():
        return steps[args.cmd](args)
    else:
        return run_nocmd(parser)
