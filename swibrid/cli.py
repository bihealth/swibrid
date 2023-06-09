"""Command line interface for swibrid.
This is the main program entry point for the ``swibrid`` executable and its sub commands.
"""

import argparse
import logging
import warnings
from logzero import logger, loglevel
from pathlib import Path
from ruamel.yaml import YAML

from .scripts import (
    demultiplex,
    filter_reads,
    process_alignments,
    construct_msa,
    get_gaps,
    construct_linkage,
    find_clusters,
    plot_clustering,
    find_variants,
    get_summary,
    get_gap_stats,
    get_switch_homology,
    get_switch_motifs,
    create_bed,
    get_synthetic_reads,
    get_alignment_pars,
    get_unique_clones_bed,
    analyze_clustering,
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
    if args.slurm:
        s_command = """#!/bin/bash\n""" + s_command + " --profile=cubi-dev\n"
        logger.info("run script content:\n" + s_command)
        run_script = Path("run_pipeline.sh")
        run_script.write_text(s_command)
        command = "sbatch -t 168:00:00 --mem=4G -n 1 -p {queue} run_pipeline.sh".format(
            queue=args.queue
        )
    else:
        command = s_command
    # run
    logger.info("executing\n" + command)
    os.system(command)


def run_setup(args):
    """setup new directory"""

    import os

    snakefile = Path(__file__).parent.joinpath("pipeline.snake")
    configfile = Path(__file__).parent.joinpath("config.yaml")

    if not os.path.isfile("pipeline.snake") or args.overwrite:
        logger.info("copying pipeline.snake")
        os.system("cp -i {0} .".format(snakefile))
    else:
        logger.warn("refusing to overwrite pipeline.snake")
    if not os.path.isfile("config.yaml") or args.overwrite:
        logger.info("copying config.yaml")
        os.system("cp -i {0} .".format(configfile))
    else:
        logger.warn("refusing to overwrite config.yaml")


def main(argv=None):
    """Main entry point before parsing command line arguments."""

    description = """
    main commands:

    swibrid setup                  copy config and snakefile to current directory
    swibrid run                    run the entire pipeline

    subcommands to run pipeline steps individually:

    swibrid demultiplex            demultiplex minION run
    swibrid filter_reads           filter reads
    swibrid process_alignments     process alignment output
    swibrid get_alignment_pars     get alignment parameters (mismatch + gap rates) from align output
    swibrid create_bed             create bed files and summary tables for insert-containing reads
    swibrid construct_msa          construct a (pseudo) MSA
    swibrid get_gaps               find gaps in MSA
    swibrid construct_linkage      construct hierarchical clustering linkage
    swibrid find_clusters          get read clustering from linkage
    swibrid find_variants          variant calling
    swibrid plot_clustering        plot the read clustering
    swibrid get_summary            get sample summary stats and plot
    swibrid get_gap_stats          get breakpoint histogram stats
    swibrid get_switch_homology    get switch sequence homology
    swibrid get_switch_motifs      get switch sequence motifs
    swibrid analyze_clustering     analyze clustering results
    swibrid collect results        collect results for multiple samples

    additional subcommands to create synthetic reads for testing or benchmarking:

    swibrid get_unique_clones_bed  get bed file with unique clones
    swibrid get_synthetic_reads    create synthetic reads from bed file
    """

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
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

    setup_parser = subparsers.add_parser("setup", help="main command: setup current directory")
    setup_parser.add_argument(
        "-f",
        "--overwrite",
        dest="overwrite",
        help="""overwrite files if present [no]""",
    )

    pipeline_parser = subparsers.add_parser("run", help="main command: run pipeline")
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

    demultiplex.setup_argparse(subparsers.add_parser("demultiplex", help="demultiplex minION run"))
    filter_reads.setup_argparse(subparsers.add_parser("filter_reads", help="filter reads"))
    process_alignments.setup_argparse(
        subparsers.add_parser("process_alignments", help=argparse.SUPPRESS)
    )
    get_alignment_pars.setup_argparse(
        subparsers.add_parser("get_alignment_pars", help=argparse.SUPPRESS)
    )
    create_bed.setup_argparse(subparsers.add_parser("create_bed", help="step: create bed file"))
    construct_msa.setup_argparse(
        subparsers.add_parser("construct_msa", help="step: construct (pseudo) MSA")
    )
    get_gaps.setup_argparse(subparsers.add_parser("get_gaps", help="step: find gaps in MSA"))
    construct_linkage.setup_argparse(
        subparsers.add_parser("construct_linkage", help="step: construct linkage")
    )
    find_clusters.setup_argparse(subparsers.add_parser("find_clusters", help="step: find clusters"))
    find_variants.setup_argparse(subparsers.add_parser("find_variants", help="step: find variants"))
    plot_clustering.setup_argparse(
        subparsers.add_parser("plot_clustering", help="step: plot clustering")
    )
    get_summary.setup_argparse(
        subparsers.add_parser("get_summary", help="step: get sample summary and plot")
    )
    get_gap_stats.setup_argparse(
        subparsers.add_parser("get_gap_stats", help="step: get gap statistics")
    )
    get_switch_homology.setup_argparse(
        subparsers.add_parser("get_switch_homology", help="step: get switch sequence homology")
    )
    get_switch_motifs.setup_argparse(
        subparsers.add_parser("get_switch_motifs", help="step: get switch sequence motifs")
    )
    analyze_clustering.setup_argparse(
        subparsers.add_parser("analyze_clustering", help="step: analyze clustering")
    )
    get_synthetic_reads.setup_argparse(
        subparsers.add_parser("get_synthetic_reads", help="step: create synthetic reads")
    )
    get_unique_clones_bed.setup_argparse(
        subparsers.add_parser("get_unique_clones_bed", help="step: get unique clones from bed")
    )
    collect_results.setup_argparse(
        subparsers.add_parser("collect_results", help="step: collect results")
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
        "filter_reads": filter_reads.run,
        "process_alignments": process_alignments.run,
        "get_alignment_pars": get_alignment_pars.run,
        "create_bed": create_bed.run,
        "construct_msa": construct_msa.run,
        "get_gaps": get_gaps.run,
        "construct_linkage": construct_linkage.run,
        "find_clusters": find_clusters.run,
        "find_variants": find_variants.run,
        "plot_clustering": plot_clustering.run,
        "get_summary": get_summary.run,
        "get_gap_stats": get_gap_stats.run,
        "get_switch_homology": get_switch_homology.run,
        "get_switch_motifs": get_switch_motifs.run,
        "analyze_clustering": analyze_clustering.run,
        "get_synthetic_reads": get_synthetic_reads.run,
        "get_unique_clones_bed": get_unique_clones_bed.run,
        "collect_results": collect_results.run,
    }

    if args.cmd == "run":
        return run_pipeline(args, extra)
    elif args.cmd in steps.keys():
        return steps[args.cmd](args)
    else:
        return run_nocmd(parser)
