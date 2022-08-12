"""Command line interface for swibrid.
This is the main program entry point for the ``swibrid`` executable and its sub commands.
"""

import argparse
import logging
import warnings
from logzero import logger, loglevel
from pathlib import Path
from ruamel.yaml import YAML

from .utils import (
    demultiplex,
    filter_reads,
    process_last_output,
    construct_msa,
    get_gaps,
    construct_graph,
    construct_linkage,
    find_clusters,
    plot_clustering,
    find_mutations,
    get_QC_plots,
    get_stats,
    create_bed,
)

# from tool import __version__


def run_nocmd(_, parser):
    """No command given, print help and ``exit(1)``."""
    parser.print_help()
    parser.exit(1)


def run_pipeline(args, snake_options):
    """run the pipeline"""

    import os

    snakefile = Path(__file__).parent.joinpath("pipeline.snake")
    Path("logs").mkdir(exist_ok=True)
    config = YAML(typ="safe").load(open(args.config))
    s_command = """export SBATCH_DEFAULTS=" --output=logs/%x-%j.log"\nsnakemake --snakefile {sfile} --configfile {cfile} -j 100 -k --rerun-incomplete --latency-wait 60 -p --retries 10 --use-conda --conda-prefix {conda}""".format(
        sfile=snakefile, cfile=args.config, conda=config["CONDA"]
    )
    s_command += " " + " ".join(snake_options)
    if args.slurm:
        s_command = """#!/bin/bash\n""" + s_command + " --profile=cubi-v1\n"
        logger.info("run script content:\n" + s_command)
        run_script = Path("run_pipeline.sh")
        run_script.write_text(s_command)
        command = "sbatch -t 168:00:00 --mem=4G -n 1 -p medium run_pipeline.sh"
    else:
        command = s_command
    # run
    logger.info("executing\n" + command)
    os.system(command)


def main(argv=None):
    """Main entry point before parsing command line arguments."""

    description = """
    main command:

    swibrid run                    run the entire pipeline

    subcommands: some steps can be run individually

    swibrid demultiplex            demultiplex minION run
    swibrid filter_reads           filter reads
    swibrid process_last_output    process output of LAST
    swibrid create_bed             create bed files and summary tables for insert-containing reads
    swibrid construct_msa          construct a (pseudo) MSA
    swibrid get_gaps               find gaps in MSA
    swibrid construct_graph        construct a read graph
    swibrid construct_linkage      construct hierarchical clustering on graph
    swibrid find_clusters          get read clustering from linkage
    swibrid find_mutations         call germline and somatic mutations
    swibrid plot_clustering        plot the read clustering
    swibrid get_QC_plots           get plots for various QC stats
    swibrid get_stats              get summary stats

    optional arguments:
     -h, --help            show this help message and exit
     --version             show version information
     --verbose             Increase verbosity
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Increase verbosity.",
    )
    # parser.add_argument("--version", action="version", version="%%(prog)s %s" % __version__)

    subparsers = parser.add_subparsers(
        dest="cmd",
        title="commands",
        description=description,
        help=argparse.SUPPRESS,
    )

    pipeline_parser = subparsers.add_parser(
        "run", help="main command: run pipeline"
    )
    pipeline_parser.add_argument(
        "-c",
        "--config",
        dest="config",
        default="config.yaml",
        help="config file [config.yaml]",
    )
    pipeline_parser.add_argument(
        "-s",
        "--slurm",
        dest="slurm",
        action="store_true",
        help="submit to slurm",
    )
    pipeline_parser.add_argument(
        "snake_options",
        nargs=argparse.REMAINDER,
        help="pass options to snakemake (...)",
    )

    demultiplex.setup_argparse(
        subparsers.add_parser("demultiplex", help="demultiplex minION run")
    )
    filter_reads.setup_argparse(
        subparsers.add_parser("filter_reads", help="filter reads")
    )
    process_last_output.setup_argparse(
        subparsers.add_parser("process_last_output", help=argparse.SUPPRESS)
    )
    create_bed.setup_argparse(
        subparsers.add_parser("create_bed", help="step: create bed file")
    )
    construct_msa.setup_argparse(
        subparsers.add_parser(
            "construct_msa", help="step: construct (pseudo) MSA"
        )
    )
    get_gaps.setup_argparse(
        subparsers.add_parser("get_gaps", help="step: find gaps in MSA")
    )
    construct_graph.setup_argparse(
        subparsers.add_parser(
            "construct_graph", help="step: construct read graph"
        )
    )
    construct_linkage.setup_argparse(
        subparsers.add_parser(
            "construct_linkage", help="step: construct linkage"
        )
    )
    find_clusters.setup_argparse(
        subparsers.add_parser("find_clusters", help="step: find clusters")
    )
    find_mutations.setup_argparse(
        subparsers.add_parser("find_mutations", help="step: find mutations")
    )
    plot_clustering.setup_argparse(
        subparsers.add_parser("plot_clustering", help="step: plot clustering")
    )
    get_QC_plots.setup_argparse(
        subparsers.add_parser("get_QC_plots", help="step: plot QC results")
    )
    get_stats.setup_argparse(
        subparsers.add_parser("get_stats", help="step: get stats")
    )

    args, extra = parser.parse_known_args(argv)

    # Setup logging verbosity.
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    loglevel(level=level)

    steps = {
        "demultiplex": demultiplex.run,
        "filter_reads": filter_reads.run,
        "process_last_output": process_last_output.run,
        "create_bed": create_bed.run,
        "construct_msa": construct_msa.run,
        "get_gaps": get_gaps.run,
        "construct_graph": construct_graph.run,
        "construct_linkage": construct_linkage.run,
        "find_clusters": find_clusters.run,
        "find_mutations": find_mutations.run,
        "plot_clustering": plot_clustering.run,
        "get_QC_plots": get_QC_plots.run,
        "get_stats": get_stats.run,
    }

    if args.cmd == "run":
        return run_pipeline(args, extra)
    elif args.cmd in steps.keys():
        return steps[args.cmd](args)
    else:
        return run_nocmd
