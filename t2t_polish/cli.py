#!/usr/bin/env python3
"""
T2T-Polish v4.0 — cli.py

Command-line interface: argparse setup and ``main()`` entry-point.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys

from t2t_polish import __version__
from t2t_polish.constants import DEFAULT_KMER_SIZE
from t2t_polish.kcov import sub_computekcov
from t2t_polish.polish import sub_polish
from t2t_polish.utils import (
    diagnostic_mode,
    log_tool_versions,
    setup_logging,
    validate_dependencies,
)

logger = logging.getLogger(__name__)


class ShowHelpOnErrorArgumentParser(argparse.ArgumentParser):
    """Prints full help (not just usage) when the user makes a mistake."""

    def error(self, message: str) -> None:  # type: ignore[override]
        self.print_help(sys.stderr)
        sys.stderr.write(f"\nerror: {message}\n")
        sys.exit(2)


def _build_parser() -> ShowHelpOnErrorArgumentParser:
    """Construct the argument parser — separated out for testability."""
    parser = ShowHelpOnErrorArgumentParser(
        description=(
            "Integrated T2T Polishing Pipeline with K-mer Coverage Optimization. "
            "Subcommands: polish, computekcov."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"T2T-Polish v{__version__}"
    )

    subparsers = parser.add_subparsers(dest="subcommand")

    # ── computekcov ──────────────────────────────────────────────────
    parser_ck = subparsers.add_parser(
        "computekcov",
        help="Compute optimal K-coverage and a fitted K-mer histogram using Jellyfish + GenomeScope2.",
    )
    parser_ck.add_argument("-r", "--reads", required=True, help="Path to reads (FASTA/FASTQ).")
    parser_ck.add_argument(
        "-o", "--prefix", default="KCoverage",
        help="Prefix for Jellyfish/GenomeScope2 outputs (default: KCoverage)",
    )
    parser_ck.add_argument(
        "-k", "--kmer_size", type=int, default=21,
        help="K-mer size (default: 21, recommended by Jellyfish / GenomeScope)",
    )
    parser_ck.add_argument("-t", "--threads", type=int, default=8, help="Number of threads (default: 8)")
    parser_ck.add_argument(
        "--ploidy", choices=["haploid", "diploid"], default="haploid",
        help="Genome ploidy for GenomeScope2. (default: haploid => ploidy = 1)",
    )
    parser_ck.add_argument(
        "--jellyfish-hash-size", type=int, default=1_000_000_000,
        help="Hash size (number of k-mers) for Jellyfish counting (default: 1e9)",
    )
    parser_ck.set_defaults(func=sub_computekcov)

    # ── polish ───────────────────────────────────────────────────────
    parser_pol = subparsers.add_parser(
        "polish",
        help=(
            "Run iterative polishing with DeepVariant-based variant calling, Merfin, etc. "
            "Use --optimized to automatically run computekcov and input the calculated values automatically."
        ),
    )
    # Required
    parser_pol.add_argument("-d", "--draft", help="Draft assembly FASTA")
    parser_pol.add_argument(
        "-r", "--reads",
        help="Reads in FASTQ or FASTA format, GZIP is currently not accepted",
    )
    parser_pol.add_argument(
        "--singularity_sif", required=True,
        help="Path to GPU-accelerated DeepVariant Singularity image.",
    )
    parser_pol.add_argument(
        "--deepseq_type",
        choices=["WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA"],
        default="PACBIO",
        help="DeepVariant model type for read input. (default: PACBIO)",
    )

    # Optional
    parser_pol.add_argument(
        "-o", "--prefix", default="AutoPolisher", help="Output prefix (default: AutoPolisher)"
    )
    parser_pol.add_argument(
        "-t", "--threads", type=int, default=32,
        help="Number of threads. Equal to number of shards in DeepVariant. (default: 32)",
    )
    parser_pol.add_argument(
        "--optimized", action="store_true",
        help="Enable automatic coverage/hist calculation (Merfin -peak / -prob).",
    )
    parser_pol.add_argument(
        "--resume", action="store_true",
        help=(
            "Resume from a partially completed run. If set, computationally expensive steps "
            "can be skipped if their output already exists."
        ),
    )
    parser_pol.add_argument(
        "--resume-from", type=int, default=0,
        help=(
            "Resume from a specific polishing step "
            "(0=all, 1=Meryl DB, 2=Winnowmap, 3=FalconC, 4=DeepVariant, 5=Merfin, 6=Consensus)"
        ),
    )
    parser_pol.add_argument(
        "-s", "--seq_type",
        help="Optionally set winnowmap read type (pb, ont). If not set, we may guess from --deepseq_type.",
    )
    parser_pol.add_argument(
        "--dv_tmp_dir", default=None,
        help="Directory to use for DeepVariant intermediate results (default: <iteration_folder>/dv_tmp)",
    )
    parser_pol.add_argument(
        "--read_corrector", choices=["hifiasm", "herro", "flye", "unknown"],
        help="Specify the read corrector used to assign default quality scores for FASTA input.",
    )
    parser_pol.add_argument(
        "-m", "--readmers",
        help="Meryl DB of read k-mers. Calculated automatically with --optimized.",
    )
    parser_pol.add_argument(
        "-k", "--kmer_size", type=int, default=DEFAULT_KMER_SIZE,
        help="K-mer size for Meryl. (default: 31, recommended by developers)",
    )
    parser_pol.add_argument(
        "-i", "--iterations", type=int, default=3,
        help="Number of polishing iterations (default: 3)",
    )
    parser_pol.add_argument(
        "--meryl-memory", type=int, default=50,
        help="Maximum memory (in GB) to hand off to Meryl (default: 50)",
    )
    parser_pol.add_argument(
        "--cleanup", action="store_true",
        help="If set, intermediate files (SAMs, uncompressed VCFs) are deleted after use to save disk space.",
    )
    parser_pol.add_argument(
        "--log-file",
        help="Path to write detailed log output (default: none—only console)",
    )
    parser_pol.add_argument(
        "--quiet", action="store_true",
        help="Suppress INFO-level logs to console; only WARNING+ERROR will print",
    )

    # Peak settings
    parser_pol.add_argument(
        "--ideal_dpeak", type=float,
        help="Coverage peak to use for haploid assemblies (was --peak_haploid).",
    )
    parser_pol.add_argument(
        "--ideal_hpeak", type=float,
        help="Coverage peak to use for diploid assemblies (was --peak_diploid).",
    )
    parser_pol.add_argument(
        "--fitted_hist", help="Fitted histogram location if not run as optimized"
    )
    parser_pol.add_argument(
        "--ploidy", choices=["haploid", "diploid"], default="haploid",
        help="Assembly ploidy (haploid or diploid) (default: haploid)",
    )

    # Config + JSON summary
    parser_pol.add_argument(
        "--config_file",
        help="Optional JSON config file with overrides to defaults, instead of command line entries.",
    )
    parser_pol.add_argument(
        "--json-summary", action="store_true",
        help="Emit a machine-readable JSON summary of each iteration's QV, completeness, and runtime.",
    )
    parser_pol.set_defaults(func=sub_polish)

    # ── diagnostics ──────────────────────────────────────────────────
    parser_diag = subparsers.add_parser(
        "diagnostics", help="Run a tool and environment diagnostic check"
    )
    parser_diag.set_defaults(func=lambda args: diagnostic_mode())

    return parser


def main() -> None:
    """Entry-point for the T2T-Polish CLI."""
    validate_dependencies()

    parser = _build_parser()
    args = parser.parse_args()

    log_file = getattr(args, "log_file", None)
    quiet = getattr(args, "quiet", False)
    setup_logging(log_file, quiet)

    # Re-bind module-level logger after root handlers are configured
    global logger
    logger = logging.getLogger(__name__)
    logger.info("Pipeline parameters: %s", vars(args))

    if hasattr(args, "threads"):
        if args.threads < 1:
            logger.error("Number of threads must be >= 1.")
            sys.exit(1)
        elif args.threads > 128:
            logger.warning(
                "You requested more than 128 threads. This may exceed available system resources."
            )

    if not args.subcommand:
        parser.print_help()
        sys.exit(1)

    # Log run parameters
    prefix = getattr(args, "prefix", None)
    if prefix:
        try:
            param_log_path = prefix + ".run_parameters.json"
            params = vars(args).copy()
            params.pop("func", None)
            with open(param_log_path, "w") as f:
                json.dump(params, f, indent=4)
            logger.info("Saved run parameters to %s", param_log_path)
        except Exception as e:
            logger.warning("Could not save run parameters: %s", e)

        log_tool_versions(prefix)

    # Validate critical inputs
    draft = getattr(args, "draft", None)
    reads = getattr(args, "reads", None)
    for path in [draft, reads]:
        if path and not os.path.exists(path):
            logger.error("Required input file not found: %s", path)
            sys.exit(1)
        elif path and os.path.getsize(path) == 0:
            logger.error("Input file %s exists but is empty.", path)
            sys.exit(1)

    # Dispatch
    try:
        args.func(args)
    except Exception:
        logger.exception("Pipeline failed with unhandled exception")
        sys.exit(1)


if __name__ == "__main__":
    main()

# T2T-Polish v4.0
# Any usage is subject to this software's license.
