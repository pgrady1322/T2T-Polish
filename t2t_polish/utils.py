#!/usr/bin/env python3
"""
T2T-Polish v4.0 — utils.py

Utility helpers: path building, dependency validation, logging setup,
tool-version logging, diagnostics, FASTA→FASTQ conversion, and readmers.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
from subprocess import PIPE

from t2t_polish.constants import (
    MERYL,
    TOOL_NAMES,
)
from t2t_polish.runner import conditional_run

logger = logging.getLogger(__name__)


# ── Path helpers ─────────────────────────────────────────────────────


def make_path(base: str, *parts: str, ext: str = "") -> str:
    """Build an output file path from *base*, dot-joined *parts*, and *ext*."""
    return f"{base}." + "_".join(parts) + ext


# ── Dependency validation ────────────────────────────────────────────


def validate_dependencies() -> None:
    """Check that every required external tool is on ``$PATH``.

    Exits with code 1 and a helpful message if any tool is missing.
    """
    missing: list[str] = []
    for tool in TOOL_NAMES:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        logger.error(
            "Required tool(s) not found in PATH: %s\n"
            "  Hint: activate the conda environment (`conda activate t2t_polish`)\n"
            "  or run inside the Docker / Singularity container.",
            ", ".join(missing),
        )
        sys.exit(1)


def needs_kcov(args: object) -> bool:
    """Return *True* when the polishing run requires automatic k-mer coverage."""
    return (
        getattr(args, "optimized", False)
        and (
            not getattr(args, "fitted_hist", None)
            or (
                getattr(args, "ploidy", "") == "haploid"
                and not getattr(args, "ideal_dpeak", None)
            )
            or (
                getattr(args, "ploidy", "") == "diploid"
                and not getattr(args, "ideal_hpeak", None)
            )
        )
    )


# ── Tool version logging ────────────────────────────────────────────


def log_tool_versions(prefix: str) -> str:
    """Write all external tool versions to ``<prefix>.tool_versions.txt``."""
    version_log = prefix + ".tool_versions.txt"
    with open(version_log, "w") as f:
        for tool in TOOL_NAMES:
            f.write(f"{tool} version:\n")
            try:
                out = subprocess.run(
                    [tool, "--version"],
                    stdout=PIPE,
                    stderr=PIPE,
                    text=True,
                    timeout=10,
                )
                f.write(out.stdout.strip() + "\n\n")
            except subprocess.TimeoutExpired:
                f.write(f"Error: {tool} --version timed out after 10 seconds\n\n")
            except Exception as e:
                f.write(f"Error retrieving version: {e}\n\n")
    logger.info("Tool versions written to %s", version_log)
    return version_log


# ── Diagnostics ──────────────────────────────────────────────────────


def diagnostic_mode() -> None:
    """Print a quick summary of tool availability and system resources."""
    logger.info("[Diagnostics] Checking required tools and system status...\n")
    for tool in TOOL_NAMES:
        path = shutil.which(tool)
        if path:
            logger.info("[OK] %s found at %s", tool, path)
        else:
            logger.info("[Missing] %s not found in PATH", tool)

    total, used, free = shutil.disk_usage(os.getcwd())
    logger.info(
        "Diagnostics: Disk space - Total: %.2f GB, Used: %.2f GB, Free: %.2f GB",
        total / 1e9,
        used / 1e9,
        free / 1e9,
    )

    for pkg in ["pysam", "json", "argparse", "subprocess", "concurrent.futures"]:
        try:
            __import__(pkg)
            logger.info("[OK] Python package %s is installed", pkg)
        except ImportError:
            logger.warning("Missing Python package: %s", pkg)

    logger.info("\n[Diagnostics] Check complete.")


# ── Logging setup ────────────────────────────────────────────────────


def setup_logging(log_file: str | None = None, quiet: bool = False) -> None:
    """Configure the root logger with console (+ optional file) handlers."""
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    fmt = "%(asctime)s %(levelname)-8s %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"
    handlers: list[logging.Handler] = []

    console_level = logging.WARNING if quiet else logging.INFO
    ch = logging.StreamHandler()
    ch.setLevel(console_level)
    ch.setFormatter(logging.Formatter(fmt, datefmt))
    handlers.append(ch)

    if log_file:
        fh = logging.FileHandler(log_file, mode="w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter(fmt, datefmt))
        handlers.append(fh)

    logging.basicConfig(level=logging.DEBUG, handlers=handlers)


# ── FASTA ↔ FASTQ ───────────────────────────────────────────────────


def fasta_to_fastq(
    fasta_file: str, fastq_file: str, quality_score: int = 40
) -> None:
    """Convert a FASTA file to FASTQ using a flat quality score."""
    qual_char = chr(quality_score + 33)
    with open(fasta_file) as fa, open(fastq_file, "w") as fq:
        name: str | None = None
        seq: list[str] = []
        for line in fa:
            line = line.strip()
            if line.startswith(">"):
                if name and seq:
                    s = "".join(seq)
                    fq.write(f"@{name}\n{s}\n+\n{qual_char * len(s)}\n")
                elif name and not seq:
                    logger.warning(
                        "FASTA record %s has no sequence data and was skipped.", name
                    )
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name and seq:
            s = "".join(seq)
            fq.write(f"@{name}\n{s}\n+\n{qual_char * len(s)}\n")
        elif name and not seq:
            logger.warning(
                "FASTA record %s has no sequence data and was skipped.", name
            )


# ── Auto-compute readmers ───────────────────────────────────────────


def compute_readmers_db(
    reads: str, output_readmers: str, k: int = 31, threads: int = 8
) -> None:
    """Build a Meryl k-mer database from *reads* if it doesn't already exist."""
    cmd = [
        MERYL,
        f"k={k}",
        f"threads={threads}",
        "count",
        reads,
        "output",
        output_readmers,
    ]
    if os.path.exists(output_readmers) and not os.path.isdir(output_readmers):
        logger.warning(
            "%s exists as a file; removing it so meryl can create a directory.",
            output_readmers,
        )
        os.remove(output_readmers)
    conditional_run(
        outfile=output_readmers,
        resume_flag=True,
        step_desc="[compute_readmers_db] Computing readmers database with Meryl",
        command=cmd,
        max_log_lines=0,
    )


# T2T-Polish v4.0
# Any usage is subject to this software's license.
