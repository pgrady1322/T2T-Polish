#!/usr/bin/env python3
"""
T2T-Polish v4.1 — kcov.py

``computekcov`` subcommand: Jellyfish count / histo → GenomeScope2 → kcov.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

from __future__ import annotations

import json
import logging
import os
import re

from t2t_polish.constants import GENOMESCOPE, JELLYFISH
from t2t_polish.exceptions import InputValidationError
from t2t_polish.runner import conditional_run

logger = logging.getLogger(__name__)


def sub_computekcov(args: object) -> tuple[str | None, str]:
    """Run Jellyfish + GenomeScope2 to determine ideal k-mer coverage.

    Returns ``(discovered_kcov, fitted_hist_path)``.
    """
    reads: str = getattr(args, "reads", "")
    if not os.path.exists(reads):
        raise InputValidationError(f"Reads file not found: {reads}")

    prefix: str = getattr(args, "prefix", "KCoverage")
    resume: bool = getattr(args, "resume", False)
    kmer_size: int = getattr(args, "kmer_size", 21)
    threads: int = getattr(args, "threads", 8)
    ploidy: str = getattr(args, "ploidy", "haploid")
    jellyfish_hash_size: int = getattr(args, "jellyfish_hash_size", 1_000_000_000)

    # ── Resume fast-path ─────────────────────────────────────────────
    hist_file = f"{prefix}.histo"
    gs_lookup = os.path.join(f"{prefix}.genomescope_out", "lookup_table.txt")
    if resume and os.path.exists(hist_file) and os.path.exists(gs_lookup):
        logger.info("compute_kcov: Resume enabled; skipping compute-kcov steps as results exist.")
        params_txt = os.path.join(f"{prefix}.genomescope_out", "params.txt")
        discovered_kcov = _parse_kcov_from_params(params_txt)
        fitted_hist = gs_lookup
        logger.info("compute_kcov: GenomeScope2 finished (resumed).")
        logger.info(
            "compute_kcov: Ideal KCov = %s",
            discovered_kcov if discovered_kcov else "Unknown",
        )
        logger.info("compute_kcov: Fitted histogram location = %s", fitted_hist)
        return discovered_kcov, fitted_hist

    # ── Full computation ─────────────────────────────────────────────
    logger.info(
        "compute_kcov: Starting computation of k-mer coverage using Jellyfish and GenomeScope2"
    )
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)

    jf_out = f"{prefix}.jf"
    jf_histo = f"{prefix}.histo"
    gs_outdir = f"{prefix}.genomescope_out"

    # Step 1 — Jellyfish count
    logger.info("compute_kcov: Step 1: Running Jellyfish K-mer counting (this may take a while)...")
    conditional_run(
        outfile=jf_out,
        resume_flag=resume,
        step_desc="[computekcov step 1] Jellyfish K-mer counting",
        command=[
            JELLYFISH,
            "count",
            "-C",
            "-m",
            str(kmer_size),
            "-s",
            str(jellyfish_hash_size),
            "-t",
            str(threads),
            "-o",
            jf_out,
            reads,
        ],
    )

    # Step 2 — Jellyfish histo
    conditional_run(
        outfile=jf_histo,
        resume_flag=(resume and os.path.exists(jf_histo)),
        step_desc="[computekcov step 2] Jellyfish histogram generation",
        command=[
            JELLYFISH,
            "histo",
            "-t",
            str(threads),
            jf_out,
            "-o",
            jf_histo,
        ],
    )

    # Step 3 — GenomeScope2
    logger.info(
        "compute_kcov: Step 3: Running GenomeScope2 for histogram fitting and coverage estimation..."
    )
    gs_cmd = [
        GENOMESCOPE,
        "-i",
        jf_histo,
        "-o",
        gs_outdir,
        "-k",
        str(kmer_size),
        "--fitted_hist",
    ]
    if ploidy == "haploid":
        gs_cmd += ["--ploidy", "1"]
    else:
        gs_cmd += ["--ploidy", "2"]

    result = conditional_run(
        outfile=os.path.join(gs_outdir, "lookup_table.txt"),
        resume_flag=resume,
        step_desc=(
            "[computekcov step 3] GenomeScope2 fitted histogram file creation "
            "and coverage calculation for pipeline"
        ),
        command=gs_cmd,
    )

    # Parse kcov from stdout first, then fall back to params.txt
    discovered_kcov = None
    if result and result.stdout:
        match = re.search(r"kcov:([0-9]+(?:\.[0-9]+)?)", result.stdout)
        if match:
            discovered_kcov = match.group(1)

    params_txt = os.path.join(gs_outdir, "params.txt")
    if discovered_kcov is None:
        discovered_kcov = _parse_kcov_from_params(params_txt)

    fitted_hist = os.path.join(gs_outdir, "lookup_table.txt")

    logger.info("compute_kcov: GenomeScope2 finished.")
    logger.info(
        "compute_kcov: Ideal KCov = %s",
        discovered_kcov if discovered_kcov else "Unknown",
    )
    logger.info("compute_kcov: Fitted histogram location = %s", fitted_hist)
    logger.info(
        "If running compute_kcov separately, pass --fitted_hist %s and --peak %s "
        "to the polish step. Otherwise, the automated pipeline will incorporate them.",
        fitted_hist,
        discovered_kcov,
    )

    # Persist metadata for resume
    meta = {"kcov": discovered_kcov, "fitted_hist": fitted_hist}
    meta_prefix = prefix[:-5] if prefix.endswith(".kcov") else prefix
    meta_file = f"{meta_prefix}.kcov.json"
    with open(meta_file, "w") as mf:
        json.dump(meta, mf)
    return discovered_kcov, fitted_hist


# ── Private helpers ──────────────────────────────────────────────────


def _parse_kcov_from_params(params_txt: str) -> str | None:
    """Extract *kcov* value from a GenomeScope2 ``params.txt`` file."""
    if not os.path.exists(params_txt):
        return None
    with open(params_txt) as f:
        for line in f:
            match = re.match(r"kcov:\s*(\S+)", line)
            if match:
                return match.group(1)
    return None


# T2T-Polish v4.1
# Any usage is subject to this software's license.
