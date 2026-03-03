#!/usr/bin/env python3
"""
T2T-Polish v4.1 — evaluation.py

QV evaluation helpers via **Merfin** (hist + completeness) and **Merqury**.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from subprocess import PIPE

from t2t_polish.constants import MERFIN, MERQURY_SH
from t2t_polish.runner import CommandResult, conditional_run, run_command

logger = logging.getLogger(__name__)


def run_merfin_eval(
    consensus: str,
    readmers: str,
    optimized: bool,
    ideal_kcov: float | str | None,
    fitted_hist_location: str | None,
) -> str:
    """Return a filtered Merfin QV / completeness summary string."""
    hist_out = f"{consensus}.merfin_hist.txt"
    comp_out = f"{consensus}.merfin_completeness.txt"

    peak_val = str(ideal_kcov) if optimized and ideal_kcov else "106.7"

    # ── histogram command ────────────────────────────────────────────
    hist_cmd = [
        MERFIN, "-hist",
        "-sequence", consensus,
        "-readmers", readmers,
        "-peak", peak_val,
        "-output", hist_out,
    ]
    if optimized and fitted_hist_location:
        hist_cmd += ["-prob", fitted_hist_location]

    # ── completeness command ─────────────────────────────────────────
    comp_cmd = [
        MERFIN, "-completeness",
        "-sequence", consensus,
        "-readmers", readmers,
        "-peak", peak_val,
        "-output", comp_out,
    ]
    if optimized and fitted_hist_location:
        comp_cmd += ["-prob", fitted_hist_location]

    # 1) Histogram — use run_command for consistent logging/timing
    hist_result: CommandResult = run_command(
        hist_cmd, description="[Merfin Eval] Histogram calculation"
    )
    if hist_result.returncode != 0:
        logger.error("[Merfin Eval] Histogram calculation failed: %s", hist_result.stderr)
        hist_text = ""
    else:
        hist_text = (hist_result.stdout or hist_result.stderr).strip()

    # 2) Completeness
    comp_result: CommandResult = run_command(
        comp_cmd, description="[Merfin Eval] Completeness calculation"
    )
    if comp_result.returncode != 0:
        logger.error("[Merfin Eval] Completeness calculation failed: %s", comp_result.stderr)
        comp_text = ""
    else:
        comp_text = (comp_result.stdout or comp_result.stderr).strip()

    hist_lines = [ln for ln in hist_text.splitlines() if re.search(r"(QV|Mean)", ln)]
    comp_lines = [ln for ln in comp_text.splitlines() if re.search(r"Completeness", ln)]
    return "\n".join(hist_lines + comp_lines)


def run_merqury_eval(
    consensus: str,
    readmers: str,
    resume: bool = False,
) -> str:
    """Return a two-line Merqury QV + completeness summary string."""
    iter_dir = os.path.dirname(os.path.abspath(consensus))
    base = Path(consensus).stem
    outprefix = os.path.join(iter_dir, f"{base}_merqury")

    cmd = ["bash", MERQURY_SH, readmers, consensus, outprefix]

    conditional_run(
        outfile=outprefix + ".qv",
        resume_flag=resume,
        step_desc="[Merqury Eval] Running merqury.sh",
        command=cmd,
    )

    if not os.path.exists(outprefix + ".qv"):
        logger.error(
            "Merqury failed: %s did not produce %s", MERQURY_SH, outprefix + ".qv"
        )
        return "Merqury QV: ERROR\nMerqury Completeness: ERROR"

    # Extract QV (column 4)
    qv_value = _read_column(outprefix + ".qv", col=3)

    # Extract completeness (column 5)
    comp_value = _read_column(outprefix + ".completeness.stats", col=4)

    lines: list[str] = []
    lines.append(f"Merqury QV: {qv_value}" if qv_value else "Merqury QV: N/A")
    lines.append(
        f"Merqury Completeness: {comp_value}" if comp_value else "Merqury Completeness: N/A"
    )
    return "\n".join(lines)


# ── Private helpers ──────────────────────────────────────────────────


def _read_column(path: str, col: int) -> str | None:
    """Read the first qualifying row from a whitespace-delimited file."""
    if not os.path.exists(path):
        logger.warning("Expected file not found: %s", path)
        return None
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) > col:
                return parts[col]
    return None


# T2T-Polish v4.1
# Any usage is subject to this software's license.
