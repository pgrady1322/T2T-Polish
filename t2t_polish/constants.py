#!/usr/bin/env python3
"""
T2T-Polish v4.1 — constants.py

Shared constants: tool names, default quality scores, and k-mer size.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

from __future__ import annotations

# External bioinformatics tool executable names
WINNOWMAP: str = "winnowmap"
FALCONC: str = "falconc"
MERYL: str = "meryl"
MERFIN: str = "merfin"
BCFTOOLS: str = "bcftools"
JELLYFISH: str = "jellyfish"
GENOMESCOPE: str = "genomescope2"
SAMTOOLS: str = "samtools"
MERQURY_SH: str = "merqury.sh"

TOOL_NAMES: list[str] = [
    WINNOWMAP,
    FALCONC,
    MERYL,
    MERFIN,
    BCFTOOLS,
    JELLYFISH,
    GENOMESCOPE,
    SAMTOOLS,
    MERQURY_SH,
]

# Default quality scores for read correctors (Phred scale)
DEFAULT_CORRECTOR_Q: dict[str, int] = {
    "hifiasm": 40,
    "herro": 35,
    "flye": 30,
    "unknown": 30,
}

DEFAULT_KMER_SIZE: int = 31

# T2T-Polish v4.0
# Any usage is subject to this software's license.
