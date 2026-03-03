#!/usr/bin/env python3
"""
T2T-Polish v4.1 — tests/test_utils.py

Unit tests for t2t_polish.utils (dependency checks, logging,
FASTA→FASTQ, readmers).
"""

from __future__ import annotations

import logging
from unittest.mock import patch

import pytest

from t2t_polish.exceptions import ToolNotFoundError
from t2t_polish.utils import (
    fasta_to_fastq,
    needs_kcov,
    setup_logging,
    validate_dependencies,
)

# ── validate_dependencies ────────────────────────────────────────────


class TestValidateDependencies:
    def test_all_found(self, all_tools_on_path):
        # Should not raise
        validate_dependencies()

    def test_missing_tool_raises(self, no_tools_on_path):
        with pytest.raises(ToolNotFoundError, match="not found in PATH"):
            validate_dependencies()

    def test_partial_missing_raises(self):
        def _fake_which(name):
            if name == "merfin":
                return None
            return "/usr/bin/" + name

        with patch("t2t_polish.utils.shutil.which", side_effect=_fake_which):
            with pytest.raises(ToolNotFoundError, match="merfin"):
                validate_dependencies()


# ── needs_kcov ───────────────────────────────────────────────────────


class TestNeedsKcov:
    def test_not_optimized(self):
        class A:
            optimized = False
            fitted_hist = None
            ploidy = "haploid"
            ideal_dpeak = None
            ideal_hpeak = None

        assert needs_kcov(A()) is False

    def test_optimized_missing_hist(self):
        class A:
            optimized = True
            fitted_hist = None
            ploidy = "haploid"
            ideal_dpeak = None
            ideal_hpeak = None

        assert needs_kcov(A()) is True

    def test_optimized_all_present_haploid(self):
        class A:
            optimized = True
            fitted_hist = "/some/path"
            ploidy = "haploid"
            ideal_dpeak = 115.0
            ideal_hpeak = None

        assert needs_kcov(A()) is False

    def test_optimized_diploid_missing_hpeak(self):
        class A:
            optimized = True
            fitted_hist = "/some/path"
            ploidy = "diploid"
            ideal_dpeak = None
            ideal_hpeak = None

        assert needs_kcov(A()) is True


# ── setup_logging ────────────────────────────────────────────────────


class TestSetupLogging:
    def test_console_only(self):
        setup_logging()
        root = logging.getLogger()
        # At least one handler (console)
        assert len(root.handlers) >= 1

    def test_with_file(self, tmp_path):
        log_file = str(tmp_path / "test.log")
        setup_logging(log_file=log_file)
        root = logging.getLogger()
        assert any(isinstance(h, logging.FileHandler) for h in root.handlers)

    def test_quiet_flag(self):
        setup_logging(quiet=True)
        root = logging.getLogger()
        stream_handlers = [
            h
            for h in root.handlers
            if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        ]
        for sh in stream_handlers:
            assert sh.level >= logging.WARNING


# ── fasta_to_fastq ──────────────────────────────────────────────────


class TestFastaToFastq:
    def test_basic_conversion(self, mock_fasta, tmp_path):
        out = str(tmp_path / "out.fastq")
        fasta_to_fastq(mock_fasta, out, quality_score=40)
        with open(out) as f:
            lines = f.read().splitlines()
        # Two records → 8 lines (4 per record)
        assert len(lines) == 8
        assert lines[0].startswith("@")
        assert lines[2] == "+"

    def test_quality_char(self, mock_fasta, tmp_path):
        out = str(tmp_path / "q30.fastq")
        fasta_to_fastq(mock_fasta, out, quality_score=30)
        with open(out) as f:
            lines = f.read().splitlines()
        expected_char = chr(30 + 33)
        assert all(c == expected_char for c in lines[3])

    def test_empty_record_skipped(self, mock_fasta_empty_record, tmp_path):
        out = str(tmp_path / "skip.fastq")
        fasta_to_fastq(mock_fasta_empty_record, out, quality_score=40)
        with open(out) as f:
            content = f.read()
        # Only seq2 should appear
        assert "@seq2" in content
        assert "@empty_seq" not in content


# T2T-Polish v4.1
# Any usage is subject to this software's license.
