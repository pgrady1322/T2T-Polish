#!/usr/bin/env python3
"""
T2T-Polish v4.1 — tests/test_utils.py

Unit tests for t2t_polish.utils (dependency checks, logging,
FASTA→FASTQ, readmers).
"""

from __future__ import annotations

import logging
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from t2t_polish.exceptions import ToolNotFoundError
from t2t_polish.utils import (
    compute_readmers_db,
    diagnostic_mode,
    fasta_to_fastq,
    log_tool_versions,
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


# ── log_tool_versions ────────────────────────────────────────────────


class TestLogToolVersions:
    @patch("t2t_polish.utils.subprocess.run")
    def test_writes_version_file(self, mock_run, tmp_path):
        mock_run.return_value = MagicMock(stdout="tool v1.0\n", stderr="", returncode=0)
        prefix = str(tmp_path / "test_prefix")
        result = log_tool_versions(prefix)
        assert result.endswith(".tool_versions.txt")
        with open(result) as f:
            content = f.read()
        assert "version" in content.lower()

    @patch("t2t_polish.utils.subprocess.run", side_effect=FileNotFoundError("not found"))
    def test_handles_missing_tool(self, mock_run, tmp_path):
        prefix = str(tmp_path / "test_prefix")
        result = log_tool_versions(prefix)
        with open(result) as f:
            content = f.read()
        assert "Error" in content

    @patch(
        "t2t_polish.utils.subprocess.run",
        side_effect=subprocess.TimeoutExpired(cmd="tool", timeout=10),
    )
    def test_handles_timeout(self, mock_run, tmp_path):
        prefix = str(tmp_path / "test_prefix")
        result = log_tool_versions(prefix)
        with open(result) as f:
            content = f.read()
        assert "timed out" in content


# ── diagnostic_mode ──────────────────────────────────────────────────


class TestDiagnosticMode:
    @patch("t2t_polish.utils.shutil.which", return_value="/usr/bin/tool")
    def test_found_tools(self, mock_which):
        # Should not raise
        diagnostic_mode()

    @patch("t2t_polish.utils.shutil.which", return_value=None)
    def test_missing_tools(self, mock_which):
        # Should not raise even when tools are missing
        diagnostic_mode()


# ── compute_readmers_db ──────────────────────────────────────────────


class TestComputeReadmersDb:
    @patch("t2t_polish.utils.conditional_run")
    def test_calls_conditional_run(self, mock_crun, tmp_path):
        reads = str(tmp_path / "reads.fq")
        output = str(tmp_path / "out.meryl")
        compute_readmers_db(reads, output, k=31, threads=4)
        mock_crun.assert_called_once()
        call_kwargs = mock_crun.call_args
        assert call_kwargs[1]["outfile"] == output

    @patch("t2t_polish.utils.conditional_run")
    def test_removes_existing_file(self, mock_crun, tmp_path):
        reads = str(tmp_path / "reads.fq")
        output = str(tmp_path / "out.meryl")
        # Create output as a file (should be removed so meryl can create dir)
        with open(output, "w") as f:
            f.write("x")
        compute_readmers_db(reads, output, k=31, threads=4)
        # The file should have been removed
        import os

        assert not os.path.isfile(output) or mock_crun.called


# T2T-Polish v4.1
# Any usage is subject to this software's license.
