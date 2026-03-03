#!/usr/bin/env python3
"""
T2T-Polish v4.1 — tests/test_kcov.py

Unit tests for t2t_polish.kcov (sub_computekcov, _parse_kcov_from_params).
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

from t2t_polish.exceptions import InputValidationError
from t2t_polish.kcov import _parse_kcov_from_params, sub_computekcov

# ── _parse_kcov_from_params ──────────────────────────────────────────


class TestParseKcovFromParams:
    def test_found(self, tmp_path):
        p = tmp_path / "params.txt"
        p.write_text("kcov: 115.3\nother: foo\n")
        assert _parse_kcov_from_params(str(p)) == "115.3"

    def test_not_found(self, tmp_path):
        p = tmp_path / "params.txt"
        p.write_text("no_kcov_here: 42\n")
        assert _parse_kcov_from_params(str(p)) is None

    def test_missing_file(self, tmp_path):
        assert _parse_kcov_from_params(str(tmp_path / "nope.txt")) is None


# ── sub_computekcov ──────────────────────────────────────────────────


class TestSubComputekcov:
    def test_missing_reads_raises(self, tmp_path):
        class Args:
            reads = str(tmp_path / "nonexistent.fq")
            prefix = str(tmp_path / "out")
            resume = False
            kmer_size = 21
            threads = 1
            ploidy = "haploid"
            jellyfish_hash_size = 100

        with pytest.raises(InputValidationError, match="not found"):
            sub_computekcov(Args())

    @patch("t2t_polish.kcov.conditional_run")
    def test_resume_fast_path(self, mock_crun, tmp_path):
        """When resume=True and outputs exist, skip computation."""
        prefix = str(tmp_path / "kc")
        # Create the files that the resume check looks for
        hist = tmp_path / "kc.histo"
        hist.write_text("1 100\n")
        gs_dir = tmp_path / "kc.genomescope_out"
        gs_dir.mkdir()
        (gs_dir / "lookup_table.txt").write_text("lookup data\n")
        params = gs_dir / "params.txt"
        params.write_text("kcov: 88\n")

        _prefix = prefix  # avoid name clash with class attr

        class Args:
            reads = str(tmp_path / "reads.fq")
            prefix = _prefix
            resume = True
            kmer_size = 21
            threads = 1
            ploidy = "haploid"
            jellyfish_hash_size = 100

        # Also create reads file so it passes existence check
        (tmp_path / "reads.fq").write_text("@r\nACGT\n+\nIIII\n")

        kcov, fitted = sub_computekcov(Args())
        assert kcov == "88"
        assert "lookup_table.txt" in fitted
        # conditional_run should NOT have been called (fast-pathed)
        mock_crun.assert_not_called()

    @patch("t2t_polish.kcov.conditional_run")
    def test_full_computation_path(self, mock_crun, tmp_path):
        """When resume=False, all three steps run and metadata is written."""
        from t2t_polish.runner import CommandResult

        prefix = str(tmp_path / "kc")
        reads = tmp_path / "reads.fq"
        reads.write_text("@r\nACGT\n+\nIIII\n")

        # Create the genomescope output dir so metadata can be written
        gs_dir = tmp_path / "kc.genomescope_out"
        gs_dir.mkdir()
        (gs_dir / "lookup_table.txt").write_text("lookup")
        (gs_dir / "params.txt").write_text("kcov: 92.5\n")

        mock_crun.return_value = CommandResult(
            stdout="some output kcov:92.5 done",
            stderr="",
            returncode=0,
        )

        _prefix = prefix

        class Args:
            reads = str(tmp_path / "reads.fq")
            prefix = _prefix
            resume = False
            kmer_size = 21
            threads = 1
            ploidy = "haploid"
            jellyfish_hash_size = 100

        kcov, fitted = sub_computekcov(Args())
        assert kcov == "92.5"
        assert "lookup_table.txt" in fitted
        # Three steps: jellyfish count, jellyfish histo, genomescope2
        assert mock_crun.call_count == 3

        # Metadata JSON should have been written
        import json

        meta_file = tmp_path / "kc.kcov.json"
        assert meta_file.exists()
        data = json.loads(meta_file.read_text())
        assert data["kcov"] == "92.5"

    @patch("t2t_polish.kcov.conditional_run")
    def test_falls_back_to_params_txt(self, mock_crun, tmp_path):
        """When stdout has no kcov, parser falls back to params.txt."""
        from t2t_polish.runner import CommandResult

        prefix = str(tmp_path / "kc")
        reads = tmp_path / "reads.fq"
        reads.write_text("@r\nACGT\n+\nIIII\n")

        gs_dir = tmp_path / "kc.genomescope_out"
        gs_dir.mkdir()
        (gs_dir / "lookup_table.txt").write_text("lookup")
        (gs_dir / "params.txt").write_text("kcov: 77.0\n")

        # Stdout has no kcov pattern
        mock_crun.return_value = CommandResult(
            stdout="no kcov data here",
            stderr="",
            returncode=0,
        )

        _prefix = prefix

        class Args:
            reads = str(tmp_path / "reads.fq")
            prefix = _prefix
            resume = False
            kmer_size = 21
            threads = 1
            ploidy = "diploid"
            jellyfish_hash_size = 100

        kcov, _ = sub_computekcov(Args())
        assert kcov == "77.0"


# T2T-Polish v4.1
# Any usage is subject to this software's license.
