#!/usr/bin/env python3
"""
T2T-Polish v4.0 — tests/test_kcov.py

Unit tests for t2t_polish.kcov (sub_computekcov, _parse_kcov_from_params).
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

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
    def test_missing_reads_exits(self, tmp_path):
        class Args:
            reads = str(tmp_path / "nonexistent.fq")
            prefix = str(tmp_path / "out")
            resume = False
            kmer_size = 21
            threads = 1
            ploidy = "haploid"
            jellyfish_hash_size = 100

        with pytest.raises(SystemExit):
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


# T2T-Polish v4.0
# Any usage is subject to this software's license.
