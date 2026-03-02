#!/usr/bin/env python3
"""
T2T-Polish v4.0 — tests/test_evaluation.py

Unit tests for t2t_polish.evaluation (Merfin and Merqury QV helpers).
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from t2t_polish.evaluation import _read_column, run_merfin_eval, run_merqury_eval

# ── _read_column ─────────────────────────────────────────────────────


class TestReadColumn:
    def test_reads_correct_column(self, tmp_path):
        f = tmp_path / "data.tsv"
        f.write_text("a\tb\tc\td\te\n")
        assert _read_column(str(f), col=3) == "d"

    def test_missing_file(self, tmp_path):
        assert _read_column(str(tmp_path / "nope"), col=0) is None

    def test_short_row(self, tmp_path):
        f = tmp_path / "short.tsv"
        f.write_text("only\n")
        assert _read_column(str(f), col=5) is None


# ── run_merfin_eval ──────────────────────────────────────────────────


class TestRunMerfinEval:
    @patch("t2t_polish.evaluation.subprocess.run")
    def test_extracts_qv_lines(self, mock_run):
        # Simulate Merfin producing QV output
        mock_run.return_value = MagicMock(
            stdout="QV: 42.3\nMean: 1.5\nOther noise\n",
            stderr="",
            returncode=0,
        )
        result = run_merfin_eval(
            consensus="fake.fa",
            readmers="fake.meryl",
            optimized=False,
            ideal_kcov=None,
            fitted_hist_location=None,
        )
        assert "QV" in result or "Mean" in result

    @patch("t2t_polish.evaluation.subprocess.run")
    def test_handles_failure(self, mock_run):
        mock_run.return_value = MagicMock(
            stdout="", stderr="error msg", returncode=1
        )
        result = run_merfin_eval(
            "fake.fa", "fake.meryl", False, None, None
        )
        # Empty on failure — should not crash
        assert isinstance(result, str)


# ── run_merqury_eval ─────────────────────────────────────────────────


class TestRunMerquryEval:
    @patch("t2t_polish.evaluation.conditional_run")
    def test_qv_and_completeness(self, mock_crun, tmp_path):
        consensus = str(tmp_path / "cons.fasta")
        # Merqury derives output prefix from stem: cons_merqury
        qv_file = tmp_path / "cons_merqury.qv"
        qv_file.write_text("all\t100\t0.5\t45.2\n")
        comp_file = tmp_path / "cons_merqury.completeness.stats"
        comp_file.write_text("set\tall\t100\t98\t99.5\n")
        with open(consensus, "w") as f:
            f.write(">seq\nACGT\n")

        result = run_merqury_eval(consensus, "readmers.meryl", resume=False)
        assert "45.2" in result
        assert "99.5" in result

    @patch("t2t_polish.evaluation.conditional_run")
    def test_missing_qv_file(self, mock_crun, tmp_path):
        consensus = str(tmp_path / "cons.fasta")
        with open(consensus, "w") as f:
            f.write(">seq\nACGT\n")
        result = run_merqury_eval(consensus, "readmers.meryl", resume=False)
        assert "ERROR" in result


# T2T-Polish v4.0
# Any usage is subject to this software's license.
