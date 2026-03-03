#!/usr/bin/env python3
"""
T2T-Polish v4.1 — tests/test_polish.py

Unit tests for t2t_polish.polish (PolishIteration, PolishConfig, and step orchestration).
"""

from __future__ import annotations

import json
import os
from unittest.mock import MagicMock, patch

import pytest

from t2t_polish.exceptions import InputValidationError
from t2t_polish.polish import (
    PolishConfig,
    PolishIteration,
    _apply_config_file,
    _apply_resume_kcov,
    _check_optimized_resume,
    _compute_kcov_if_needed,
    _convert_fasta_reads,
    _resolve_chosen_peak,
    _resolve_seq_type,
    _validate_inputs,
    _write_summaries,
)

# ── PolishIteration ──────────────────────────────────────────────────


class TestPolishIteration:
    def test_paths_populated(self, tmp_path):
        folder = str(tmp_path / "iter_01")
        os.makedirs(folder, exist_ok=True)
        pi = PolishIteration(folder=folder, prefix="iter_1")
        assert "iter_1" in pi.bam
        assert "falconc" in pi.bam
        assert os.path.isdir(pi.dv_outdir)

    def test_consensus_path(self, tmp_path):
        folder = str(tmp_path / "iter_03")
        os.makedirs(folder, exist_ok=True)
        pi = PolishIteration(folder=folder, prefix="iter_3")
        assert pi.consensus.endswith("iter_3.consensus.fasta")


# ── PolishConfig ─────────────────────────────────────────────────────


class TestPolishConfig:
    def test_from_args_basic(self):
        class FakeArgs:
            draft = "/tmp/draft.fa"
            reads = "/tmp/reads.fq"
            prefix = "TestPrefix"
            threads = 16
            iterations = 5
            kmer_size = 21
            meryl_memory = 64
            singularity_sif = "/tmp/dv.sif"
            deepseq_type = "ONT_R104"
            seq_type = "ont"
            ploidy = "diploid"
            optimized = True
            resume = False
            resume_from = 0
            cleanup = True
            dv_tmp_dir = None
            read_corrector = "herro"
            readmers = "/tmp/readmers.meryl"
            ideal_dpeak = None
            ideal_hpeak = 55.0
            fitted_hist = "/tmp/hist.txt"
            config_file = None
            jellyfish_hash_size = 2_000_000_000
            json_summary = True
            log_file = None
            quiet = False

        cfg = PolishConfig.from_args(FakeArgs())
        assert cfg.draft == "/tmp/draft.fa"
        assert cfg.threads == 16
        assert cfg.ploidy == "diploid"
        assert cfg.ideal_hpeak == 55.0
        assert cfg.json_summary is True
        assert cfg.jellyfish_hash_size == 2_000_000_000

    def test_from_args_defaults(self):
        """Missing attributes on the namespace should fall back to defaults."""

        class MinimalArgs:
            draft = "d.fa"
            reads = "r.fq"

        cfg = PolishConfig.from_args(MinimalArgs())
        assert cfg.prefix == "AutoPolisher"
        assert cfg.threads == 32
        assert cfg.iterations == 3
        assert cfg.optimized is False

    def test_dataclass_is_mutable(self):
        cfg = PolishConfig(draft="a", reads="b")
        cfg.reads = "c"
        assert cfg.reads == "c"


# ── run_polish_iteration (mocked steps) ──────────────────────────────


class TestRunPolishIteration:
    @patch("t2t_polish.polish.run_merfin_eval", return_value="QV: 40")
    @patch("t2t_polish.polish.run_merqury_eval", return_value="Merqury QV: 40")
    @patch("t2t_polish.polish.conditional_run")
    def test_step_sequence(self, mock_crun, mock_merq, mock_merf, tmp_path):
        """Ensure conditional_run is invoked multiple times for pipeline steps."""
        from t2t_polish.polish import run_polish_iteration

        def _side_effect(outfile, **kwargs):
            """Simulate conditional_run creating the expected output file."""
            if outfile and not os.path.exists(outfile):
                os.makedirs(os.path.dirname(outfile), exist_ok=True)
                # Create as dir if it looks like a meryl DB, else as file
                if outfile.endswith("meryl") or outfile.endswith("merylDB"):
                    os.makedirs(outfile, exist_ok=True)
                else:
                    with open(outfile, "w") as f:
                        f.write("mock")
            return MagicMock(returncode=0, stdout="", stderr="")

        mock_crun.side_effect = _side_effect

        iter_folder = str(tmp_path / "iter_01")
        draft = str(tmp_path / "draft.fa")
        reads = str(tmp_path / "reads.fq")
        readmers = str(tmp_path / "readmers.meryl")
        sif = str(tmp_path / "dv.sif")

        # Create dummy input files
        for p in (draft, reads, sif):
            with open(p, "w") as f:
                f.write("dummy")

        run_polish_iteration(
            iter_folder=iter_folder,
            iteration_id=1,
            draft_fasta=draft,
            reads=reads,
            readmers=readmers,
            k_mer_size=31,
            num_threads=4,
            winnowmap_seq_type="pb",
            deepseq_type="PACBIO",
            singularity_sif=sif,
        )
        # conditional_run handles individual pipeline steps
        assert mock_crun.call_count > 0


# ── Polish helper functions ──────────────────────────────────────────


class TestValidateInputs:
    def test_missing_draft_raises(self, tmp_path):
        reads = tmp_path / "reads.fq"
        reads.write_text("@r\nACGT\n+\nIIII\n")
        cfg = PolishConfig(draft="", reads=str(reads))
        with pytest.raises(InputValidationError, match="draft"):
            _validate_inputs(cfg)

    def test_missing_reads_raises(self, tmp_path):
        draft = tmp_path / "draft.fa"
        draft.write_text(">seq\nACGT\n")
        cfg = PolishConfig(draft=str(draft), reads="/nonexistent")
        with pytest.raises(InputValidationError, match="reads"):
            _validate_inputs(cfg)

    def test_valid_inputs_passes(self, tmp_path):
        draft = tmp_path / "draft.fa"
        reads = tmp_path / "reads.fq"
        draft.write_text(">seq\nACGT\n")
        reads.write_text("@r\nACGT\n+\nIIII\n")
        cfg = PolishConfig(draft=str(draft), reads=str(reads))
        _validate_inputs(cfg)  # Should not raise


class TestResolveSeqType:
    def test_already_set(self):
        cfg = PolishConfig(draft="d", reads="r", seq_type="ont")
        result = _resolve_seq_type(cfg)
        assert result.seq_type == "ont"

    def test_pacbio_auto(self):
        cfg = PolishConfig(draft="d", reads="r", deepseq_type="PACBIO")
        result = _resolve_seq_type(cfg)
        assert result.seq_type == "pb"

    def test_ont_auto(self):
        cfg = PolishConfig(draft="d", reads="r", deepseq_type="ONT_R104")
        result = _resolve_seq_type(cfg)
        assert result.seq_type == "ont"

    def test_unknown_raises(self):
        cfg = PolishConfig(draft="d", reads="r", deepseq_type="UNKNOWN")
        with pytest.raises(InputValidationError, match="isn't automatically mapped"):
            _resolve_seq_type(cfg)


class TestConvertFastaReads:
    def test_non_fasta_is_noop(self, tmp_path):
        cfg = PolishConfig(draft="d", reads=str(tmp_path / "reads.fq"))
        result = _convert_fasta_reads(cfg)
        assert result.reads == cfg.reads

    def test_fasta_converts(self, tmp_path):
        fasta = tmp_path / "reads.fasta"
        fasta.write_text(">r\nACGT\n")
        cfg = PolishConfig(
            draft="d",
            reads=str(fasta),
            prefix=str(tmp_path / "out"),
            read_corrector="hifiasm",
        )
        result = _convert_fasta_reads(cfg)
        assert result.reads.endswith(".converted_reads.fastq")

    def test_fasta_no_corrector_defaults(self, tmp_path):
        fasta = tmp_path / "reads.fa"
        fasta.write_text(">r\nACGT\n")
        cfg = PolishConfig(
            draft="d",
            reads=str(fasta),
            prefix=str(tmp_path / "out"),
        )
        result = _convert_fasta_reads(cfg)
        assert result.read_corrector == "unknown"


class TestApplyConfigFile:
    def test_no_config_is_noop(self):
        cfg = PolishConfig(draft="d", reads="r")
        result = _apply_config_file(cfg)
        assert result.reads == "r"

    def test_loads_config_values(self, tmp_path):
        config = tmp_path / "config.json"
        config.write_text(json.dumps({"reads": "/new/reads.fq", "readmers": "/new/readmers.meryl"}))
        cfg = PolishConfig(draft="d", reads="", config_file=str(config))
        result = _apply_config_file(cfg)
        assert result.reads == "/new/reads.fq"
        assert result.readmers == "/new/readmers.meryl"

    def test_config_optimized_peaks(self, tmp_path):
        config = tmp_path / "config.json"
        config.write_text(
            json.dumps(
                {
                    "peak_haploid": 115.0,
                    "fitted_hist_location": "/path/hist.txt",
                }
            )
        )
        cfg = PolishConfig(draft="d", reads="r", optimized=True, config_file=str(config))
        result = _apply_config_file(cfg)
        assert result.ideal_dpeak == 115.0
        assert result.fitted_hist == "/path/hist.txt"


class TestResolveChosenPeak:
    def test_no_resume_returns_none(self):
        cfg = PolishConfig(draft="d", reads="r")
        assert _resolve_chosen_peak(cfg) is None

    def test_resume_haploid_with_peak(self):
        cfg = PolishConfig(
            draft="d",
            reads="r",
            resume=True,
            fitted_hist="/path",
            ploidy="haploid",
            ideal_dpeak=100.0,
        )
        result = _resolve_chosen_peak(cfg)
        assert result == 100.0

    def test_resume_diploid_with_peak(self):
        cfg = PolishConfig(
            draft="d",
            reads="r",
            resume=True,
            fitted_hist="/path",
            ploidy="diploid",
            ideal_hpeak=55.0,
        )
        result = _resolve_chosen_peak(cfg)
        assert result == 55.0


class TestApplyResumeKcov:
    def test_no_resume_is_noop(self, tmp_path):
        cfg = PolishConfig(draft="d", reads="r", prefix=str(tmp_path / "out"))
        result = _apply_resume_kcov(cfg)
        assert result.fitted_hist is None

    def test_loads_kcov_meta(self, tmp_path):
        prefix = str(tmp_path / "out")
        meta_file = tmp_path / "out.kcov.json"
        meta_file.write_text(json.dumps({"kcov": "88.5", "fitted_hist": "/path/hist.txt"}))
        cfg = PolishConfig(draft="d", reads="r", prefix=prefix, resume=True)
        result = _apply_resume_kcov(cfg)
        assert result.fitted_hist == "/path/hist.txt"
        assert result.ideal_dpeak == "88.5"


class TestComputeKcovIfNeeded:
    def test_not_optimized_is_noop(self):
        cfg = PolishConfig(draft="d", reads="r", optimized=False)
        peak, result_cfg = _compute_kcov_if_needed(cfg, None)
        assert peak is None

    def test_optimized_haploid_with_peak(self):
        cfg = PolishConfig(
            draft="d",
            reads="r",
            optimized=True,
            ploidy="haploid",
            ideal_dpeak=100.0,
            fitted_hist="/path/hist.txt",
        )
        peak, _ = _compute_kcov_if_needed(cfg, None)
        assert peak == 100.0

    def test_optimized_diploid_with_peak(self):
        cfg = PolishConfig(
            draft="d",
            reads="r",
            optimized=True,
            ploidy="diploid",
            ideal_hpeak=55.0,
            fitted_hist="/path/hist.txt",
        )
        peak, _ = _compute_kcov_if_needed(cfg, None)
        assert peak == 55.0


class TestCheckOptimizedResume:
    def test_not_optimized_is_noop(self):
        cfg = PolishConfig(draft="d", reads="r")
        _check_optimized_resume(cfg)  # Should not raise

    def test_optimized_resume_missing_meta_raises(self, tmp_path):
        cfg = PolishConfig(
            draft="d",
            reads="r",
            prefix=str(tmp_path / "out"),
            optimized=True,
            resume=True,
        )
        with pytest.raises(InputValidationError, match="Missing kcov metadata"):
            _check_optimized_resume(cfg)

    def test_optimized_resume_with_meta_passes(self, tmp_path):
        prefix = str(tmp_path / "out")
        meta = tmp_path / "out.kcov.json"
        meta.write_text(json.dumps({"kcov": "88", "fitted_hist": "/x"}))
        cfg = PolishConfig(
            draft="d",
            reads="r",
            prefix=prefix,
            optimized=True,
            resume=True,
        )
        _check_optimized_resume(cfg)  # Should not raise


class TestWriteSummaries:
    def test_writes_text_summary(self, tmp_path):
        prefix = str(tmp_path / "out")
        cfg = PolishConfig(draft="d", reads="r", prefix=prefix, iterations=2)
        records = [(1, "Merfin QV: 40", "Merqury QV: 40"), (2, "Merfin QV: 42", "Merqury QV: 42")]
        _write_summaries(cfg, records)
        summary = tmp_path / "out.QV_Completeness_summary.txt"
        assert summary.exists()
        content = summary.read_text()
        assert "Iteration 1" in content
        assert "Iteration 2" in content

    def test_writes_json_summary(self, tmp_path):
        prefix = str(tmp_path / "out")
        cfg = PolishConfig(
            draft="d",
            reads="r",
            prefix=prefix,
            iterations=1,
            json_summary=True,
        )
        records = [(1, "Merfin QV: 40", "Merqury QV: 40")]
        _write_summaries(cfg, records)
        json_file = tmp_path / "out.summary.json"
        assert json_file.exists()
        data = json.loads(json_file.read_text())
        assert data["pipeline"] == "T2T-Polish"
        assert data["version"] == "4.1.2"
        assert len(data["records"]) == 1


# T2T-Polish v4.1
# Any usage is subject to this software's license.
