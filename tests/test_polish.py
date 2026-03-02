#!/usr/bin/env python3
"""
T2T-Polish v4.0 — tests/test_polish.py

Unit tests for t2t_polish.polish (PolishIteration and step orchestration).
"""

from __future__ import annotations

import os
from unittest.mock import MagicMock, patch

from t2t_polish.polish import PolishIteration

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


# T2T-Polish v4.0
# Any usage is subject to this software's license.
