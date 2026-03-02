#!/usr/bin/env python3
"""
T2T-Polish v4.0 — tests/test_cli.py

Unit tests for t2t_polish.cli (argument parsing and --version).
"""

from __future__ import annotations

import pytest

from t2t_polish.cli import _build_parser

# ── parser structure ─────────────────────────────────────────────────


class TestBuildParser:
    def test_parser_has_subcommands(self):
        parser = _build_parser()
        assert parser is not None

    def test_version_flag(self, capsys):
        parser = _build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["--version"])
        assert exc_info.value.code == 0
        captured = capsys.readouterr()
        assert "4.0.0" in captured.out


# ── polish subcommand ────────────────────────────────────────────────


class TestPolishSubcommand:
    def test_valid_args(self, tmp_path):
        parser = _build_parser()
        draft = str(tmp_path / "asm.fa")
        reads = str(tmp_path / "reads.fq")
        sif = str(tmp_path / "dv.sif")
        for p in (draft, reads, sif):
            with open(p, "w") as f:
                f.write("dummy")
        args = parser.parse_args([
            "polish",
            "-d", draft,
            "-r", reads,
            "--singularity_sif", sif,
            "-t", "8",
            "-i", "2",
        ])
        assert args.draft == draft
        assert args.threads == 8
        assert args.iterations == 2

    def test_missing_singularity_sif(self):
        parser = _build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["polish", "-d", "asm.fa", "-r", "reads.fq"])
        assert exc_info.value.code != 0


# ── computekcov subcommand ───────────────────────────────────────────


class TestComputeKcovSubcommand:
    def test_parses_reads(self, tmp_path):
        parser = _build_parser()
        reads = str(tmp_path / "reads.fq")
        with open(reads, "w") as f:
            f.write("dummy")
        args = parser.parse_args([
            "computekcov",
            "-r", reads,
            "-t", "4",
        ])
        assert args.reads == reads
        assert args.threads == 4


# ── diagnostics subcommand ───────────────────────────────────────────


class TestDiagnosticsSubcommand:
    def test_diagnostics_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["diagnostics"])
        assert args.subcommand == "diagnostics"


# T2T-Polish v4.0
# Any usage is subject to this software's license.
