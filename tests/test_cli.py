#!/usr/bin/env python3
"""
T2T-Polish v4.1 — tests/test_cli.py

Unit tests for t2t_polish.cli (argument parsing, --version, required args).
"""

from __future__ import annotations

from unittest.mock import patch

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
        assert "4.1.2" in captured.out


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
        args = parser.parse_args(
            [
                "polish",
                "-d",
                draft,
                "-r",
                reads,
                "--singularity_sif",
                sif,
                "-t",
                "8",
                "-i",
                "2",
            ]
        )
        assert args.draft == draft
        assert args.threads == 8
        assert args.iterations == 2

    def test_missing_singularity_sif(self):
        parser = _build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["polish", "-d", "asm.fa", "-r", "reads.fq"])
        assert exc_info.value.code != 0

    def test_missing_draft_exits(self):
        """--draft is now required; omitting it should trigger argparse exit."""
        parser = _build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["polish", "-r", "reads.fq", "--singularity_sif", "dv.sif"])
        assert exc_info.value.code != 0

    def test_missing_reads_exits(self):
        """--reads is now required; omitting it should trigger argparse exit."""
        parser = _build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["polish", "-d", "asm.fa", "--singularity_sif", "dv.sif"])
        assert exc_info.value.code != 0


# ── computekcov subcommand ───────────────────────────────────────────


class TestComputeKcovSubcommand:
    def test_parses_reads(self, tmp_path):
        parser = _build_parser()
        reads = str(tmp_path / "reads.fq")
        with open(reads, "w") as f:
            f.write("dummy")
        args = parser.parse_args(
            [
                "computekcov",
                "-r",
                reads,
                "-t",
                "4",
            ]
        )
        assert args.reads == reads
        assert args.threads == 4


# ── diagnostics subcommand ───────────────────────────────────────────


class TestDiagnosticsSubcommand:
    def test_diagnostics_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["diagnostics"])
        assert args.subcommand == "diagnostics"


# ── main() integration ───────────────────────────────────────────────


class TestMain:
    def test_no_subcommand_exits(self):
        """main() with no subcommand prints help and exits 1."""
        with patch("sys.argv", ["t2t-polish"]):
            with pytest.raises(SystemExit) as exc_info:
                from t2t_polish.cli import main

                main()
            assert exc_info.value.code == 1

    def test_thread_count_zero_exits(self, tmp_path):
        """main() exits when threads < 1."""
        draft = tmp_path / "d.fa"
        reads = tmp_path / "r.fq"
        sif = tmp_path / "dv.sif"
        for p in (draft, reads, sif):
            p.write_text("x")
        with patch(
            "sys.argv",
            [
                "t2t-polish",
                "polish",
                "-d",
                str(draft),
                "-r",
                str(reads),
                "--singularity_sif",
                str(sif),
                "-t",
                "0",
            ],
        ):
            with pytest.raises(SystemExit) as exc_info:
                from t2t_polish.cli import main

                main()
            assert exc_info.value.code == 1

    def test_high_threads_warns(self, tmp_path):
        """main() warns but does not exit when threads > 128."""
        draft = tmp_path / "d.fa"
        reads = tmp_path / "r.fq"
        sif = tmp_path / "dv.sif"
        for p in (draft, reads, sif):
            p.write_text("x")
        with (
            patch(
                "sys.argv",
                [
                    "t2t-polish",
                    "polish",
                    "-d",
                    str(draft),
                    "-r",
                    str(reads),
                    "--singularity_sif",
                    str(sif),
                    "-t",
                    "200",
                ],
            ),
            patch("t2t_polish.cli.validate_dependencies"),
            patch("t2t_polish.cli.log_tool_versions"),
            patch("t2t_polish.cli.sub_polish"),
        ):
            from t2t_polish.cli import main

            main()  # Should not raise

    def test_missing_input_file_exits(self, tmp_path):
        """main() exits when --draft file does not exist."""
        with (
            patch(
                "sys.argv",
                [
                    "t2t-polish",
                    "polish",
                    "-d",
                    str(tmp_path / "nonexistent.fa"),
                    "-r",
                    str(tmp_path / "also_missing.fq"),
                    "--singularity_sif",
                    "dv.sif",
                ],
            ),
            patch("t2t_polish.cli.validate_dependencies"),
        ):
            with pytest.raises(SystemExit) as exc_info:
                from t2t_polish.cli import main

                main()
            assert exc_info.value.code == 1

    def test_diagnostics_dispatch(self):
        """main() dispatches to diagnostic_mode for diagnostics subcommand."""
        with (
            patch("sys.argv", ["t2t-polish", "diagnostics"]),
            patch("t2t_polish.cli.diagnostic_mode") as mock_diag,
        ):
            from t2t_polish.cli import main

            main()
            mock_diag.assert_called_once()

    def test_empty_input_file_exits(self, tmp_path):
        """main() exits when --draft file exists but is empty."""
        draft = tmp_path / "empty.fa"
        reads = tmp_path / "reads.fq"
        sif = tmp_path / "dv.sif"
        draft.write_text("")
        reads.write_text("data")
        sif.write_text("data")
        with (
            patch(
                "sys.argv",
                [
                    "t2t-polish",
                    "polish",
                    "-d",
                    str(draft),
                    "-r",
                    str(reads),
                    "--singularity_sif",
                    str(sif),
                ],
            ),
            patch("t2t_polish.cli.validate_dependencies"),
            patch("t2t_polish.cli.log_tool_versions"),
        ):
            with pytest.raises(SystemExit) as exc_info:
                from t2t_polish.cli import main

                main()
            assert exc_info.value.code == 1


# T2T-Polish v4.1
# Any usage is subject to this software's license.
