#!/usr/bin/env python3
"""
T2T-Polish v4.0 — tests/test_runner.py

Unit tests for t2t_polish.runner (CommandResult, run_command, conditional_run).
"""

from __future__ import annotations

import os
import subprocess
from unittest.mock import MagicMock, patch

from t2t_polish.runner import CommandResult, conditional_run, run_command

# ── CommandResult ────────────────────────────────────────────────────


class TestCommandResult:
    def test_fields(self):
        cr = CommandResult(stdout="ok", stderr="", returncode=0, elapsed=1.5)
        assert cr.stdout == "ok"
        assert cr.returncode == 0
        assert cr.elapsed == 1.5

    def test_defaults(self):
        cr = CommandResult(stdout="", stderr="", returncode=0)
        assert cr.elapsed == 0.0


# ── run_command ──────────────────────────────────────────────────────


class TestRunCommand:
    @patch("t2t_polish.runner.subprocess.run")
    def test_success(self, mock_run):
        mock_run.return_value = MagicMock(
            stdout="hello\n", stderr="", returncode=0
        )
        result = run_command(["echo", "hello"], description="echo test")
        assert result.returncode == 0
        assert "hello" in result.stdout

    @patch("t2t_polish.runner.subprocess.run")
    def test_failure_returns_result(self, mock_run):
        err = subprocess.CalledProcessError(1, "fail")
        err.stdout = "out"
        err.stderr = "err"
        mock_run.side_effect = err
        result = run_command(["fail"], description="fail test")
        assert result.returncode == 1

    @patch("t2t_polish.runner.subprocess.run")
    def test_logfile_written(self, mock_run, tmp_path):
        mock_run.return_value = MagicMock(
            stdout="data\n", stderr="", returncode=0
        )
        logfile = str(tmp_path / "cmd.log")
        run_command(["echo"], description="log test", logfile=logfile)
        assert os.path.exists(logfile)
        with open(logfile) as f:
            assert "log test" in f.read()


# ── conditional_run ──────────────────────────────────────────────────


class TestConditionalRun:
    def test_skips_when_output_exists(self, tmp_path):
        outfile = str(tmp_path / "existing.txt")
        with open(outfile, "w") as f:
            f.write("data")
        result = conditional_run(
            outfile=outfile,
            resume_flag=True,
            step_desc="skip test",
            command=["echo", "should not run"],
        )
        assert result is None

    @patch("t2t_polish.runner.run_command")
    def test_runs_when_output_missing(self, mock_rc, tmp_path):
        outfile = str(tmp_path / "missing.txt")
        mock_rc.return_value = CommandResult("ok", "", 0, 0.1)
        result = conditional_run(
            outfile=outfile,
            resume_flag=True,
            step_desc="run test",
            command=["echo"],
        )
        mock_rc.assert_called_once()
        assert result.returncode == 0

    @patch("t2t_polish.runner.subprocess.Popen")
    def test_stream_mode(self, mock_popen, tmp_path):
        outfile = str(tmp_path / "stream.txt")
        proc = MagicMock()
        proc.communicate.return_value = (None, "")
        proc.returncode = 0
        mock_popen.return_value = proc
        conditional_run(
            outfile=outfile,
            resume_flag=False,
            step_desc="stream test",
            command=["cat"],
            stream=True,
        )
        assert os.path.exists(outfile)


# T2T-Polish v4.0
# Any usage is subject to this software's license.
