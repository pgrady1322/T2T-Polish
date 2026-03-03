#!/usr/bin/env python3
"""
T2T-Polish v4.1 — tests/test_runner.py

Unit tests for t2t_polish.runner (CommandResult, run_command, conditional_run).
"""

from __future__ import annotations

import os
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from t2t_polish.exceptions import PipelineStepError
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

    @patch("t2t_polish.runner.subprocess.run")
    def test_timeout_raises_pipeline_step_error(self, mock_run):
        mock_run.side_effect = subprocess.TimeoutExpired(cmd="slow", timeout=5)
        with pytest.raises(PipelineStepError, match="timed out"):
            run_command(["slow"], description="timeout test", timeout=5)

    @patch("t2t_polish.runner.subprocess.run")
    def test_timeout_passed_to_subprocess(self, mock_run):
        mock_run.return_value = MagicMock(
            stdout="ok\n", stderr="", returncode=0
        )
        run_command(["echo"], description="timeout passthrough", timeout=60)
        _, kwargs = mock_run.call_args
        assert kwargs["timeout"] == 60


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

    @patch("t2t_polish.runner.run_command")
    def test_nonzero_exit_raises_pipeline_step_error(self, mock_rc, tmp_path):
        outfile = str(tmp_path / "fail.txt")
        mock_rc.return_value = CommandResult("", "oops", 1, 0.5)
        with pytest.raises(PipelineStepError, match="failed with exit code 1"):
            conditional_run(
                outfile=outfile,
                resume_flag=False,
                step_desc="fail step",
                command=["false"],
            )

    @patch("t2t_polish.runner.subprocess.Popen")
    def test_stream_nonzero_raises(self, mock_popen, tmp_path):
        outfile = str(tmp_path / "stream_fail.txt")
        proc = MagicMock()
        proc.communicate.return_value = (None, "stream error")
        proc.returncode = 2
        mock_popen.return_value = proc
        with pytest.raises(PipelineStepError, match="exit code 2"):
            conditional_run(
                outfile=outfile,
                resume_flag=False,
                step_desc="stream fail",
                command=["bad"],
                stream=True,
            )

    @patch("t2t_polish.runner.subprocess.Popen")
    def test_stream_timeout_raises(self, mock_popen, tmp_path):
        outfile = str(tmp_path / "stream_timeout.txt")
        proc = MagicMock()
        proc.communicate.side_effect = subprocess.TimeoutExpired(cmd="slow", timeout=3)
        proc.kill.return_value = None
        # After kill, communicate succeeds
        def _comm_after_kill(*a, **kw):
            return (None, "")
        proc.communicate.side_effect = [
            subprocess.TimeoutExpired(cmd="slow", timeout=3),
            (None, ""),
        ]
        mock_popen.return_value = proc
        with pytest.raises(PipelineStepError, match="timed out"):
            conditional_run(
                outfile=outfile,
                resume_flag=False,
                step_desc="stream timeout",
                command=["slow"],
                stream=True,
                timeout=3,
            )


# T2T-Polish v4.1
# Any usage is subject to this software's license.
