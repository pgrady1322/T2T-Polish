#!/usr/bin/env python3
"""
T2T-Polish v4.1 — runner.py

Subprocess execution helpers: ``CommandResult``, ``run_command``, and
the resume-aware ``conditional_run`` wrapper.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

from __future__ import annotations

import logging
import os
import subprocess
import time
from dataclasses import dataclass, field

from t2t_polish.exceptions import PipelineStepError

logger = logging.getLogger(__name__)


@dataclass
class CommandResult:
    """Stores the output of a single subprocess invocation."""

    stdout: str
    stderr: str
    returncode: int
    elapsed: float = field(default=0.0)


def run_command(
    cmd: list[str] | str,
    description: str | None = None,
    use_shell: bool = False,
    flush: bool = True,
    logfile: str | None = None,
    max_log_lines: int = 100,
    timeout: int | None = None,
) -> CommandResult:
    """Run *cmd* as a subprocess and return a :class:`CommandResult`.

    Parameters
    ----------
    timeout:
        Maximum wall-clock seconds to wait before killing the process.
        ``None`` (the default) means wait indefinitely.

    On non-zero exit the result is still returned (``returncode != 0``) so that
    callers can decide how to proceed.
    """

    def _truncate(text: str, limit: int | None) -> str:
        if limit is None:
            return text
        lines = text.splitlines()
        if len(lines) <= limit:
            return text
        return "\n".join(lines[:limit] + [f"... ({len(lines) - limit} more lines truncated)"])

    if description:
        logger.info("Running: %s", description)
    logger.debug("Command: %s", cmd if isinstance(cmd, str) else " ".join(cmd))

    start = time.time()
    try:
        result = subprocess.run(
            cmd,
            shell=use_shell,
            check=True,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        elapsed = time.time() - start

        stdout = result.stdout.strip() or "No output"
        stderr = result.stderr.strip() or "No errors"

        logger.info("[STDOUT]\n%s", _truncate(stdout, max_log_lines))
        logger.info("[STDERR]\n%s", _truncate(stderr, max_log_lines))

        if logfile:
            with open(logfile, "a") as f:
                f.write(
                    f"[{description}]\n"
                    f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}\n"
                    f"STDOUT:\n{stdout}\nSTDERR:\n{stderr}\n"
                    f"Elapsed: {elapsed:.2f} seconds\n\n"
                )

        return CommandResult(stdout, stderr, result.returncode, elapsed)

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        logger.error(
            "[Timeout] %s exceeded %s s (elapsed: %.2f s)",
            description if description else cmd,
            timeout,
            elapsed,
        )
        raise PipelineStepError(
            step=description or str(cmd),
            returncode=-9,
            stderr=f"Process timed out after {timeout} seconds",
        ) from None

    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start
        stdout = e.stdout.strip() if e.stdout else "No output"
        stderr = e.stderr.strip() if e.stderr else "No errors"

        logger.error(
            "[Fatal Error] %s (elapsed: %.2f s)",
            description if description else cmd,
            elapsed,
        )
        logger.error("[STDOUT]\n%s", _truncate(stdout, max_log_lines))
        logger.error("[STDERR]\n%s", _truncate(stderr, max_log_lines))

        return CommandResult(stdout, stderr, e.returncode, elapsed)


def conditional_run(
    outfile: str,
    resume_flag: bool,
    step_desc: str,
    command: list[str],
    stream: bool = False,
    timeout: int | None = None,
    use_shell: bool = False,
    flush: bool = True,
    logfile: str | None = None,
    max_log_lines: int = 100,
) -> CommandResult | None:
    """Run *command* only when *outfile* is missing (resume-aware).

    If *stream* is ``True`` the command's stdout is piped directly into
    *outfile* via :class:`subprocess.Popen`.

    Parameters
    ----------
    timeout:
        Maximum wall-clock seconds to wait.  Applies to both stream and
        non-stream modes.

    Raises :class:`PipelineStepError` on non-zero exit or timeout.
    """
    if resume_flag and os.path.exists(outfile):
        logger.info("Resume: skipping %s (exists: %s)", step_desc, outfile)
        return None

    if stream:
        logger.info("Running: %s (streaming)", step_desc)
        with open(outfile, "w") as fh:
            proc = subprocess.Popen(
                command,
                stdout=fh,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            try:
                _, stderr = proc.communicate(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.communicate()
                raise PipelineStepError(
                    step=step_desc,
                    returncode=-9,
                    stderr=f"Process timed out after {timeout} seconds",
                ) from None
            if proc.returncode != 0:
                logger.error("[STDERR]\n%s", stderr.strip() or "No stderr")
                raise PipelineStepError(
                    step=step_desc,
                    returncode=proc.returncode,
                    stderr=stderr.strip() if stderr else "",
                )
            logger.info("Completed: %s (streamed to %s)", step_desc, outfile)
        return None

    result = run_command(
        command,
        description=step_desc,
        timeout=timeout,
        use_shell=use_shell,
        flush=flush,
        logfile=logfile,
        max_log_lines=max_log_lines,
    )
    if result.returncode != 0:
        raise PipelineStepError(
            step=step_desc,
            returncode=result.returncode,
            stderr=result.stderr,
        )
    return result


# T2T-Polish v4.1
# Any usage is subject to this software's license.
