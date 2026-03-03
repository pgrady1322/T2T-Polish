#!/usr/bin/env python3
"""
T2T-Polish v4.1 — exceptions.py

Custom exception hierarchy.  Only ``cli.main()`` should call ``sys.exit()``;
all library code raises one of the exceptions defined here.

Author: Patrick Grady
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

from __future__ import annotations


class T2TPolishError(Exception):
    """Base exception for all T2T-Polish errors."""


class ToolNotFoundError(T2TPolishError):
    """One or more required external tools are missing from ``$PATH``."""


class PipelineStepError(T2TPolishError):
    """A pipeline step (subprocess) exited with a non-zero return code."""

    def __init__(self, step: str, returncode: int, stderr: str = "") -> None:
        self.step = step
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(
            f"Step '{step}' failed with exit code {returncode}"
            + (f": {stderr[:300]}" if stderr else "")
        )


class InputValidationError(T2TPolishError):
    """A required input file is missing, empty, or invalid."""


# T2T-Polish v4.1
# Any usage is subject to this software's license.
