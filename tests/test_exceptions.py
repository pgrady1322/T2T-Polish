#!/usr/bin/env python3
"""
T2T-Polish v4.1 — tests/test_exceptions.py

Unit tests for the custom exception hierarchy in t2t_polish.exceptions.
"""

from __future__ import annotations

import pytest

from t2t_polish.exceptions import (
    InputValidationError,
    PipelineStepError,
    T2TPolishError,
    ToolNotFoundError,
)


class TestExceptionHierarchy:
    def test_base_is_exception(self):
        assert issubclass(T2TPolishError, Exception)

    def test_tool_not_found_inherits(self):
        assert issubclass(ToolNotFoundError, T2TPolishError)

    def test_pipeline_step_inherits(self):
        assert issubclass(PipelineStepError, T2TPolishError)

    def test_input_validation_inherits(self):
        assert issubclass(InputValidationError, T2TPolishError)


class TestPipelineStepError:
    def test_fields(self):
        err = PipelineStepError(step="Step 1", returncode=42, stderr="bad stuff")
        assert err.step == "Step 1"
        assert err.returncode == 42
        assert err.stderr == "bad stuff"

    def test_str_includes_step(self):
        err = PipelineStepError(step="Meryl", returncode=1)
        assert "Meryl" in str(err)
        assert "exit code 1" in str(err)

    def test_str_includes_stderr(self):
        err = PipelineStepError(step="x", returncode=2, stderr="disk full")
        assert "disk full" in str(err)

    def test_catch_as_base(self):
        with pytest.raises(T2TPolishError):
            raise PipelineStepError(step="test", returncode=1)


class TestToolNotFoundError:
    def test_message(self):
        err = ToolNotFoundError("merfin not found")
        assert "merfin" in str(err)

    def test_catch_as_base(self):
        with pytest.raises(T2TPolishError):
            raise ToolNotFoundError("missing tool")


class TestInputValidationError:
    def test_message(self):
        err = InputValidationError("bad input file")
        assert "bad input" in str(err)


# T2T-Polish v4.1
# Any usage is subject to this software's license.
