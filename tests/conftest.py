#!/usr/bin/env python3
"""
T2T-Polish v4.0 — tests/conftest.py

Shared pytest fixtures for the T2T-Polish test suite.
"""

from __future__ import annotations

import textwrap
from unittest.mock import patch

import pytest


@pytest.fixture()
def tmp_dir(tmp_path):
    """Return a temporary directory *Path* that is cleaned up automatically."""
    return tmp_path


@pytest.fixture()
def mock_fasta(tmp_path):
    """Write a small two-record FASTA and return its path."""
    fasta = tmp_path / "test.fasta"
    fasta.write_text(
        textwrap.dedent("""\
        >seq1
        ACGTACGTACGTACGT
        >seq2
        TTTTAAAACCCCGGGG
        """)
    )
    return str(fasta)


@pytest.fixture()
def mock_fasta_empty_record(tmp_path):
    """FASTA with one record that has no sequence lines."""
    fasta = tmp_path / "empty_record.fasta"
    fasta.write_text(">empty_seq\n>seq2\nACGT\n")
    return str(fasta)


@pytest.fixture()
def mock_fastq(tmp_path):
    """Write a small FASTQ and return its path."""
    fq = tmp_path / "test.fastq"
    fq.write_text(
        textwrap.dedent("""\
        @read1
        ACGTACGT
        +
        IIIIIIII
        """)
    )
    return str(fq)


@pytest.fixture()
def all_tools_on_path():
    """Patch ``shutil.which`` so every tool appears to be installed."""
    with patch("shutil.which", return_value="/usr/local/bin/fake_tool"):
        yield


@pytest.fixture()
def no_tools_on_path():
    """Patch ``shutil.which`` so every tool appears missing."""
    with patch("shutil.which", return_value=None):
        yield


# T2T-Polish v4.1
# Any usage is subject to this software's license.
