#!/usr/bin/env python3
"""
T2T-Polish v4.1
Fully automatic, iterative K-mer based polishing of genome assemblies.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

__version__ = "4.1.2"

from t2t_polish.constants import (  # noqa: F401
    DEFAULT_CORRECTOR_Q,
    DEFAULT_KMER_SIZE,
    TOOL_NAMES,
)
from t2t_polish.exceptions import (  # noqa: F401
    InputValidationError,
    PipelineStepError,
    T2TPolishError,
    ToolNotFoundError,
)

# T2T-Polish v4.1
# Any usage is subject to this software's license.
