# Changelog

All notable changes to T2T-Polish will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **CI/CD workflows** for automated releases on `v*` tags:
  - **GitHub Release** — auto-generates a changelog-based Release page.
  - **PyPI Publish** — builds sdist + wheel and uploads via OIDC trusted publishing.
  - **GHCR Docker Publish** — builds and pushes the Docker image to
    `ghcr.io/pgrady1322/t2t-polish` with semver tags and layer caching.
- **Concurrency controls** on CI — in-flight runs are cancelled when new
  commits are pushed to the same branch/PR.

## [4.1.0] — 2025-07-18

### Added

- **`PolishConfig` dataclass** — typed, mutable configuration object replaces
  the old pattern of mutating an `argparse.Namespace` with `type: ignore`
  annotations.  All `# noqa: C901` and `# type: ignore[attr-defined]` lines
  in `sub_polish` have been eliminated.
- **Custom exception hierarchy** (`t2t_polish/exceptions.py`):
  `T2TPolishError`, `ToolNotFoundError`, `PipelineStepError`,
  `InputValidationError`.  Library code no longer calls `sys.exit()`.
- **Subprocess timeout support** — `run_command()` and `conditional_run()`
  accept an optional `timeout` parameter (seconds).  On expiry a
  `PipelineStepError` is raised.
- **`--json-summary`** flag on the `polish` subcommand emits a
  machine-readable JSON file with per-iteration QV, completeness, and runtime.
- **`diagnostics` subcommand** checks tool availability, disk space, and
  Python package status without requiring external tools on `$PATH`.
- **`--draft` and `--reads` are now `required=True`** in argparse, so
  missing arguments produce a clear error message instead of a `NoneType`
  traceback deep in the pipeline.
- **CI hardening** — GitHub Actions matrix (Python 3.10 / 3.11 / 3.12),
  ruff format check, mypy enforced, `--cov-fail-under=70`.
- `.pre-commit-config.yaml`, `.github/dependabot.yml`, `.dockerignore`.
- `CONTRIBUTING.md` with development setup, style, and PR guidelines.
- Comprehensive test suite additions: `test_exceptions.py`, timeout tests,
  `PipelineStepError` assertion tests, `PolishConfig` tests, required-arg
  tests for the CLI parser.

### Changed

- **`sub_polish` refactored** into small, testable helper functions
  (`_apply_resume_kcov`, `_apply_config_file`, `_resolve_chosen_peak`,
  `_convert_fasta_reads`, `_validate_inputs`, `_ensure_readmers`,
  `_resolve_seq_type`, `_compute_kcov_if_needed`, `_check_optimized_resume`,
  `_run_iterations`, `_write_summaries`).
- Removed all `redirect_stdout(devnull)` / `redirect_stderr(devnull)`
  wrappers in `polish.py` — Meryl output is already suppressed via
  `max_log_lines=0` in `conditional_run`.
- `evaluation.py` now uses `run_command()` instead of raw `subprocess.run`,
  so all subprocess calls are routed through the runner infrastructure.
- `Dockerfile` and `APv4_Singularity.def` updated to Python 3.12, properly
  copy the `t2t_polish/` package, and `pip install .`.
- JSON summary now reads `__version__` at runtime instead of a hard-coded
  string.

### Fixed

- `--help` and `--version` no longer blocked by `validate_dependencies()`.
- `Dockerfile` and `Singularity` definitions were broken after the v4
  modular refactor (only copied the legacy `APv4.py`).
- `diagnostic_mode()` now checks `pysam` and `tqdm` availability.
- Stale README content: CI badge URL, diagnostics syntax, Python version
  requirement, tool path references.

### Removed

- Dead `make_path()` function from `utils.py`.
- `CKArgs` inner class hack in `sub_polish`.
- All `sys.exit()` calls outside of `cli.py`.

## [4.0.0] — 2025-06-01

### Added

- Full modular rewrite: `t2t_polish/` package with `cli.py`, `runner.py`,
  `polish.py`, `evaluation.py`, `kcov.py`, `utils.py`, `constants.py`.
- PEP 517/518 packaging via `pyproject.toml`.
- `--resume` and `--resume-from` for checkpoint-based restarts.
- `--optimized` auto-computes k-mer coverage via Jellyfish + GenomeScope2.
- Parallel QV evaluation: Merfin + Merqury via `ThreadPoolExecutor`.
- Conda environment file `APv4.yml` with fully pinned dependencies.
- GitHub Actions CI pipeline.
- Initial test suite.

[4.1.0]: https://github.com/arangrhie/T2T-Polish/compare/v4.0.0...v4.1.0
[4.0.0]: https://github.com/arangrhie/T2T-Polish/releases/tag/v4.0.0
