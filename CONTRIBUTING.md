# Contributing to T2T-Polish

Thank you for your interest in contributing to T2T-Polish!  This guide
covers the development setup, code style, testing, and pull request process.

## Development Setup

```bash
# Clone the repository
git clone https://github.com/arangrhie/T2T-Polish.git
cd T2T-Polish

# Create a Python virtual environment (3.10+)
python3 -m venv .venv
source .venv/bin/activate

# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# (Optional) Install pre-commit hooks
pip install pre-commit
pre-commit install
```

## Code Style

- **Formatter / linter:** [Ruff](https://docs.astral.sh/ruff/) — configured
  in `pyproject.toml`.
- **Type checking:** [mypy](https://mypy-lang.org/) with
  `ignore_missing_imports = true`.
- **Line length:** 100 characters.
- **Imports:** sorted by `isort` rules (handled by Ruff `I` rules).
- **Docstrings:** use triple-double quotes.  Public functions should have at
  least a one-line summary.

Run the full lint & format check locally:

```bash
make lint        # ruff check + mypy
make format      # ruff format --check
```

Or use the combined CI target:

```bash
make ci          # lint + format + test (matches GitHub Actions)
```

## Testing

Tests live in `tests/` and use [pytest](https://docs.pytest.org/):

```bash
# Run all tests with verbose output
make test

# Run with coverage report
pytest --cov=t2t_polish --cov-report=term-missing tests/

# Run a single test file
pytest tests/test_runner.py -v
```

### Writing Tests

- Place new tests in `tests/test_<module>.py`.
- Use `pytest.raises` with the appropriate **custom exception** from
  `t2t_polish.exceptions` — library code does **not** call `sys.exit()`.
- Mock external tools via `unittest.mock.patch` on the runner functions
  (`run_command`, `conditional_run`) rather than patching `subprocess.run`
  directly in each module.
- Use the shared fixtures in `conftest.py` (`tmp_dir`, `mock_fasta`,
  `all_tools_on_path`, etc.).

## Pull Request Process

1. **Fork** the repository and create a feature branch:
   `git checkout -b feature/my-change`
2. Make your changes with clear, atomic commits.
3. Run `make ci` locally to ensure lint, format, and tests pass.
4. Open a Pull Request against the `main` branch.
5. Fill in the PR template — describe *what* changed and *why*.
6. A maintainer will review and may request changes.

### Commit Messages

Use conventional-style messages when possible:

```
feat: add --timeout flag to runner
fix: correct step numbering in kcov
docs: update README quick-start section
test: add PipelineStepError assertion tests
```

## Reporting Issues

Open an issue on GitHub with:

- A clear title and description.
- Steps to reproduce (command, inputs, expected vs. actual output).
- Environment details: OS, Python version, conda env, container image.

## License

T2T-Polish is a **U.S. Government Work** (NHGRI/NIH) and is in the
**public domain** within the United States.  By contributing, you agree that
your contributions are also released into the public domain.
