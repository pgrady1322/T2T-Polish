.PHONY: help install dev lint format typecheck test test-cov ci docker clean

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

# ── Setup ────────────────────────────────────────────────────────────

install:  ## Install package
	pip install -e .

dev:  ## Install with dev dependencies
	pip install -e ".[dev]"

# ── Quality ──────────────────────────────────────────────────────────

lint:  ## Run ruff linter
	ruff check t2t_polish/ tests/

format:  ## Auto-format code
	ruff format t2t_polish/ tests/
	ruff check --fix t2t_polish/ tests/

typecheck:  ## Run mypy type checker
	mypy t2t_polish/ --ignore-missing-imports

# ── Testing ──────────────────────────────────────────────────────────

test:  ## Run tests
	pytest tests/ -v --tb=short

test-cov:  ## Run tests with coverage
	pytest tests/ -v --tb=short --cov=t2t_polish --cov-report=term-missing --cov-fail-under=70

# ── CI ───────────────────────────────────────────────────────────────

ci: lint typecheck test-cov  ## Run lint + typecheck + tests (mirrors CI pipeline)
	@ruff format --check t2t_polish/ tests/

# ── Docker ───────────────────────────────────────────────────────────

docker:  ## Build Docker image
	docker build -t t2t-polish .

# ── Cleanup ──────────────────────────────────────────────────────────

clean:  ## Remove caches and artifacts
	rm -rf __pycache__ .pytest_cache .mypy_cache .ruff_cache
	rm -rf t2t_polish/__pycache__ tests/__pycache__
	find . -name '*.pyc' -delete
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	rm -rf *.egg-info dist build htmlcov
