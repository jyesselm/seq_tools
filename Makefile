.PHONY: help lint format test install clean build publish publish-test type-check coverage

help: ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  %-15s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install package in editable mode with dev dependencies
	pip install -e .
	pip install pre-commit
	pre-commit install

lint: ## Run linting checks (ruff)
	ruff check seq_tools/ test/

lint-fix: ## Run linting checks and auto-fix issues
	ruff check --fix seq_tools/ test/

format: ## Format code with ruff
	ruff format seq_tools/ test/

type-check: ## Run type checking with mypy
	mypy seq_tools/

test: ## Run tests
	pytest test/ -v --tb=short

test-fast: ## Run tests quickly (no verbose output)
	pytest test/ --tb=short

coverage: ## Run tests with coverage (minimum 90% required)
	pytest --cov=seq_tools --cov-report=term-missing --cov-report=html --cov-fail-under=90
	@echo "Coverage report generated in htmlcov/index.html"

clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache
	rm -rf .mypy_cache
	rm -rf .ruff_cache
	rm -rf htmlcov/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -r {} +
	find . -type f -name "*.pyc" -delete

check: format lint type-check coverage ## Run all checks (format, lint, type-check, coverage)

ci: format check ## Format code and run all checks (for CI)

build: clean ## Build distribution packages
	python -m pip install --upgrade build
	python -m build

publish-test: build ## Publish to TestPyPI (for testing)
	python -m pip install --upgrade twine
	python -m twine upload --repository testpypi dist/*

publish: build ## Publish to PyPI (production)
	python -m pip install --upgrade twine
	python -m twine upload dist/*

