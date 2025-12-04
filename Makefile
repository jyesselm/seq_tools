.PHONY: help lint format test install clean build publish publish-test

help: ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  %-15s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install package in editable mode
	pip install -e .

lint: ## Run linting checks (black check)
	black --check seq_tools/ test/

format: ## Format code with black
	black seq_tools/ test/

test: ## Run tests
	pytest test/ -v --tb=short

test-fast: ## Run tests quickly (no verbose output)
	pytest test/ --tb=short

clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -r {} +
	find . -type f -name "*.pyc" -delete

check: lint test ## Run all checks (lint + test)

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

