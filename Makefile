.PHONY: help lint format test install clean

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

