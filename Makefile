# Makefile for SyndrumNET pipeline

.PHONY: all data pipeline evaluate figures test lint api clean appendix help

# Default target
all: data pipeline evaluate figures

# Build all data
data:
	python scripts/build_all_data.py --config configs/default.yaml

# Run pipeline
pipeline:
	python scripts/run_pipeline.py --config configs/default.yaml

# Evaluate predictions
evaluate:
	python scripts/evaluate.py --config configs/default.yaml

# Generate figures
figures:
	python scripts/make_figures.py --config configs/default.yaml

# Run tests
test:
	pytest tests/ -v --cov=syndrumnet --cov-report=html

# Lint (same gate CI runs)
lint:
	ruff check .

# Regenerate the API reference from source
api:
	python scripts/gen_api_docs.py > docs/API.md

# Clean generated files
clean:
	rm -rf data/raw/* data/interim/* data/processed/*
	rm -rf reports/figures/* reports/tables/*
	rm -rf logs/*
	rm -rf htmlcov/ .pytest_cache/ .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

# Generate code appendix (requires pandoc)
appendix:
	@echo "Generating code appendix..."
	@find src/ -name "*.py" -exec echo "# {}" \; -exec cat {} \; > /tmp/syndrumnet_code.txt
	@pandoc /tmp/syndrumnet_code.txt -o Appendix_Code.pdf --pdf-engine=xelatex
	@echo "Generated Appendix_Code.pdf"

# Help
help:
	@echo "SyndrumNET Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all        - Run complete pipeline (data + pipeline + evaluate + figures)"
	@echo "  data       - Download and preprocess all data"
	@echo "  pipeline   - Run synergy prediction pipeline"
	@echo "  evaluate   - Evaluate predictions against known synergies"
	@echo "  figures    - Generate all visualization figures"
	@echo "  test       - Run test suite"
	@echo "  lint       - Run ruff (same gate as CI)"
	@echo "  api        - Regenerate docs/API.md from source"
	@echo "  clean      - Remove generated files"
	@echo "  appendix   - Generate code appendix PDF"
	@echo "  help       - Show this help message"