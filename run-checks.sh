#!/bin/bash
# Lint and format
poetry run black smarca5/ tests/
poetry run pydocstyle smarca5/
poetry run flake8 smarca5/ tests/

# Type check
poetry run mypy --disallow-untyped-defs smarca5/
poetry run mypy --check-untyped-defs tests/

# Run unit tests and code coverage.
poetry run coverage run -m pytest tests/
poetry run coverage report -m