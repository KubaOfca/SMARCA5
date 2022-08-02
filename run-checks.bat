:: Lint and format
call poetry run black smarca5/ tests/
call poetry run pydocstyle seshat/
call poetry run flake8 smarca5/ tests/

:: Type check
call poetry run mypy --disallow-untyped-defs smarca5/
call poetry run mypy --check-untyped-defs tests/

:: Run unit tests and code coverage.
call poetry run coverage run -m pytest tests/
call poetry run coverage report -m
