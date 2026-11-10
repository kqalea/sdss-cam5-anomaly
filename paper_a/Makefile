PY ?= $(if $(wildcard ../.venv/bin/python),../.venv/bin/python,python3)

.PHONY: all check test figures calibrate

all: test check calibrate figures

check:
	$(PY) code/analyse.py --check

calibrate:
	$(PY) code/calibrate.py --check

test:
	$(PY) tests/run_tests.py

figures:
	$(PY) code/figures.py
