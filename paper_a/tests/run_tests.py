#!/usr/bin/env python
"""Run the analysis unit tests: plain asserts, no pytest and no data files.

Usage (from the repository root):
  .venv/bin/python analysis/tests/run_tests.py
  .venv/bin/python analysis/tests/run_tests.py test_trend
"""
from __future__ import annotations

import importlib
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULES = ("test_spin_core", "test_significance", "test_trend", "test_sdss_more")


def main(argv: list[str]) -> int:
    names = argv[1:] or list(MODULES)
    sys.path.insert(0, str(HERE))
    sys.path.insert(0, str(HERE.parent / "code"))
    n_pass = n_fail = 0
    failures = []
    t0 = time.perf_counter()
    for name in names:
        module = importlib.import_module(name)
        tests = [(n, f) for n, f in vars(module).items()
                 if n.startswith("test_") and callable(f)]
        for test_name, fn in tests:
            try:
                fn()
            except Exception:
                n_fail += 1
                failures.append((f"{name}.{test_name}", traceback.format_exc()))
                print(f"FAIL {name}.{test_name}")
            else:
                n_pass += 1
                print(f"PASS {name}.{test_name}")
    for label, tb in failures:
        print(f"\n=== {label} ===")
        print(tb.rstrip())
    print(f"\n{n_pass} passed, {n_fail} failed in {time.perf_counter() - t0:.2f} s")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
