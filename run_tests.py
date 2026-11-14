#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Find and run the test suite.

Exits with a non-zero status when any test fails, so CI and pre-commit hooks
can rely on the return code. Discovery is anchored to this file, so the script
works from any working directory.
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT / "bruges"


def run_tests() -> unittest.TestSuite:
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    test_dirs = sorted(p for p in PACKAGE.rglob("test") if p.is_dir())
    if not test_dirs:
        raise SystemExit(f"no test directories found under {PACKAGE}")
    for test_dir in test_dirs:
        suite.addTests(loader.discover(str(test_dir), pattern="*test.py"))
    return suite


def main() -> int:
    result = unittest.TextTestRunner(verbosity=2).run(run_tests())
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
