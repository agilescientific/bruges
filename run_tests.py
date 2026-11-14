#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run the test suite.

Exits with a non-zero status when any test fails, so CI and pre-commit hooks
can rely on the return code. Uses pytest, which is already a test dependency
and which discovers the `*test.py` modules on every supported Python version
(3.9-3.13); `unittest discover` fails on those directory names before 3.12.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TESTS = ROOT / "bruges"


def main() -> int:
    try:
        import pytest
    except ImportError:  # pragma: no cover
        raise SystemExit("pytest is required: pip install .[test]")
    return int(pytest.main(["-q", str(TESTS)]))


if __name__ == "__main__":
    sys.exit(main())
