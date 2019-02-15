#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Find and run tests.
"""
import unittest
import os


def run_tests():
    tests = []
    for root, _, _ in os.walk('bruges'):
        if 'test' in root:
            tests.append(
                unittest.TestLoader().discover(root,
                                               pattern='*test.py'))
    suite = unittest.TestSuite(tests)
    return(suite)


if __name__ == '__main__':
    suite = run_tests()
    unittest.TextTestRunner(verbosity=2).run(suite)
