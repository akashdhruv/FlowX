"""Run the test suite."""

import sys
import unittest

tests = ['grid.grid', 'poisson.jacobi', 'poisson.SOR','poisson.gauss']

suite = unittest.TestSuite()

for test in tests:
    suite.addTest(unittest.defaultTestLoader.loadTestsFromName(test))

res = unittest.TextTestRunner().run(suite).wasSuccessful()
rc = (0 if res else 1)
sys.exit(rc)
