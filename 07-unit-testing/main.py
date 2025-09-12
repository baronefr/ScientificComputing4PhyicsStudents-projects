import numpy as np

import unittest

# target routine to test: DAXPY (Double-precision AX Plus Y)
def daxpy(a, x, y):
    return a * x + y

# building python unittest
class tests_DAXPY(unittest.TestCase):
    def test_daxpy_small(self):
        # test on a small, controlled example
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([4.0, 5.0, 6.0])
        res = daxpy(2.0, x, y) # 2.0 * x + y
        self.assertTrue( np.allclose(res, [6.0, 9.0, 12.0]) )

    def test_daxpy_zero(self):
        # test with a = 0, should return y
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([4.0, 5.0, 6.0])
        res = daxpy(0.0, x, y)
        self.assertTrue(np.allclose(res, y))

    def test_daxpy_rand(self):
        # test with random data
        test_size = int(1e5)
        rng = np.random.default_rng(42) # fixed seed for reproducibility
        x = rng.standard_normal(test_size)
        y = rng.standard_normal(test_size)
        a = rng.standard_normal()
        res = daxpy(a, x, y)
        self.assertTrue(np.allclose(res, a * x + y))

# NOTE: to run the tests, use the command line:
# python -m unittest main.py
