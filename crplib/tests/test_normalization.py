# coding=utf-8

"""
Test method for quantile normalization
"""

import unittest as unittest
import importlib as importlib
import numpy as np
import numpy.testing as nptest


class TestNormalization(unittest.TestCase):
    def setUp(self):
        norm_mod = importlib.import_module('crplib.numalg.normalization')
        self.qnorm = norm_mod.nonzero_qnorm
        self.sample_a = np.array([5, 2, 3, 4])
        self.sample_b = np.array([4, 1, 4, 2])
        self.sample_c = np.array([3, 4, 6, 8])
        self.sample_mat = np.array([self.sample_a, self.sample_b, self.sample_c])

    def tearDown(self):
        pass

    def test_quant_norm(self):

        expect = [np.array([5.66667, 2., 3., 4.66667], dtype=np.float64),
                  np.array([4.66667, 2., 4.66667, 3.], dtype=np.float64),
                  np.array([2., 3., 4.66667, 5.66667], dtype=np.float64)]

        got = self.qnorm(self.sample_mat)

        # the tolerance is set according to expected output for this test case
        # and not based on possible precision for dtype=np.float64
        nptest.assert_allclose(got, expect, rtol=0.0001)
        return
