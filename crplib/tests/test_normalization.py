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
        self.nz_qnorm = norm_mod.nonzero_qnorm
        self.simple_qnorm = norm_mod.simple_quantile_normalization
        self.ppc_qnorm = norm_mod.preprocess_core_qnorm

        # Example from Wikipedia, 3 samples with 4 observations each
        self.wiki_ex_fxs = np.array([5, 4, 3,
                                     2, 1, 4,
                                     3, 4, 6,
                                     4, 2, 8], dtype=np.float64).reshape(4, 3)
        # Change to sample X feature layout
        self.wiki_ex_sxf = np.transpose(self.wiki_ex_fxs)

        # Expected results for Wikipedia example
        self.wiki_res_fxs = np.array([5.666667, 4.666667, 2,
                                      2, 2, 3,
                                      3, 4.666667, 4.666667,
                                      4.666667, 3, 5.666667], dtype=np.float64).reshape(4, 3)
        self.wiki_res_sxf = np.transpose(self.wiki_res_fxs)

        # Example by Rafael Irizarry, 4 samples with 5 observations each
        self.irizarry_ex_fxs = np.array([2, 4, 4, 5,
                                         5, 14, 4, 7,
                                         4, 8, 6, 9,
                                         3, 8, 5, 8,
                                         3, 9, 3, 5], dtype=np.float64).reshape(5, 4)
        # Change to sample X feature layout
        self.irizarry_ex_sxf = np.transpose(self.irizarry_ex_fxs)

        # Expected results for Irizarry examples
        # For explanation, see crplib.numalg.normalization module
        self.irizarry_code_res_fxs = np.array([3.5, 3.5, 5.25, 4.25,
                                               8.5, 8.5, 5.25, 5.5,
                                               6.5, 5.25, 8.5, 8.5,
                                               5.25, 5.25, 6.5, 6.5,
                                               5.25, 6.5, 3.5, 4.25], dtype=np.float64).reshape(5, 4)
        self.irizarry_code_res_sxf = np.transpose(self.irizarry_code_res_fxs)

        self.irizarry_desc_res_fxs = np.array([3.5, 3.5, 5., 5.,
                                               8.5, 8.5, 5.5, 5.5,
                                               6.5, 5., 8.5, 8.5,
                                               5., 5.5, 6.5, 6.5,
                                               5.5, 6.5, 3.5, 3.5], dtype=np.float64).reshape(5, 4)
        self.irizarry_desc_res_sxf = np.transpose(self.irizarry_desc_res_fxs)

        # Just for robustness testing, some artificial datasets

        self.zero_mat = np.zeros((5, 6), dtype=np.int32)

        self.one_mat = np.ones((5, 6), dtype=np.int32)

    def tearDown(self):
        pass

    def test_simple_qnorm_wiki(self):
        """
        - simple q-norm expects data in "Feature x Sample" layout
        - uses numpy.argsort for ranking, i.e., "first" for handling ties
        :return:
        """
        wiki_res_fxs = self.simple_qnorm(self.wiki_ex_fxs)
        res_shape = wiki_res_fxs.shape
        true_shape = self.wiki_res_fxs.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(wiki_res_fxs, self.wiki_res_fxs, rtol=0.0001)
        return

    def test_simple_qnorm_irizarry_desc(self):
        """
        - simple q-norm expects data in "Feature x Sample" layout
        - uses numpy.argsort for ranking, i.e., "first" for handling ties
        :return:
        """
        irizarry_res_fxs = self.simple_qnorm(self.irizarry_ex_fxs)
        res_shape = irizarry_res_fxs.shape
        true_shape = self.irizarry_desc_res_fxs.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(irizarry_res_fxs, self.irizarry_desc_res_fxs, rtol=0.0001)
        return

    def test_simple_qnorm_irizarry_code(self):
        """
        - simple q-norm expects data in "Feature x Sample" layout
        - uses numpy.argsort for ranking, i.e., "first" for handling ties
        :return:
        """
        irizarry_res_fxs = self.simple_qnorm(self.irizarry_ex_fxs)
        res_shape = irizarry_res_fxs.shape
        true_shape = self.irizarry_code_res_fxs.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(irizarry_res_fxs, self.irizarry_code_res_fxs, rtol=0.0001)
        return

    def test_ppc_qnorm_wiki(self):
        """
        - hard-coded results as produced with R's preprocessCore package
        - expects data in "Feature x Sample" layout
        - presumably, uses "average" for assigning ranks to ties
        :return:
        """
        wiki_res_fxs = self.ppc_qnorm(self.wiki_ex_fxs)
        res_shape = wiki_res_fxs.shape
        true_shape = self.wiki_res_fxs.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(wiki_res_fxs, self.wiki_res_fxs, rtol=0.0001)
        return

    def test_ppc_qnorm_irizarry_desc(self):
        """
        - simple q-norm expects data in "Feature x Sample" layout
        - uses numpy.argsort for ranking, i.e., "first" for handling ties
        :return:
        """
        irizarry_res_fxs = self.ppc_qnorm(self.irizarry_ex_fxs)
        res_shape = irizarry_res_fxs.shape
        true_shape = self.irizarry_desc_res_fxs.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(irizarry_res_fxs, self.irizarry_desc_res_fxs, rtol=0.0001)
        return

    def test_ppc_qnorm_irizarry_code(self):
        """
        - simple q-norm expects data in "Feature x Sample" layout
        - uses numpy.argsort for ranking, i.e., "first" for handling ties
        :return:
        """
        irizarry_res_fxs = self.ppc_qnorm(self.irizarry_ex_fxs)
        res_shape = irizarry_res_fxs.shape
        true_shape = self.irizarry_code_res_fxs.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(irizarry_res_fxs, self.irizarry_code_res_fxs, rtol=0.0001)
        return

    def test_nz_qnorm_wiki(self):
        """
        - implementation based on "min" ranking that ignores all-zero features (columns)
        - expects "Sample x Feature" layout
        :return:
        """
        wiki_res_sxf = self.nz_qnorm(self.wiki_ex_sxf)
        res_shape = wiki_res_sxf.shape
        true_shape = self.wiki_res_sxf.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(wiki_res_sxf, self.wiki_res_sxf, rtol=0.0001)
        return

    def test_nz_qnorm_irizarry_desc(self):
        """
        - implementation based on "min" ranking that ignores all-zero features (columns)
        - expects "Sample x Feature" layout
        :return:
        """
        irizarry_res_sxf = self.nz_qnorm(self.irizarry_ex_sxf)
        res_shape = irizarry_res_sxf.shape
        true_shape = self.irizarry_desc_res_sxf.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(irizarry_res_sxf, self.irizarry_desc_res_sxf, rtol=0.0001)
        return

    def test_nz_qnorm_irizarry_code(self):
        """
        - implementation based on "min" ranking that ignores all-zero features (columns)
        - expects "Sample x Feature" layout
        :return:
        """
        irizarry_res_sxf = self.nz_qnorm(self.irizarry_ex_sxf)
        res_shape = irizarry_res_sxf.shape
        true_shape = self.irizarry_code_res_sxf.shape
        self.assertEqual(res_shape, true_shape,
                         'Shape mismatch: {} vs {}'.format(res_shape, true_shape))
        nptest.assert_allclose(irizarry_res_sxf, self.irizarry_code_res_sxf, rtol=0.0001)
        return
