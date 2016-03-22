# coding=utf-8
"""
Unit tests module
"""

import sys as sys
import unittest as unittest
import logging as logging

from crplib.tests.import_tests import TestModuleImports


def _load_test_suites(logger):
    """
    :param logger:
    :return:
    """
    logger.debug('Loading test suites')
    all_test_suites = []
    import_test_suite = unittest.TestLoader().loadTestsFromTestCase(TestModuleImports)
    all_test_suites.append(import_test_suite)
    logger.debug('Test suites loaded')
    return all_test_suites


def run_tests():
    """
    :return:
    """
    logger = logging.getLogger(__name__)
    all_test_suites = _load_test_suites(logger)
    logger.info('Running tests...')
    test_runner = unittest.TextTestRunner(stream=sys.stdout, verbosity=2)
    list(map(test_runner.run, all_test_suites))
    logger.debug('Tests done')
    return 0