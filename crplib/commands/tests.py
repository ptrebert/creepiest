# coding=utf-8
"""
Unit tests module
"""

import os as os
import sys as sys
import unittest as unittest


def identify_toplevel_dir(path):
    """
    :param path:
    :return:
    """
    src = path
    while not path.endswith('/creepiest'):
        if path == '/':
            raise EnvironmentError('Could not identify library path starting from {}'.format(src))
        path, _ = os.path.split(path)
    assert path in sys.path, 'CREEPIEST path not in PYTHONPATH: {}'.format(sys.path)
    return path


def run_tests(args, modpath):
    """ To ease automatic test discovery, receives full path to module, e.g.
    /root/other/path/creepiest/crplib/commands/tests.py
    :param args:
    :param modpath: full path to tests module (in commands)
    :return:
    """
    logger = args.module_logger
    path = identify_toplevel_dir(modpath)
    logger.debug('Starting test discovery, top level: {}'.format(path))
    loader = unittest.TestLoader()
    testsuite = loader.discover(start_dir=os.path.join(path, 'crplib', 'tests'),
                                pattern='test_*.py', top_level_dir=path)
    logger.debug('Discovery finished, found {} tests'.format(testsuite.countTestCases()))
    test_runner = unittest.TextTestRunner(stream=sys.stderr, verbosity=2)
    test_runner.run(testsuite)
    logger.debug('Tests done')
    return 0
