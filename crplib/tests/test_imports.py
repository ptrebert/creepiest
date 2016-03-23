# coding=utf-8

"""
Hopefully testing for all necessary imports
"""

import sys as sys
import os as os
import unittest as unittest
import importlib as importlib


class TestModuleImports(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def _crplib_path_heuristic(self):
        """ Try to guess full path to library to collect modules
        :return:
        """
        for p in sys.path:
            if p.endswith('creepiest'):
                return p
        return ''

    def test_python_imports(self):
        test_modules = ['socket', 'multiprocessing', 'queue', 'random', 'pandas',
                        'matplotlib', 'sklearn', 'numpy', 'scipy', 'resource', 'pdb',
                        'pstats', 'statsmodels', 'threading', 'bz2', 'json', 'select',
                        'twobitreader']
        for module in test_modules:
            _ = importlib.import_module(module)
        return

    def test_python_version(self):
        self.assertTrue(sys.version_info >= (3, 4))
        return

    def test_creepiest_imports(self):
        def ferr(err):
            raise err
        import_path = 'crplib.commands.'
        load_path = os.path.join(self._crplib_path_heuristic(), 'crplib', 'commands')
        all_modules = os.listdir(load_path)
        all_modules = [mod for mod in all_modules if mod != '__init__.py' or mod != '__pycache__']
        imp_modules = [import_path + mod.strip('.py') for mod in all_modules]
        for module in imp_modules:
            _ = importlib.import_module(module)
        return
