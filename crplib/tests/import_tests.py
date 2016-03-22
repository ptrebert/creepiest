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

    def _crp_path_heuristic(self):
        """ Try to guess full path to library to collect modules
        :return:
        """
        for p in sys.path:
            if 'crplib' in p:
                return p
        return ''

    def test_builtin_imports(self):
        test_modules = ['socket', 'multiprocessing', 'queue', 'time', 'random',
                        'configparser', 'traceback', 'matplotlib', 'datetime',
                        'sklearn', 'numpy', 'scipy', 'resource', 'pdb', 'cProfile', 'pstats',
                        'statsmodels', 'threading', 'bz2', 'json', 'select']
        for module in test_modules:
            _ = importlib.import_module(module)
        return

    def test_python_version(self):
        sys = importlib.import_module('sys')
        self.assertTrue(sys.version_info >= (3, 4))
        return

    def test_custom_imports(self):
        def ferr(err):
            raise err
        lib_path = os.path.join(self._crp_path_heuristic(), 'lib')
        all_modules = ['lib.crp_config']
        for root, dirs, files in os.walk(lib_path, topdown=True, onerror=ferr, followlinks=False):
            if '.svn' in root or '__pycache__' in root:
                continue
            if set(dirs) <= {'__pycache__', '.svn'}:
                if 'lib/datahandlers' in root:
                    all_modules.extend(['lib.datahandlers.' + f[:-3] for f in files if f != '__init__.py'])
                elif 'lib/datastructs' in root:
                    all_modules.extend(['lib.datastructs.' + f[:-3] for f in files if f != '__init__.py'])
                elif 'lib/mllib' in root:
                    all_modules.extend(['lib.mllib.' + f[:-3] for f in files if f != '__init__.py'])
                elif 'lib/modules' in root:
                    all_modules.extend(['lib.modules.' + f[:-3] for f in files if f != '__init__.py'])
                elif 'lib/commands' in root:
                    all_modules.extend(['lib.commands.' + f[:-3] for f in files if f != '__init__.py'])
                elif 'lib/numlib' in root:
                    all_modules.extend(['lib.numlib.' + f[:-3] for f in files if f != '__init__.py'])
                else:
                    pass
        for module in all_modules:
            _ = importlib.import_module(module)
        return
