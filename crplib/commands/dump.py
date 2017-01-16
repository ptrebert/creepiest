# coding=utf-8

"""
Small module to dump data from HDF format files into text formats (BED-like)
"""

import numpy as np
import pandas as pd


def dump_index():
    """
    :return:
    """


def dump_signal():
    """
    :return:
    """
    stat = {'mean': np.mean, 'max': np.max, 'min': np.min,
            'median': np.median, 'sum': np.sum, 'product': np.product}



def run_dump_data(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    logger.debug('Opening HDF file {}'.format(args.inputfile))
    with pd.HDFStore(args.inputfile, 'r') as hdf:
        print('=== General HDF structure')
        print(hdf)
        print('<<<')
        print('=== Groups in file')
        print(list(sorted(hdf.keys())))
        print('<<<')
        if '/metadata' in hdf.keys():
            print('=== Metadata in file')
            print(hdf['metadata'])
            print('<<<')
    logger.debug('Info done')
    return 0
