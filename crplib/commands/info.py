# coding=utf-8

"""
Small module to inspect HDF format files
"""

import pandas as pd


def run_print_info(args):
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
