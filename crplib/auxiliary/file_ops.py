# coding=utf-8

"""
Convenience module for file operations
"""

import os as os
import gzip as gz
import bz2 as bz
import pandas as pd


def text_file_mode(fpath):
    """ Naive determination of file type and appropriate selection
    of opening function and mode
    :param fpath:
    :return:
    """
    assert os.path.isfile(fpath), 'Invalid path to file: {}'.format(fpath)
    ext = fpath.split('.')[-1].lower()
    # TODO should read magic numbers instead...
    if ext in ['gz', 'gzip']:
        f, m = gz.open, 'rt'  # Python 3.4+
    elif ext in ['bz', 'bz2', 'bzip', 'bzip2']:
        f, m = bz.open, 'rt'
    else:
        f, m = open, 'r'
    return f, m


def get_valid_hdf5_groups(filepath, prefix):
    """
    :param filepath:
    :param prefix:
    :return:
    """
    if not prefix.startswith('/'):
        prefix = '/' + prefix
    groups = []
    with pd.HDFStore(filepath, 'r', complevel=9, complib='blosc') as hdf:
        groups = [grp for grp in hdf.keys() if grp.startswith(prefix)]

    return groups
