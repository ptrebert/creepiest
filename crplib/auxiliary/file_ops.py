# coding=utf-8

"""
Convenience module for file operations
"""

import os as os
import numpy as np
import gzip as gz
import bz2 as bz
import pandas as pd


from crplib.auxiliary.text_parsers import get_chain_iterator


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
    # this is necessary since hdf.keys() returns the groups
    # always starting with a '/' though access-wise,
    # '/root' and 'root' are equivalent
    if not prefix.startswith('/'):
        prefix = '/' + prefix
    with pd.HDFStore(filepath, 'r') as hdf:
        groups = [grp for grp in hdf.keys() if grp.startswith(prefix)]
    return groups


def build_conservation_mask(chainfile, chrom, csize):
    """ Build a mask that indicates
    1: is not conserved (= is masked)
    0: is conserved (= is not masked)
    :param chainfile:
    :param chrom:
    :param csize:
    :return:
    """
    mask = np.ones(csize, dtype=np.bool)
    opn, mode = text_file_mode(chainfile)
    with opn(chainfile, mode) as cf:
        chainit = get_chain_iterator(cf, select=chrom)
        for aln in chainit:
            mask[aln[1]:aln[2] + 1] = 0
    return mask


def load_masked_sigtrack(hdfsig, chainfile, group, chrom, csize):
    """
    :return:
    """
    mask = build_conservation_mask(chainfile, chrom, csize)
    with pd.HDFStore(hdfsig, 'r') as hdf:
        load_group = os.path.join(group, chrom)
        signal = np.ma.array(hdf[load_group].values, mask=mask)
    return signal
