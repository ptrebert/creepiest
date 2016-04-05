# coding=utf-8

"""
Convenience module for some operations on HDF files
"""

import os as os
import numpy as np
import pandas as pd

from crplib.auxiliary.file_ops import text_file_mode
from crplib.auxiliary.text_parsers import get_chain_iterator, chromsize_from_chain


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


def build_conservation_mask(chainfile, chrom, csize=None):
    """ Build a mask that indicates
    1: is not conserved (= is masked)
    0: is conserved (= is not masked)
    :param chainfile:
    :param chrom:
    :param csize:
    :return:
    """
    if csize is not None:
        mask = np.ones(csize, dtype=np.bool)
    else:
        chromsize = chromsize_from_chain(chainfile, chrom)
        mask = np.ones(chromsize, dtype=np.bool)
    opn, mode = text_file_mode(chainfile)
    with opn(chainfile, mode) as cf:
        chainit = get_chain_iterator(cf, select=chrom)
        for aln in chainit:
            mask[aln[1]-1:aln[2]] = 0
    return mask


def load_masked_sigtrack(hdfsig, chainfile, group, chrom, csize=None):
    """
    :return:
    """
    mask = build_conservation_mask(chainfile, chrom, csize)
    with pd.HDFStore(hdfsig, 'r') as hdf:
        load_group = os.path.join(group, chrom)
        signal = np.ma.array(hdf[load_group].values, mask=mask)
    return signal
