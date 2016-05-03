# coding=utf-8

"""
Convenience module for some operations on HDF files
"""

import os as os
import numpy as np
import pandas as pd
import collections as col

from crplib.auxiliary.file_ops import text_file_mode
from crplib.auxiliary.text_parsers import get_chain_iterator, chromsize_from_chain

from crplib.auxiliary.constants import TRGIDX_MASK, TRGIDX_SPLITS, TRGIDX_SELECT


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


def get_valid_chrom_group(filepath, chrom):
    """
    :param filepath:
    :param chrom:
    :return:
    """
    all_groups = get_valid_hdf5_groups(filepath, '')
    single_group = list(filter(lambda g: g.endswith(chrom), all_groups))
    assert len(single_group) == 1, 'Could not identify matching group in file {} for {}'.format(filepath, chrom)
    return single_group[0]


def get_trgindex_groups(fpath, grproot):
    """
    :param fpath:
    :param grproot:
    :return:
    """
    groups = get_valid_hdf5_groups(fpath, grproot)
    assert groups, 'No data found in target index {} with group {}'.format(fpath, grproot)
    infos = col.defaultdict(dict)
    for g in groups:
        chrom = os.path.split(g)[1]
        if TRGIDX_MASK in g:
            infos[chrom]['mask'] = g
        elif TRGIDX_SPLITS in g:
            infos[chrom]['splits'] = g
        elif TRGIDX_SELECT in g:
            infos[chrom]['select'] = g
        else:
            raise ValueError('Unexpected group in target index {}: {}'.format(fpath, g))
    return infos


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
    num_aln = 0
    with opn(chainfile, mode) as cf:
        chainit = get_chain_iterator(cf, tselect=chrom)
        for aln in chainit:
            mask[aln[1]:aln[2]] = 0
            num_aln += 1
    return mask, num_aln


def load_masked_sigtrack(hdfsig, chainfile, group, chrom, csize=None, mask=None):
    """
    :param hdfsig:
    :param chainfile:
    :param group:
    :param chrom:
    :param csize:
    :param mask:
    :return:
    """
    if mask is None:
        mask, num_aln = build_conservation_mask(chainfile, chrom, csize)
    with pd.HDFStore(hdfsig, 'r') as hdf:
        if group.endswith(chrom):
            load_group = group
        else:
            load_group = os.path.join(group, chrom)
        signal = np.ma.array(hdf[load_group].values, mask=mask)
    return signal
