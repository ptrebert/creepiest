# coding=utf-8

"""
Convenience module for some operations on HDF files
"""

import os as os
import re as re
import numpy as np
import pandas as pd
import collections as col

from crplib.auxiliary.file_ops import text_file_mode
from crplib.auxiliary.text_parsers import get_chain_iterator, chromsize_from_chain

from crplib.auxiliary.constants import MAPIDX_MASK, MAPIDX_SPLITS, MAPIDX_SELECT, MAPIDX_ORDER


def get_valid_hdf5_groups(filepath, prefix, exclmd=True):
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
    if exclmd:
        groups = [grp for grp in groups if grp != '/metadata']
    return sorted(groups)


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


def get_default_group(filepath):
    """
    :param filepath:
    :return:
    """
    group_root = set()
    with pd.HDFStore(filepath, 'r') as hdf:
        try:
            mdf = hdf['metadata']
        except KeyError:
            # no metadata in file, maybe it's left over from a failed run
            return ''
        for row in mdf.itertuples(index=False):
            group_root.add(os.path.split(row.group)[0])
    assert len(group_root) == 1,\
        'Cannot identify default group in file {} - several groups per chromosome: {}'.format(os.path.basename(filepath), group_root)
    return group_root.pop()


def check_path_infos(filepath, arggroup=None):
    """
    :param filepath:
    :return:
    """
    parts = filepath.split(':')
    label, group, fp = '', '', ''
    if len(parts) == 2:
        label = None
        group = parts[0].strip()
        fp = parts[1].strip()
    elif len(parts) == 3:
        label = parts[0].strip()
        group = parts[1].strip()
        fp = parts[2].strip()
    else:
        assert filepath.count(':') < 3, 'Cannot decode information prepended to file path: {}'.format(filepath)
    # clearly distinguish between no information and a file path that
    # erroneously still contains a colon, e.g., :/path/to/file
    lab = label if label else None
    grp = group if group else None
    fp = fp if fp else filepath
    if arggroup is not None:
        assert grp == arggroup or not grp or not arggroup,\
            'Group mismatch between argument ({}) and path ({})'.format(arggroup, grp)
        if arggroup and not grp:
            grp = arggroup
    if grp is None or grp.lower() in ['default', 'auto']:
        if not os.path.isfile(fp):
            # this can happen for output files that
            # do not yet exist
            grp = ''
        else:
            grp = get_default_group(fp)
    return lab, grp, fp


def get_chrom_list(filepath, verify=False):
    """
    :param filepath:
    :param verify:
    :return:
    """
    with pd.HDFStore(filepath, 'r') as hdf:
        chroms = hdf['metadata'].chrom.tolist()
        if verify:
            for grp in hdf.keys():
                if grp == '/metadata':
                    continue
                root, chrom = os.path.split(grp)
                assert chrom.strip('/') in chroms, \
                    'Could not find data group for chrom {} in file {}'.format(chrom, os.path.basename(filepath))
    return chroms


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
        if MAPIDX_MASK in g:
            infos[chrom]['mask'] = g
        elif MAPIDX_SPLITS in g:
            infos[chrom]['splits'] = g
        elif MAPIDX_SELECT in g:
            infos[chrom]['select'] = g
        elif MAPIDX_ORDER in g:
            infos[chrom]['order'] = g
        else:
            raise ValueError('Unexpected group in target index {}: {}'.format(fpath, g))
    return infos


def load_data_group(filepath, group, chrom='', allow_none=False):
    """
    :param filepath:
    :param group:
    :param chrom:
    :return:
    """
    if chrom:
        if not group.endswith(chrom):
            group = os.path.join(group, chrom)
    with pd.HDFStore(filepath, 'r') as hdf:
        try:
            data_group = hdf[group]
        except KeyError as ke:
            if allow_none:
                data_group = None
            else:
                raise ke
    return data_group


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
    with opn(chainfile, mode=mode, encoding='ascii') as cf:
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


def extract_chromsizes_from_map(mapfile, which, select):
    """
    :param mapfile:
    :param which:
    :param select:
    :return:
    """
    chromselect = re.compile(select)
    chromosomes = dict()
    with pd.HDFStore(mapfile, 'r') as hdf:
        md = hdf['metadata']
        assembly = md.loc[md.key == which, 'value'].str.cat()
        assm_chroms = md.loc[md.key.str.contains(assembly)]
        assert not assm_chroms.empty, 'No chromosomes selected from metadata for assembly {}'.format(assembly)
        for row in assm_chroms.itertuples():
            a, c = os.path.split(row.key)
            assert a == assembly, 'Unexpected key for assembly {}: {} ({})'.format(assembly, a, which)
            if chromselect.match(c) is not None:
                chromosomes[c] = int(row.value)
    assert chromosomes, 'No chromosomes sizes extracted from file: {} - {} - {}'.format(os.path.basename(mapfile),
                                                                                        which,
                                                                                        select)
    return chromosomes
