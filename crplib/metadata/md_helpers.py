# coding=utf-8

"""
Some small helper functions for dealing with metadata
"""

import os as os


def normalize_group_path(group, suffix=None):
    """
    :param group:
    :param suffix:
    :return:
    """
    group = os.path.join('/', group)
    if suffix is not None:
        if not group.endswith(suffix):
            group = os.path.join(group, suffix.rstrip('/'))
    return group


def normalize_chrom_name(chrom):
    """
    Add 'chr' prefix if not present
    Remove all trailing or inner white spaces
    :param chrom:
    :return:
    """
    if not chrom.startswith('chr'):
        chrom = 'chr' + chrom
    chrom = chrom.strip().replace(' ', '')
    return chrom


def update_metadata_index(mdframe, group):
    """
    :param mdframe:
    :param group:
    :return:
    """
    assert 'group' in mdframe.columns, 'Non-standard metadata dataframe - no "group" column'
    update_index = None
    if group in mdframe.group.values:
        tmp = mdframe.where(mdframe.group == group).dropna().index
        assert len(tmp) == 1, 'Group {} multiple times in metadata'.format(group)
        update_index = tmp[0]
    return update_index


def flaterator(iterable):
    """ This iterator is only intended for small iterables
    (due to the simple buffer construct)
    :param iterable:
    :return:
    """
    buffer = []
    for thing in iterable:
        if isinstance(thing, str):
            yield thing
        elif not hasattr(thing, '__iter__'):
            yield thing
        else:
            buffer.extend(list(flaterator(thing)))
    for item in buffer:
        yield item
    return
