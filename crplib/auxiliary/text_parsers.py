# coding=utf-8

"""
Some helper functions to parse text-based files
with more or less no well-defined format
"""

import re as re

from crplib.auxiliary.file_ops import text_file_mode


def read_chromosome_sizes(fpath, keep='.+'):
    """
    :param fpath:
    :param keep:
    :return:
    """
    chroms = dict()
    keeper = re.compile(keep)
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        for line in infile:
            if not line.strip():
                continue
            cols = line.strip().split()
            cname = cols[0].strip()
            m = keeper.match(cname)
            if m is not None:
                csize = int(cols[1])
                chroms[cname] = csize
    assert chroms, 'No chromosomes from file {} selected with pattern {}'.format(fpath, keep)
    return chroms
