# coding=utf-8

"""
Module to encapsulate (currently) the import of the
twobitreader module; may change in the future
"""

import twobitreader as tbr


def get_twobit_seq(fpath, chrom):
    """ Important to remember here that the return value
    is a TwoBitSequence object, allows sliced access
    :param fpath:
    :param chrom:
    :return:
     :rtype: TwoBitSequence
    """
    seqfile = tbr.TwoBitFile(fpath)
    return seqfile[chrom]


def add_seq_regions(regions, seqfile, chrom):
    """ One mild assumption is that the number of regions
    is between a few hundred and a couple of thousand
    :param regions:
    :param seqfile:
    :param chrom:
    :return:
    """
    seqobj = get_twobit_seq(seqfile, chrom)
    for reg in regions:
        reg['seq'] = seqobj[reg['start']:reg['end']]
    return regions
