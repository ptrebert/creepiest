#!/usr/bin/env python3
# coding=utf-8

"""
Stand-alone script to parse the raw ENCODE metadata and create
reasonable sample annotation files
"""

import sys as sys
import operator as op

from pitlib.readers.rowdata import read_row_data


ENCODE_GET = ['File accession', 'Assembly', 'Biosample term name', 'Experiment target', 'Biosample sex', 'Biosample term id', 'File download URL', 'md5sum']

ENCODE_SAMPLES = {'hg19': ['GM12878', 'H1-hESC', 'K562', 'HepG2', 'HUVEC'],
                  'mm9': ['CH12.LX', 'ES-Bruce4', 'ES-E14', 'MEL']}




if __name__ == '__main__':
    mdfile = sys.argv[1]
    outfile = sys.argv[2]
    rows, fd = read_row_data(mdfile)
    itemgetter = op.itemgetter(*tuple(ENCODE_GET))