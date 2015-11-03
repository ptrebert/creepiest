#!/usr/bin/env python3
# coding=utf-8

"""
A stand-alone Python script to rename and cleanup the downloaded
ENCODE files
"""

import sys as sys
import os as os
import operator as op
import collections as col

from pitlib.readers.rowdata import rows_to_buffer

DLFOLDER = '/TL/deep/fhgfs/projects/pebert/crppipe/encode/rawDL'
TRGFOLDER = '/TL/deep/fhgfs/projects/pebert/crppipe/encode/cleanFiles'

MDHG19 = '/home/pebert/work/code/mpggit/creepiest/datasrc/encode/hg19_ENC_select.tsv'
MDMM9 = '/home/pebert/work/code/mpggit/creepiest/datasrc/encode/mm9_ENC_select.tsv'


if __name__ == '__main__':
    allfiles = os.listdir(DLFOLDER)
    hg19info = rows_to_buffer(MDHG19)
    mm9info = rows_to_buffer(MDMM9)
    joined_info = dict()
    select = op.itemgetter(*tuple(['File accession', 'Biosample term name', 'Biosample term id',
                                   'Experiment target', 'Assembly', 'Biological replicate(s)']))
    reorder = op.itemgetter(*tuple([4, 1, 3, 2, 5]))
    for md in [hg19info, mm9info]:
        for entry in md:
            infos = select(entry)
            assert infos[0] not in joined_info, 'Duplicate ENCODE ID: {}'.format(infos[0])
            joined_info[infos[0]] = reorder(infos)
    # now associate files with infos
    # and check for biol replicates
    repcheck = col.defaultdict(list)
    for f in allfiles:
        fn, ext = f.split('.', 1)
        # this must find something
        info = joined_info[fn]
        key = tuple(info[:4] + (ext,))
        val = (fn, info[4])
        repcheck[key].append(val)
    
    for k, v in repcheck.items():
        if len(v) > 1:
            print(k, v)

    sys.exit(0)