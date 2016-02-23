# coding=utf-8

import os as os
import datetime as dt
import pandas as pd
import numpy as np

from crplib.auxiliary.constants import DIV_B_TO_MB

_PCT_LEVELS = [0, 5, 25, 50, 75, 95, 100]

_PCT_LABELS = ['pct' + str(p).zfill(3) for p in _PCT_LEVELS]

MD_SIGNAL_COLDEFS = ['group', 'chrom', 'ctime', 'mtime', 'size_mb', 'length', 'srcfile'] + _PCT_LABELS


def gen_obj_and_md(mdframe, rootpath, chrom, srcfiles, datavals):
    """
    :param mdframe:
    :param rootpath:
    :param chrom:
    :param srcfiles:
    :param datavals:
    :return:
    """
    grp = rootpath + '/' + chrom
    pct_scores = np.percentile(datavals, q=_PCT_LEVELS)
    tmpsrc = ','.join([os.path.basename(f) for f in srcfiles])
    ctime = dt.datetime.now()
    dataobj = pd.Series(data=datavals, dtype='float64')
    size_mem = (dataobj.values.nbytes + dataobj.index.nbytes) / DIV_B_TO_MB
    datalen = datavals.size
    mtime = dt.datetime.now()
    entries = [grp, chrom, ctime, mtime, int(size_mem), datalen, tmpsrc] + pct_scores
    tmp = pd.DataFrame([entries, ], columns=MD_SIGNAL_COLDEFS)
    mdframe = mdframe.append(tmp, ignore_index=True)
    return grp, dataobj, mdframe
