# coding=utf-8

import os as os
import datetime as dt
import pandas as pd
import numpy as np

from crplib.auxiliary.constants import DIV_B_TO_MB
from crplib.metadata.md_helpers import normalize_group_path, update_metadata_index

_PCT_LEVELS = [0, 5, 10, 25, 50, 75, 90, 95, 100]

_PCT_LABELS = ['pct' + str(p).zfill(3) for p in _PCT_LEVELS]

MD_SIGNAL_COLDEFS = ['group', 'chrom', 'mtime', 'size_mb', 'length', 'srcfile'] + _PCT_LABELS


def gen_obj_and_md(mdframe, group, chrom, srcfiles, datavals):
    """
    :param mdframe:
    :param group:
    :param chrom:
    :param srcfiles:
    :param datavals:
    :return:
    """
    group = normalize_group_path(group, chrom)
    pct_scores = np.percentile(datavals, q=_PCT_LEVELS)
    if isinstance(srcfiles, (list, tuple)):
        tmpsrc = ','.join([os.path.basename(f) for f in srcfiles])
    elif isinstance(srcfiles, str):
        tmpsrc = os.path.basename(srcfiles)
    else:
        raise TypeError('Cannot handle source file information: {}'.format(srcfiles))
    mtime = dt.datetime.now()
    if isinstance(datavals, np.ndarray) and datavals.dtype == np.float64:
        dataobj = datavals
    else:
        dataobj = np.array(datavals, dtype=np.float64)
    size_mem = dataobj.nbytes / DIV_B_TO_MB
    datalen = dataobj.size
    entries = [group, chrom, mtime, int(size_mem), datalen, tmpsrc]
    entries.extend(list(pct_scores))
    upd_idx = update_metadata_index(mdframe, group)
    if upd_idx is not None:
        mdframe.iloc[upd_idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_SIGNAL_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return group, dataobj, mdframe
