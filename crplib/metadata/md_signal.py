# coding=utf-8

import os as os
import datetime as dt
import pandas as pd
import numpy as np

from crplib.auxiliary.constants import DIV_B_TO_MB
from crplib.metadata.md_helpers import normalize_group_path, \
    update_metadata_index, normalize_chrom_name

MD_SIGNAL_COLDEFS = ['group', 'chrom', 'mtime', 'size_mb', 'length', 'signaltype', 'srcfile']


def gen_obj_and_md(mdframe, group, chrom, srcfiles, datavals):
    """
    :param mdframe:
    :param group:
    :param chrom:
    :param srcfiles:
    :param datavals:
    :return:
    """
    chrom = normalize_chrom_name(chrom)
    group = normalize_group_path(group, chrom)
    if isinstance(srcfiles, (list, tuple)):
        tmpsrc = ','.join([os.path.basename(f) for f in srcfiles])
    elif isinstance(srcfiles, str):
        tmpsrc = os.path.basename(srcfiles)
    else:
        raise TypeError('Cannot handle source file information: {}'.format(srcfiles))
    mtime = dt.datetime.now()
    if np.issubdtype(datavals.dtype, np.float):
        dataobj = pd.Series(datavals, dtype=np.float64)
        sigtype = 'values'
    elif np.issubdtype(datavals.dtype, np.integer):
        dataobj = pd.Series(datavals, dtype=np.int32)
        sigtype = 'pctranks'
    else:
        raise ValueError('Cannot handle data type of converted object: {}'.format(datavals.dtype))
    size_mem = dataobj.nbytes / DIV_B_TO_MB
    datalen = dataobj.size
    entries = [group, chrom, mtime, int(size_mem), int(datalen), sigtype, tmpsrc]
    upd_idx = update_metadata_index(mdframe, group)
    if upd_idx is not None:
        mdframe.iloc[upd_idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_SIGNAL_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return group, dataobj, mdframe
