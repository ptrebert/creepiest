# coding=utf-8

import os as os
import datetime as dt
import pandas as pd

from crplib.auxiliary.constants import DIV_B_TO_MB
from crplib.metadata.md_helpers import normalize_group_path, update_metadata_index, flaterator


MD_REGION_COLDEFS = ['group', 'chrom', 'mtime', 'size_mb', 'numreg', 'covbp', 'srcfile']


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
    srcfiles = list(flaterator(srcfiles))
    srcfiles = ','.join([os.path.basename(fp) for fp in srcfiles])
    if isinstance(datavals, list):
        numreg = len(datavals)
        datavals = [(reg[1], reg[2], reg[3]) for reg in datavals]
        dataobj = pd.DataFrame(data=datavals, columns=['start', 'end', 'name'])
    else:
        dataobj = datavals  # assume it is already a DataFrame
        numreg = dataobj.shape[0]
    mtime = dt.datetime.now()
    covbp = (dataobj.end - dataobj.start)
    assert (covbp > 0).all(), 'Some regions malformed; region length smaller than 1'
    covbp = covbp.sum()
    size_mem = (dataobj.values.nbytes + dataobj.index.nbytes) / DIV_B_TO_MB
    entries = [group, chrom, mtime, int(size_mem), numreg, covbp, srcfiles]
    upd_idx = update_metadata_index(mdframe, group)
    if upd_idx is not None:
        mdframe.iloc[upd_idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_REGION_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return group, dataobj, mdframe
