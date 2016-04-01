# coding=utf-8

import os as os
import datetime as dt
import pandas as pd

from crplib.auxiliary.constants import DIV_B_TO_MB

MD_REGION_COLDEFS = ['group', 'chrom', 'mtime', 'size_mb', 'numreg', 'covbp', 'srcfile']


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
    if isinstance(srcfiles, (list, tuple)):
        tmpsrc = ','.join([os.path.basename(f) for f in srcfiles])
    elif isinstance(srcfiles, str):
        tmpsrc = os.path.basename(srcfiles)
    else:
        raise TypeError('Cannot handle source file references: {}'.format(srcfiles))
    numreg = len(datavals)
    datavals = [(reg[1], reg[2], reg[3]) for reg in datavals]
    dataobj = pd.DataFrame(data=datavals, columns=['start', 'end', 'name'])
    mtime = dt.datetime.now()
    covbp = (dataobj.end - dataobj.start).sum()
    size_mem = (dataobj.values.nbytes + dataobj.index.nbytes) / DIV_B_TO_MB
    entries = [grp, chrom, mtime, int(size_mem), numreg, covbp, tmpsrc]
    if grp in mdframe.group.values:
        tmp = mdframe.where(mdframe.group == 'grp').dropna().index
        assert len(tmp) == 1, 'Group {} multiple times in metadata'.format(grp)
        idx = tmp[0]
        mdframe.iloc[idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_REGION_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return grp, dataobj, mdframe
