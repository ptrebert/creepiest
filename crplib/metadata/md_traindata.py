# coding=utf-8

import os as os
import datetime as dt
import pandas as pd

from crplib.auxiliary.constants import DIV_B_TO_MB

MD_TRAINDATA_COLDEFS = ['group', 'chrom', 'ctime', 'mtime', 'size_mb', 'numsamples', 'sigfile', 'chainfile']


def gen_obj_and_md(mdframe, rootpath, chrom, dtype, sigfile, chainfile, datavals):
    """
    :param mdframe:
    :param rootpath:
    :param chrom:
    :param dtype:
    :param sigfile:
    :param chainfile:
    :param datavals:
    :return:
    """
    grp = os.path.join(rootpath, dtype, chrom)
    ctime = dt.datetime.now()
    numreg = len(datavals)
    dataobj = pd.DataFrame.from_dict(datavals)
    mtime = dt.datetime.now()
    numsamples = len(datavals)
    size_mem = (dataobj.values.nbytes + dataobj.index.nbytes) / DIV_B_TO_MB
    entries = [grp, chrom, ctime, mtime, int(size_mem), numsamples,
               os.path.basename(sigfile), os.path.basename(chainfile)]
    if grp in mdframe.group.values:
        tmp = mdframe.where(mdframe.group == 'grp').dropna().index
        assert len(tmp) == 1, 'Group {} multiple times in metadata'.format(grp)
        idx = tmp[0]
        mdframe.iloc[idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_TRAINDATA_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return grp, dataobj, mdframe
