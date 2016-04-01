# coding=utf-8

import os as os
import datetime as dt
import pandas as pd

from crplib.auxiliary.constants import DIV_B_TO_MB

MD_TRAINDATA_COLDEFS = ['group', 'chrom', 'mtime', 'size_mb', 'numsamples',
                        'resolution', 'features', 'kmers', 'sigfile', 'chainfile']


def gen_obj_and_md(mdframe, args, chrom, datavals):
    """
    :param mdframe:
    :param args:
    :param chrom:
    :param datavals:
    :return:
    """
    grp = os.path.join(args.grouproot, chrom)
    mtime = dt.datetime.now()
    dataobj = pd.DataFrame.from_dict(datavals)
    numsamples = len(datavals)
    resolution = args.resolution
    features = ','.join(args.features)
    kmers = ','.join(map(str, args.kmers))
    size_mem = (dataobj.values.nbytes + dataobj.index.nbytes) / DIV_B_TO_MB
    entries = [grp, chrom, mtime, int(size_mem), numsamples, resolution, features,
               kmers, os.path.basename(args.inputfile), os.path.basename(args.chainfile)]
    if grp in mdframe.group.values:
        tmp = mdframe.where(mdframe.group == 'grp').dropna().index
        assert len(tmp) == 1, 'Group {} multiple times in metadata'.format(grp)
        idx = tmp[0]
        mdframe.iloc[idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_TRAINDATA_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return grp, dataobj, mdframe
