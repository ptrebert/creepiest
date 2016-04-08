# coding=utf-8

import os as os
import datetime as dt
import pandas as pd

from crplib.auxiliary.constants import DIV_B_TO_MB
from crplib.metadata.md_helpers import normalize_group_path, update_metadata_index

MD_TRAINDATA_COLDEFS = ['group', 'chrom', 'mtime', 'size_mb', 'numsamples',
                        'resolution', 'features', 'kmers', 'srcfile', 'chainfile']


def gen_obj_and_md(mdframe, group, chrom, args, datavals):
    """
    :param mdframe:
    :param args:
    :param chrom:
    :param datavals:
    :return:
    """
    group = normalize_group_path(group, chrom)
    mtime = dt.datetime.now()
    dataobj = pd.DataFrame.from_dict(datavals)
    numsamples = len(datavals)
    resolution = args.resolution
    features = ','.join(args.features)
    kmers = ','.join(map(str, args.kmers))
    size_mem = (dataobj.values.nbytes + dataobj.index.nbytes) / DIV_B_TO_MB
    entries = [group, chrom, mtime, int(size_mem), numsamples, resolution, features,
               kmers, os.path.basename(args.inputfile), os.path.basename(args.chainfile)]
    upd_idx = update_metadata_index(mdframe, group)
    if upd_idx is not None:
        mdframe.iloc[upd_idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_TRAINDATA_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return group, dataobj, mdframe
