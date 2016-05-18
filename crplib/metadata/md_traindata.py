# coding=utf-8

import os as os
import datetime as dt
import pandas as pd

from crplib.auxiliary.constants import DIV_B_TO_MB
from crplib.metadata.md_helpers import normalize_group_path, update_metadata_index

MD_TRAINDATA_COLDEFS = ['group', 'chrom', 'mtime', 'size_mb', 'numsamples',
                        'resolution', 'features', 'kmers', 'srcfile', 'indexfile']


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
    if isinstance(datavals, list):
        dataobj = pd.DataFrame.from_dict(datavals)
        numsamples = len(datavals)
    else:
        dataobj = datavals
        numsamples = dataobj.shape[0]
    if args.task == 'regsig':
        resolution = args.resolution
    elif args.task == 'clsreg':
        resolution = 'N/A'
    else:
        raise ValueError('Cannot create metadata for unknown task: {}'.format(args.task))
    features = ','.join(args.features)
    kmers = ','.join(map(str, args.kmers))
    size_mem = (dataobj.values.nbytes + dataobj.index.nbytes) / DIV_B_TO_MB
    if 'msig' in args.features:
        srcfiles = os.path.basename(args.inputfile) + ',' + os.path.basename(args.signalfile)
    else:
        srcfiles = os.path.basename(args.inputfile)
    entries = [group, chrom, mtime, int(size_mem), numsamples, resolution, features,
               kmers, srcfiles, os.path.basename(args.targetindex)]
    upd_idx = update_metadata_index(mdframe, group)
    if upd_idx is not None:
        mdframe.iloc[upd_idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_TRAINDATA_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return group, dataobj, mdframe
