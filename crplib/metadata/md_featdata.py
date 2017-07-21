# coding=utf-8

import os as os
import datetime as dt
import pandas as pd

from crplib.auxiliary.constants import DIV_B_TO_MB
from crplib.metadata.md_helpers import normalize_group_path, update_metadata_index, flaterator

MD_FEATDATA_COLDEFS = ['group', 'chrom', 'mtime', 'size_mb', 'numsamples',
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
    if args.task == 'regress':
        resolution = args.resolution
    elif args.task == 'matched':
        resolution = 'N/A'
    elif args.task == 'classify':
        resolution = args.window
    else:
        raise ValueError('Cannot create metadata for unknown task: {}'.format(args.task))
    features = ','.join(args.features)
    kmers = ','.join(map(str, args.kmers))
    size_mem = (dataobj.values.nbytes + dataobj.index.nbytes) / DIV_B_TO_MB
    srcfiles = []
    for item in flaterator(args.inputfile):
        srcfiles.append(item)
    srcfiles = ','.join(srcfiles)
    if 'msig' in args.features:
        srcfiles += ',' + ','.join([os.path.basename(sf) for sf in args.sigfile])
    if 'roi' in args.features:
        srcfiles += ',' + ','.join([os.path.basename(sf) for sf in args.roifile])
    if 'tfm' in args.features:
        srcfiles += ',' + os.path.basename(args.tfmotifs)
    entries = [group, chrom, mtime, int(size_mem), numsamples, resolution, features,
               kmers, srcfiles, 'n/a' if not args.mapfile else os.path.basename(args.mapfile)]
    upd_idx = update_metadata_index(mdframe, group)
    if upd_idx is not None:
        mdframe.iloc[upd_idx, ] = entries
    else:
        tmp = pd.DataFrame([entries, ], columns=MD_FEATDATA_COLDEFS)
        mdframe = mdframe.append(tmp, ignore_index=True)
    return group, dataobj, mdframe
