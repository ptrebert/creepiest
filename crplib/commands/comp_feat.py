# coding=utf-8

"""
Command to compute features for regions
"""

import os as os
import pandas as pd
import multiprocessing as mp

from crplib.auxiliary.file_ops import get_valid_hdf5_groups
from crplib.auxiliary.seq_parsers import add_seq_regions


def print_feature_description():
    """
    :return:
    """


def compute_features(params):
    """
    :param params:
    :return:
    """
    mypid = mp.current_process().pid
    with pd.HDFStore(params['inputfile'], 'r', complevel=9, complib='blosc') as hdf:
        regions = hdf[params['group']].to_dict('records')
    if os.path.isfile(params['addseq']):
        regions = add_seq_regions(regions, params['addseq'], params['chrom'])

    if os.path.isdir(params['tfmotifdb']) and 'tfmotif' in params['features']:
        pass
    if os.path.isfile(params['signal']) and 'covsig' in params['features']:
        pass
    if os.path.isfile(params['pwaln']) and 'mapsig' in params['features']:
        pass

    return mypid, params['group'], regions


def assemble_worker_args(args, datagroups):
    """
    :param args:
    :param datagroups:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfile
    commons['addseq'] = args.addseq
    commons['features'] = args.features
    commons['tfmotifdb'] = args.tfmotifdb
    commons['signal'] = args.signal
    commons['pwaln'] = args.pwaln
    commons['kmers'] = args.kmers
    for grp in datagroups:
        tmp = commons
        tmp['assembly'] = grp.strip('/').split('/')[0]
        tmp['chrom'] = grp.strip('/').split('/')[-1]
        tmp['group'] = grp
        arglist.append(tmp)
    return arglist


def run_compute_features(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    if args.descfeat:
        print_feature_description()
        return 0
    datagroups = get_valid_hdf5_groups(args.inputfile, args.groups)
    assert datagroups, 'No valid groups found in file {} with prefix {}'.format(args.inputfile, args.groups)
    arglist = assemble_worker_args(args, datagroups)
    logger.debug('Start feature computation...')
    with pd.HDFStore(args.outputfile, 'a', complevel=9, complib='blosc') as hdfout:
        with mp.Pool(args.workers) as pool:
            mapres = pool.map_async(compute_features, arglist, chunksize=1)
            for pid, grp, obj in mapres.get():
                pass

    return 0
