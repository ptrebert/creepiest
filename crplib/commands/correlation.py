# coding=utf-8

"""
Module to compute pairwise correlation between signal tracks
"""

import os as os
import multiprocessing as mp
import importlib as imp
import numpy as np
import pandas as pd
import json as json

from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups, load_masked_sigtrack, build_conservation_mask
from crplib.numalg.iterators import iter_consecutive_blocks


def assemble_worker_params(args):
    """
    :param args:
    :return:
    """
    groups_a = get_valid_hdf5_groups(args.inputfilea, args.inputgroupa)
    groups_b = get_valid_hdf5_groups(args.inputfileb, args.inputgroupb)

    chroms_a = set([os.path.split(g)[1] for g in groups_a])
    chroms_b = set([os.path.split(g)[1] for g in groups_b])

    chrom_union = chroms_a.intersection(chroms_b)

    assert chrom_union, 'No shared chromosomes between input files/groups'
    commons = {'inputfilea': args.inputfilea, 'inputfileb': args.inputfileb,
               'inputgroupa': args.inputgroupa, 'inputgroupb': args.inputgroupb,
               'chainfile': args.chainfile, 'corrtype': args.corrtype}
    arglist = []
    for chrom in chrom_union:
        tmp = dict(commons)
        tmp['chrom'] = chrom
        arglist.append(tmp)
    return arglist


def get_corr_fun(corrtype, masked):
    """
    :param corrtype:
    :param masked:
    :return:
    """
    if masked:
        module = imp.import_module('scipy.stats.mstats')
    else:
        module = imp.import_module('scipy.stats')
    funcs = {'pearson': module.pearsonr,
             'spearman': module.spearmanr}
    return funcs[corrtype]


def compute_corr_cons(params):
    """ Compute
    :param params:
    :return:
    """
    mypid = mp.current_process().pid
    cons_mask, num_aln = build_conservation_mask(params['chainfile'], params['chrom'])
    dataset1 = load_masked_sigtrack(params['inputfilea'], params['chainfile'],
                                    params['inputgroupa'], params['chrom'], mask=cons_mask)
    dataset2 = load_masked_sigtrack(params['inputfileb'], params['chainfile'],
                                    params['inputgroupb'], params['chrom'], mask=cons_mask)

    data1_avg = np.zeros(num_aln, dtype=np.float64)
    data2_avg = np.zeros(num_aln, dtype=np.float64)
    unmask_pos = np.arange(len(dataset1), dtype=np.int32)[~dataset1.mask]
    for idx, (start, end) in enumerate(iter_consecutive_blocks(unmask_pos)):
        data1_avg[idx] = np.average(dataset1[start:end])
        data2_avg[idx] = np.average(dataset2[start:end])
    corr_fun = get_corr_fun(params['corrtype'], masked=False)
    corr, pv = corr_fun(data1_avg, data2_avg)
    return mypid, params['chrom'], corr, pv


def compute_corr_active(params):
    """ Compute correlation only for active positions, i.e. where
    at least one of the two signal tracks is non-zero
    :param params:
    :return:
    """
    mypid = mp.current_process().pid
    with pd.HDFStore(params['inputfilea'], 'r') as hdf:
        load_group = os.path.join(params['inputgroupa'], params['chrom'])
        data = hdf[load_group].values
        dataset1 = np.ma.masked_where(data > 0, data)
    with pd.HDFStore(params['inputfileb'], 'r') as hdf:
        load_group = os.path.join(params['inputgroupb'], params['chrom'])
        data = hdf[load_group].values
        dataset2 = np.ma.masked_where(data > 0, data)
    comb_mask = np.ma.getmask(dataset1) & np.ma.getmask(dataset2)
    dataset1 = np.ma.array(dataset1.data, mask=comb_mask)
    dataset2 = np.ma.array(dataset2.data, mask=comb_mask)
    corr_fun = get_corr_fun(params['corrtype'], masked=True)
    corr, pv = corr_fun(dataset1, dataset2)
    return mypid, params['chrom'], corr, pv.data


def compute_corr_full(params):
    """
    :param params:
    :return:
    """
    mypid = mp.current_process().pid
    with pd.HDFStore(params['inputfilea'], 'r') as hdf:
        load_group = os.path.join(params['inputgroupa'], params['chrom'])
        dataset1 = hdf[load_group].values
    with pd.HDFStore(params['inputfileb'], 'r') as hdf:
        load_group = os.path.join(params['inputgroupb'], params['chrom'])
        dataset2 = hdf[load_group].values
    corr_fun = get_corr_fun(params['corrtype'], masked=False)
    corr, pv = corr_fun(dataset1, dataset2)
    return mypid, params['chrom'], corr, pv


def run_compute_correlation(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    if args.task == 'cons':
        assert os.path.isfile(args.chainfile), 'No chain file specified for task "cons"'
    run_funcs = {'cons': compute_corr_cons, 'full': compute_corr_full, 'active': compute_corr_active}
    exec_fun = run_funcs[args.task]
    arglist = assemble_worker_params(args)
    output = {'file_A': os.path.basename(args.inputfilea),
              'file_B': os.path.basename(args.inputfileb),
              'chainfile': 'None' if not args.chainfile else os.path.basename(args.chainfile),
              'group_A': args.inputgroupa, 'group_B': args.inputgroupb,
              'task': args.task, 'corrtype': args.corrtype,
              'correlations': []}
    logger.debug('Initializing worker pool')
    with mp.Pool(args.workers) as pool:
        mapres = pool.map_async(exec_fun, arglist, chunksize=1)
        for pid, chrom, corr, pv in mapres.get():
            logger.debug('Process {} finished correlation for chromosome {}'.format(pid, chrom))
            output['correlations'].append((chrom, corr, pv))
    logger.debug('Finished computation')
    with open(args.outputfile, 'w') as outfile:
        json.dump(output, outfile, indent=1, sort_keys=True)
    logger.debug('Output written to file {}'.format(args.outputfile))
    return 0

