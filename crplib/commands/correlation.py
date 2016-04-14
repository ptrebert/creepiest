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

from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups, get_trgindex_groups, load_masked_sigtrack, build_conservation_mask
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
               'targetindex': args.targetindex, 'measure': args.measure}
    index_groups = get_trgindex_groups(args.targetindex, args.indexgroup)
    arglist = []
    for chrom in chrom_union:
        tmp = dict(commons)
        tmp['chrom'] = chrom
        tmp['loadgroupa'] = os.path.join(args.inputgroupa, chrom)
        tmp['loadgroupb'] = os.path.join(args.inputgroupb, chrom)
        tmp.update(index_groups[chrom])
        arglist.append(tmp)
    return arglist


def get_corr_fun(measure, masked):
    """
    :param measure:
    :param masked:
    :return:
    """
    if masked:
        module = imp.import_module('scipy.stats.mstats')
    else:
        module = imp.import_module('scipy.stats')
    funcs = {'pearson': module.pearsonr,
             'spearman': module.spearmanr}
    try:
        corrfun = funcs[measure]
    except KeyError:
        module = imp.import_module('sklearn.metrics')
        corrfun = module.r2_score
    return corrfun


def compute_corr_cons(params):
    """ Compute
    :param params:
    :return:
    """
    mypid = mp.current_process().pid
    # loading index data
    fun_avg = np.vectorize(np.average, otypes=[np.float64])
    with pd.HDFStore(params['targetindex'], 'r') as idx:
        splits = idx[params['splits']].values
        select = idx[params['select']].values
    with pd.HDFStore(params['inputfilea'], 'r') as hdf1:
        dataset1 = hdf1[params['loadgroupa']].values
        data1_avg = fun_avg(np.compress(select, np.split(dataset1, splits)))
    with pd.HDFStore(params['inputfileb'], 'r') as hdf2:
        dataset2 = hdf2[params['loadgroupb']].values
        data2_avg = fun_avg(np.compress(select, np.split(dataset2, splits)))
    results = dict()
    for ms in params['measure']:
        corr_fun = get_corr_fun(ms, masked=False)
        res = corr_fun(data1_avg, data2_avg)
        try:
            corr, pv = res
        except (ValueError, TypeError):
            corr, pv = res, -1
        infos = {'stat': corr, 'pv': pv}
        results[ms] = infos
    return mypid, params['chrom'], results


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
    results = dict()
    for ms in params['measure']:
        corr_fun = get_corr_fun(params[ms], masked=True)
        res = corr_fun(dataset1, dataset2)
        try:
            corr, pv = res
        except (ValueError, TypeError):
            corr, pv = res, -1
        infos = {'stat': corr, 'pv': pv.data}
        results[ms] = infos
    return mypid, params['chrom'], results


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
    results = dict()
    for ms in params['measure']:
        corr_fun = get_corr_fun(params[ms], masked=False)
        res = corr_fun(dataset1, dataset2)
        try:
            corr, pv = res
        except (ValueError, TypeError):
            corr, pv = res, -1
        infos = {'stat': corr, 'pv': pv}
        results[ms] = infos
    return mypid, params['chrom'], results


def run_compute_correlation(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    if args.task == 'cons':
        assert os.path.isfile(args.targetindex), 'No target index specified for task "cons"'
    run_funcs = {'cons': compute_corr_cons, 'full': compute_corr_full, 'active': compute_corr_active}
    exec_fun = run_funcs[args.task]
    logger.debug('Statistics to compute: {}'.format(args.measure))
    arglist = assemble_worker_params(args)
    output = {'file_A': os.path.basename(args.inputfilea),
              'file_B': os.path.basename(args.inputfileb),
              'targetindex': 'None' if not args.targetindex else os.path.basename(args.targetindex),
              'indexgroup': 'None' if not args.indexgroup else args.indexgroup,
              'group_A': args.inputgroupa, 'group_B': args.inputgroupb,
              'task': args.task, 'measure': args.measure,
              'correlations': []}
    logger.debug('Initializing worker pool')
    with mp.Pool(args.workers) as pool:
        mapres = pool.map_async(exec_fun, arglist, chunksize=1)
        for pid, chrom, results in mapres.get():
            logger.debug('Process {} finished correlation for chromosome {}'.format(pid, chrom))
            output['correlations'].append((chrom, results))
    logger.debug('Finished computation')
    with open(args.outputfile, 'w') as outfile:
        json.dump(output, outfile, indent=1, sort_keys=True)
    logger.debug('Output written to file {}'.format(args.outputfile))
    return 0

