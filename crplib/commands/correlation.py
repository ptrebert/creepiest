# coding=utf-8

"""
Module to compute pairwise correlation between signal tracks
"""

import os as os
import multiprocessing as mp
import importlib as imp
import numpy as np
import numpy.ma as msk
import pandas as pd
import json as json
import itertools as itt

from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups, get_trgindex_groups, get_default_group
from crplib.auxiliary.file_ops import create_filepath, text_file_mode


def assemble_worker_params(args):
    """
    :param args:
    :return:
    """
    groups_a = get_valid_hdf5_groups(args.inputfilea, args.inputgroupa)
    groups_b = get_valid_hdf5_groups(args.inputfileb, args.inputgroupb)
    ingrp_a = args.inputgroupa
    if not ingrp_a:
        ingrp_a = get_default_group(args.inputfilea)
    ingrp_b = args.inputgroupb
    if not ingrp_b:
        ingrp_b = get_default_group(args.inputfileb)

    chroms_a = set([os.path.split(g)[1] for g in groups_a])
    chroms_b = set([os.path.split(g)[1] for g in groups_b])
    chrom_union = chroms_a.intersection(chroms_b)
    assert chrom_union, 'No shared chromosomes between input files/groups'
    commons = {'inputfilea': args.inputfilea, 'inputfileb': args.inputfileb,
               'inputgroupa': ingrp_a, 'inputgroupb': ingrp_b,
               'targetindex': args.targetindex, 'measure': args.measure,
               'roifile': args.roifile, 'roilimit': args.roilimit,
               'skipsize': args.skipsize}
    index_groups = get_trgindex_groups(args.targetindex, '')
    arglist = []
    for chrom in chrom_union:
        tmp = dict(commons)
        tmp['chrom'] = chrom
        tmp['loadgroupa'] = os.path.join(ingrp_a, chrom)
        tmp['loadgroupb'] = os.path.join(ingrp_b, chrom)
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
    corrfun = funcs[measure]
    return corrfun


def adapt_mask_to_roi(consmask, roifp, chrom):
    """
    :param consmask:
    :param roifp:
    :param chrom:
    :return:
    """
    roimask = np.ones(consmask.size, dtype=np.bool)
    opn, mode = text_file_mode(roifp)
    with opn(roifp, mode) as rois:
        for line in rois:
            if not line or line.startswith('#'):
                continue
            cols = line.strip().split()
            if cols[0] != chrom:
                continue
            roimask[int(cols[1]):int(cols[2])] = 0
    comb_mask = consmask | roimask  # 0 | 0 = 0
    pos = np.arange(consmask.size)
    pos = msk.array(pos, mask=comb_mask)
    select = np.array([k for k, g in itt.groupby(~comb_mask)], dtype=np.bool)
    splits = np.array(list(itt.chain.from_iterable((item.start, item.stop) for item in msk.clump_unmasked(pos))), dtype=np.int32)
    return splits, select


def compute_corr_cons(params):
    """ Compute
    :param params:
    :return:
    """
    chrom = params['chrom']
    skipsize = params['skipsize']
    # loading index data
    fun_avg = np.vectorize(np.average, otypes=[np.float64])
    with pd.HDFStore(params['targetindex'], 'r') as idx:
        splits = idx[params['splits']].values
        select = idx[params['select']].values
        if params['roilimit']:
            splits, select = adapt_mask_to_roi(idx[params['mask']].values, params['roifile'], chrom)
    with pd.HDFStore(params['inputfilea'], 'r') as hdf1:
        dataset1 = hdf1[params['loadgroupa']].values
        dataset1 = np.compress(select, np.split(dataset1, splits))
        data1_avg = fun_avg([a for a in dataset1 if a.size >= skipsize])
    with pd.HDFStore(params['inputfileb'], 'r') as hdf2:
        dataset2 = hdf2[params['loadgroupb']].values
        dataset2 = np.compress(select, np.split(dataset2, splits))
        data2_avg = fun_avg([a for a in dataset2 if a.size >= skipsize])
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
    results['num_regions'] = data1_avg.size
    return chrom, results


def compute_corr_active(params):
    """ Compute correlation only for active positions, i.e. where
    at least one of the two signal tracks is non-zero
    :param params:
    :return:
    """
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
    return params['chrom'], results


def compute_corr_full(params):
    """
    :param params:
    :return:
    """
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
    return params['chrom'], results


def compute_corr_roi(params):
    """
    :param params:
    :return:
    """
    pass


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
              'roi_file': os.path.basename(args.roifile),
              'roi_limit': args.roilimit,
              'targetindex': 'None' if not args.targetindex else os.path.basename(args.targetindex),
              'group_A': args.inputgroupa, 'group_B': args.inputgroupb,
              'task': args.task, 'measure': args.measure,
              'correlations': []}
    logger.debug('Initializing worker pool')
    with mp.Pool(args.workers) as pool:
        resit = pool.imap_unordered(exec_fun, arglist, chunksize=1)
        for chrom, results in resit:
            logger.debug('Computed correlation for chromosome {}'.format(chrom))
            output['correlations'].append((chrom, results))
    logger.debug('Finished computation')
    fpath = create_filepath(args.outputfile, logger)
    with open(fpath, 'w') as outfile:
        json.dump(output, outfile, indent=1, sort_keys=True)
    logger.debug('Output written to file {}'.format(args.outputfile))
    return 0

