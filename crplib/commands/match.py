# coding=utf-8

"""
Module to search for matching background regions
Create negative set for learning algorithm
"""

import os as os
import pandas as pd
import multiprocessing as mp
import operator as op
import numpy as np
import numpy.random as rng
import threading as thd
import scipy.spatial as spat

from crplib.auxiliary.seq_parsers import add_seq_regions, get_twobit_seq
from crplib.mlfeat.featdef import get_online_version
from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups
from crplib.auxiliary.constants import CHROMOSOME_BOUNDARY
from crplib.metadata.md_regions import gen_obj_and_md, MD_REGION_COLDEFS


def get_feature_selector(region):
    """
    :param region:
    :return:
    """
    feat = []
    for entry in region.keys():
        if entry.startswith('ft'):
            feat.append(entry)
    feat = sorted(feat)
    return feat, len(feat)


def start_relaxed_search(nntree, chromseq, mask, compfeat, featselect, limits, params):
    """
    :param nntree:
    :param chromseq:
    :param compfeat:
    :param limits:
    :param params:
    :return:
    """
    minlen, maxlen, numsamples = limits
    chrom = params['chrom']
    found = set()
    matched = []
    relax = params['relax_init']
    knn = 1
    timeout = thd.Event()
    if params['timeout'] > 0:
        timer = thd.Timer(params['timeout'] * 60, timeout.set)
        timer.start()
    while not timeout.is_set() and len(found) != numsamples:
        free_pos = np.arange(len(mask))[~mask]
        num_candidates = int((numsamples - len(found)) * 1.2)
        rand_pos = rng.choice(free_pos, num_candidates, replace=False)
        for pos in rand_pos:
            halflen = rng.randint(minlen, maxlen)
            s, e = pos - halflen, pos + halflen
            if s < 0 or e > len(chromseq):  # might happen if looking for matches for large peaks
                continue
            if mask[s:e].any():
                continue
            candseq = chromseq[s:e]
            if candseq.count('N') > (e - s) * 0.05:
                continue
            cand = {'start': s, 'end': e, 'seq': candseq}
            cand = compfeat(cand)
            neighbors = nntree.query(featselect(cand), k=knn, eps=0, p=np.inf)
            if knn > 1:
                for dist, index in zip(neighbors[0], neighbors[1]):
                    if dist <= relax and index not in found:
                        found.add(index)
                        mask[s:e] = 1
                        cand['index'] = index
                        cand['name'] = 'match_' + chrom + '_idx' + str(index)
                        matched.append(cand)
                        break
            else:
                if neighbors[0] <= relax and neighbors[1] not in found:
                    found.add(neighbors[1])
                    mask[s:e] = 1
                    cand['index'] = neighbors[1]
                    cand['name'] = 'match_' + chrom + '_idx' + str(neighbors[1])
                    matched.append(cand)
        knn += 1
        relax += params['relax_step']
        if relax > params['relax_limit']:
            knn = 2
            relax = params['relax_init']
    return matched


def find_background_regions(params):
    """
    :param params:
    :return:
    """
    mult_hi = 1 + params['relax_limit'] / 100.
    mult_lo = 1 - params['relax_limit'] / 100.
    with pd.HDFStore(params['inputfile'], 'r') as hdfin:
        fg = hdfin[params['group_fg']]
        assert not fg.empty, 'Loaded empty group {} from file {}'.format(params['group_fg'], params['inputfile'])
        maxlen = (fg.end - fg.start).max()
        yardstick = int(maxlen * mult_hi)
        minlen = (fg.end - fg.start).min() * mult_lo
        idx = fg.index
        fg = fg.to_dict(orient='record')
        for i, rec in zip(idx, fg):
            rec['index'] = i
        maxlen = int(maxlen / 2.)
        minlen = int(minlen / 2.)
    fg = add_seq_regions(fg, params['seqfile'], params['chrom'])
    compfeat = get_online_version(params['features'], params['kmers'], yardstick)
    fg = map(compfeat, fg)
    fg = sorted(fg, key=lambda d: d['index'])
    featnames, numfeat = get_feature_selector(fg[0])
    getvals = op.itemgetter(*tuple(featnames))
    featmatrix = np.array([list(getvals(d)) for d in fg], dtype=np.float64)
    kdtree = spat.KDTree(featmatrix)
    chromseq = get_twobit_seq(params['seqfile'], params['chrom'])
    mask = np.zeros(len(chromseq), dtype=np.bool)
    for rec in fg:
        mask[rec['start']:rec['end']] = 1
    mask[0:CHROMOSOME_BOUNDARY] = 1
    limits = minlen, maxlen, len(fg)
    matched = start_relaxed_search(kdtree, chromseq, mask, compfeat, getvals, limits, params)
    if not matched and not params['allownomatch']:
        raise AssertionError('No matches found for file {} / group {} (query size: {})'.format(params['inputfile'], params['group_fg'], len(fg)))
    return params['chrom'], fg, params['group_fg'], matched, params['group_bg']


def assemble_worker_params(args):
    """
    :param args:
    :return:
    """
    commons = dict()
    commons['inputfile'] = args.inputfile
    commons['seqfile'] = args.seqfile
    commons['features'] = args.features
    commons['kmers'] = args.kmers
    commons['timeout'] = args.timeout
    commons['relax_init'] = args.relaxinit
    commons['relax_step'] = args.relaxstep
    commons['relax_limit'] = args.relaxlimit
    commons['allownomatch'] = args.allownomatch
    ingroups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    arglist = []
    for ig in ingroups:
        tmp = dict(commons)
        _, chrom = os.path.split(ig)
        tmp['chrom'] = chrom
        tmp['group_fg'] = ig
        tmp['group_bg'] = os.path.join(args.outputgroup, chrom)
        arglist.append(tmp)
    return arglist


def regions_to_dataframe(regions):
    """
    :param regions:
    :return:
    """
    featnames, numfeat = get_feature_selector(regions[0])
    columns = ['start', 'end', 'name'] + featnames + ['seq']
    getvals = op.itemgetter(*tuple(columns))
    regions = sorted(regions, key=lambda d: d['index'])
    indices = [reg['index'] for reg in regions]
    df = pd.DataFrame([list(getvals(reg)) for reg in regions], columns=columns, index=indices)
    return df


def run_background_match(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    if 5 < len(args.features) + len(args.kmers) <= 10:
        logger.warning('Large number of features specified; the more features are to be used to find'
                       ' similar regions, the more inefficient the nearest neighbor search is'
                       ' (i.e. the high dimensional space is sparsely populated).')
    if len(args.features) + len(args.kmers) > 10:
        logger.error('Too many features specified - aborting search...')
        raise AssertionError('Too many features for nearest neighbor search.')
    arglist = assemble_worker_params(args)
    logger.debug('Compiled argument list of size {} to process'.format(len(arglist)))
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdfout:
        logger.debug('Initializing worker pool of size {}'.format(args.workers))
        if 'metadata' in hdfout:
            metadata = hdfout['metadata']
        else:
            metadata = pd.DataFrame(columns=MD_REGION_COLDEFS)
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(find_background_regions, arglist, chunksize=1)
            logger.debug('Waiting for results...')
            for chrom, fgreg, fggrp, bgreg, bggrp in resit:
                logger.debug('Found {} matches (of max {}) for {}'.format(len(bgreg), len(fgreg), chrom))
                if not bgreg:
                    logger.debug('No matches for group {}, discarding input set'.format(fgreg))
                    continue
                df_fg = regions_to_dataframe(fgreg)
                grp, df_fg, metadata = gen_obj_and_md(metadata, fggrp, chrom, os.path.basename(args.inputfile), df_fg)
                hdfout.put(grp, df_fg, format='fixed')
                hdfout.flush()
                df_bg = regions_to_dataframe(bgreg)
                grp, df_bg, metadata = gen_obj_and_md(metadata, bggrp, chrom, os.path.basename(args.inputfile), df_bg)
                hdfout.put(grp, df_bg, format='fixed')
                hdfout.flush()
                logger.debug('Saved matched regions')
            hdfout.put('/metadata', metadata, format='table')
            hdfout.flush()
            logger.debug('Stored metadata in HDF')
    logger.debug('Processing done')
    return 0
