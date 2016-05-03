# coding=utf-8

"""
Module to apply a previously trained model to estimate the epigenome
for a specific cell type in a different species
"""

import pandas as pd
import numpy as np
import operator as op
import multiprocessing as mp
import json as json
import pickle as pck

from scipy.interpolate import LSQUnivariateSpline as kspline

from crplib.auxiliary.seq_parsers import get_twobit_seq
from crplib.auxiliary.hdf_ops import load_masked_sigtrack, get_valid_hdf5_groups, get_trgindex_groups
from crplib.auxiliary.file_ops import create_filepath
from crplib.mlfeat.featdef import feat_mapsig, get_online_version
from crplib.auxiliary.constants import CHROMOSOME_BOUNDARY
from crplib.metadata.md_signal import MD_SIGNAL_COLDEFS, gen_obj_and_md


def smooth_signal_estimate(signal, res):
    """
    :param signal:
    :param res:
    :return:
    """
    # for default resolution of 25, spans roughly one nucleosome
    window = res * 10
    for pos in range(CHROMOSOME_BOUNDARY, len(signal), window):
        smoother = kspline(range(pos, pos+window),  # x-axis
                           signal[pos:pos+window],  # y-axis
                           t=[j+res for j in range(pos, pos + window - res, res)],  # knots
                           k=3)  # degree of polynomial, cubic is default
        signal[pos:pos+window] = smoother(range(pos, pos+window))
    signal = np.clip(signal, 0, signal.max())
    return signal


def make_signal_estimate(params):
    """
    :param params:
    :return:
    """
    chrom = params['chrom']
    mypid = mp.current_process().pid
    model = pck.load(open(params['modelfile'], 'rb'))
    seq = get_twobit_seq(params['seqfile'], params['chrom'])
    chromlen = len(seq)
    res = params['resolution']
    index_groups = get_trgindex_groups(params['targetindex'], '')
    with pd.HDFStore(params['targetindex'], 'r') as idx:
        mask = idx[index_groups[chrom]['mask']]
    map_sig = load_masked_sigtrack(params['inputfile'], '',
                                   params['inputgroup'], chrom, chromlen, mask=mask)
    mapfeat = feat_mapsig
    est_sig = np.zeros(chromlen, dtype=np.float64)
    lolim = CHROMOSOME_BOUNDARY
    hilim = int(chromlen // res * res)
    comp_seqfeat = get_online_version(params['features'], params['kmers'])
    get_values = op.itemgetter(*tuple(params['feature_order']))
    chunks = []
    positions = []
    for pos in range(lolim, hilim, res):
        chunk = {'seq': seq[pos:pos+res]}
        positions.append(pos)
        chunk.update(mapfeat(map_sig[pos:pos+res]))
        chunk = comp_seqfeat(chunk)
        chunks.append(get_values(chunk))
        if len(chunks) >= 10000:
            y_hat = model.predict(np.array(chunks))
            for idx, val in zip(positions, y_hat):
                est_sig[idx:idx+res] = val
            chunks = []
            positions = []
    if len(chunks) > 0:
        y_hat = model.predict(np.array(chunks))
        for idx, val in zip(positions, y_hat):
            est_sig[idx:idx+res] = val
    if not params['nosmooth']:
        est_sig = smooth_signal_estimate(est_sig, res)
    return mypid, params['chrom'], est_sig


def assemble_params_estsig(args):
    """
    :param args:
    :return:
    """
    all_groups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    if not args.modelmetadata:
        fpath_md = args.modelfile.rsplit('.', 1)[0] + '.json'
    else:
        fpath_md = args.modelmetadata
    model_md = json.load(open(fpath_md, 'r'))
    commons = {'modelfile': args.modelfile, 'resolution': int(model_md['resolution']),
               'seqfile': args.seqfile, 'targetindex': args.targetindex, 'inputfile': args.inputfile,
               'inputgroup': args.inputgroup, 'features': model_md['features'],
               'kmers': model_md['kmers'], 'feature_order': model_md['feature_order'],
               'nosmooth': args.nosmooth}
    arglist = []
    for g in all_groups:
        chrom = g.rsplit('/', 1)[1]
        tmp = dict(commons)
        tmp['chrom'] = chrom
        arglist.append(tmp)
    return arglist


def run_estimate_signal(logger, args):
    """
    :param logger:
    :param args:
    :return:
    """
    logger.debug('Assembling worker parameters')
    arglist = assemble_params_estsig(args)
    with pd.HDFStore(args.outputfile, 'a', complevel=9, complib='blosc') as hdfout:
        if 'metadata' in hdfout:
            metadata = hdfout['metadata']
        else:
            metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
        with mp.Pool(args.workers) as pool:
            logger.debug('Start processing...')
            mapres = pool.imap_unordered(make_signal_estimate, arglist, chunksize=1)
            for chrom, valobj in mapres:
                logger.debug('Processed chromosome {}'.format(chrom))
                group, valobj, metadata = gen_obj_and_md(metadata, args.outputgroup, chrom, args.inputfile, valobj)
                hdfout.put(group, valobj, format='fixed')
                hdfout.flush()
                logger.debug('Estimated signal data saved')
        logger.debug('Saving metadata')
        hdfout.put('metadata', metadata)
        hdfout.flush()
    logger.debug('All chromosomes processed')
    return 0


def run_apply_model(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    rv = 0
    _ = create_filepath(args.outputfile, logger)
    if args.task == 'estsig':
        logger.debug('Running task: estimate signal')
        rv = run_estimate_signal(logger, args)
    else:
        raise ValueError('Unknown task: {}'.format(args.task))
    return rv
