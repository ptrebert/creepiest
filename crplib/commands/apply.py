# coding=utf-8

"""
Module to apply a previously trained model to estimate the epigenome
for a specific cell type in a different species
"""

import os as os
import pandas as pd
import numpy as np
import operator as op
import multiprocessing as mp
import json as json
import pickle as pck

from scipy.interpolate import LSQUnivariateSpline as kspline

from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

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
        try:
            smoother = kspline(range(pos, pos+window),  # x-axis
                               signal[pos:pos+window],  # y-axis
                               t=[j+res for j in range(pos, pos + window - res, res)],  # knots
                               k=3)  # degree of polynomial, cubic is default
            signal[pos:pos+window] = smoother(range(pos, pos+window))
        except Exception as err:
            if pos + window > len(signal):
                break
            else:
                raise err
    signal = np.clip(signal, 0, signal.max())
    return signal


def make_signal_estimate(params):
    """
    :param params:
    :return:
    """
    chrom = params['chrom']
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
    return chrom, est_sig


def assemble_params_estsig(args):
    """
    :param args:
    :return:
    """
    all_groups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    model_md = load_model_metadata(args)
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
    logger.debug('Created parameter list of size {}'.format(len(arglist)))
    with pd.HDFStore(args.outputfile, 'w', complevel=9, complib='blosc', encoding='utf-8') as hdfout:
        metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
        with mp.Pool(args.workers) as pool:
            logger.debug('Start processing...')
            resit = pool.imap_unordered(make_signal_estimate, arglist, chunksize=1)
            for chrom, valobj in resit:
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


def load_model_metadata(args):
    """
    :param args:
    :return:
    """
    if not args.modelmetadata:
        fpath_md = args.modelfile.rsplit('.', 1)[0] + '.json'
    else:
        fpath_md = args.modelmetadata
    model_md = json.load(open(fpath_md, 'r'))
    return model_md


def load_region_data(fpath, groups, features, labelcol):
    """
    :param fpath:
    :param groups:
    :param features:
    :param labelcol:
    :return:
    """

    with pd.HDFStore(fpath, 'r') as hdf:
        full_dataset = pd.concat([hdf[grp] for grp in groups], ignore_index=True)
        if labelcol in full_dataset.columns:
            classlabels = full_dataset.loc[:, labelcol].astype(np.int32, copy=False)
        else:
            classlabels = None
        dataset = full_dataset.loc[:, features]
    return dataset, classlabels


def run_classify_regions(logger, args):
    """
    :param logger:
    :param args:
    :return:
    """
    logger.debug('Loading model metadata')
    model_md = load_model_metadata(args)
    logger.debug('Loading model')
    model = pck.load(open(args.modelfile, 'rb'))
    feat_order = model_md['feature_order']
    load_groups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    logger.debug('Loading region dataset')
    dataset, class_true = load_region_data(args.inputfile, load_groups, feat_order, args.classlabels)
    logger.debug('Loaded dataset of size {}'.format(dataset.shape))
    y_pred = model.predict(dataset)
    y_prob = model.predict_proba(dataset)
    class_order = list(map(int, model.classes_))
    if class_true is not None:
        # to serialize this to JSON, need to cast
        # all numeric values to Python types (from numpy)
        if len(class_order) == 2:  # binary classification
            auc = roc_auc_score(class_true, y_prob[:, 1])
            f1 = f1_score(class_true, y_pred)
            acc = accuracy_score(class_true, y_pred)
            scores = {'f1': float(f1), 'accuracy': float(acc), 'roc_auc': float(auc)}
        else:  # multiclass classification
            acc = -1
            auc = -1
            f1 = f1_score(class_true, y_pred, average='weighted', pos_label=None)
            scores = {'f1': float(f1), 'accuracy': acc, 'roc_auc': auc}
        true_labels = list(map(int, class_true))
    else:
        scores = {'f1': -1, 'accuracy': -1, 'roc_auc': -1}
        true_labels = []
    class_pred = list(map(int, y_pred))
    class_probs = [list(map(float, entry)) for entry in y_prob]
    dump = {'scores': scores, 'class_true': true_labels, 'class_pred': class_pred,
            'class_probs': class_probs, 'class_order': class_order,
            'inputfile': os.path.basename(args.inputfile), 'modelfile': os.path.basename(args.modelfile),
            'n_samples': dataset.shape[0], 'n_features': dataset.shape[1]}
    logger.debug('Dumping prediction to output file')
    with open(args.outputfile, 'w') as outfile:
        _ = json.dump(dump, outfile, indent=1)
    logger.debug('Done')
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
    elif args.task == 'clsreg':
        logger.debug('Running task: classify regions')
        rv = run_classify_regions(logger, args)
    else:
        raise ValueError('Unknown task: {}'.format(args.task))
    return rv
