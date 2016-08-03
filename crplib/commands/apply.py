# coding=utf-8

"""
Module to apply a previously trained model to estimate the epigenome
for a specific cell type in a different species
"""

import os as os
import pandas as pd
import numpy as np
import collections as col
import numpy.random as rng
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
from crplib.metadata.md_signal import MD_SIGNAL_COLDEFS
from crplib.metadata.md_signal import gen_obj_and_md as gen_sigobj
from crplib.metadata.md_regions import MD_REGION_COLDEFS
from crplib.metadata.md_regions import gen_obj_and_md as genregobj


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
                group, valobj, metadata = gen_sigobj(metadata, args.outputgroup, chrom, args.inputfile, valobj)
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


def load_region_data(fpath, groups, features, labelcol, labeltype='class'):
    """
    :param fpath:
    :param groups:
    :param features:
    :param labelcol:
    :return:
    """
    if labeltype == 'class':
        coltype = np.int32
    elif labeltype == 'value':
        coltype = np.float64
    else:
        raise ValueError('Unknown label type: {}'.format(labeltype))
    with pd.HDFStore(fpath, 'r') as hdf:
        if isinstance(groups, (list, tuple)):
            full_dataset = pd.concat([hdf[grp] for grp in sorted(groups)], ignore_index=True)
        else:
            full_dataset = hdf[groups]
        if labelcol in full_dataset.columns:
            classlabels = full_dataset.loc[:, labelcol].astype(coltype, copy=False)
        else:
            classlabels = None
        region_names = [cn for cn in full_dataset.columns if cn in ['name', 'source']]
        if len(region_names) > 0:
            region_names = full_dataset.loc[:, region_names[0]].values
        else:
            region_names = np.array([])
        dataset = full_dataset.loc[:, features]
    return dataset, region_names, classlabels


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
    feat_order = model_md['feat_order']
    load_groups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    logger.debug('Loading region dataset')
    dataset, names, class_true = load_region_data(args.inputfile, load_groups, feat_order, args.classlabels, args.labeltype)
    logger.debug('Loaded dataset of size {}'.format(dataset.shape))
    y_pred = model.predict(dataset)
    class_order = list(map(int, model.classes_))
    y_prob = pd.DataFrame(model.predict_proba(dataset), columns=class_order)
    pred_class_prob = list(map(float, y_prob.lookup(np.arange(y_prob.shape[0]), y_pred)))
    names = list(map(str, names.tolist()))
    if class_true is not None:
        # to serialize this to JSON, need to cast
        # all numeric values to Python types (from numpy)
        true_class_prob = list(map(float, y_prob.lookup(np.arange(y_prob.shape[0]), class_true)))
        if len(class_order) == 2:  # binary classification
            auc = roc_auc_score(class_true, y_prob[:, 1])
            f1 = f1_score(class_true, y_pred)
            acc = accuracy_score(class_true, y_pred)
            scores = {'f1': float(f1), 'accuracy': float(acc), 'roc_auc': float(auc)}
        else:  # multiclass classification
            auc = -1
            acc = accuracy_score(class_true, y_pred)
            f1 = f1_score(class_true, y_pred, average='weighted', pos_label=None)
            scores = {'f1': float(f1), 'accuracy': float(acc), 'roc_auc': auc}
        priors = col.Counter(class_true)
        asc_classes = sorted(class_order)
        num_samples = len(class_true)
        priors = np.array([priors[cnt] / num_samples for cnt in asc_classes], dtype=np.float64)
        runif_f1 = []
        runif_acc = []
        rprior_f1 = []
        rprior_acc = []
        for _ in range(1000):
            runif_f1.append(f1_score(class_true, rng.choice(asc_classes, num_samples), average='weighted', pos_label=None))
            runif_acc.append(accuracy_score(class_true, rng.choice(asc_classes, num_samples)))
            rprior_f1.append(f1_score(class_true, rng.choice(asc_classes, num_samples, p=priors), average='weighted', pos_label=None))
            rprior_acc.append(accuracy_score(class_true, rng.choice(asc_classes, num_samples, p=priors)))
        rand_scores = {'f1_runif_mean': float(np.mean(runif_f1)), 'f1_runif_max': float(np.max(runif_f1)),
                       'f1_runif_95pct': float(np.percentile(runif_f1, 95)),
                       'f1_rprior_mean': float(np.mean(rprior_f1)), 'f1_rprior_max': float(np.max(rprior_f1)),
                       'f1_rprior_95pct': float(np.percentile(rprior_f1, 95)),
                       'acc_runif_mean': float(np.mean(runif_acc)), 'acc_runif_max': float(np.max(runif_acc)),
                       'acc_runif_95pct': float(np.percentile(runif_acc, 95)),
                       'acc_rprior_mean': float(np.mean(rprior_acc)), 'acc_rprior_max': float(np.max(rprior_acc)),
                       'acc_rprior_95pct': float(np.percentile(rprior_acc, 95))}
        scores.update(rand_scores)
        true_labels = list(map(int, class_true))
    else:
        scores = {'f1': -1, 'accuracy': -1, 'roc_auc': -1}
        true_labels = []
        true_class_prob = []
    class_pred = list(map(int, y_pred))
    class_probs = [list(map(float, row)) for row in y_prob.as_matrix().tolist()]
    dump = {'scores': scores, 'class_true': true_labels, 'class_pred': class_pred,
            'class_probs': class_probs, 'class_order': class_order,
            'inputfile': os.path.basename(args.inputfile), 'modelfile': os.path.basename(args.modelfile),
            'n_samples': dataset.shape[0], 'n_features': dataset.shape[1],
            'names': names, 'class_pred_prob': pred_class_prob, 'class_true_prob': true_class_prob}
    logger.debug('Dumping prediction to output file')
    with open(args.outputfile, 'w') as outfile:
        _ = json.dump(dump, outfile, indent=1)
    logger.debug('Done')
    return 0


def build_output_dataframe(dset, pred, probs, classes, merge, reduce):
    """
    :param dset:
    :param pred:
    :param probs:
    :param classes:
    :param merge:
    :return:
    """
    dset = dset.assign(class_pred=pred)
    class_cols = ['class_prob_' + cls for cls in map(str, list(map(int, classes)))]
    dset = pd.concat([dset, pd.DataFrame(probs, columns=class_cols)], axis='columns')
    if reduce:
        dset = dset.loc[dset.class_pred.isin(reduce), :]
    if merge:
        assert len(reduce) == 1, 'Merging overlapping regions in dataset w/o reducing to single class not supported'
        dset.drop([col for col in dset.columns if col.startswith('ft')], axis='columns', inplace=True)
        dset.sort_values(by=['start', 'end'], axis='index', ascending=True, inplace=True)
        dset.index = np.arange(dset.shape[0])
        # TODO
        # in spare time, find out if there is a more
        # native way in Pandas to merge overlapping intervals...
        new_rows = []
        cur_start = dset.loc[0, 'start']
        cur_end = dset.loc[0, 'end']
        cur_probs = []
        cur_names = set()
        get_name = 'name' if 'name' in dset.columns else 'source'
        assert get_name in dset.columns, 'No naming column exists for regions'
        class_prob = 'class_prob_' + str(reduce[0])
        for row in dset.itertuples(index=False):
            if row.start <= cur_end:
                cur_end = row.end
                cur_probs.append(row.__getattribute__(class_prob))
                cur_names.add(row.__getattribute__(get_name))
            else:
                regname = cur_names.pop()  # return any or unique
                new_rows.append([cur_start, cur_end, regname, np.average(cur_probs)])
                cur_names = set()
                cur_probs = []
                cur_start = row.start
                cur_end = row.end
                cur_probs.append(row.__getattribute__(class_prob))
                cur_names.add(row.__getattribute__(get_name))
        new_cols = ['start', 'end', 'name', 'class_prob']
        dset = pd.DataFrame(new_rows, columns=new_cols)
    return dset


def model_scan_regions(params):
    """
    :param params:
    :return:
    """
    model_md = json.load(open(params['modelmetadata'], 'r'))
    model = pck.load(open(params['modelfile'], 'rb'))
    feat_order = model_md['feature_order']
    featdata, _, _ = load_region_data(params['inputfile'], params['inputgroup'], feat_order, params['classlabels'])
    y_pred = model.predict(featdata)
    y_prob = model.predict_proba(featdata)
    class_order = model.classes_
    with pd.HDFStore(params['inputfile'], 'r') as hdf:
        full_dataset = hdf[params['inputgroup']]
    df = build_output_dataframe(full_dataset, y_pred, y_prob, class_order, params['merge'], params['reduce'])
    return params['chrom'], df


def assemble_params_scnreg(args, logger):
    """
    :return:
    """
    all_groups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    logger.debug('Identified {} valid groups in input file'.format(len(all_groups)))
    merge_regions = False
    with pd.HDFStore(args.inputfile, 'r') as hdf:
        md = hdf['metadata']
        res = md.loc[0, 'resolution']
        if res > 0:
            merge_regions = True
    logger.debug('Detected - merge regions: {}'.format(merge_regions))
    if not args.modelmetadata:
        fpath_md = args.modelfile.rsplit('.', 1)[0] + '.json'
    else:
        fpath_md = args.modelmetadata
    commons = vars(args)
    del commons['module_logger']
    del commons['execute']
    commons['modelmetadata'] = fpath_md
    arglist = []
    for g in all_groups:
        tmp = dict(commons)
        tmp['inputgroup'] = g
        _, tmp['chrom'] = os.path.split(g)
        tmp['merge'] = merge_regions
        arglist.append(tmp)
    logger.debug('Build argument list of size {} to process'.format(len(arglist)))
    return arglist


def run_scan_regions(logger, args):
    """
    :param logger:
    :param args:
    :return:
    """
    arglist = assemble_params_scnreg(args, logger)
    with pd.HDFStore(args.outputfile, 'w') as hdf:
        metadata = pd.DataFrame(columns=MD_REGION_COLDEFS)
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(model_scan_regions, arglist, chunksize=1)
            for chrom, dataobj in resit:
                logger.debug('Received data for chromosome {}'.format(chrom))
                grp, dataobj, metadata = genregobj(metadata, args.outputgroup, chrom, [args.inputfile, args.modelfile], dataobj)
                hdf.put(grp, dataobj, format='fixed')
                hdf.flush()
                logger.debug('Flushed data to file')
        hdf.put('metadata', metadata, format='table')
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
    elif args.task == 'scnreg':
        logger.debug('Running task: scan regions')
        rv = run_scan_regions(logger, args)
    else:
        raise ValueError('Unknown task: {}'.format(args.task))
    return rv
