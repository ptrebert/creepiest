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

from sklearn.cross_validation import permutation_test_score as permtest
from sklearn.cross_validation import StratifiedKFold as skf

from crplib.auxiliary.seq_parsers import get_twobit_seq
from crplib.auxiliary.hdf_ops import load_masked_sigtrack, get_valid_hdf5_groups, get_trgindex_groups
from crplib.auxiliary.file_ops import create_filepath
from crplib.auxiliary.modeling import select_dataset_subset, load_model, \
    load_model_metadata, get_scorer
from crplib.mlfeat.featdef import feat_mapsig, get_online_version
from crplib.auxiliary.constants import CHROMOSOME_BOUNDARY
from crplib.metadata.md_signal import MD_SIGNAL_COLDEFS
from crplib.metadata.md_signal import gen_obj_and_md as gen_sigobj
from crplib.metadata.md_regions import MD_REGION_COLDEFS
from crplib.metadata.md_regions import gen_obj_and_md as genregobj


def load_dataset(fpath, groups, features, subset='', ycol='', ytype=None):
    """
    :param fpath:
    :param groups:
    :param features:
    :param subset:
    :param ycol:
    :param ytype:
    :return:
    """
    with pd.HDFStore(fpath, 'r') as hdf:
        if isinstance(groups, (list, tuple)):
            dataset = pd.concat([hdf[grp] for grp in sorted(groups)], ignore_index=True)
        else:
            dataset = hdf[groups]
        dataset = select_dataset_subset(dataset, subset)
        y_depvar = None
        if ytype is not None:
            y_depvar = dataset.loc[:, ycol].astype(ytype, copy=True)
        name_col = [cn for cn in dataset.columns if cn in ['name', 'source']]
        sample_names = []
        if name_col:
            sample_names = dataset.loc[:, name_col[0]].tolist()
        predictors = dataset.loc[:, features]
    return predictors, sample_names, y_depvar


def run_permutation_test(data, output, model, cvperm, numperm, workers, scorer):
    """
    :param data:
    :param output:
    :param model:
    :param cvperm:
    :param numperm:
    :param workers:
    :param scorer:
    :return:
    """
    true_score, perm_scores, pval = permtest(model, data.as_matrix(), output.values, cv=cvperm,
                                             n_permutations=numperm, n_jobs=workers,
                                             scoring=scorer, verbose=0)
    # based on: Ojala et al., J. ML Res. 2010
    # Permutation Tests for Studying Classifier Performance
    prec = float(1 / (2 * np.sqrt(numperm)))
    perm_params = {'perm_num': numperm, 'perm_cvfolds': cvperm, 'bound_test_prec': prec}
    perm_stats = {'perm_mean': float(np.mean(perm_scores)), 'perm_max': float(np.max(perm_scores)),
                  'perm_median': float(np.median(perm_scores)), 'perm_std': float(np.std(perm_scores)),
                  'perm_pval': float(pval)}
    perm_scores = perm_scores.tolist()  # for JSONify
    return true_score, perm_scores, perm_stats, perm_params


def run_randomization_test(output, classes, numrand, scorer):
    """
    :param output:
    :param classes:
    :param numrand:
    :param scorer:
    :return:
    """
    class_count = col.Counter(output)
    asc_classes = sorted(classes)
    propensities = np.array([class_count[cnt] / len(output) for cnt in asc_classes], dtype=np.float64)
    rand_scores = []
    for _ in range(numrand):
        rand_scores.append(scorer(output, rng.choice(asc_classes, len(output), p=propensities)))
    rand_stats = {'rand_mean': float(np.mean(rand_scores)),
                  'rand_max': float(np.max(rand_scores)), 'rand_median': float(np.median(rand_scores)),
                  'rand_95pct': float(np.percentile(rand_scores, 95))}
    return rand_stats


############################


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


def run_classification_testdata(args, model, model_md, dataset, y_true, logger):
    """
    :param args:
    :param model:
    :param model_md:
    :param dataset:
    :param y_true:
    :return:
    """
    logger.debug('Making prediction for test dataset')
    scorer = get_scorer(model_md['scoring'])
    y_pred = model.predict(dataset)
    model_perf = scorer(model, dataset, y_true)
    class_order = list(map(int, model.classes_))
    logger.debug('Getting class probabilities')
    y_prob = pd.DataFrame(model.predict_proba(dataset), columns=class_order)
    pred_class_prob = (y_prob.lookup(y_prob.index, y_pred)).tolist()
    true_class_prob = (y_prob.lookup(y_prob.index, y_true)).tolist()
    out_metadata = {}
    if not args.noperm:
        logger.debug('Running permutation test')
        perf_score, permscores, permstats, permparams = run_permutation_test(dataset, y_true, model, args.cvperm,
                                                                             args.numperm, args.workers, scorer)
        if not np.isclose(model_perf, perf_score, rtol=1e-05, atol=1e-05):
            # I do not think this can happen - is deterministic,
            # numerical instabilities are ignored
            logger.warning('Performance scores not close: {} vs {}'.format(model_perf, perf_score))
        out_metadata['permutation_test'] = {'perm_scores': permscores, 'perm_params': permparams,
                                            'perm_stats': permstats}
        logger.debug('Running randomization test')
        randstats = run_randomization_test(y_true, model.classes_, args.numperm, scorer)
        out_metadata['randomization_test'] = {'rand_num': args.numperm, 'rand_stats': randstats}
    logger.debug('Collecting metadata')
    out_metadata['scoring'] = model_md['scoring']
    out_metadata['labels'] = {'true': list(map(int, y_true)),
                              'pred': list(map(int, y_pred)),
                              'order': class_order}
    out_metadata['probabilities'] = {'true': true_class_prob,
                                     'pred': pred_class_prob,
                                     'all': [list(map(float, row)) for row in y_prob.as_matrix().tolist()]}
    out_metadata['model_perf'] = float(model_perf)
    logger.debug('Classification of test data finished')
    return out_metadata


def run_classification_newdata(model, dataset, logger):
    """
    :param model:
    :param dataset:
    :return:
    """
    y_pred = model.predict(dataset)
    class_order = list(map(int, model.classes_))
    y_prob = pd.DataFrame(model.predict_proba(dataset), columns=class_order)
    pred_class_prob = list(map(float, y_prob.lookup(y_prob.index, y_pred)))
    out_metadata = {}
    logger.debug('Collecting metadata')
    out_metadata['labels'] = {'pred': list(map(int, y_pred)),
                              'order': class_order}
    out_metadata['probabilities'] = {'pred': pred_class_prob,
                                     'all': [list(map(float, row)) for row in y_prob.as_matrix().tolist()]}
    logger.debug('Classification of test data finished')
    return out_metadata


def run_classification(args, model, modelmd, loadgroups, logger):
    """
    :param args:
    :param modelmd:
    :param logger:
    :return:
    """
    _ = create_filepath(args.outputfile, logger)
    logger.debug('Loading dataset')
    if args.task == 'test':
        dataset, samplenames, output = load_dataset(args.inputfile, loadgroups, modelmd['feat_order'],
                                                    subset=args.subset, ycol=args.outvariable, ytype=np.int32)
        out_md = run_classification_testdata(args, model, modelmd, dataset, output, logger)
    elif args.task == 'pred':
        dataset, samplenames, output = load_dataset(args.inputfile, loadgroups, modelmd['feat_order'],
                                                    subset=args.subset)
        assert output is None, 'Loaded sample outputs from dataset for prediction task'
        out_md = run_classification_newdata(model, dataset, logger)
    else:
        raise ValueError('Unknown task for classification: {}'.format(args.task))
    out_md['model_task'] = args.task
    #out_md['model_name'] = modelmd['model_name']
    out_md['model_type'] = modelmd['model_type']
    out_md['num_samples'] = dataset.shape[0]
    out_md['num_features'] = dataset.shape[1]
    out_md['sample_names'] = list(map(str, samplenames))
    logger.debug('Writing metadata of run...')
    with open(args.outputfile, 'w') as outfile:
        _ = json.dump(out_md, outfile, indent=1, sort_keys=True)
    logger.debug('Metadata saved')
    return 0


def run_apply_model(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    logger.debug('Loading model and metadata...')
    model_md = load_model_metadata(args)
    model = load_model(args.modelfile)
    model_type = model_md['model_type']
    logger.debug('Loading groups from input data file')
    load_groups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    if model_type == 'classifier':
        logger.debug('Applying classification model {} on task: {}'.format(model_md['model'], args.task))
        _ = run_classification(args, model, model_md, load_groups, logger)
    elif model_type == 'regressor':
        logger.debug('Applying regression model {} on task: {}'.format(model_md['model'], args.task))
        #_ = run_regression()
    else:
        raise NotImplementedError('No support for model of type: {} '
                                  '(just classifier or regressor are supported)'.format(model_type))
    return 0
