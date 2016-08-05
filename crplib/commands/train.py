# coding=utf-8

"""
Command to train models on predefined training datasets
"""

import os as os
import json as json
import numpy as np
import pickle as pck
import pandas as pd

from sklearn.grid_search import GridSearchCV
from sklearn.base import clone
from sklearn.cross_validation import StratifiedKFold

from crplib.auxiliary.modeling import load_model, get_scorer
from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups
from crplib.auxiliary.file_ops import create_filepath, text_file_mode
from crplib.mlfeat.featdef import get_prefix_list, get_classes_from_names


def train_nocv(model, params, traindata, outputs, sampleweights):
    """
    :param model:
    :param params:
    :param traindata:
    :param outputs:
    :return:
    """
    model = model.set_params(**params)
    if sampleweights is not None:
        if isinstance(sampleweights, list) and len(sampleweights) == 0:
            sampleweights = None
        model = model.fit(traindata, outputs, sample_weight=sampleweights)
    else:
        model = model.fit(traindata, outputs, sample_weight=None)
    return model


def train_gridcv(model, params, traindata, outputs, folds, njobs, sampleweights):
    """
    :param model:
    :param params:
    :param traindata:
    :param outputs:
    :param folds:
    :param njobs:
    :return:
    """
    scorer = get_scorer(params['scoring'])
    param_grid = dict(params)
    del param_grid['scoring']
    fit_params = None
    # TODO this needs to be simplified
    if sampleweights is not None:
        if isinstance(sampleweights, list) and len(sampleweights) == 0:
            fit_params = None
        else:
            fit_params = {'sample_weight': sampleweights}
    tune_model = GridSearchCV(model, param_grid, scoring=scorer, pre_dispatch=njobs,
                              cv=folds, n_jobs=njobs, refit=True, fit_params=fit_params)
    tune_model = tune_model.fit(traindata, outputs)
    return tune_model


def load_training_data(filepath, prefix, onlyfeatures, depvar, sampleweights):
    """
    :param filepath:
    :param prefix:
    :param onlyfeatures:
    :param depvar:
    :param weights:
    :return:
    """
    all_groups = get_valid_hdf5_groups(filepath, prefix)
    with pd.HDFStore(filepath, 'r') as hdf:
        full_dset = pd.concat([hdf[g] for g in all_groups], ignore_index=True)
        feat_order = sorted([ft for ft in full_dset.columns if ft.startswith('ft')])
        assert feat_order, 'No features selected from column names: {}'.format(full_dset.columns)
        mdf = hdf['metadata']
        load_group = all_groups[0]
        if 'features' in mdf.columns:
            # regular training dataset
            ft_classes = mdf.where(mdf.group == load_group).dropna()['features']
            ft_classes = ft_classes.values[0].split(',')
            ft_kmers = mdf.where(mdf.group == load_group).dropna()['kmers']
            ft_kmers = list(map(int, ft_kmers.values[0].split(',')))
            res = mdf.where(mdf.group == load_group).dropna()['resolution']
            res = res.values[0]
        else:
            # apparently, just regions
            ft_classes = get_classes_from_names(feat_order)
            ft_kmers = []
            res = 0
    weights = []
    outputs = pd.Series(full_dset.loc[:, depvar])
    names = full_dset.loc[:, 'name'].tolist()
    if sampleweights is not None:
        if 'weight' in full_dset.columns:
            full_dset.drop('weight', axis='columns', inplace=True)
        # this merge ensures that no matter the sorting of the sample weights in the
        # annotation file, all samples are assigned the correct weight
        full_dset = full_dset.merge(sampleweights, how='outer', on='name', suffixes=('', ''), copy=False)
        weights = pd.Series(full_dset.loc[:, 'weight']).values
    if onlyfeatures:
        prefixes = get_prefix_list(onlyfeatures)
        feat_order = list(filter(lambda x: any([x.startswith(p) for p in prefixes]), feat_order))
        assert feat_order, 'No features left after filtering for: {}'.format(onlyfeatures)
        ft_classes = onlyfeatures
    traindata = full_dset.loc[:, feat_order]
    md_info = {'features': ft_classes, 'kmers': ft_kmers, 'feat_order': feat_order,
               'resolution': res, 'names': names, 'weights': list(map(float, weights))}
    return md_info, traindata, outputs, weights


def simplify_cv_scores(cvfolds):
    """
    :param cvfolds:
     :type: list of sklearn.grid_search._CVScoreTuple
    :return:
     :rtype: list of dict
    """
    grid = []
    for index, entry in enumerate(cvfolds):
        this_comb = {'grid_index': index}
        values = np.array(entry.cv_validation_scores)
        mean = np.mean(values)
        std = np.std(values)
        min_score = np.min(values)
        max_score = np.max(values)
        this_comb['scores'] = {'mean': mean, 'std': std,
                               'max': max_score, 'min': min_score}
        this_comb['params'] = entry.parameters
        grid.append(this_comb)
    return grid


def calc_sample_weights(model, traindata, outputs, names, k):
    """
    :param model:
    :param traindata:
    :param outputs:
    :return:
    """
    cv_indices = StratifiedKFold(outputs, n_folds=k)
    ra_names = pd.Series(names)
    # clone model since fit returns self
    partial_model = clone(model)
    weights = []
    for train_index, test_index in cv_indices:
        partial_model = partial_model.fit(traindata.iloc[train_index, :], outputs.iloc[train_index])
        y_probs = pd.DataFrame(partial_model.predict_proba(traindata.iloc[test_index, :]), columns=list(model.classes_))
        true_class_prob = y_probs.lookup(y_probs.index, outputs.iloc[test_index])
        for n, p in zip(ra_names.iloc[test_index], true_class_prob):
            weights.append([n, 1 - p])
    return weights


def load_sample_weights(fpath):
    """
    :param fpath:
    :return:
    """
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode) as mdfile:
        annotation = json.load(mdfile)
        smpwt = pd.DataFrame(annotation['sample_weights'], columns=['name', 'weight'])
        assert not smpwt.empty, 'Entry sample_weights in annotation file {} is empty'.format(fpath)
    return smpwt


def run_train_model(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    _ = create_filepath(args.modelout, logger)
    logger.debug('Loading model specification from {}'.format(args.modelspec))
    model_params = json.load(open(args.modelspec))
    model = load_model(model_params['module_path'], model_params['model_name'])
    sample_weights = None
    if args.sampleweights:
        logger.debug('Loading user requested sample weights from file {}'.format(args.sampleweights))
        sample_weights = load_sample_weights(args.sampleweights)
    logger.debug('Loading training data')
    md_info, traindata, outputs, weights = load_training_data(args.traindata, args.traingroup, args.onlyfeatures, args.depvar, sample_weights)
    logger.debug('Training model')
    if args.notuning:
        params = model_params['default']
        model = train_nocv(model, params, traindata, outputs, weights)
        metadata = {'final_params': params}
    else:
        params = model_params['cvtune']
        tune_info = train_gridcv(model, params, traindata, outputs, args.cvfolds, args.workers, weights)
        model = tune_info.best_estimator_
        metadata = {'final_params': tune_info.best_params_}
        metadata['grid_scores'] = simplify_cv_scores(tune_info.grid_scores_)
        metadata['best_score'] = tune_info.best_score_
        metadata['scoring'] = params['scoring']
    logger.debug('Training finished')
    if 'store_attributes' in model_params:
        logger.debug('Storing user requested model attributes')
        for attr in model_params['store_attributes']:
            if hasattr(model, attr):
                metadata[attr] = list(getattr(model, attr))
            else:
                logger.debug('Skipping attribute {} - does not exist'.format(attr))
    if args.calcweights:
        logger.debug('Calculating sample weights...')
        assert model_params['model_type'] == 'classifier',\
            'Calculating sample weights not implemented for regression models'
        wt = calc_sample_weights(model, traindata, outputs, md_info['names'], args.cvfolds)
        metadata['sample_weights'] = wt
        logger.debug('Sample weights calculated')
    logger.debug('Saving model and metadata')
    metadata['model_spec'] = os.path.basename(args.modelspec)
    metadata['init_params'] = params
    metadata['model'] = model_params['model_name']
    metadata['model_type'] = model_params['model_type']
    metadata['model_name'] = model_params['model_name']
    metadata['traindata'] = os.path.basename(args.traindata)
    metadata['traingroup'] = args.traingroup
    metadata['traindata_size'] = list(traindata.shape)
    metadata['modelfile'] = os.path.basename(args.modelout)
    metadata['depvar'] = args.depvar
    metadata.update(md_info)

    logger.debug('Writing model file...')
    with open(args.modelout, 'wb') as outfile:
        pck.dump(model, outfile)

    if not args.metadataout:
        mdout = args.modelout.rsplit('.', 1)[0] + '.json'
    else:
        mdout = args.metadataout
    _ = create_filepath(mdout, logger)
    logger.debug('Writing model metadata...')
    with open(mdout, 'w') as outfile:
        _ = json.dump(metadata, outfile, indent=0)

    return 0
