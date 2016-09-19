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

from crplib.auxiliary.modeling import load_model, get_scorer, load_ml_dataset, extract_model_attributes
from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups
from crplib.auxiliary.file_ops import create_filepath


def train_nocv(model, params, traindata, outputs, sampleweights):
    """
    :param model:
    :param params:
    :param traindata:
    :param outputs:
    :return:
    """
    model = model.set_params(**params)
    model = model.fit(traindata, outputs, sample_weight=sampleweights)
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
    scorer = get_scorer(params['scoring'], sampleweights)
    param_grid = dict(params)
    del param_grid['scoring']
    fit_params = None
    if sampleweights is not None:
        fit_params = {'sample_weight': sampleweights}
    tune_model = GridSearchCV(model, param_grid, scoring=scorer, pre_dispatch=njobs,
                              cv=folds, n_jobs=njobs, refit=True, fit_params=fit_params)
    tune_model = tune_model.fit(traindata, outputs)
    return tune_model


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


def run_train_model(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    _ = create_filepath(args.modelout, logger)
    logger.debug('Loading model specification from {}'.format(args.modelspec))
    model_spec = json.load(open(args.modelspec))
    model = load_model(model_spec['module_path'], model_spec['model_name'])
    load_groups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    traindata, targets, dtinfo, sminfo, ftinfo = load_ml_dataset(args.inputfile, load_groups, None, args, logger)
    assert traindata.shape[0] > 1, 'No samples (rows) in training data'
    assert traindata.shape[1] > 1, 'No features (columns) in training data'
    if targets is not None:
        assert targets.size == traindata.shape[0], 'Mismatch num targets {} and num samples {}'.format(targets.size, traindata.shape[0])
    run_metadata = {'dataset_info': dtinfo, 'sample_info': sminfo,
                    'feature_info': ftinfo, 'model_info': dict()}
    logger.debug('Training model')
    if args.notuning:
        params = model_spec['default']
        model = train_nocv(model, params, traindata, targets, sminfo['weights'])
        run_metadata['model_info']['params'] = params
        run_metadata['model_info']['tuned'] = False
    else:
        params = model_spec['cvtune']
        tune_info = train_gridcv(model, params, traindata, targets, args.cvfolds, args.workers, sminfo['weights'])
        model = tune_info.best_estimator_
        run_metadata['model_info']['params'] = tune_info.best_params_
        run_metadata['model_info']['tuned'] = True
        run_metadata['training_info'] = dict()
        run_metadata['training_info']['cv_scores'] = simplify_cv_scores(tune_info.grid_scores_)
        run_metadata['training_info']['best_score'] = tune_info.best_score_
        run_metadata['training_info']['scoring'] = params['scoring']
    run_metadata['model_info']['name'] = model_spec['model_name']
    run_metadata['model_info']['type'] = model_spec['model_type']
    logger.debug('Training finished')
    if 'store_attributes' in model_spec:
        logger.debug('Storing user requested model attributes')
        attribs = extract_model_attributes(model, model_spec['store_attributes'], logger)
        run_metadata['attribute_info'] = attribs
    if args.calcweights:
        raise NotImplementedError('Currently not functional')

    logger.debug('Saving model and metadata')
    run_metadata['run_info'] = dict()
    run_metadata['run_info']['model_spec'] = os.path.basename(args.modelspec)
    run_metadata['run_info']['model_file'] = os.path.basename(args.modelout)
    run_metadata['run_info']['train_data'] = os.path.basename(args.inputfile)
    run_metadata['run_info']['train_group'] = args.inputgroup

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
        _ = json.dump(run_metadata, outfile, indent=1, sort_keys=True)
    logger.debug('Done')
    return 0
