# coding=utf-8

"""
Command to train models on predefined training datasets
"""

import os as os
import json as json
import importlib as imp
import numpy as np
import pickle as pck
import pandas as pd

from sklearn.grid_search import GridSearchCV

from crplib.auxiliary.file_ops import get_valid_hdf5_groups


def train_nocv(model, params, traindata, outputs):
    """
    :param model:
    :param params:
    :param traindata:
    :param outputs:
    :return:
    """
    model = model.set_params(**params)
    model = model.fit(traindata, outputs)
    return model


def train_gridcv(model, params, traindata, outputs, folds, njobs):
    """
    :param model:
    :param params:
    :param traindata:
    :param outputs:
    :param folds:
    :param njobs:
    :return:
    """
    tune_model = GridSearchCV(model, params, scoring='r2', pre_dispatch=njobs,
                              cv=folds, n_jobs=njobs, refit=True)
    tune_model = tune_model.fit(traindata, outputs)
    return tune_model


def load_training_data(filepath, prefix):
    """
    :param filepath:
    :param prefix:
    :return:
    """
    all_groups = get_valid_hdf5_groups(filepath, prefix)
    with pd.HDFStore(filepath, 'r') as hdf:
        full_dset = pd.concat([hdf[g] for g in all_groups], ignore_index=True)
    outputs = full_dset.loc[:, 'y_depvar']
    feat_order = sorted([ft for ft in full_dset.columns if ft.startswith('ft')])
    traindata = full_dset.loc[:, feat_order]
    return feat_order, traindata, outputs


def load_model(modname, modpath):
    """
    :param modname:
    :param modpath:
    :return:
    """
    module = imp.import_module(modpath)
    model = module.__dict__[modname]()
    return model


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
        this_comb['scores'] = {'mean': np.round(mean, 2), 'std': np.round(std, 2),
                               'max': np.round(max_score, 2), 'min': np.round(min_score, 2)}
        this_comb['params'] = entry.parameters
        grid.append(this_comb)
    return grid


def run_train_model(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    logger.debug('Loading model specification from {}'.format(args.modelspec))
    model_params = json.load(open(args.modelspec))
    model = load_model(model_params['model_name'], model_params['module_path'])
    feat_order, traindata, outputs = load_training_data(args.traindata, args.traingroup)
    if args.notuning:
        params = model_params['default']
        model = train_nocv(model, params, traindata, outputs)
        metadata = {'final_params': params}
    else:
        params = model_params['cvtune']
        tune_info = train_gridcv(model, params, traindata, outputs, args.cvfolds, args.workers)
        model = tune_info.best_estimator_
        metadata = {'final_params': tune_info.best_params_}
        metadata['grid_scores'] = simplify_cv_scores(tune_info.grid_scores_)
        metadata['best_score'] = tune_info.best_score_
        metadata['scoring'] = 'r2_score'

    metadata['model_spec'] = os.path.basename(args.modelspec)
    metadata['init_params'] = params
    metadata['model'] = model_params['model_name']
    metadata['traindata'] = os.path.basename(args.traindata)
    metadata['traingroup'] = args.traingroup
    metadata['features'] = feat_order
    metadata['modelfile'] = os.path.basename(args.modelout)

    with open(args.modelout, 'wb') as outfile:
        pck.dump(model, outfile)

    if not args.metadataout:
        mdout = args.modelout.rsplit('.', 1)[0] + '.json'
    else:
        mdout = args.metadataout
    with open(mdout, 'w') as outfile:
        _ = json.dump(metadata, outfile, indent=0)

    return 0
