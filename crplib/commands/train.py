# coding=utf-8

"""
Command to train models on predefined training datasets
"""

import json as json
import importlib as imp
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
    model.fit(traindata, outputs)
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
    tune_model.fit(traindata, outputs)
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
    model = module.__dict__[modname]
    return model


def run_train_model(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    logger.debug('Loading model parameters from file {}'.format(args.modelfile))
    model_params = json.load(open(args.modelfile))
    model = load_model(model_params['model_name'], model_params['module_path'])
    feat_order, traindata, outputs = load_training_data(args.inputfile, args.inputgroup)
    if args.nocv:
        params = model_params['default']
        model = train_nocv(model, params, traindata, outputs)
    else:
        params = model_params['cvtune']
        model = train_gridcv(model, params, traindata, outputs, args.folds, args.workers)

    return 0
