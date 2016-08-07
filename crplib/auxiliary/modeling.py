# coding=utf-8

import os as os
import json as json
import importlib as imp
import pickle as pck

import sklearn.metrics as sklmet

from crplib.auxiliary.text_parsers import text_file_mode


def get_scorer(name, smpwt=None):
    """
    :param name:
    :param smpwt:
    :return:
    """
    if name == 'f1_score' or name == 'f1':
        scorer = sklmet.make_scorer(sklmet.f1_score, average='weighted', sample_weight=smpwt, pos_label=None)
    elif name == 'accuracy_score' or name == 'accuracy':
        scorer = sklmet.make_scorer(sklmet.accuracy_score, normalize=True, sample_weight=smpwt)
    elif name == 'roc_auc_score' or name == 'roc_auc':
        scorer = sklmet.make_scorer(sklmet.roc_auc_score, average='weighted', sample_weight=smpwt)
    elif name == 'mse' or name == 'mean_squared_error':
        scorer = sklmet.make_scorer(sklmet.mean_squared_error, sample_weight=smpwt, greater_is_better=False)
    elif name == 'r2_score' or name == 'r2':
        scorer = sklmet.make_scorer(sklmet.r2_score, sample_weight=smpwt, multioutput='uniform_average')
    else:
        try:
            scorer = sklmet.make_scorer(getattr(sklmet, name), sample_weight=smpwt)
        except AttributeError:
            raise AttributeError('Requested scoring function {} does not exist in sklearn.metrics'.format(name))
    return scorer


def load_subset_names(fpath):
    """
    :param fpath:
    :return:
    """
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode) as infile:
        try:
            content = json.load(infile)['subset']
        except json.JSONDecodeError:
            content = infile.read().split()
    content = set(content)
    assert len(content) > 1, 'Subset selecting names list loaded from file {} has length <=1: {}'.format(fpath, content)
    return content


def select_dataset_subset(dataset, subset):
    """
    :param dataset:
    :param subset:
    :return:
    """
    rows, cols = dataset.shape
    if not subset:
        pass
    elif os.path.isfile(subset):
        subset_names = load_subset_names(subset)
        dataset = dataset[dataset.name.isin(subset_names)]
    else:
        subset = subset.strip('"')
        dataset = dataset.query(subset)
    assert not dataset.empty, 'Dataset empty after subsetting with {}; initial size {} x {}'.format(subset, rows, cols)
    return dataset


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


def load_model(path, name=None):
    """
    :param path:
    :param name:
    :return:
    """
    if os.path.isfile(path):
        with open(path, 'rb') as modelfile:
            model = pck.load(modelfile)
    else:
        assert name is not None, 'Need model name to load model from module {}'.format(path)
        module = imp.import_module(path)
        model = getattr(module, name)()
    return model
