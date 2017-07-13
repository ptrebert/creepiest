# coding=utf-8

import os as os
import json as json
import importlib as imp
import pickle as pck
import numpy as np
import functools as fnt

import pandas as pd
import sklearn.metrics as sklmet
import scipy.stats as stats

from crplib.auxiliary.text_parsers import text_file_mode
from crplib.mlfeat.featdef import get_prefix_list, get_classes_from_names


def kendall_tau_scorer(x, y):
    """
    :param x:
    :param y:
    :return:
    """
    statistic, p_value = stats.kendalltau(x, y, nan_policy='raise')
    return statistic


def get_scorer(name, smpwt=None):
    """
    :param name:
    :param smpwt:
    :return:
    """
    if smpwt is not None:
        assert hasattr(smpwt, '__iter__'), 'Sample weights not supplied as iterable: {}'.format(smpwt)
        assert len(smpwt) > 0, 'Sample weights supplied as empty iterable - has to be None or not empty'
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
    elif name == 'kendall_tau' or 'kendall_tau_score':
        scorer = sklmet.make_scorer(kendall_tau_scorer, greater_is_better=True)
    else:
        try:
            scorer = sklmet.make_scorer(getattr(sklmet, name), sample_weight=smpwt)
        except AttributeError:
            raise AttributeError('Requested scoring function {} does not exist in sklearn.metrics'.format(name))
    return scorer


def determine_scoring_method(args, metadata, logger):
    """
    :param args:
    :param metadata:
    :param logger:
    :return:
    """
    scm = ''
    if hasattr(args, 'scoring'):
        scm = args.scoring
    if scm and 'scoring' in metadata['model_info']:
        logger.warning('Storing scoring method as part of the model_info is deprecated')
        model_scm = metadata['model_info']['scoring']
        if scm != model_scm:
            logger.warning('Scoring parameter set to {}, but loaded model metadata say {}.'
                           ' Could that be an error, human? Using scoring method: {}'.format(scm, model_scm, scm))
    elif scm and 'scoring' in metadata['training_info']:
        training_scm = metadata['training_info']['scoring']
        if scm != training_scm:
            logger.warning('Scoring parameter set to {}, but {} was used for model training.'
                           ' Could that be an error, human? Using scoring method: {}'.format(scm, training_scm, scm))
    elif scm:
        pass
    elif not scm:
        if 'scoring' in metadata['model_info']:
            scm = metadata['model_info']['scoring']
        else:
            try:
                scm = metadata['training_info']['scoring']
            except KeyError:
                logger.error('Could not find scoring method information in metadata')
                raise
    else:
        raise ValueError('Logic error in determine_scoring_method')
    return scm


def load_ml_dataset(fpath, groups, modelfeat, args, logger):
    """
    :param fpath:
    :param groups:
    :param modelfeat:
    :param args:
    :param logger:
    :return:
    """
    trgname, wtname, usedtype = None, None, 'n/a'
    dataset_info = dict()
    with pd.HDFStore(fpath, 'r') as hdf:
        metadata = hdf['metadata']
        if isinstance(groups, (list, tuple)):
            dataset = pd.concat([hdf[grp] for grp in sorted(groups)], ignore_index=True)
            anygroup = groups[0]
        else:
            dataset = hdf[groups]
            anygroup = groups
    logger.debug('Raw dataset loaded')
    dataset_info['rows'] = dataset.shape[0]
    dataset_info['columns'] = dataset.shape[1]
    if hasattr(args, 'task') and args.task == 'est':
        pass
    else:
        logger.debug('Determining target variable for dataset')
        dataset, trgname, usedtype = augment_with_target(dataset, args.targetvar, args.derivetarget)
    datrows = dataset.shape[0]
    datcols = dataset.shape[1]
    assert datrows > 1, 'No samples (rows) remain in dataset after determining the target variable'
    assert datcols > 1, 'No features (columns) remain in dataset after determining target variable'
    if hasattr(args, 'loadweights') and args.loadweights:
        logger.debug('Loading sample weights from file')
        dataset, wtname = augment_with_weights_file(dataset, args.loadweights)
    dataset_info['subset'] = args.subset
    logger.debug('Selecting subset of data (if applicable)')
    dataset = select_dataset_subset(dataset, args.subset, args.crpmetadata)
    datrows = dataset.shape[0]
    datcols = dataset.shape[1]
    assert datrows > 1, 'No samples (rows) remain in dataset after subset selection'
    assert datcols > 1, 'No features (columns) remain in dataset after subset selection'
    if trgname is not None:
        targets = dataset.loc[:, trgname].values
        assert np.unique(targets).size >= 2, 'Less than 2 unique targets remain - maybe, you selected' \
                                             ' a pointless subset of the data, human?'
        store_targets = targets.tolist()
    else:
        targets = None
        store_targets = None
    if wtname is not None:
        weights = dataset.loc[:, wtname].values
    else:
        weights = None
    sample_names = extract_sample_names(dataset)
    if modelfeat is not None:
        logger.debug('Selecting features based on model information')
        datafeat = [c for c in dataset.columns if c.startswith('ft')]
        missing = set(modelfeat) - set(datafeat)
        assert len(missing) == 0, 'Missing features in dataset: {}'.format(missing)
        predictors = dataset.loc[:, modelfeat]
        assert predictors.shape[1] > 1, 'No features selected using model information: {}'.format(modelfeat)
        feat_info = {'order': modelfeat}
    else:
        logger.debug('Extracting feature information')
        feat_info = extract_feature_information(dataset.columns, metadata, anygroup)
        if hasattr(args, 'usefeatures') and args.usefeatures:
            prefixes = get_prefix_list(args.usefeatures)
            feat_info['order'] = list(filter(lambda x: any([x.startswith(p) for p in prefixes]), feat_info['order']))
            feat_info['classes'] = args.usefeatures
        predictors = dataset.loc[:, feat_info['order']]
        assert predictors.shape[1] > 1, 'No features selected from dataset columns'
    assert predictors.notnull().all(axis=1).all(), 'Detected invalid NULL values in predictor matrix'
    logger.debug('Assembling metadata')
    dataset_info['num_samples'] = predictors.shape[0]
    dataset_info['num_features'] = predictors.shape[1]
    dataset_info['target_var'] = trgname
    dataset_info['derive_target'] = args.derivetarget
    dataset_info['target_type'] = usedtype
    sample_info = {'weights': weights, 'names': sample_names, 'targets': store_targets}
    logger.debug('Dataset loading done')
    return predictors, targets, dataset_info, sample_info, feat_info


def extract_feature_information(datacols, mdframe, onegroup):
    """
    :param datacols:
    :param mdframe:
    :return:
    """
    feat_order = sorted([ft for ft in datacols if ft.startswith('ft')])
    assert len(feat_order) > 0, 'No features selected from column names: {}'.format(datacols)
    if 'features' in mdframe.columns:
        grp_index = mdframe.where(mdframe.group == onegroup).dropna().index[0]
        # regular training dataset
        ft_classes = mdframe.loc[grp_index, 'features']
        ft_classes = ft_classes.values[0].split(',')
        ft_kmers = mdframe.loc[grp_index, 'kmers']
        ft_kmers = list(map(int, ft_kmers.values[0].split(',')))
        res = mdframe.loc[grp_index, 'resolution']
        res = res.values[0]
    else:
        # apparently, just regions
        ft_classes = get_classes_from_names(feat_order)
        ft_kmers = []
        res = 0
    feat_info = {'order': feat_order, 'classes': ft_classes,
                 'kmers': ft_kmers, 'resolution': res}
    return feat_info


def augment_with_target(data, trgname, expr):
    """
    :param data:
    :param trgname:
    :param expr:
    :return:
    """
    if trgname:
        assert trgname in data.columns, 'Supplied name of target variable {} is not a column name'.format(trgname)
        trg_dtype = data[trgname].values.dtype
        if np.issubdtype(bool, trg_dtype):
            usedtype = 'bool/int8'
        elif np.issubdtype(float, trg_dtype):
            usedtype = 'float/float64'
        elif np.issubdtype(int, trg_dtype):
            usedtype = 'int/int64'
        else:
            raise TypeError(
                'Cannot infer type ({}) of derived target (expr.: {})'.format(data.target.values.dtype, expr))
    else:
        assert 'data' in expr, 'Could not find "data" keyword in expression to evaluate: {}'.format(expr)
        data = data.assign(target=pd.eval(expr.strip('"')))
        # numpy.issubdtype(arg1, arg2)
        # Returns True if first argument is a typecode lower/equal in type hierarchy.
        if np.issubdtype(bool, data.target.values.dtype):
            data.target = data.target.astype(np.int8, copy=False)
            usedtype = 'bool/int8'
        elif np.issubdtype(float, data.target.values.dtype):
            data.target = data.target.astype(np.float64, copy=False)
            usedtype = 'float/float64'
        elif np.issubdtype(int, data.target.values.dtype):
            data.target = data.target.astype(np.int64, copy=False)
            usedtype = 'int/int64'
        else:
            raise TypeError('Cannot infer type ({}) of derived target (expr.: {})'.format(data.target.values.dtype, expr))
        trgname = 'target'
    assert data[trgname].unique().size >= 2, 'Number of unique targets smaller than 2 for name {} / expr {}'.format(trgname, expr)
    return data, trgname, usedtype


def augment_with_weights_file(data, wtfile):
    """
    :param data:
    :param wtfile:
    :return:
    """
    loaded_weights = load_sample_weights(wtfile)
    if 'weight' in data.columns:
        data.drop('weight', axis='columns', inplace=True)
    data = data.merge(loaded_weights, how='outer', on='name', suffixes=('', ''), copy=False)
    return data, 'weight'


def extract_sample_names(dataset):
    """
    :param data:
    :return:
    """
    name_col = [cn for cn in dataset.columns if cn in ['name', 'source']]
    sample_names = []
    if name_col:
        sample_names = dataset.loc[:, name_col[0]].tolist()
    return sample_names


def load_sample_weights(fpath):
    """
    :param fpath:
    :return:
    """
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode) as mdfile:
        try:
            annotation = json.load(mdfile)
            smpwt = pd.DataFrame(annotation['sample_weights'], columns=['name', 'weight'])
        except json.JSONDecodeError:
            smpwt = pd.read_csv(mdfile, sep='\t', header=0, names=['name', 'weight'])
    assert not smpwt.empty, 'Entry sample_weights in annotation file {} is empty'.format(fpath)
    return smpwt


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


def load_names_from_metadata(fpath, expr):
    """
    :param fpath:
    :param expr:
    :return:
    """
    with open(fpath, 'r') as infile:
        try:
            md = json.load(infile)
            assert 'md' in expr, 'Missing metdata keyword "md" in expression: {}'.format(expr)
            sample_names = md['sample_info']['names']
            selector = eval(expr)
            assert len(sample_names) == len(selector), 'Evaluating the expression > {} < resulted in {} selectors,' \
                                                       ' but there are {} sample names' \
                                                       ' in the metadata file {}'.format(expr, len(selector), len(sample_names), fpath)
            subset_names = [n for s, n in zip(selector, sample_names) if s]
        except json.JSONDecodeError:
            raise TypeError('Expected CREEPIEST metadata file (JSON) - cannot decode this file: {}'.format(fpath))
    return subset_names


def select_dataset_subset(data, subset, crpmd=None):
    """
    :param data:
    :param subset:
    :param crpmd:
    :return:
    """
    rows, cols = data.shape
    if not subset:
        pass
    elif os.path.isfile(subset):
        if crpmd is not None:
            subset_names = load_names_from_metadata(subset, crpmd)
        else:
            subset_names = load_subset_names(subset)
        data = data[data.name.isin(subset_names)]
    else:
        subset = subset.strip('"')
        data = data.query(subset)
    assert not data.empty, 'Dataset empty after subsetting with {}; initial size {} x {}'.format(subset, rows, cols)
    return data


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


def extract_model_attributes(model, attributes, logger=None):
    """
    :param model:
    :param attributes:
    :return:
    """
    attribs = dict()
    for attr in attributes:
        if hasattr(model, attr):
            attr_obj = getattr(model, attr)
            if hasattr(attr_obj, 'tolist'):  # numpy array
                attr_obj = attr_obj.tolist()
            elif np.issubdtype(attr_obj, bool):
                attr_obj = bool(attr_obj)  # that should not be unnecessary
            elif np.issubdtype(attr_obj, int):
                attr_obj = int(attr_obj)
            elif np.issubdtype(attr_obj, float):
                attr_obj = float(attr_obj)
            else:
                if logger is not None:
                    logger.warning('Cannot cast type {} of model attribute {}'
                                   ' to Python primitive - skipping'.format(type(attr_obj), attr))
            attribs[attr] = attr_obj
        else:
            if logger is not None:
                logger.debug('Skipping attribute {} - does not exist'.format(attr))
    if logger is not None:
        logger.debug('Extracted {} attributes of {} requested ones'.format(len(attribs), len(attributes)))
    return attribs
