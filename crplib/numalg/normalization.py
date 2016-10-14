# coding=utf-8

import pandas as pd
import numpy as np
import numpy.ma as msk
import scipy.stats as stats


def quant_norm_deprecated(*args, statistic='mean', nzonly=False):
    """ Perform quantile normalization as proposed by Bolstad et al.
    Since the handling of ties during the ranking is not specified in the original paper,
    this method uses dense ranking (rank always increases by 1 between groups)
    Hence, the result of this function is slightly different to, e.g., the R function
    normalize.quantiles from the preprocessCore package

    :param args:
     :type: numpy.ndarray, tuple of numpy.ndarray
    :param statistic: which statistic to use; default: mean
     :type: str
    :param nzonly: reduce data to columns where at least one sample has a non-zero entry
    :return: the quantile normalized data
     :rtype: pandas.core.frame.DataFrame
    """
    #
    # This method consumes too much memory (Pandas overhead???)
    # Use only for small data sizes (num data points < 10e6)
    #
    # dataframe with samples (rows) X values (cols)
    # avoid: "can use starred expression only as assignment target"
    # df = pd.DataFrame([*args], dtype='float64')
    df = pd.DataFrame([arg for arg in args], dtype='float64')
    if nzonly:
        raise NotImplementedError('Using only non-zero columns not yet implemented')
    # compute ranks for each sample
    rk = df.rank(axis=1, method='dense')
    # sort values in each sample (row)
    df = df.apply(np.sort, axis=1, raw=True)
    # generate series with data (statistic of values, for each col)
    # and index (= ranking of statistic)
    if statistic == 'mean':
        col_stats = pd.Series(df.mean(axis=0).values, index=df.mean(axis=0).rank(method='dense'))
    elif statistic == 'median':
        col_stats = pd.Series(df.median(axis=0).values, index=df.median(axis=0).rank(method='dense'))
    else:
        raise ValueError('Unknown statistic: {}'.format(statistic))
    # sort for faster access
    # ...could be premature optimization here...
    col_stats.sort_index(inplace=True)
    # replace rank values with averages
    rk.replace(to_replace=col_stats.index, value=col_stats.values, inplace=True)
    return rk


def quant_norm(*args):
    """ Perform quantile normalization as proposed by Bolstad et al.
    Since the handling of ties during the ranking is not specified in the original paper,
    this method uses dense ranking (rank always increases by 1 between groups)
    Hence, the result of this function is slightly different to, e.g., the R function
    normalize.quantiles from the preprocessCore package

    :param args: arbitrary number of numpy arrays with dtype float64
    :return: quantile normalized arrays, same number as input
    """
    ranks = list()
    for arg in args:
        rk = stats.rankdata(arg, method='dense')
        ranks.append(rk)
        arg.sort()
    col_stat = np.mean([arg for arg in args], axis=0, dtype='float64')
    col_stat = np.unique(col_stat)
    col_ranks = stats.rankdata(col_stat, method='dense')
    for idx, rarr in enumerate(ranks):
        indices = np.digitize(rarr, col_ranks, right=True)
        rarr = col_stat[indices]
        ranks[idx] = rarr
    return ranks


def merge_1d_datasets(*args, mergestat, qnorm):
    """
    :param args:
    :param mergestat:
    :param qnorm:
    :return:
    """
    if qnorm:
        args = quant_norm(*args)
    mergers = {'mean': np.mean, 'median': np.median, 'max': np.max, 'min': np.min}
    col_mrg = mergers[mergestat]([arg for arg in args], axis=0)
    return col_mrg


def transform_to_pct_ranks(signal):
    """
    :param signal:
    :return:
    """
    m = msk.masked_where(signal == 0., signal)
    percentiles = np.percentile(m.compressed(), [10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    ranks = np.digitize(signal, percentiles, right=True)
    # this is necessary to set all zero values to rank 0 again
    ranks += 1
    ranks[m.mask] = 0
    ranks = ranks.astype(np.int32, copy=False)
    assert ranks.min() == 0, 'Minimum rank is not 0: {}'.format(ranks.min())
    assert ranks.max() == 10, 'Maximum rank is not 10: {}'.format(ranks.max())
    return ranks
