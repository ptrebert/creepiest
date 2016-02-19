# coding=utf-8

import pandas as pd
import numpy as np


def quant_norm(*args, statistic='mean'):
    """ Perform quantile normalization as proposed by Bolstad et al.
    Specifically, this method mimicks the behavior of the R library preprocessCore,
    i.e. of the function normalize.quantiles (v1.32) [NB here: samples X values]
    Ties when ranking the data are broken as follows during the replacement step:
    value(rk=N.f) = (value(rk=N) + value(rk=N+1)) / 2
    More conservative normalization is possible using the median.
    This function does not handle missing data.

    :param args:
     :type: numpy.ndarray, list
    :param statistic: which statistic to use; default: mean
     :type: str
    :return: the quantile normalized data
     :rtype: pandas.core.frame.DataFrame
    """
    # dataframe with samples (rows) X values (cols)
    df = pd.DataFrame([*args], dtype='float64')
    # compute ranks for each sample
    # use average to identify ties later
    rk = df.rank(axis=1, method='average')
    # sort values in each sample (row)
    df = df.apply(np.sort, axis=1)
    # generate series with data (statistic of values, for each col)
    # and index (= ranking of statistic)
    if statistic == 'mean':
        col_stats = pd.Series(df.mean(axis=0).values, index=df.mean(axis=0).rank(method='dense'))
    elif statistic == 'median':
        col_stats = pd.Series(df.median(axis=0).values, index=df.median(axis=0).rank(method='dense'))
    else:
        raise ValueError('Unknown statistic: {}'.format(statistic))
    # drop duplicates to get single element when accessing col_means
    col_stats.drop_duplicates(keep='first', inplace=True)
    # sort for faster access
    col_stats.sort_index(inplace=True)
    # get indices of those elements that have an averaged rank (e.g. 3.5)
    idx_ties = list(rk[~rk.isin(col_stats.index)].stack().index)
    # replace all non-average ranks with corresponding means
    rk.replace(to_replace=col_stats.index, value=col_stats.values, inplace=True)
    for row, col in idx_ties:
        # replace values with an averaged rank with the average of the two means (floor and ceiling positon)
        rk.iat[row, col] = (col_stats.loc[int(rk.iat[row, col])] + col_stats.loc[int(rk.iat[row, col]) + 1]) / 2.
    return rk