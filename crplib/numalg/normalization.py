# coding=utf-8

import numpy as np
import numpy.random as rng
import numpy.ma as msk
import scipy.stats as stats


def nonzero_qnorm(mat):
    """
    Perform quantile normalization as below, but restrict to columns in mat
    that have at least one non-zero entry, i.e., exclude all all-zero columns
    Intended for, e.g., chromosome-wide normalizations of histone data that
    are mostly zero
    :param mat:
    :return:
    """
    add_zero_columns = False
    col_idx = mat.sum(axis=0) > 0
    if np.sum(col_idx) != mat.shape[1]:
        # at least one all-zero column
        add_zero_columns = True
        mat = mat[:, col_idx]
    ranks = np.zeros_like(mat)
    for row_idx in range(mat.shape[0]):
        ranks[row_idx, :] = stats.rankdata(mat[row_idx, :], method='dense')
    mat.sort(axis=1)
    col_means = np.unique(mat.mean(axis=0))
    mean_ranks = stats.rankdata(col_means, method='dense')
    for row_idx in range(ranks.shape[0]):
        indices = np.digitize(ranks[row_idx, :], mean_ranks, right=True)
        ranks[row_idx, :] = col_means[indices]
    norm_mat = ranks
    if add_zero_columns:
        add_zeros = np.zeros((ranks.shape[0], col_idx.size))
        col_indices = np.arange(col_idx.size)[col_idx]
        add_zeros[:, col_indices] = ranks[:]
        norm_mat = add_zeros
    return norm_mat


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


def transform_to_dec_ranks(signal):
    """
    :param signal:
    :return:
    """
    m = msk.masked_where(signal == 0., signal)
    try:
        percentiles = np.percentile(m.compressed(), [10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        ranks = np.digitize(signal, percentiles, right=True)
    except ValueError as ve:
        if isinstance(ve.args[0], str) and ve.args[0].startswith('bins must be monotonically'):
            # very "flat" signal, add a little bit of noise
            # Observed for E015_mm9_ESE14_H3K4me1
            # the following raises again if it's still to flat
            m += (rng.random(m.size) / 100.)
            percentiles = np.percentile(m.compressed(), [10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
            ranks = np.digitize(signal, percentiles, right=True)
        else:
            raise ve
    # this is necessary to set all zero values to rank 0 again
    ranks += 1
    ranks[m.mask] = 0
    ranks = ranks.astype(np.int32, copy=False)
    assert ranks.min() == 0, 'Minimum rank is not 0: {}'.format(ranks.min())
    assert ranks.max() == 10, 'Maximum rank is not 10: {}'.format(ranks.max())
    return ranks
