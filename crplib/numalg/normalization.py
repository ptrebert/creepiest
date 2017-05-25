# coding=utf-8

import numpy as np
import numpy.random as rng
import numpy.ma as msk
import multiprocessing as mp
import scipy.stats as stats

import time as ti


_shm_input_matrix = None
_shm_rank_matrix = None


def nonzero_qnorm(mat, workers=1):
    """ Perform quantile normalization as proposed by Bolstad et al., but restrict to columns
    in matrix that have at least one non-zero entry, i.e., exclude columns with no signal.
    Since the handling of ties during the ranking is not specified in the original paper,
    this method uses dense ranking (rank always increases by 1 between groups)
    Hence, the result of this function is slightly different to, e.g., the R function
    normalize.quantiles from the preprocessCore package
    Intended for chromosome-wide normalizations of, e.g., histone data that commonly
    have a substantial amount of all-zero positions in the respective signal tracks

    :param mat: [rows] samples X genes [columns]
    :param workers: number of CPU core to use, default 1
    :return:
    """
    if workers <= 1:
        return nonzero_qnorm_sequential(mat)
    else:
        return nonzero_qnorm_parallel(mat, workers)


def nonzero_qnorm_sequential(mat):
    """
    :param mat:
    :return:
    """
    add_zero_columns = False
    col_idx = mat.sum(axis=0) > 0
    if np.sum(col_idx) != mat.shape[1]:
        # at least one all-zero column
        add_zero_columns = True
        mat = mat[:, col_idx]
    ranks = np.zeros_like(mat, dtype=np.float64)
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


def _rank_sort_matrix(params):
    """
    :param params:
    :return:
    """
    global _shm_input_matrix
    global _shm_rank_matrix
    row_idx, num_cols = params
    start, end = row_idx * num_cols, row_idx * num_cols + num_cols
    _shm_rank_matrix[start:end] = stats.rankdata(_shm_input_matrix[start:end], method='dense')
    _shm_input_matrix[start:end] = np.sort(_shm_input_matrix[start:end])
    return


def _mean_replace_matrix(params):
    """
    :param params:
    :return:
    """
    row_idx, num_cols, ranks, means = params
    global _shm_rank_matrix
    start, end = row_idx * num_cols, row_idx * num_cols + num_cols
    indices = np.digitize(_shm_rank_matrix[start:end], ranks, right=True)
    _shm_rank_matrix[start:end] = means[indices]
    return


def nonzero_qnorm_parallel(mat, workers):
    """ Current issues: multiprocessing.Array is just a wrapper around
    the Python Array module, not around numpy array (quite logically so, actually...).
    It follows that multidimensional arrays are not supported.
    All the copying and reshaping here consumes more
    memory and make the code slower than the sequential version. In order
    to avoid that, need to rewrite rest of the code or get the numpy.ctypeslib
    to work (currently throws an error about strided arrays not being supported)

    :param mat:
    :param workers:
    :return:
    """
    raise RuntimeError('Abort: parallel version of q-norm too inefficient to be used')
    global _shm_input_matrix
    global _shm_rank_matrix
    add_zero_columns = False
    num_rows = mat.shape[0]
    col_idx = mat.sum(axis=0) > 0
    if np.sum(col_idx) != mat.shape[1]:
        # at least one all-zero column
        add_zero_columns = True
        mat = mat[:, col_idx]
        value_cols = mat.shape[1]
    else:
        value_cols = mat.shape[1]
    num_items = num_rows * value_cols
    _shm_input_matrix = mp.Array('d', num_items, lock=False)
    _shm_input_matrix[:] = mat.ravel()
    _shm_rank_matrix = mp.Array('d', num_items, lock=False)
    with mp.Pool(workers) as pool:
        _ = pool.map(_rank_sort_matrix, [(idx, value_cols) for idx in range(num_rows)])
        col_means = np.unique(np.array(_shm_input_matrix[:]).reshape(num_rows, value_cols).mean(axis=0))
        mean_ranks = stats.rankdata(col_means, method='dense')
        _ = pool.map(_mean_replace_matrix, [(idx, value_cols, mean_ranks, col_means) for idx in range(num_rows)])
    if add_zero_columns:
        add_zeros = np.zeros((num_rows, col_idx.size))
        col_indices = np.arange(col_idx.size)[col_idx]
        add_zeros[:, col_indices] = np.array(_shm_rank_matrix[:]).reshape(num_rows, value_cols)
        return add_zeros
    else:
        return np.array(_shm_rank_matrix[:]).reshape(num_rows, value_cols)


def merge_1d_datasets(*args, mergestat, qnorm):
    """
    :param args:
    :param mergestat:
    :param qnorm:
    :return:
    """
    ncols = len(args[0])
    data_matrix = np.array([arg for arg in args], dtype=args[0].dtype)
    mat_rows = data_matrix.shape[0]
    mat_cols = data_matrix.shape[1]
    assert mat_cols == ncols,\
        'Building data matrix from samples failed: {} resulted in {}'.format(ncols, data_matrix.shape)
    # I guess the following should generally be the case for real datasets...
    assert mat_rows < mat_cols,\
        'Data matrix has more columns (genes) than rows (samples): {}'.format(data_matrix.shape)
    # if the merge statistic is mean, than running q-norm first is unnecessary
    # maybe, this should be checked here
    if qnorm:
        data_matrix = nonzero_qnorm(data_matrix)
    mergers = {'mean': np.mean, 'median': np.median, 'max': np.max, 'min': np.min}
    col_mrg = mergers[mergestat](data_matrix, axis=0)
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
