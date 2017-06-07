# coding=utf-8

import numpy as np
import numpy.random as rng
import numpy.ma as msk
import multiprocessing as mp
import scipy.stats as stats

import pandas as pd

_shm_input_matrix = None
_shm_rank_matrix = None


def nonzero_qnorm(mat, workers=1):
    """ Perform quantile normalization as proposed by Bolstad et al., but restrict to columns
    in matrix that have at least one non-zero entry, i.e., exclude columns with no signal.
    Since the handling of ties during the ranking is not specified in the original paper,
    this method uses min ranking.
    Hence, the result of this function is different to, e.g., the R function
    normalize.quantiles from the preprocessCore package
    Intended for chromosome-wide normalizations of, e.g., histone data that commonly
    have a substantial amount of all-zero positions in the respective signal tracks

    :param mat: [rows] samples X features [columns]
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
    assert np.isfinite(mat).all(), 'Matrix contains invalid (NaN, INF) values'
    mat = mat.round(decimals=3)
    add_zero_columns = False
    col_idx = np.array(mat.sum(axis=0) > 0, dtype=np.bool)
    non_zero_cols = np.sum(col_idx)
    if non_zero_cols == 0:
        return mat
    orig_shape = mat.shape
    if non_zero_cols != mat.shape[1]:
        # at least one all-zero column
        add_zero_columns = True
        mat = mat[:, col_idx]
    # determine data ranks
    data_ranks = np.zeros(mat.shape, dtype=np.int32)
    for row_idx in range(mat.shape[0]):
        data_ranks[row_idx, :] = stats.rankdata(mat[row_idx, :], method='min')
    data_ranks = data_ranks.astype(np.int32)
    uniq_data_ranks = np.unique(data_ranks)

    # get column means
    mat.sort(axis=1)
    col_means = mat.mean(axis=0)
    col_means = col_means.round(decimals=3).astype(np.float32)

    # ranks for column means
    col_mean_ranks = stats.rankdata(col_means, method='min')
    col_mean_ranks = col_mean_ranks.astype(np.int32)

    # lookup table for rank -> mean replacement
    replace_means = pd.Series(col_means, index=col_mean_ranks,
                              dtype=np.float32)
    # reduce to unique entries (i.e., group by index value = rank of mean value)
    replace_means = replace_means.groupby(by=replace_means.index).first()

    # Check if there are ranks missing
    # There are cases where data_ranks contains values that do
    # not appear in col_mean_ranks - apparently, the ranking of
    # floating point numbers is the issue here
    missing_ranks = set(uniq_data_ranks).difference(set(replace_means.index))
    if len(missing_ranks) > 0:
        missing_means = []
        missing_positions = []
        for mr in sorted(missing_ranks):
            # this "cannot" raise a ValueError since this would
            # imply that missing ranks contains rank 1; cannot happen
            # if no other error causes this situation
            lo = replace_means[replace_means.index < mr].idxmax()
            lo_val = replace_means[lo]
            try:
                # this happened indeed once for hg19/H3K36me3
                hi = replace_means[replace_means.index > mr].idxmin()
                hi_val = replace_means[hi]
            except ValueError:
                # catches: argmin of empty sequence
                hi_val = replace_means.max()
            assert lo_val <= hi_val,\
                'Invalid combination of low and high replacement means: {} and {}'.format(lo_val, hi_val)
            new_val = np.round((lo_val + hi_val) / 2., 3)
            missing_means.append(new_val)
            missing_positions.append(mr)

        additional_means = pd.Series(missing_means, index=missing_positions, dtype=np.float32)
        replace_means = replace_means.append(additional_means)

    for row_idx in range(mat.shape[0]):
        # by constructions, this must give a valid replacement vector
        mat[row_idx, :] = replace_means[data_ranks[row_idx, :]]

    if add_zero_columns:
        add_zeros = np.zeros(orig_shape, dtype=np.float32)
        col_indices = np.arange(col_idx.size)[col_idx]
        add_zeros[:, col_indices] = mat[:]
        mat = add_zeros
    # Note that this is ok since all zero matrices are returned
    # immediately in the beginning
    assert mat.max() > 0, 'Normalization failed: ends with all zero data matrix'
    assert np.isfinite(mat).all(), 'Result matrix contains invalid (NaN, INF) values'
    return mat


def simple_quantile_normalization(mat):
    """
    This code: simple q-norm implementation ignoring issues like breaking ties,
    normalizing against target distribution etc.

    Code found at:
    http://lists.open-bio.org/pipermail/biopython/2010-March/012505.html

    Presumed author: Vincent Davis

    Original docstring:
    "anarray with samples in the columns and probes across the rows"

    My changes to improve readability:
    data := A
    new_data := AA
    ranks := I

    :param mat:
    :return:
    """
    new_data = np.zeros_like(mat, dtype=mat.dtype)
    ranks = np.argsort(mat, axis=0)
    new_data[ranks, np.arange(mat.shape[1])] = np.mean(mat[ranks, np.arange(mat.shape[1])], axis=1)[:, np.newaxis]
    return new_data


def preprocess_core_qnorm(mat):
    """
    This method is just a mock-up to keep a "reference" to R's preprocessCore
    package w/o introducing an additional dependency

    R version 3.1.2 (2014-10-31) -- "Pumpkin Helmet"
    preprocessCore_1.28.0

    # Example from Wikipedia page
    # note that the Wikipedia example uses "min" as tie breaking
    # rule whereas the R default is "average"
    # Presumably, this default also applies to preprocessCore ?

    > mat <- matrix(c(5,2,3,4,4,1,4,2,3,4,6,8), nrow=4, ncol=3)
    > normalize.quantiles(mat)
         [,1]     [,2]     [,3]
    [1,] 5.666667 5.166667 2.000000
    [2,] 2.000000 2.000000 3.000000
    [3,] 3.000000 5.166667 4.666667
    [4,] 4.666667 3.000000 5.666667

    # Example by Rafael Irizarry posted on Twitter
    # http://t.co/lCyy6YyNB8

    > mat <- matrix(c(2,5,4,3,3,4,14,8,8,9,4,4,6,5,3,5,7,9,8,5), nrow=5, ncol=4)
    > normalize.quantiles(mat)
        [,1] [,2] [,3] [,4]
    [1,] 3.50 3.50 5.25 4.25
    [2,] 8.50 8.50 5.25 5.50
    [3,] 6.50 5.25 8.50 8.50
    [4,] 5.25 5.25 6.50 6.50
    [5,] 5.25 6.50 3.50 4.25

    :param mat: data matrix with layout samples (col) X probes (row)
    :return:
    """
    wiki_example = np.array([5, 4, 3,
                             2, 1, 4,
                             3, 4, 6,
                             4, 2, 8], dtype=np.float32).reshape(4, 3)
    irizarry_example = np.array([2, 4, 4, 5,
                                 5, 14, 4, 7,
                                 4, 8, 6, 9,
                                 3, 8, 5, 8,
                                 3, 9, 3, 5], dtype=np.float32).reshape(5, 4)
    if wiki_example.shape == mat.shape and np.allclose(mat, wiki_example):
        wiki_ex_result = np.array([5.666667, 5.166667, 2,
                                   2, 2, 3,
                                   3, 5.166667, 4.666667,
                                   4.666667, 3, 5.666667], dtype=np.float32).reshape(4, 3)
        return wiki_ex_result
    elif irizarry_example.shape == mat.shape and np.allclose(mat, irizarry_example):
        # This result matrix is the result when applying the R function
        # as indicated in the docstring
        irizarry_ex_result = np.array([3.5, 3.5, 5.25, 4.25,
                                       8.5, 8.5, 5.25, 5.5,
                                       6.5, 5.25, 8.5, 8.5,
                                       5.25, 5.25, 6.5, 6.5,
                                       5.25, 6.5, 3.5, 4.25], dtype=np.float32).reshape(5, 4)
        # This result is the one given in Irizarry's Twitter post
        # Why this is different is unclear, could be again a change
        # in the tie breaking strategy
        # irizarry_ex_result2 = np.array([3.5, 3.5, 5., 5.,
        #                                8.5, 8.5, 5.5, 5.5,
        #                                6.5, 5., 8.5, 8.5,
        #                                5., 5.5, 6.5, 6.5,
        #                                5.5, 6.5, 3.5, 3.5], dtype=np.float64).reshape(5, 4)
        return irizarry_ex_result
    else:
        return np.zeros_like(mat)


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
