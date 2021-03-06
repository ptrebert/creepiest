# coding=utf-8

"""
Small command module to merge existing datasets or extend them with additional data
"""

import os as os
import pandas as pd
import multiprocessing as mp
from string import ascii_uppercase as asciiup

from crplib.auxiliary.hdf_ops import get_default_group, get_chrom_list,\
    load_data_group, check_path_infos
from crplib.auxiliary.file_ops import create_filepath

from crplib.metadata.md_regions import MD_REGION_COLDEFS
from crplib.metadata.md_regions import gen_obj_and_md as region_generate


def assemble_worker_params(args):
    """
    :param args:
    :return:
    """
    commons = dict(vars(args))
    del commons['module_logger']
    del commons['execute']
    commons['inputfile'] = []
    inputlabels, inputgroups = [], []
    # Note for zip:
    # "The iterator stops when the shortest input iterable is exhausted"
    # no problem when there is only one file
    for deflab, ipf in zip(list(asciiup), args.inputfiles):
        lab, grp, fp = check_path_infos(ipf)
        if lab is None:
            lab = deflab
        commons['inputfile'].append(fp)
        inputlabels.append(lab)
        inputgroups.append(grp)
    commons['inputlabel'] = inputlabels
    commons['inputgroup'] = inputgroups
    req_label = len(args.valfile) > 1
    commons['valfile'] = []
    commons['vallabel'] = []
    commons['valgroup'] = []
    for deflab, vfp in zip(['V' + c for c in list(asciiup)], args.valfile):
        lab, grp, fp = check_path_infos(vfp)
        if lab is None and req_label:
            lab = deflab
        commons['valfile'].append(fp)
        commons['vallabel'].append(lab)
        commons['valgroup'].append(grp)
    arglist = []
    for chrom in get_chrom_list(commons['inputfile'][0]):
        tmp = dict(commons)
        tmp['chrom'] = chrom
        arglist.append(tmp)
    return arglist


def merge_datasets(params):
    """
    :param params:
    :return:
    """
    chrom = params['chrom']
    fp_1, grp_1, lab_1 = params['inputfile'][0], params['inputgroup'][0], params['inputlabel'][0]
    merge_cols = params['mergeon']
    data_labels = params['inputlabel']
    merged = load_data_group(fp_1, grp_1, chrom)
    for fp, grp, lab in zip(params['inputfile'][1:],
                            params['inputgroup'][1:],
                            params['inputlabel'][1:]):
        right_df = load_data_group(fp, grp, chrom)
        for colname in merge_cols:
            shape_left, shape_right = merged[colname].shape, right_df[colname].shape
            assert shape_left == shape_right,\
                'Different number of data entries ({} vs {}) for column {}'.format(shape_left, shape_right, colname)
            assert merged[colname].isin(right_df[colname]).all(),\
                'Some entries are not shared between datasets' \
                ' ({} and {}) for column {}'.format(os.path.basename(fp_1), os.path.basename(fp), colname)
        if not lab_1:
            suffix_left = ''
        else:
            suffix_left = '_' + lab_1
        suffix_right = '_' + lab
        columns_left = set(merged.columns)
        columns_right = set(right_df.columns)
        merged = merged.merge(right_df, how='outer', on=merge_cols, suffixes=(suffix_left, suffix_right), copy=False)
        new_columns = []
        for col in merged.columns:
            if col in merge_cols:
                new_columns.append(col)
            elif any([col.endswith(x) for x in data_labels]):
                new_columns.append(col)
            # the following two take care of column names unique
            # to one dataframe - pandas suffix replacement only
            # works on overlapping column names
            elif col in columns_left:
                ncol = col + suffix_left
                new_columns.append(ncol)
            elif col in columns_right:
                ncol = col + suffix_right
                new_columns.append(ncol)
            else:
                raise AssertionError('Cannot assign column {} to left or right dataframe'.format(col))
        merged.columns = new_columns
        lab_1 = ''
    start_cols = [c for c in merged.columns if c.startswith('start')]
    end_cols = [c for c in merged.columns if c.startswith('end')]
    merged = merged.assign(start=lambda x: x[start_cols].min(axis=1))
    merged = merged.assign(end=lambda x: x[end_cols].max(axis=1))
    return merged


def extend_datasets(dataset, params):
    """
    :param dataset:
    :param params:
    :return:
    """
    chrom = params['chrom']
    merge_cols = params['mergeon']
    select_cols = params['mergeon'] + params['valcolumn']
    for valfile, group, label in zip(params['valfile'], params['valgroup'], params['vallabel']):
        suffix_val = '_' + label if label else ''
        try:
            data_val = load_data_group(valfile, group, chrom)
        except KeyError as ke:
            if params['wgdataset']:
                # data come in form of whole-genome dataset, i.e., not split by chromosome
                assert len(merge_cols) == 1, 'Only single column to merge on supported: {}'.format(merge_cols)
                merge_name = merge_cols[0]
                data_val = load_data_group(valfile, group)
                subset_idx = data_val[merge_name].isin(dataset[merge_name])
                subset = data_val.loc[subset_idx, :]
                data_val = subset.copy()
            else:
                raise ke
        for colname in merge_cols:
            shape_a, shape_b = dataset[colname].shape, data_val[colname].shape
            try:
                assert shape_a == shape_b,\
                    'Different number of data entries ({} vs {}) for column: "{}"'.format(shape_a, shape_b, colname)
                assert dataset[colname].isin(data_val[colname]).all(),\
                    'Some entries are not shared between dataset' \
                    ' and value file {} for column {}'.format(os.path.basename(valfile), colname)
            except AssertionError as ae:
                if params['ignoremissing']:
                    pass
                else:
                    raise ae
        data_val = data_val[select_cols]
        dataset = dataset.merge(data_val, how='outer', on=merge_cols, suffixes=('', suffix_val))
        all_rows = dataset.shape[0]
        if params['ignoremissing']:
            dataset = dataset.dropna(axis=0, how='any', inplace=False)
            if dataset.isnull().any(axis=0).any():
                assert dataset.shape[0] < all_rows, 'No rows dropped for chrom: {}'.format(chrom)
            assert not dataset.empty, 'Dataset empty after dropping N/A rows'
    return dataset


def merge_extend_datasets(params):
    """
    :param params:
    :return:
    """
    chrom = params['chrom']
    if len(params['inputfile']) > 1:
        chrom_data = merge_datasets(params)
    else:
        fp = params['inputfile'][0]
        chrom_data = load_data_group(fp, params['inputgroup'][0], chrom)
    if len(params['valfile']) > 0:
        chrom_data = extend_datasets(chrom_data, params)
    assert not chrom_data.empty, 'Merged and extended dataset is empty'
    return chrom, chrom_data


def stack_datasets(params):
    """
    :param params:
    :return:
    """
    dsets = []
    chrom = params['chrom']
    for fp, grp, lab in zip(params['inputfile'],
                            params['inputgroup'],
                            params['inputlabel']):
        data = load_data_group(fp, grp, chrom)
        if params['indicator']:
            data.insert(0, 'indicator', pd.Series([lab] * data.shape[0], dtype='object'))
        dsets.append(data)
    stacked = pd.concat(dsets, axis=0, ignore_index=True)
    stacked = stacked.dropna(axis=1, how='any')
    subset = [c for c in stacked.columns if c != 'indicator']
    stacked = stacked.drop_duplicates(subset=subset, keep='first', inplace=False)
    stacked = stacked.reset_index(drop=True)
    return chrom, stacked


def run_merge_datasets(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    if len(args.inputfiles) > 1:
        assert len(args.inputfiles) < 27, 'Merging of more than 26 datasets not supported at the moment'
        logger.debug('Preparing merge of datasets...')
    if len(args.valfile) > 0:
        logger.debug('Extending the final dataset with infos'
                     ' from {} additional file(s)'.format(len(args.valfile)))
    if len(args.inputfiles) < 2 and len(args.valfile) == 0:
        logger.warning('No merging possible and no additional datasets specified.'
                       ' What do you want me to do, human?')
    else:
        arglist = assemble_worker_params(args)
        logger.debug('Assembled a parameter list of size {} to process'.format(len(arglist)))
        _ = create_filepath(args.outputfile, logger)
        with pd.HDFStore(args.outputfile, args.filemode, complevel=9, complib='blosc') as hdf:
            if 'metadata' in hdf:
                metadata = hdf['metadata']
            else:
                metadata = pd.DataFrame(columns=MD_REGION_COLDEFS)
            with mp.Pool(args.workers) as pool:
                logger.debug('Initialized {} worker process(es)'.format(args.workers))
                if args.stack:
                    resit = pool.imap_unordered(stack_datasets, arglist, chunksize=1)
                else:
                    resit = pool.imap_unordered(merge_extend_datasets, arglist, chunksize=1)
                for chrom, dataobj in resit:
                    logger.debug('Received data for chromosome {}'.format(chrom))
                    grp, dataobj, metadata = region_generate(metadata, args.outputgroup, chrom, [args.inputfiles, args.valfile], dataobj)
                    hdf.put(grp, dataobj, format='fixed')
                    hdf.flush()
                    logger.debug('Flushed data to file')
            hdf.put('metadata', metadata, format='table')
            hdf.flush()
    logger.debug('Merging complete')
    return 0
