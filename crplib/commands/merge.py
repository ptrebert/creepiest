# coding=utf-8

"""
Small command module to merge existing datasets or extend them with additional data
"""

import os as os
import pandas as pd
import multiprocessing as mp

from crplib.auxiliary.hdf_ops import get_default_group, get_chrom_list, load_data_group
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
    inputlabels, inputgroups = dict(), dict()
    # Note for zip:
    # "The iterator stops when the shortest input iterable is exhausted"
    # no problem when there is only one file
    for deflab, ipf in zip(['A', 'B'], args.inputfile):
        if os.path.isfile(ipf):
            label, group, fp = deflab, get_default_group(ipf), ipf
        else:
            assert ipf.count(':') == 2, 'Please specify the input as LABEL:GROUP:FILEPATH' \
                                        ' (LABEL or GROUP may be empty, but exactly two ":" are needed)' \
                                        ' You specified as: {}'.format(ipf)
            label, group, fp = ipf.split(':')
            assert os.path.isfile(fp), 'Invalid path to file: {}'.format(fp)
            if not label:
                label = deflab
            if not group or group in ['auto', 'default']:
                group = get_default_group(fp)
        commons['inputfile'].append(fp)
        inputlabels[fp] = label
        inputgroups[fp] = group
    commons['inputlabel'] = inputlabels
    commons['inputgroup'] = inputgroups
    req_label = len(args.valfile) > 1
    commons['valfile'] = []
    commons['vallabel'] = dict()
    commons['valgroup'] = dict()
    for vfp in args.valfile:
        if req_label or not os.path.isfile(vfp):
            assert vfp.count(':') == 2, 'Please specify the value inputs as LABEL:GROUP:FILEPATH' \
                                        ' (GROUP may be empty, but exactly two ":" are needed)' \
                                        ' You specified as: {}'.format(vfp)
        if os.path.isfile(vfp):
            label, group, fp = '', get_default_group(vfp), vfp
        else:
            label, group, fp = vfp.split(':')
            assert os.path.isfile(fp), 'Invalid path to file: {}'.format(fp)
            if not group or group in ['auto', 'default']:
                group = get_default_group(fp)
        commons['valfile'].append(fp)
        commons['vallabel'][fp] = label
        commons['valgroup'][fp] = group
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
    fp_a = params['inputfile'][0]
    data_a = load_data_group(fp_a, params['inputgroup'][fp_a], chrom)
    fp_b = params['inputfile'][1]
    data_b = load_data_group(fp_b, params['inputgroup'][fp_b], chrom)
    merge_cols = params['mergeon']
    for colname in merge_cols:
        shape_a, shape_b = data_a[colname].shape, data_b[colname].shape
        assert shape_a == shape_b,\
            'Different number of data entries ({} vs {}) for column {}'.format(shape_a, shape_b, colname)
        assert data_a[colname].isin(data_b[colname]).all(),\
            'Some entries are not shared between datasets' \
            ' ({} and {}) for column {}'.format(os.path.basename(fp_a), os.path.basename(fp_b), colname)
    label_a = params['inputlabel'][fp_a]
    suffix_a = '_' + label_a
    label_b = params['inputlabel'][fp_b]
    suffix_b = '_' + label_b
    merged_data = data_a.merge(data_b, how='outer', on=merge_cols, suffixes=(suffix_a, suffix_b))
    merged_data = merged_data.assign(start=lambda x: x[['start' + suffix_a, 'start' + suffix_b]].min(axis=1))
    merged_data = merged_data.assign(end=lambda x: x[['end' + suffix_a, 'end' + suffix_b]].max(axis=1))
    return merged_data


def extend_datasets(dataset, params):
    """
    :param dataset:
    :param params:
    :return:
    """
    chrom = params['chrom']
    merge_cols = params['mergeon']
    select_cols = params['mergeon'] + params['valcolumns']
    for valfile in params['valfile']:
        label_val = params['vallabel'][valfile]
        suffix_val = '_' + label_val if label_val else ''
        group_val = params['valgroup'][valfile]
        data_val = load_data_group(valfile, group_val, chrom)
        for colname in merge_cols:
            shape_a, shape_b = dataset[colname].shape, data_val[colname].shape
            assert shape_a == shape_b,\
                'Different number of data entries ({} vs {}) for column {}'.format(shape_a, shape_b, colname)
            assert dataset[colname].isin(data_val[colname]).all(),\
                'Some entries are not shared between dataset' \
                ' and value file {} for column {}'.format(os.path.basename(valfile), colname)
        data_val = data_val[select_cols]
        dataset = dataset.merge(data_val, how='outer', on=merge_cols, suffixes=('', suffix_val))
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
        chrom_data = load_data_group(fp, params['inputgroup'][fp], chrom)
    if len(params['valfile']) > 0:
        chrom_data = extend_datasets(chrom_data, params)
    assert not chrom_data.empty, 'Merged and extended dataset is empty'

    return chrom, chrom_data


def run_merge_datasets(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    if len(args.inputfile) > 1:
        assert len(args.inputfile) == 2,\
            'Merging more than two datasets in a single run not supported: {} != 2'.format(len(args.inputfile))
        logger.debug('Preparing merge of datasets...')
    if len(args.valfile) > 0:
        logger.debug('Extending the final dataset with infos'
                     ' from {} additional file(s)'.format(len(args.valfile)))
    if len(args.inputfile) < 2 and len(args.valfiles) == 0:
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
                resit = pool.imap_unordered(merge_extend_datasets, arglist, chunksize=1)
                for chrom, dataobj in resit:
                    logger.debug('Received data for chromosome {}'.format(chrom))
                    grp, dataobj, metadata = region_generate(metadata, args.outputgroup, chrom, [args.inputfile, args.valfile], dataobj)
                    hdf.put(grp, dataobj, format='fixed')
                    hdf.flush()
                    logger.debug('Flushed data to file')
            hdf.put('metadata', metadata, format='table')
            hdf.flush()
    return 0
