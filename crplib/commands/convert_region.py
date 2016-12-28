# coding=utf-8

"""
Module to handle conversion of region files into HDF5 format
"""

import pandas as pd
import multiprocessing as mp
import collections as col
import numpy as np
import re as re
import scipy.stats as stats

from crplib.auxiliary.file_ops import text_file_mode, create_filepath
from crplib.metadata.md_regions import gen_obj_and_md, MD_REGION_COLDEFS
from crplib.auxiliary.text_parsers import determine_text_table_type


def assemble_worker_args(args, logger):
    """
    :param args:
    :return:
    """
    arglist = []
    assert args.nameidx not in {0, 1, 2},\
        'Specified name index field as {} - first three fields need to be "chrom - start - end"'.format(args.nameidx)
    assert args.scoreidx not in {0, 1, 2},\
        'Specified score index field as {} - first three fields need to be "chrom - start - end"'.format(args.scoreidx)
    for fp in args.inputfiles:
        skip, delim, colnames = determine_text_table_type(fp, args.useheader, logger)
        if colnames:
            colnames[0] = 'chrom'
            colnames[1] = 'start'
            colnames[2] = 'end'
            if args.nameidx != -1:
                colnames[args.nameidx] = 'name'
        else:
            colnames = ['chrom', 'start', 'end']
            # note if both are set to -1, then max(-1,-1) + 1 = 0
            # so range(3,0,1) results in empty list
            for idx in range(3, max(args.nameidx, args.scoreidx) + 1, 1):
                if idx == args.nameidx:
                    colnames.append('name')
                elif idx == args.scoreidx:
                    colnames.append('score')
                else:
                    colnames.append('col' + str(idx))
        colnames = validate_column_names(colnames, args.nameidx, logger)
        commons = dict()
        commons['skip'] = skip
        commons['delimiter'] = delim
        commons['colnames'] = colnames
        commons['inputfile'] = fp
        commons['selectchroms'] = args.selectchroms
        commons['keeptop'] = args.keeptop
        commons['scoreidx'] = args.scoreidx
        commons['useheader'] = args.useheader
        commons['filtersize'] = args.filtersize
        arglist.append(commons)
    return arglist


def validate_column_names(colnames, nameidx, logger):
    """
    :param colnames:
    :param nameidx:
    :return:
    """
    if any([cn.startswith('ft') for cn in colnames]):
        logger.error('Detected column names starting with "ft" - conversion of such fields'
                     ' not supported: {}'.format(' | '.join([cn for cn in colnames if cn.startswith('ft')])))
        raise AssertionError('Prefix "ft" in column names detected')
    if len(colnames) != len(set(colnames)):  # duplicates
        logger.debug('Detected duplicates in column names')
        fix_indices = {0, 1, 2}
        if nameidx != -1:
            fix_indices.add(nameidx)
        duplicates = set()
        seen = set()
        for cn in colnames:
            if cn in seen:
                duplicates.add(cn)
            else:
                seen.add(cn)
        for idx, cn in enumerate(colnames):
            # need to check for the fix indices here since
            # the user specified name column can be anywhere
            if cn in duplicates and idx not in fix_indices:
                logger.debug('Duplicate name {} is suffixed with _dup{}'.format(cn, idx))
                colnames[idx] = cn + '_dup' + str(idx)
    return colnames


def merge_overlapping_regions(allregions, allchroms):
    """
    :param allregions:
    :param allchroms:
    :return:
    """
    allregions.sort_values(['chrom', 'start', 'end'], axis='index', inplace=True)
    merged_regions = None
    for chrom in allchroms:
        chrom_regions = allregions[allregions.chrom == chrom]
        chrom_regions = col.deque(chrom_regions.to_dict('record'))
        chrom_merged = []
        this = chrom_regions.popleft()
        while 1:
            try:
                step = chrom_regions.popleft()
            except IndexError:  # deque empty
                break
            if this['end'] < step['start']:
                chrom_merged.append(this)
                this = step
                continue
            else:
                # since all regions are sorted, must overlap now,
                # or be at least book-ended
                this['start'] = min(this['start'], step['start'])
                this['end'] = max(this['end'], step['end'])
                if 'name' in this:
                    this['name'] = this['name'] + '-' + step['name']
        if merged_regions is None:
            merged_regions = pd.DataFrame.from_dict(chrom_merged, orient='columns')
        else:
            merged_regions = pd.concat([merged_regions, pd.DataFrame.from_dict(chrom_merged, orient='columns')],
                                       ignore_index=True, axis=0)
    merged_regions.sort_values(['chrom', 'start', 'end'], axis='index', inplace=True)
    merged_regions.index = np.arange(merged_regions.shape[0])
    assert not merged_regions.empty, 'Something went wrong - no regions left after merging'
    return merged_regions


def process_regions(params):
    """
    :param params:
    :return:
    """
    fpath = params['inputfile']
    chr_match = re.compile(params['selectchroms'])
    score_col_idx = params['scoreidx']
    if score_col_idx != -1:
        score_col_name = params['colnames'][score_col_idx]
        datatypes = {'start': np.int32, 'end': np.int32, score_col_name: np.float64}
    else:
        datatypes = {'start': np.int32, 'end': np.int32}

    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        regions = pd.read_csv(infile, sep=params['delimiter'], names=params['colnames'],
                              index_col=False, dtype=datatypes, header=0,
                              skipinitialspace=True, skiprows=params['skip'], skip_blank_lines=True,
                              encoding='utf-8', comment='#', usecols=params['colnames'])
    chroms_in_file = regions.chrom.drop_duplicates().tolist()
    remove_chroms = set(filter(lambda x: chr_match.match(x) is None, chroms_in_file))
    drop_columns = ['filter_for_chrom']
    regions = regions.assign(filter_for_chrom=lambda x: x.chrom.isin(remove_chroms))
    regions.drop(regions.index[regions.filter_for_chrom], inplace=True, axis='index')
    if params['filtersize'] > 0:
        drop_columns.append('filter_for_length')
        regions = regions.assign(filter_for_length=lambda x: x.end - x.start)
        regions.drop(regions[regions.filter_for_length < params['filtersize']].index, inplace=True, axis='index')
    if score_col_idx != -1 and params['keeptop'] < 100:
        drop_columns.append('filter_for_score')
        # heuristic to check if score column seems to be reasonable
        assert regions[score_col_name].var() > 0, \
            'Scores have 0 variance in file {} for selected column {}'.format(fpath, score_col_idx)
        lower_threshold = stats.scoreatpercentile(regions[score_col_name].values, 100 - params['keeptop'])
        regions = regions.assign(filter_for_score=lambda x: x[score_col_name] < lower_threshold)
        regions.drop(regions.index[regions.filter_for_score], inplace=True, axis='index')
    if not params['useheader']:
        for col in regions.columns:
            if col in ['chrom', 'start', 'end', 'name']:
                continue
            drop_columns.append(col)
    regions.drop(drop_columns, axis='columns', inplace=True)
    reordered_columns = reorder_columns(regions.columns.tolist())
    regions = regions[reordered_columns]
    regions.sort_values(['chrom', 'start', 'end'], axis='index', inplace=True)
    regions.index = np.arange(regions.shape[0])
    assert not regions.empty, 'No regions read from file {} (or are left after filtering)'.format(fpath)
    return regions, set(regions.chrom.drop_duplicates().tolist())


def reorder_columns(colnames):
    """
    :param colnames:
    :return:
    """
    reordered = []
    idx = 5
    for col in colnames:
        if col == 'chrom':
            reordered.append((0, col))
        elif col == 'start':
            reordered.append((1, col))
        elif col == 'end':
            reordered.append((2, col))
        elif col == 'name':
            reordered.append((3, col))
        elif col == 'score':
            reordered.append((4, col))
        else:
            reordered.append((idx, col))
            idx += 1
    return [t[1] for t in sorted(reordered)]


def run_region_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    if args.useheader:
        assert len(args.inputfiles) == 1, 'Too many input files. Cannot use header information when merging' \
                                         ' several input files (since merging only works on overlapping' \
                                         ' intervals defined by start and end coordinate).'
    arglist = assemble_worker_args(args, logger)
    logger.debug('Start processing {} region file(s)'.format(len(args.inputfiles)))
    _ = create_filepath(args.outputfile, logger)
    with pd.HDFStore(args.outputfile, args.filemode, complevel=9, complib='blosc') as hdfout:
        if 'metadata' in hdfout:
            metadata = hdfout['metadata']
        else:
            metadata = pd.DataFrame(columns=MD_REGION_COLDEFS)
        with mp.Pool(args.workers) as pool:
            all_chroms = set()
            all_regions = None
            logger.debug('Iterating results')
            resit = pool.imap_unordered(process_regions, arglist, chunksize=1)
            for regobj, chroms in resit:
                logger.debug('Received {} regions'.format(regobj.shape[0]))
                # collect all chromosomes in dataset(s)
                all_chroms |= chroms
                if all_regions is None:
                    all_regions = regobj
                else:
                    # note to self: concat does not accept aliases 'index' and 'columns' for parameter axis
                    all_regions = pd.concat([all_regions, regobj], axis=0, ignore_index=True, join='inner')
            logger.debug('All files processed...')
            if len(args.inputfiles) > 1:
                logger.debug('Merging {} files...'.format(len(args.inputfiles)))
                all_regions = merge_overlapping_regions(all_regions, all_chroms)
                logger.debug('Merging resulted in {} regions'.format(all_regions.shape[0]))
            # note here that if the file contains a "name" field in the header
            # the user does not need to specify name-idx
            if args.nameidx == -1 and 'name' not in all_regions.columns:
                all_regions = all_regions.assign(name=lambda x: ['region_' + str(idx) for idx in x.index])
            logger.debug('Identified {} chromosomes in dataset(s)'.format(len(all_chroms)))
            for chrom in sorted(all_chroms):
                chrom_regions = all_regions[all_regions.chrom == chrom]
                chrom_regions = chrom_regions.drop(['chrom'], axis='columns', inplace=False)
                if chrom_regions.empty:
                    continue
                grp, valobj, metadata = gen_obj_and_md(metadata, args.outputgroup, chrom, args.inputfiles, chrom_regions)
                hdfout.put(grp, valobj, format='fixed')
                hdfout.flush()
                logger.debug('Processed chromosome {}'.format(chrom))
        hdfout.put('metadata', metadata, format='table')
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
