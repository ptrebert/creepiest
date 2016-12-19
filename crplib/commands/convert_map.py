# coding=utf-8

"""
Module to handle conversion of map files into HDF5 format (target index)
"""

import os as os
import re as re
import pandas as pd
import numpy as np
import itertools as itt
import operator as op

from crplib.auxiliary.file_ops import get_checksum, shm_file_object, create_filepath
from crplib.auxiliary.text_parsers import iter_map_blocks, read_chromosome_sizes
from crplib.metadata.md_helpers import normalize_group_path

from crplib.auxiliary.constants import MAPIDX_MASK, MAPIDX_SELECT, \
    MAPIDX_SPLITS, MAPIDX_BLOCKS, DIV_B_TO_GB


def assemble_worker_args(args, target, trgchroms, query, qrychroms):
    """
    :param args:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfiles[0]

    # first two for loops: construct fast index structures
    # for correlation computations by saving the conserved
    # positions for each chromosome in form of a bool mask
    # and respective access patterns
    query_match = '(' + '|'.join([qc for qc in qrychroms.keys()]) + ')$'
    for chrom, size in trgchroms.items():
        tmp = dict(commons)
        tmp['chrom'] = chrom
        tmp['size'] = size
        tmp['role'] = 'target'
        tmp['match'] = query_match
        tmp['assembly'] = target
        arglist.append(tmp)

    target_match = '(' + '|'.join([tc for tc in trgchroms.keys()]) + ')$'
    for chrom, size in qrychroms.items():
        tmp = dict(commons)
        tmp['chrom'] = chrom
        tmp['size'] = size
        tmp['role'] = 'query'
        tmp['match'] = target_match
        tmp['assembly'] = query
        arglist.append(tmp)

    # this for loop: store the actual mapping in form of a dataframe
    # to avoid reading the text file/jumping around in the text file
    for tchrom, _ in trgchroms.items():
        for qchrom, _ in qrychroms.items():
            tmp = dict(commons)
            tmp['tchrom'] = tchrom
            tmp['qchrom'] = qchrom
            tmp['size'] = 0
            tmp['role'] = 'both'
            tmp['assembly'] = '{}_from_{}'.format(query, target)
            arglist.append(tmp)
    assert arglist, 'No arguments created for workers'
    return arglist


def build_fast_index_structures(mapit, csize):
    """
    :param mapit:
    :param csize:
    :return:
    """
    # Default: all positions masked (1), all positions not conserved
    # Set position to 0 (unmask) if the position is conserved
    consmask = np.ones(csize, dtype=np.bool)
    for block in mapit:
        consmask[block.start:block.end] = 0
    # this creates a boolean select pattern to select
    # all conserved regions after splitting a signal track
    # using the split indices based on the alignment blocks
    pos = np.arange(consmask.size)
    pos = np.ma.array(pos, mask=consmask)
    splits = np.array(list(itt.chain.from_iterable((item.start, item.stop) for item in np.ma.clump_unmasked(pos))), dtype=np.int32)
    select = np.array([k for k, g in itt.groupby(~consmask)], dtype=np.bool)
    assert consmask.sum() < csize, 'No regions masked'
    return consmask, splits, select


def build_map_index_structure(mapit):
    """
    :param mapit:
    :return:
    """
    index = []
    maps = []
    for block in mapit:
        index.append(block.name)
        tstrand = 1 if block.tstrand == '+' else 0
        assert tstrand == 1, 'Negative target strand detected - invalid line in map: {}'.format(block)
        qstrand = 1 if block.qstrand == '+' else 0
        maps.append((block.tstart, block.tend, tstrand, block.qstart, block.qend, qstrand))
    assert len(index) == len(maps), 'Size mismatch between index and maps: {} vs {}'.format(len(index), len(maps))
    assert maps or len(index) == 0, 'No mapping blocks extracted'
    return index, maps


def process_maps(params):
    """
    :param params:
    :return:
    """
    fpath = params['inputfile']
    opener, ftype = shm_file_object(fpath)
    global _shm_mapfile
    assert _shm_mapfile.tell() == 0, 'Shm map not at position 0: {}'.format(_shm_mapfile.tell())
    if ftype == 'raw':
        with _shm_mapfile as blockdata:
            infos, values = iterate_blockdata(params, blockdata)
    elif ftype == 'gzip':
        with opener(fileobj=_shm_mapfile, mode='r') as blockdata:
            infos, values = iterate_blockdata(params, blockdata)
    elif ftype == 'bzip':
        with opener(filename=_shm_mapfile, mode='r') as blockdata:
            infos, values = iterate_blockdata(params, blockdata)
    else:
        raise ValueError('Unexpected file type for shm object')
    _shm_mapfile.seek(0)
    return infos, values


def iterate_blockdata(params, blockdata):
    """
    :param params:
    :param blockdata:
    :return:
    """
    if params['role'] == 'both':
        # build full map index
        tchrom = params['tchrom']
        tchrom_select = re.compile(tchrom + '$')
        qchrom = params['qchrom']
        qchrom_select = re.compile(qchrom + '$')
        mapit = iter_map_blocks(blockdata, tselect=tchrom_select, qselect=qchrom_select, mode='both')
        index, maps = build_map_index_structure(mapit)
        infos, values = (params['role'], params['assembly'], tchrom, qchrom), (index, maps)
    else:
        # build fast index structures
        # for correlations
        ref_chrom = params['chrom']
        ref_chrom_select = re.compile(ref_chrom + '$')
        csize = params['size']
        match_chroms = re.compile(params['match'])
        if params['role'] == 'target':
            mapit = iter_map_blocks(blockdata, tselect=ref_chrom_select, qselect=match_chroms, mode='target')
            mask, splits, select = build_fast_index_structures(mapit, csize)
            infos, values = (params['role'], params['assembly'], ref_chrom), (mask, splits, select)
        elif params['role'] == 'query':
            mapit = iter_map_blocks(blockdata, tselect=match_chroms, qselect=ref_chrom_select, mode='query')
            mask, splits, select = build_fast_index_structures(mapit, csize)
            infos, values = (params['role'], params['assembly'], ref_chrom), (mask, splits, select)
        else:
            raise ValueError('Unexpected role in parameter set: {}'.format(params['role']))
    return infos, values

# ==================================
# everything above likely deprecated
# ==================================


def determine_assembly_names(tname, tfile, qname, qfile):
    """
    :param tname:
    :param qname:
    :param tfile:
    :param qfile:
    :return:
    """
    assm_re = re.compile('(?P<ASSM>[a-zA-Z0-9]+)')
    if tname:
        target = tname
    else:
        fn = os.path.basename(tfile)
        mobj = assm_re.match(fn)
        assert mobj is not None, 'Could not identify target assembly name in filename {}'.format(fn)
        target = mobj.group('ASSM')
    if qname:
        query = qname
    else:
        fn = os.path.basename(qfile)
        mobj = assm_re.match(fn)
        assert mobj is not None, 'Could not identify query assembly name in filename {}'.format(fn)
        query = mobj.group('ASSM')
    return target, query


def built_index_metadata(trg, qry, htype, check, trgsizes, qrysizes, minmap):
    """
    :param trg:
    :param qry:
    :param htype:
    :param check:
    :param trgsizes:
    :param qrysizes:
    :param minmap:
    :return:
    """
    rows = [('target', trg), ('query', qry), ('hashtype', htype), ('hash', check)]
    for chrom, size in trgsizes.items():
        rows.append((os.path.join(trg, chrom), str(size)))
    for chrom, size in qrysizes.items():
        rows.append((os.path.join(qry, chrom), str(size)))
    if minmap:
        rows.append(('maptype', 'minimal'))
    else:
        rows.append(('maptype', 'complete'))
    md = pd.DataFrame(rows, columns=['key', 'value'])
    return md


def filter_hdf_paths(chrompairs, select, which):
    """
    :param chrompairs:
    :param select:
    :param which:
    :return:
    """
    if which == 'target':
        load_paths = [t[2] for t in chrompairs if t[1] == select]
    elif which == 'query':
        load_paths = [t[2] for t in chrompairs if t[0] == select]
    else:
        raise ValueError('Reference has to be specified as target or query, not as {}'.format(which))
    return load_paths


def save_conservation_masks(hdfobj, chroms, chrompairs, minmap, which, logger):
    """
    :param hdfobj:
    :param chroms:
    :param chrompairs:
    :param minmap:
    :param which:
    :param logger:
    :return:
    """
    select_entries = {'target': op.itemgetter(*(0, 1)),
                      'query': op.itemgetter(*(2, 3))}
    get_coords = select_entries[which]
    logger.debug('Building conservation masks for {}'.format(which))
    for chrom, size in chroms.items():
        loadpaths = filter_hdf_paths(chrompairs, chrom, which)
        consmask = np.ones(size, dtype=np.bool)
        unmasked = 0
        for lp in loadpaths:
            mapdf = hdfobj[lp]
            indices = mapdf.apply(lambda x: slice(*get_coords(x)), axis=1, raw=True)
            for idx in indices:
                assert np.all(consmask[idx] == 1), '{} / {} overlap detected at: {}'.format(which, chrom, idx)
                consmask[idx] = 0
                unmasked += idx.stop - idx.start
        masked = consmask.sum()
        assert masked < size, 'No regions masked for {} chrom: {}'.format(which, chrom)
        # just paranoid, but if blocks would not be disjoint (= overlap)
        # this should raise
        delta = size - masked - unmasked
        assert delta == 0, '{} / {} mismatch - masked {} -' \
                           ' unmasked {} - size {} - delta {}'.format(which, chrom, masked, unmasked, size, delta)
        grp_mask = normalize_group_path(os.path.join(which, MAPIDX_MASK), suffix=chrom)
        hdfobj.put(grp_mask, pd.Series(consmask, dtype=np.bool), format='fixed')
        if not minmap:
            pos = np.arange(consmask.size)
            pos = np.ma.array(pos, mask=consmask)
            splits = np.array(list(itt.chain.from_iterable((item.start, item.stop) for item in np.ma.clump_unmasked(pos))), dtype=np.int32)
            select = np.array([k for k, g in itt.groupby(~consmask)], dtype=np.bool)
            grp_splits = normalize_group_path(os.path.join(which, MAPIDX_SPLITS), suffix=chrom)
            hdfobj.put(grp_splits, pd.Series(splits, dtype=np.int32), format='fixed')
            grp_select = normalize_group_path(os.path.join(which, MAPIDX_SELECT), suffix=chrom)
            hdfobj.put(grp_select, pd.Series(select, dtype=np.bool), format='fixed')
            hdfobj.flush()
        hdfobj.flush()
        logger.debug('Done for {} {}'.format(which, chrom))
    return


def save_splits(hdfobj, trgchroms, qrychroms, mapdf, logger):
    """
    :param hdfobj:
    :param trgchroms:
    :param qrychroms:
    :param mapdf:
    :param logger:
    :return:
    """
    logger.debug('Writing individual chromosome pairs...')
    subset_names = ['tstart', 'tend', 'qstart', 'qend', 'qstrand']
    subset_dtypes = {'tstart': np.int32, 'tend': np.int32, 'qstart': np.int32,
                     'qend': np.int32, 'qstrand': np.int8}
    chrom_pairs = []
    for qc, qsize in qrychroms.items():
        for tc, tsize in trgchroms.items():
            subset = mapdf.loc[(mapdf.tchrom == tc) & (mapdf.qchrom == qc), subset_names]
            if subset.empty:
                logger.debug('No combination target {} - query {}'.format(tc, qc))
                continue
            subset = subset.astype(subset_dtypes, copy=True)
            subset = subset.sort_index(axis='index')
            hdf_path = os.path.join('qt', MAPIDX_BLOCKS, qc, tc)
            grp_mapidx = normalize_group_path(hdf_path)
            hdfobj.put(grp_mapidx, subset, format='fixed')
            hdfobj.flush()
            chrom_pairs.append((qc, tc, grp_mapidx))
    logger.debug('All data saved')
    return chrom_pairs


def read_split_map(args, trgchroms, qrychroms, logger):
    """
    :param args:
    :return:
    """
    names = ['tchrom', 'tstart', 'tend', 'tstrand',
             'qchrom', 'qstart', 'qend', 'qstrand']
    names.insert(args.indexcol, 'index')
    datatypes = {'tchrom': str, 'tstart': np.int32, 'tend': np.int32, 'tstrand': str,
                 'qchrom': str, 'qstart': np.int32, 'qend': np.int32, 'qstrand': str,
                 'index': np.int32}
    logger.debug('Reading map file...')
    mapdf = pd.read_csv(args.inputfiles[0], sep='\t', names=names, index_col=args.indexcol,
                        low_memory=True, dtype=datatypes, compression='infer', encoding='utf-8')
    size_in_mem = mapdf.values.nbytes
    logger.debug('Reading done - full map size: ~{}GiB with {} rows'.format(round(size_in_mem / DIV_B_TO_GB, 2), mapdf.shape[0]))
    mapdf.replace({'qstrand': {'+': 1, '-': -1}}, inplace=True)
    tchroms = dict((k, trgchroms[k]) for k in mapdf.tchrom.unique())
    logger.debug('Identified {} chromosomes for target assembly in map file'.format(len(tchroms)))
    qchroms = dict((k, qrychroms[k]) for k in mapdf.qchrom.unique())
    logger.debug('Identified {} chromosomes for query assembly in map file'.format(len(qchroms)))
    _ = create_filepath(args.outputfile, logger)
    with pd.HDFStore(args.outputfile, args.filemode, complevel=9, complib='blosc', encoding='utf-8') as hdfout:
        chrompairs = save_splits(hdfout, tchroms, qchroms, mapdf, logger)
        save_conservation_masks(hdfout, tchroms, chrompairs, args.minmap, 'target', logger)
        save_conservation_masks(hdfout, qchroms, chrompairs, args.minmap, 'query', logger)
    return tchroms, qchroms, chrompairs


def run_map_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    logger.debug('Determining checksum for map file: {}'.format(args.inputfiles[0]))
    hashtype, checksum = get_checksum(args.inputfiles[0])
    logger.debug('{} checksum: {}'.format(hashtype, checksum))
    target, query = determine_assembly_names(args.target, args.targetchrom,
                                             args.query, args.querychrom)
    logger.debug('Assembly names identified: from TARGET {} >>> to QUERY {}'.format(target, query))
    trgchroms = read_chromosome_sizes(args.targetchrom, args.selectchroms)
    qrychroms = read_chromosome_sizes(args.querychrom, args.selectchroms)
    logger.debug('Files with chromosome sizes read')
    size_on_disk = os.path.getsize(args.inputfiles[0])
    if (size_on_disk / DIV_B_TO_GB) > 2:
        logger.warning('The map file seems to be very large - are you sure this is correct, human?')
    trgchroms, qrychroms, chrompairs = read_split_map(args, trgchroms, qrychroms, logger)
    with pd.HDFStore(args.outputfile, 'a', complevel=9, complib='blosc', encoding='utf-8') as hdfout:
        logger.debug('Creating metadata')
        mdf = built_index_metadata(target, query, hashtype, checksum, trgchroms, qrychroms, args.minmap)
        hdfout.put('/metadata', mdf, format='table')
        hdfout.flush()
        logger.debug('Full index created, storing metadata')
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
