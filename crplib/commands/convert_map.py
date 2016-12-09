# coding=utf-8

"""
Module to handle conversion of map files into HDF5 format (target index)
"""

import os as os
import re as re
import pandas as pd
import multiprocessing as mp
import operator as op
import numpy as np
import itertools as itt

from crplib.auxiliary.file_ops import text_file_mode, get_checksum
from crplib.auxiliary.text_parsers import iter_reduced_blocks, read_chromosome_sizes
from crplib.metadata.md_helpers import normalize_group_path

from crplib.auxiliary.constants import MAPIDX_MASK, MAPIDX_SELECT, MAPIDX_SPLITS, MAPIDX_ORDER


def assemble_worker_args(args, target, trgchroms, query, qrychroms):
    """
    :param args:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfile[0]

    # prepare arguments for target-centric processing
    # match all possible query chromosomes for one
    # specific target chromosome
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
    assert arglist, 'No arguments created for workers'
    return arglist


def build_index_structures(mapit, csize):
    """
    :param mapit:
    :param csize:
    :return:
    """
    # need to buffer all blocks since the input map file
    # is expected to be sorted by chain/block, not by
    # genomic coordinate. It follows that the map number
    # is only meaningful in that specific condition;
    # need to sort everything by target coordinate
    # before saving the positions

    # Default: all positions masked (1), all positions not conserved
    # Set position to 0 (unmask) if the position is conserved
    consmask = np.ones(csize, dtype=np.bool)
    block_buffer = []
    for block in mapit:
        consmask[block.start:block.end] = 0
        block_buffer.append(block)
    block_buffer = sorted(block_buffer, key=lambda x: (x.start, x.end))
    # this creates a boolean select pattern to select
    # all conserved regions after splitting a signal track
    # using the split indices based on the alignment blocks
    pos = np.arange(consmask.size)
    pos = np.ma.array(pos, mask=consmask)
    splits = np.array(list(itt.chain.from_iterable((item.start, item.stop) for item in np.ma.clump_unmasked(pos))), dtype=np.int32)
    select = np.array([k for k, g in itt.groupby(~consmask)], dtype=np.bool)
    order = [block.name for block in block_buffer]
    return consmask, splits, select, order


def process_maps(params):
    """
    :param params:
    :return:
    """
    fpath = params['inputfile']
    ref_chrom = params['chrom']
    ref_chrom_select = re.compile(ref_chrom + '$')
    csize = params['size']
    match_chroms = re.compile(params['match'])
    opn, mode = text_file_mode(fpath)
    mask, splits, select, order = None, None, None, None
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        if params['role'] == 'target':
            mapit = iter_reduced_blocks(infile, tselect=ref_chrom_select, qselect=match_chroms, which='target')
        elif params['role'] == 'query':
            mapit = iter_reduced_blocks(infile, tselect=match_chroms, qselect=ref_chrom_select, which='query')
        else:
            raise ValueError('Unknown role specified: {}'.format(params['role']))
        mask, splits, select, order = build_index_structures(mapit, csize)
    return (params['role'], params['assembly'], ref_chrom), mask, splits, select, order


def determine_assembly_names(tname, qname, tfile, qfile):
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


def run_map_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    logger.debug('Determining checksum for map file: {}'.format(args.inputfile[0]))
    hashtype, checksum = get_checksum(args.inputfile[0])
    logger.debug('{} checksum: {}'.format(hashtype, checksum))
    target, query = determine_assembly_names(args.target, args.targetchrom,
                                             args.query, args.querychrom)
    logger.debug('Assembly names identified: from TARGET {} >-> to QUERY {}'.format(target, query))
    trgchroms = read_chromosome_sizes(args.targetchrom, args.selectchroms)
    qrychroms = read_chromosome_sizes(args.querychrom, args.selectchroms)
    logger.debug('Files with chromosome sizes read, assembling worker parameters')
    arglist = assemble_worker_args(args, target, trgchroms, query, qrychroms)
    logger.debug('Argument list of size {} created'
                 ' - start processing map file {}'.format(len(arglist), args.inputfile[0]))
    with pd.HDFStore(args.outputfile, args.filemode, complevel=9, complib='blosc', encoding='utf-8') as hdfout:
        with mp.Pool(args.workers) as pool:
            logger.debug('Iterating results')
            resit = pool.imap_unordered(process_maps, arglist)
            for infos, mask, splits, select, relpos in resit:
                role, assembly, chrom = infos
                logger.debug('Finished {} for {} {}'.format(chrom, role, assembly))
                grp_mask = normalize_group_path(og_mask, suffix=chrom)
                hdfout.put(grp_mask, pd.Series(mask, dtype=np.bool), format='fixed')
                grp_splits = normalize_group_path(og_splits, suffix=chrom)
                hdfout.put(grp_splits, pd.Series(splits, dtype=np.int32), format='fixed')
                grp_select = normalize_group_path(og_select, suffix=chrom)
                hdfout.put(grp_select, pd.Series(select, dtype=np.bool), format='fixed')
                grp_order = normalize_group_path(og_order, suffix=chrom)
                # the order is the block name,
                hdfout.put(grp_order, pd.Series(relpos, dtype=np.int32), format='fixed')
                hdfout.flush()
                logger.debug('Data saved')
        hdfout.flush()
        logger.debug('Full index created, storing metadata')
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
