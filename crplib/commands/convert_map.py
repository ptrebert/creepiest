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

from crplib.auxiliary.file_ops import text_file_mode
from crplib.auxiliary.text_parsers import get_map_iterator, read_chromosome_sizes
from crplib.metadata.md_helpers import normalize_group_path

from crplib.auxiliary.constants import TRGIDX_MASK, TRGIDX_SPLITS, TRGIDX_SELECT, TRGIDX_RELPOS


def assemble_worker_args(args, csizes):
    """
    :param args:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfile[0]
    commons['matchchroms'] = args.keepchroms
    commons['query'] = args.query
    for chrom, size in csizes.items():
        tmp = dict(commons)
        tmp['chrom'] = chrom
        tmp['size'] = size
        arglist.append(tmp)
    return arglist


def build_index_structures(mapit, csize, swap):
    """
    :param mapit:
    :param csize:
    :param swap:
    :return:
    """
    # need to buffer all blocks since the input map file
    # is expected to be sorted by chain/block, not by
    # genomic coordinate. It follows that the map number
    # is only meaningful in that specific condition;
    # need to sort everything by target coordinate
    # before saving the positions
    block_buffer = []
    # store minimal information necessary
    if swap:
        get_items = op.itemgetter(*(6, 7, 9))
    else:
        get_items = op.itemgetter(*(1, 2, 9))
    for block in mapit:
        block_buffer.append(get_items(block))
    # Default: all positions masked (1), all positions not conserved
    # Set position to 0 (unmask) if the position is conserved
    consmask = np.ones(csize, dtype=np.bool)
    block_buffer = sorted(block_buffer)
    for s, e, n in block_buffer:
        consmask[s:e] = 0
    # this creates a boolean select pattern to select
    # all conserved regions after splitting a signal track
    # using the split indices based on the alignment blocks
    pos = np.arange(consmask.size)
    pos = np.ma.array(pos, mask=consmask)
    splits = np.array(list(itt.chain.from_iterable((item.start, item.stop) for item in np.ma.clump_unmasked(pos))), dtype=np.int32)
    select = np.array([k for k, g in itt.groupby(~consmask)], dtype=np.bool)
    relpos = [n for s, e, n in block_buffer]
    return consmask, splits, select, relpos


def process_maps(params):
    """
    :param params:
    :return:
    """
    fpath = params['inputfile']
    select_chrom = params['chrom']
    selector = re.compile(select_chrom + '$')
    csize = params['size']
    match_chroms = re.compile(params['matchchroms'])
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        if params['query']:
            mapit = get_map_iterator(infile, tselect=match_chroms, qselect=selector)
            mask, splits, select, relpos = build_index_structures(mapit, csize, True)
        else:
            mapit = get_map_iterator(infile, tselect=selector, qselect=match_chroms)
            mask, splits, select, relpos = build_index_structures(mapit, csize, False)
    return select_chrom, mask, splits, select, relpos


def run_map_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    csizes = read_chromosome_sizes(args.chromsizes, args.keepchroms)
    arglist = assemble_worker_args(args, csizes)
    logger.debug('Start processing chain file {}'.format(args.inputfile[0]))
    og_mask = os.path.join(args.outputgroup, TRGIDX_MASK)
    og_splits = os.path.join(args.outputgroup, TRGIDX_SPLITS)
    og_select = os.path.join(args.outputgroup, TRGIDX_SELECT)
    og_relpos = os.path.join(args.outputgroup, TRGIDX_RELPOS)
    with pd.HDFStore(args.outputfile, args.filemode, complevel=9, complib='blosc', encoding='utf-8') as hdfout:
        with mp.Pool(args.workers) as pool:
            logger.debug('Iterating results')
            resit = pool.imap_unordered(process_maps, arglist, chunksize=1)
            for chrom, mask, splits, select, relpos in resit:
                logger.debug('Processed chromosome {}'.format(chrom))
                # collect all chromosomes in dataset(s)
                grp_mask = normalize_group_path(og_mask, suffix=chrom)
                hdfout.put(grp_mask, pd.Series(mask, dtype=np.bool), format='fixed')
                grp_splits = normalize_group_path(og_splits, suffix=chrom)
                hdfout.put(grp_splits, pd.Series(splits, dtype=np.int64), format='fixed')
                grp_select = normalize_group_path(og_select, suffix=chrom)
                hdfout.put(grp_select, pd.Series(select, dtype=np.bool), format='fixed')
                grp_relpos = normalize_group_path(og_relpos, suffix=chrom)
                hdfout.put(grp_relpos, pd.Series(relpos, dtype=np.int32), format='fixed')
                hdfout.flush()
                logger.debug('Saved chromosome {}'.format(chrom))
        hdfout.flush()
        # three groups per chromosome should be created
        assert len(hdfout.keys()) == len(arglist) * 4, 'Incomplete HDF file created: {}'.format(hdfout.keys())
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
