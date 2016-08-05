# coding=utf-8

"""
Module to handle conversion of region files into HDF5 format
"""

import os as os
import re as re
import pandas as pd
import multiprocessing as mp
import numpy as np
import itertools as itt

from crplib.auxiliary.file_ops import text_file_mode
from crplib.auxiliary.text_parsers import get_chain_iterator, read_chromosome_sizes
from crplib.metadata.md_helpers import normalize_group_path

from crplib.auxiliary.constants import TRGIDX_MASK, TRGIDX_SPLITS, TRGIDX_SELECT


def assemble_worker_args(args, csizes):
    """
    :param args:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfile[0]
    commons['qcheck'] = args.qcheck
    for chrom, size in csizes.items():
        tmp = dict(commons)
        tmp['chrom'] = chrom
        tmp['size'] = size
        arglist.append(tmp)
    return arglist


def build_index_structures(chainit, csize):
    """
    :param chainit:
    :param csize:
    :return:
    """
    # Default: all positions masked (1), all positions not conserved
    # Set position to 0 (unmask) if the position is conserved
    consmask = np.ones(csize, dtype=np.bool)
    splits = []
    last_end = 0
    for block in chainit:
        consmask[block[1]:block[2]] = 0
        # chainSort sorts according to number, not genomic coordinate
        # but chains are disjunct, so should not be a problem not to check
        if block[1] == 0:
            last_end = block[2]
        else:
            if last_end == block[1]:
                # we have a continuous block
                last_end = block[2]
            else:
                if last_end == 0:
                    splits.append(block[1])
                else:
                    splits.extend([last_end, block[1]])
                last_end = block[2]
    splits.append(last_end)
    # this creates a boolean select pattern to select
    # all conserved regions after splitting a signal track
    # using the split indices based on the alignment blocks
    select = np.array([k for k, g in itt.groupby(~consmask)], dtype=np.bool)
    splits = np.array(sorted(splits), dtype=np.int64)
    return consmask, splits, select


def process_chains(params):
    """
    :param params:
    :return:
    """
    fpath = params['inputfile']
    chrom = params['chrom']
    csize = params['size']
    qchroms = re.compile(params['qcheck'])
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        chainit = get_chain_iterator(infile, tselect=chrom, qcheck=qchroms)
        mask, splits, select = build_index_structures(chainit, csize)
    return chrom, mask, splits, select


def run_chain_conversion(args, logger):
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
    with pd.HDFStore(args.outputfile, 'w', complevel=9, complib='blosc', encoding='utf-8') as hdfout:
        with mp.Pool(args.workers) as pool:
            logger.debug('Iterating results')
            resit = pool.imap_unordered(process_chains, arglist, chunksize=1)
            for chrom, mask, splits, select in resit:
                logger.debug('Processed chromosome {}'.format(chrom))
                # collect all chromosomes in dataset(s)
                grp_mask = normalize_group_path(og_mask, suffix=chrom)
                hdfout.put(grp_mask, pd.Series(mask, dtype=np.bool), format='fixed')
                grp_splits = normalize_group_path(og_splits, suffix=chrom)
                hdfout.put(grp_splits, pd.Series(splits, dtype=np.int64), format='fixed')
                grp_select = normalize_group_path(og_select, suffix=chrom)
                hdfout.put(grp_select, pd.Series(select, dtype=np.bool), format='fixed')
                hdfout.flush()
                logger.debug('Saved chromosome {}'.format(chrom))
        hdfout.flush()
        # three groups per chromosome should be created
        assert len(hdfout.keys()) == len(arglist) * 3, 'Incomplete HDF file created: {}'.format(hdfout.keys())
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
