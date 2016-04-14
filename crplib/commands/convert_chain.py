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


def process_chains(params):
    """
    :param params:
    :return:
    """
    mypid = mp.current_process().pid
    fpath = params['inputfile']
    chrom = params['chrom']
    csize = params['size']
    qchroms = re.compile(params['qcheck'])
    consmask = np.ones(csize, dtype=np.bool)
    splits = []
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        chainit = get_chain_iterator(infile, tselect=chrom, qcheck=qchroms)
        for block in chainit:
            consmask[block[1]:block[2]] = 0
            if block[1] == 0:
                # chainSort sorts according to number, not genomic coordinate
                # but chains are disjunct, so should not be a problem not to check
                #assert not splits, \
                #    'Split indices not empty for 0 start in chain file {} - chrom {} and block {}'.format(fpath, chrom, block)
                splits.append(block[2])
            else:
                splits.extend([block[1], block[2]])
    # this creates a boolean select pattern to select
    # all conserved regions after splitting a signal track
    # using the split indices based on the alignment blocks
    select = np.array([k for k, g in itt.groupby(~consmask)], dtype=np.bool)
    splits = np.array(sorted(splits), dtype=np.int32)
    return mypid, chrom, consmask, splits, select


def run_chain_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    csizes = read_chromosome_sizes(args.chromsizes, args.keepchroms)
    arglist = assemble_worker_args(args, csizes)
    logger.debug('Start processing chain file {}'.format(args.inputfile[0]))
    og_mask = os.path.join(args.outputgroup, 'cons', 'mask')
    og_splits = os.path.join(args.outputgroup, 'cons', 'splits')
    og_select = os.path.join(args.outputgroup, 'cons', 'select')
    with pd.HDFStore(args.outputfile, 'a', complevel=9, complib='blosc') as hdfout:
        with mp.Pool(args.workers) as pool:
            logger.debug('Iterating results')
            mapres = pool.map_async(process_chains, arglist, chunksize=1)
            for pid, chrom, mask, splits, select in mapres.get():
                logger.debug('Worker (PID {}) completed chromosome {}'.format(pid, chrom))
                # collect all chromosomes in dataset(s)
                grp_mask = normalize_group_path(og_mask, suffix=chrom)
                hdfout.put(grp_mask, pd.Series(mask, dtype=np.bool), format='fixed')
                grp_splits = normalize_group_path(og_splits, suffix=chrom)
                hdfout.put(grp_splits, pd.Series(splits, dtype=np.int32), format='fixed')
                grp_select = normalize_group_path(og_select, suffix=chrom)
                hdfout.put(grp_select, pd.Series(select, dtype=np.bool), format='fixed')
                hdfout.flush()
                logger.debug('Saved chromosome {}'.format(chrom))
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
