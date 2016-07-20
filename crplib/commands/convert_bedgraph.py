# coding=utf-8

"""
Module to handle conversion of bedGraph signal tracks into HDF5 format
"""

import pandas as pd
import multiprocessing as mp
import numpy as np
import psutil as psu
import scipy.stats as stats
import itertools as itt

from crplib.auxiliary.text_parsers import read_chromosome_sizes
from crplib.auxiliary.file_ops import text_file_mode, create_filepath
from crplib.metadata.md_signal import gen_obj_and_md, MD_SIGNAL_COLDEFS
from crplib.numalg.normalization import merge_1d_datasets

from crplib.auxiliary.constants import DIV_B_TO_GB


def assemble_worker_args(chroms, args):
    """
    :param chroms:
    :param args:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfile
    commons['mergestat'] = args.mergestat
    commons['noqnorm'] = args.noqnorm
    commons['clip'] = args.clip
    for name, size in chroms.items():
        tmp = dict(commons)
        tmp['chrom'] = name
        tmp['size'] = size
        arglist.append(tmp)
    return arglist


def process_signal(params):
    """
    :param params:
    :return:
    """
    all_data = tuple()
    chrom = params['chrom']
    for fp in params['inputfile']:
        opn, mode = text_file_mode(fp)
        values = np.zeros(params['size'], dtype='float64')
        with opn(fp, mode=mode, encoding='ascii') as infile:
            it = itt.dropwhile(lambda x: x.split()[0] != chrom, infile)
            for line in it:
                c, s, e, v = line.split()
                if c != chrom:
                    break
                values[int(s):int(e)] = float(v)
        if params['clip'] < 100.:
            new_max = stats.scoreatpercentile(values, params['clip'])
            values = np.clip(values, 0., new_max)
        all_data += values,
    if len(all_data) > 1 and not params['noqnorm']:
        retvals = merge_1d_datasets(*all_data, mergestat=params['mergestat'], qnorm=True)
    elif len(all_data) > 1 and params['noqnorm']:  # being explicit...
        retvals = merge_1d_datasets(*all_data, mergestat=params['mergestat'], qnorm=False)
    else:
        retvals = all_data[0]
    return chrom, retvals


def run_bedgraph_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    csizes = read_chromosome_sizes(args.chromsizes, args.keepchroms)
    logger.debug('Processing {} chromosome(s)'.format(len(csizes)))
    arglist = assemble_worker_args(csizes, args)
    meminfo = round(psu.virtual_memory().available / DIV_B_TO_GB, 2)
    logger.debug('Start processing, available memory: {}GB'.format(meminfo))
    create_filepath(args.outputfile, logger)
    with pd.HDFStore(args.outputfile, args.filemode, complevel=9, complib='blosc') as hdfout:
        with mp.Pool(args.workers) as pool:
            if 'metadata' in hdfout:
                metadata = hdfout['metadata']
            else:
                metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
            resit = pool.imap_unordered(process_signal, arglist, chunksize=1)
            logger.debug('Start processing chromosomes...')
            for chrom, valobj in resit:
                logger.debug('Chromosome {} completed'.format(chrom))
                grp, valobj, metadata = gen_obj_and_md(metadata, args.outputgroup, chrom, args.inputfile, valobj)
                hdfout.put(grp, valobj, format='fixed')
                hdfout.flush()
                meminfo = round(psu.virtual_memory().available / DIV_B_TO_GB, 2)
                logger.debug('Processed chromosome {} - available memory: {}'.format(chrom, meminfo))
        hdfout.put('metadata', metadata, format='table')
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    meminfo = round(psu.virtual_memory().available / DIV_B_TO_GB, 2)
    logger.debug('Available memory: {}'.format(meminfo))
    return 0
