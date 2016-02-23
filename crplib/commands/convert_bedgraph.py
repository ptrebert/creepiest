# coding=utf-8

"""
Module to handle conversion of bedGraph signal tracks into HDF5 format
"""

import pandas as pd
import multiprocessing as mp
import numpy as np
import psutil as psu

from crplib.auxiliary.text_parsers import read_chromosome_sizes
from crplib.auxiliary.file_ops import text_file_mode
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
    for name, size in chroms.items():
        tmp = dict(commons)
        tmp['chrom'] = name
        tmp['size'] = size
        arglist.append(tmp)
    return arglist


def define_chromosome_indices(csizes):
    """
    :param csizes:
    :return:
    """
    curr_start = 0
    indices = dict()
    for name, size in csizes.items():
        indices[name] = curr_start, curr_start + size
        curr_start += size
    return curr_start, indices


def _parse_bg_line(line):
    """
    :param line:
    :return:
    """
    _, start, end, val = line.strip().split()
    return int(start), int(end), float(val)


def process_signal(**kwargs):
    """
    :param kwargs:
    :return:
    """
    all_data = tuple()
    mypid = mp.current_process().pid
    chrom = kwargs['chrom']
    for fp in kwargs['inputfile']:
        opn, mode = text_file_mode(fp)
        rdl = _parse_bg_line
        values = np.zeros(kwargs['size'], dtype='float64')
        with opn(fp, mode=mode, encoding='ascii') as infile:
            for line in infile:
                if line.startswith(chrom):
                    start, end, val = rdl(line)
                    values[start:end] = val
        all_data += values,
    if len(all_data) > 1 and not kwargs['noqnorm']:
        retvals = merge_1d_datasets(all_data, mergestat=kwargs['mergestat'], qnorm=True)
    elif len(all_data) > 1 and kwargs['noqnorm']:  # being explicit...
        retvals = merge_1d_datasets(all_data, mergestat=kwargs['mergestat'], qnorm=False)
    else:
        retvals = all_data[0]
    return mypid, chrom, retvals


def run_bedgraph_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    csizes = read_chromosome_sizes(args.chromsizes, args.keepchroms)
    logger.debug('Processing {} chromosomes'.format(len(csizes)))
    arglist = assemble_worker_args(csizes, args)
    wgsize, chromidx = define_chromosome_indices(csizes)
    # allocate enough memory for whole genome signal
    # to compute the genome-wide metadata at the end
    wgsignal = np.zeros(wgsize, dtype='float64')
    meminfo = round(psu.virtual_memory().available / DIV_B_TO_GB, 2)
    logger.debug('Whole-genome array allocated; available memory: {}GB'.format(meminfo))
    with pd.HDFStore(args.outputfile, 'a', complevel=9, complib='blosc') as hdfout:
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(process_signal, arglist, chunksize=1)
            metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
            for pid, chrom, valarray in resit.next():
                logger.debug('Worker (PID {}) completed chromosome {}'.format(pid, chrom))
                grp, obj, metadata = gen_obj_and_md(metadata, args.grouproot, chrom, args.inputfile, valarray)
                hdfout.put(grp, obj, format='fixed')
                hdfout.flush()
                del obj
                start, end = chromidx[chrom]
                wgsignal[start:end] = valarray
                meminfo = round(psu.virtual_memory().available / DIV_B_TO_GB, 2)
                logger.debug('Processed chromosome {} - available memory: {}'.format(chrom, meminfo))
        _, _, metadata = gen_obj_and_md(metadata, args.grouproot, 'wg', args.inputfile, wgsignal)
        hdfout.put('metadata', metadata, format='table')
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
