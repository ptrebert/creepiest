# coding=utf-8

"""
Map a signal track from one species/assembly to another
"""

import os as os
import numpy as np
import pandas as pd
import re as re
import collections as col
import multiprocessing as mp
import psutil as psu

from crplib.metadata.md_signal import gen_obj_and_md, MD_SIGNAL_COLDEFS
from crplib.auxiliary.text_parsers import read_chromosome_sizes, get_chain_iterator
from crplib.auxiliary.file_ops import text_file_mode, create_filepath
from crplib.auxiliary.hdf_ops import get_chrom_list,\
    check_path_infos, get_default_group
from crplib.auxiliary.constants import DIV_B_TO_GB


_shm_carr = dict()


def assemble_worker_params(inputfile, inputgroup, numalloc, chromsizes, allblocks, carrays):
    """
    :param inputfile:
    :param inputgroup:
    :param numalloc:
    :param chromsizes:
    :param allblocks:
    :return:
    """
    this_run = list(chromsizes.keys())[:numalloc]
    if not inputgroup:
        inputgroup = get_default_group(inputfile)
    for c in this_run:
        # NB: lock=False is possible since all writes to the chromosome
        # arrays happen to disjoint slices
        carrays[c] = mp.Array('d', np.zeros(chromsizes[c], dtype=np.float64), lock=False)
    params = []
    for c in this_run:
        for blocks in allblocks[c].values():
            params.append([inputfile, inputgroup] + blocks)
    for c in this_run:
        del chromsizes[c]
        del allblocks[c]
    return params


def get_alignment_blocks(chainfile, tchroms, qchroms, logger):
    """
    :param chainfile:
    :param tchroms:
    :param qchroms:
    :param logger:
    :return:
    """
    opn, mode = text_file_mode(chainfile)
    tselect = re.compile('(' + '|'.join([c + '$' for c in tchroms]) + ')')
    qselect = re.compile('(' + '|'.join([c + '$' for c in qchroms]) + ')')
    logger.debug('Chain selectors compiled')
    blocks = col.defaultdict(dict)
    read_targets = set()
    read_queries = set()
    bc = 0
    with opn(chainfile, mode=mode, encoding='ascii') as cf:
        for alnblock in get_chain_iterator(cf, tselect, qselect):
            bc += 1
            tc = alnblock[0]
            qc = alnblock[4]
            read_targets.add(tc)
            read_queries.add(qc)
            try:
                blocks[qc][tc].append(alnblock)
            except KeyError:
                blocks[qc][tc] = [alnblock]
    # returns sorted by target/reference chromosome - start - end...
    logger.debug('Read target chromosomes: {}'.format(read_targets))
    logger.debug('Read query chromosomes: {}'.format(read_queries))
    logger.debug('Read total alignment blocks: {}'.format(bc))
    return blocks


def map_signal_data(params):
    """
    :param params:
    :return:
    """
    inputfile = params[0]
    inputgroup = params[1]
    alnblocks = params[2:]
    tchrom = alnblocks[0][0]
    qchrom = alnblocks[0][4]
    global _shm_carr
    carr = _shm_carr[qchrom]
    with pd.HDFStore(inputfile, 'r') as hdf:
        load_group = os.path.join(inputgroup, tchrom)
        sig_to_map = hdf[load_group].values
        for block in alnblocks:
            assert block[0] == tchrom, 'Target chromosome mismatch: {}'.format(block)
            assert block[4] == qchrom, 'Query chromosome mismatch: {}'.format(block)
            if block[7] == '-':
                carr[block[5]:block[6]] = sig_to_map[block[1]:block[2]][::-1]
            else:
                carr[block[5]:block[6]] = sig_to_map[block[1]:block[2]]
    return True


def run_map_signal(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    args.__dict__['keepchroms'] = args.keepchroms.strip('"')
    logger.debug('Chromosome select pattern for query [map to]: {}'.format(args.keepchroms))
    _, ingroup, infile = check_path_infos(args.inputfile, args.inputgroup)
    _, outgroup, outfile = check_path_infos(args.outputfile, args.outputgroup)
    qchroms = read_chromosome_sizes(args.querychroms, args.keepchroms)
    num_qchroms = len(qchroms)
    tchroms = get_chrom_list(infile, verify=True)
    logger.debug('Chromosomes in target data file [map from]: {}'.format(tchroms))
    meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB, 2)
    logger.debug('Before reading chains; active memory: {}GB'.format(meminfo))
    alnblocks = get_alignment_blocks(args.chainfile, tchroms, list(qchroms.keys()), logger)
    meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB, 2)
    logger.debug('All chains read; active memory: {}GB'.format(meminfo))
    _ = create_filepath(args.outputfile, logger)
    logger.debug('Processing {} query chromosomes at a time'.format(args.allocate))
    meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB, 2)
    logger.debug('Start processing, active memory: {}GB'.format(meminfo))
    with pd.HDFStore(outfile, args.filemode, complib='blosc', complevel=9, encoding='utf-8') as hdf:
        if '/metadata' in hdf:
            metadata = hdf['metadata']
        else:
            metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
        while len(qchroms) > 0:
            logger.debug('Query chromosomes left: {}'.format(len(qchroms)))
            global _shm_carr
            blocklists = assemble_worker_params(infile, ingroup, args.allocate, qchroms, alnblocks, _shm_carr)
            logger.debug('Processing query chromosomes: {}'.format(sorted(_shm_carr.keys())))
            meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB, 2)
            logger.debug('Parameter list assembled, active memory: {}GB'.format(meminfo))
            # the following only works on fork platforms (no Windows support)
            with mp.Pool(args.workers) as pool:
                _ = [res for res in pool.imap_unordered(map_signal_data, blocklists, chunksize=1)]
                logger.debug('Worker pool finished')
            for chrom, valobj in _shm_carr.items():
                grp, valobj, metadata = gen_obj_and_md(metadata, outgroup, chrom,
                                                       infile, pd.Series(valobj[:], dtype=np.float64))
                hdf.put(grp, valobj, format='fixed')
                hdf.flush()
                logger.debug('Stored data for chromosome {}'.format(chrom))
                meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB, 2)
                logger.debug('Active memory: {}GB'.format(meminfo))
            for k in list(_shm_carr.keys()):
                del _shm_carr[k]
        assert len(hdf.keys()) >= num_qchroms, 'Signal mapping incomplete: {}'.format(hdf.keys())
        hdf.put('metadata', metadata, format='table')
        hdf.flush()
        logger.debug('Metadata saved')
    return 0
