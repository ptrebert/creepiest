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
from crplib.auxiliary.text_parsers import read_chromosome_sizes, get_chain_iterator,\
    get_chain_positions, get_map_positions, get_map_iterator
from crplib.auxiliary.file_ops import text_file_mode, create_filepath
from crplib.auxiliary.hdf_ops import get_chrom_list,\
    check_path_infos, get_default_group
from crplib.auxiliary.constants import DIV_B_TO_GB


_shm_carr = dict()


def assemble_worker_params(inputfile, inputgroup, mapfile, numalloc, chromsizes, mappos, carrays):
    """
    :param inputfile:
    :param inputgroup:
    :param mapfile:
    :param numalloc:
    :param chromsizes:
    :param mappos:
    :param carrays:
    :return:
    """
    this_run = list(chromsizes.keys())[:numalloc]
    if not inputgroup:
        inputgroup = get_default_group(inputfile)
    for qchrom in this_run:
        # NB: lock=False is possible since all writes to the chromosome
        # arrays happen to disjoint slices
        carrays[qchrom] = mp.Array('d', np.zeros(chromsizes[qchrom], dtype=np.float64), lock=False)
    params = []
    commons = {'inputfile': inputfile, 'inputgroup': inputgroup, 'mapfile': mapfile}
    for qchrom in this_run:
        for tchrom, positions in mappos[qchrom].items():
            tmp = dict(commons)
            tmp['qchrom'] = qchrom
            tmp['tchrom'] = tchrom
            tmp['start_pos'] = min(positions)
            tmp['num_blocks'] = len(positions)
            params.append(tmp)
    for qchrom in this_run:
        del chromsizes[qchrom]
        del mappos[qchrom]
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


def build_chainfile_index(chainfile, tchroms, qchroms, logger):
    """
    :param chainfile:
    :param tchroms:
    :param qchroms:
    :param logger:
    :return:
    """
    tselect = re.compile('(' + '|'.join([c + '$' for c in tchroms]) + ')')
    qselect = re.compile('(' + '|'.join([c + '$' for c in qchroms]) + ')')
    logger.debug('Chain selectors compiled')
    chain_pos = get_chain_positions(chainfile, tselect, qselect, tref=False)
    logger.debug('Chain file index created')
    return chain_pos


def build_mapfile_index(mapfile, tchroms, qchroms, logger):
    """
    :param mapfile:
    :param tchroms:
    :param qchroms:
    :param logger:
    :return:
    """
    tselect = re.compile('(' + '|'.join([c + '$' for c in tchroms]) + ')')
    qselect = re.compile('(' + '|'.join([c + '$' for c in qchroms]) + ')')
    logger.debug('Map selectors compiled')
    map_pos = get_map_positions(mapfile, tselect, qselect, tref=False)
    logger.debug('Map file index created')
    return map_pos


def map_signal_data(params):
    """
    :param params:
    :return:
    """
    with pd.HDFStore(params['inputfile'], 'r') as hdf:
        load_group = os.path.join(params['inputgroup'], params['tchrom'])
        sig_to_map = hdf[load_group].values

    tchrom = params['tchrom']
    qchrom = params['qchrom']
    tselect = re.compile(tchrom + '$')
    qselect = re.compile(qchrom + '$')
    opn, mode = text_file_mode(params['mapfile'])
    global _shm_carr
    carr = _shm_carr[qchrom]
    with opn(params['mapfile'], mode) as maps:
        maps.seek(params['start_pos'])
        mapit = get_map_iterator(maps, tselect, qselect, read_num=params['num_blocks'])
        buffer = []
        for block in mapit:
            assert block[0] == tchrom, 'Target chromosome mismatch: {}'.format(block)
            assert block[5] == qchrom, 'Query chromosome mismatch: {}'.format(block)
            buffer.append(block)
            if len(buffer) >= 10000:
                for b in buffer:
                    if b[8] == '-':
                        carr[b[6]:b[7]] = sig_to_map[b[1]:b[2]][::-1]
                    else:
                        carr[b[6]:b[7]] = sig_to_map[b[1]:b[2]]
                buffer = []
        for b in buffer:
            if b[8] == '-':
                carr[b[6]:b[7]] = sig_to_map[b[1]:b[2]][::-1]
            else:
                carr[b[6]:b[7]] = sig_to_map[b[1]:b[2]]
    return True


def run_map_signal(args):
    """
    :param args:
    :return:
    """
    baseline_mem = round(psu.virtual_memory().active / DIV_B_TO_GB, 2)
    logger = args.module_logger
    args.__dict__['keepchroms'] = args.keepchroms.strip('"')
    logger.debug('Chromosome select pattern for query [map to]: {}'.format(args.keepchroms))
    _, ingroup, infile = check_path_infos(args.inputfile, args.inputgroup)
    _, outgroup, outfile = check_path_infos(args.outputfile, args.outputgroup)
    qchroms = read_chromosome_sizes(args.querychroms, args.keepchroms)
    num_qchroms = len(qchroms)
    tchroms = get_chrom_list(infile, verify=True)
    logger.debug('Chromosomes in target data file [map from]: {}'.format(tchroms))
    meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
    logger.debug('Before building map index; approx. used memory: {}GB'.format(meminfo))
    chainpos = build_mapfile_index(args.mapfile, tchroms, list(qchroms.keys()), logger)
    meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
    logger.debug('Map index built; active memory: {}GB'.format(meminfo))
    _ = create_filepath(args.outputfile, logger)
    logger.debug('Processing {} query chromosomes at a time'.format(args.allocate))
    meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
    logger.debug('Start processing, active memory: {}GB'.format(meminfo))
    with pd.HDFStore(outfile, args.filemode, complib='blosc', complevel=9, encoding='utf-8') as hdf:
        if '/metadata' in hdf:
            metadata = hdf['metadata']
        else:
            metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
        while len(qchroms) > 0:
            logger.debug('Query chromosomes left: {}'.format(len(qchroms)))
            global _shm_carr
            indexlists = assemble_worker_params(infile, ingroup, args.mapfile, args.allocate,
                                                qchroms, chainpos, _shm_carr)
            logger.debug('Processing query chromosomes: {}'.format(sorted(_shm_carr.keys())))
            logger.debug('Parameter list of size {} created'.format(len(indexlists)))
            meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
            logger.debug('Parameter list assembled, active memory: {}GB'.format(meminfo))
            # the following only works on fork platforms (no Windows support)
            with mp.Pool(args.workers) as pool:
                _ = [res for res in pool.imap_unordered(map_signal_data, indexlists, chunksize=1)]
                logger.debug('Worker pool finished')
            for chrom, valobj in _shm_carr.items():
                grp, valobj, metadata = gen_obj_and_md(metadata, outgroup, chrom,
                                                       infile, pd.Series(valobj[:], dtype=np.float64))
                hdf.put(grp, valobj, format='fixed')
                hdf.flush()
                logger.debug('Stored data for chromosome {}'.format(chrom))
                meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
                logger.debug('Active memory: {}GB'.format(meminfo))
            for k in list(_shm_carr.keys()):
                del _shm_carr[k]
        assert len(hdf.keys()) >= num_qchroms, 'Signal mapping incomplete: {}'.format(hdf.keys())
        hdf.put('metadata', metadata, format='table')
        hdf.flush()
        logger.debug('Metadata saved')
    return 0
