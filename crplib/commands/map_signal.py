# coding=utf-8

"""
Map a signal track from one species/assembly to another
"""

import os as os
import numpy as np
import pandas as pd
import collections as col
import multiprocessing as mp
import psutil as psu

from crplib.metadata.md_signal import gen_obj_and_md, MD_SIGNAL_COLDEFS
from crplib.auxiliary.file_ops import create_filepath
from crplib.auxiliary.hdf_ops import get_chrom_list,\
    check_path_infos, get_default_group, extract_chromsizes_from_map
from crplib.auxiliary.constants import DIV_B_TO_GB, MAPIDX_BLOCKS


_shm_carr = dict()


def assemble_worker_params(inputfile, inputgroup, mapfile, numalloc, tchroms, qchroms, carrays):
    """
    :param inputfile:
    :param inputgroup:
    :param mapfile:
    :param numalloc:
    :param tchroms:
    :param qchroms:
    :param carrays:
    :return:
    """
    this_run = list(qchroms.keys())[:numalloc]
    if not inputgroup:
        inputgroup = get_default_group(inputfile)
    for qchrom in this_run:
        # NB: lock=False is possible since all writes to the chromosome
        # arrays happen to disjoint slices
        carrays[qchrom] = mp.Array('d', np.zeros(qchroms[qchrom], dtype=np.float64), lock=False)
    params = []
    commons = {'inputfile': inputfile, 'inputgroup': inputgroup,
               'mapfile': mapfile}
    base_path = os.path.join('/qt', MAPIDX_BLOCKS, '')
    trg_to_load = col.defaultdict(list)
    with pd.HDFStore(mapfile, 'r') as hdf:
        chrom_maps = list(hdf.keys())
        for cm in chrom_maps:
            if cm.startswith(base_path):
                query, target = os.path.split(cm.replace(base_path, ''))
                if query in this_run and target in tchroms:
                    trg_to_load[target].append((query, cm))
    for target, to_load in trg_to_load.items():
        tmp = dict(commons)
        tmp['tchrom'] = target
        tmp['blocks'] = to_load
        params.append(tmp)
    for qchrom in this_run:
        del qchroms[qchrom]
    assert params, 'No parameter list for workers created'
    return params


def map_row(row, target, query):
    """
    Row contains: tstart - tend - qstart - qend - qstrand
    :param row:
    :param target:
    :param query:
    :return:
    """
    query[row[2]:row[3]] = target[row[0]:row[1]][::row[4]]
    # this was for debugging
    # return np.sum(target[row[0]:row[1]] > 0)
    return row[1] - row[0]


def map_signal_data(params):
    """
    :param params:
    :return:
    """
    results = []
    tchrom = params['tchrom']
    with pd.HDFStore(params['inputfile'], 'r') as hdf:
        load_group = os.path.join(params['inputgroup'], params['tchrom'])
        sig_to_map = hdf[load_group].values
    global _shm_carr
    for qchrom, bpath in params['blocks']:
        assert qchrom in bpath and tchrom in bpath, \
            'Wrong block path {}: q {} - t {}'.format(bpath, qchrom, tchrom)
        with pd.HDFStore(params['mapfile'], 'r') as hdf:
            # load Pandas DataFrame from map file that holds all blocks
            # describing a mapping between query and target for one
            # particular chromosome combination
            mapblocks = hdf[bpath]
        carr = _shm_carr[qchrom]
        cov = mapblocks.apply(map_row, axis=1, raw=True, args=(sig_to_map, carr))
        results.append((qchrom, tchrom, cov.sum()))
    assert results, 'No data processed for parameter set: {}'.format(params)
    return results


def run_map_signal(args):
    """
    :param args:
    :return:
    """
    baseline_mem = round(psu.virtual_memory().active / DIV_B_TO_GB, 2)
    logger = args.module_logger
    setattr(args, 'selectchroms', args.selectchroms.strip('"'))
    logger.debug('Chromosome select pattern for query [map to]: {}'.format(args.selectchroms))
    _, ingroup, infile = check_path_infos(args.inputfile, args.inputgroup)
    _, outgroup, outfile = check_path_infos(args.outputfile, args.outputgroup)
    qchroms = extract_chromsizes_from_map(args.mapfile, 'query', args.selectchroms)
    num_qchroms = len(qchroms)
    tchroms = get_chrom_list(infile, verify=True)
    logger.debug('Chromosomes in target data file [map from]: {}'.format(tchroms))
    meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
    logger.debug('Occupied RAM: {}GB'.format(meminfo))
    _ = create_filepath(args.outputfile, logger)
    logger.debug('Processing {} query chromosomes at a time'.format(args.allocate))
    meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
    logger.debug('Start processing - occupied RAM: {}GB'.format(meminfo))
    global _shm_carr
    with pd.HDFStore(outfile, args.filemode, complib='blosc', complevel=9, encoding='utf-8') as hdf:
        if '/metadata' in hdf:
            metadata = hdf['metadata']
        else:
            metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
        while len(qchroms) > 0:
            logger.debug('Query chromosomes left: {}'.format(len(qchroms)))
            indexlists = assemble_worker_params(infile, ingroup, args.mapfile,
                                                args.allocate, tchroms, qchroms, _shm_carr)
            logger.debug('Processing query chromosomes: {}'.format(sorted(_shm_carr.keys())))
            logger.debug('Parameter list of size {} created'.format(len(indexlists)))
            meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
            logger.debug('Parameter list assembled - occupied RAM: {}GB'.format(meminfo))
            check_coverage = col.defaultdict(int)
            # the following only works on fork platforms (no Windows support)
            with mp.Pool(args.workers) as pool:
                resit = pool.imap_unordered(map_signal_data, indexlists, chunksize=1)
                for res in resit:
                    for item in res:
                        qchrom, tchrom, qcov = item
                        check_coverage[qchrom] += qcov
                logger.debug('Worker pool finished')
            for chrom, valobj in _shm_carr.items():
                grp, valobj, metadata = gen_obj_and_md(metadata, outgroup, chrom,
                                                       infile, pd.Series(valobj[:], dtype=np.float64))
                # valobj is a pandas.Series at this point
                actual_cov = (valobj > 0).sum()
                logger.debug('Non-zero coverage in chrom {}: {} (mapped pos.: {})'.format(chrom, actual_cov, check_coverage[chrom]))
                hdf.put(grp, valobj, format='fixed')
                hdf.flush()
                logger.debug('Stored data for chromosome {}'.format(chrom))
                meminfo = round(psu.virtual_memory().active / DIV_B_TO_GB - baseline_mem, 2)
                logger.debug('Occupied RAM: {}GB'.format(meminfo))
            for k in list(_shm_carr.keys()):
                del _shm_carr[k]
        assert len(hdf.keys()) >= num_qchroms, 'Signal mapping incomplete: {}'.format(hdf.keys())
        hdf.put('metadata', metadata, format='table')
        hdf.flush()
        logger.debug('Metadata saved')
    return 0
