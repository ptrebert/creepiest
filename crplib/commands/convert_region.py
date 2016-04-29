# coding=utf-8

"""
Module to handle conversion of region files into HDF5 format
"""

import pandas as pd
import multiprocessing as mp
import collections as col
import numpy as np
import re as re
import operator as op
import scipy.stats as stats

from crplib.auxiliary.file_ops import text_file_mode
from crplib.metadata.md_regions import gen_obj_and_md, MD_REGION_COLDEFS


def assemble_worker_args(args):
    """
    :param args:
    :return:
    """
    cols = (0, 1, 2)
    if args.nameidx != -1:
        cols += args.nameidx,
    if args.scoreidx != -1:
        cols += args.scoreidx,
    arglist = []
    for fp in args.inputfile:
        commons = dict()
        commons['inputfile'] = fp
        commons['keepchroms'] = args.keepchroms
        commons['keeptop'] = args.keeptop
        commons['scoreidx'] = args.scoreidx
        commons['columns'] = cols
        commons['filtersize'] = args.filtersize
        arglist.append(commons)
    return arglist


def merge_overlapping_regions(allregions):
    """
    :param allregions:
    :return:
    """
    merged = []
    this = allregions.popleft()
    while 1:
        try:
            step = allregions.popleft()
        except IndexError:  # deque empty
            break
        if this[0] == step[0]:  # check chroms are identical
            # if this_end < next_start
            if this[2] < step[1]:
                merged.append(this)
                this = step
                continue
            else:
                # since all regions are sorted, must overlap now,
                # or be at least book-ended
                if len(this) == 3:
                    this = (this[0], min(this[1], step[1]), max(this[2], step[2]))
                elif len(this) == 4:
                    this = (this[0], min(this[1], step[1]), max(this[2], step[2]), this[3] + '-' + step[3])
                else:
                    raise ValueError('Unexpected number of components for region: {}'.format(this))
        else:
            merged.append(this)
            this = step
            continue
    merged.append(this)
    return merged


def process_regions(params):
    """
    :param params:
    :return:
    """
    mypid = mp.current_process().pid
    fpath = params['inputfile']
    chr_match = re.compile(params['keepchroms'])
    getvals = op.itemgetter(*params['columns'])
    filter_size = params['filtersize']
    regions = []
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        for line in infile:
            if not line or chr_match.match(line) is None:
                continue
            reg = getvals(line.split())
            if int(reg[2]) - int(reg[1]) < filter_size:
                continue
            regions.append(reg)
    assert len(regions) > 0,\
        'No regions selected for file {} and pattern {}'.format(fpath, params['keepchroms'])
    if params['scoreidx'] != -1 and params['keeptop'] < 100.:
        # convention here: score is always last index by construction
        # and assume ranking where highest score is first one
        scores = np.array([float(reg[-1]) for reg in regions])
        # this is heuristic to check if the selected score column makes sense
        assert np.var(scores) > 0,\
            'Scores have 0 variance for file {} and column {}'.format(fpath, params['columns'][-1])
        thres = stats.scoreatpercentile(scores, 100 - params['keeptop'])
        regions = [reg[:-1] for reg in regions if float(reg[-1]) > thres]
    if len(regions[0]) == 3:
        regions = [(reg[0], int(reg[1]), int(reg[2])) for reg in regions]
    elif len(regions[0]) == 4:
        # means there was a name column specified
        regions = [(reg[0], int(reg[1]), int(reg[2]), reg[3]) for reg in regions]
    else:
        raise ValueError('Unexpected number of components for region: {}'.format(regions[0]))
    return mypid, regions


def add_names(allregions):
    """
    :param allregions:
    :return:
    """
    out = []
    for idx, reg in enumerate(allregions, start=1):
        name = 'region_' + reg[0] + '_' + str(idx)
        out.append(reg + (name,))
    return out


def run_region_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    arglist = assemble_worker_args(args)
    logger.debug('Start processing {} region file(s)'.format(len(args.inputfile)))
    with pd.HDFStore(args.outputfile, 'a', complevel=9, complib='blosc') as hdfout:
        with mp.Pool(args.workers) as pool:
            if 'metadata' in hdfout:
                metadata = hdfout['metadata']
            else:
                metadata = pd.DataFrame(columns=MD_REGION_COLDEFS)
            all_regions = list()
            chroms = set()
            logger.debug('Iterating results')
            mapres = pool.map_async(process_regions, arglist)
            for pid, regobj in mapres.get():
                logger.debug('Worker (PID {}) completed, returned {} regions'.format(pid, len(regobj)))
                # collect all chromosomes in dataset(s)
                [chroms.add(reg[0]) for reg in regobj]
                all_regions.extend(regobj)
            # TODO
            # below here: looks like it could be simplified...
            logger.debug('All files processed, sorting {} regions...'.format(len(all_regions)))
            assert len(all_regions) > 0, 'No regions selected by worker processes: {}'.format(args.inputfile)
            all_regions = sorted(all_regions)
            if len(args.inputfile) == 1:
                pass  # nothing more to do ?
            else:
                all_regions = merge_overlapping_regions(col.deque(all_regions))
                logger.debug('After merging {} regions left'.format(len(all_regions)))
            if args.nameidx == -1:
                all_regions = add_names(all_regions)
            logger.debug('Identified {} chromosomes in dataset(s)'.format(len(chroms)))
            for chrom in sorted(chroms):
                grp, valobj, metadata = gen_obj_and_md(metadata, args.outputgroup, chrom, args.inputfile,
                                                       [reg for reg in all_regions if reg[0] == chrom])
                hdfout.put(grp, valobj, format='fixed')  # not sure here... usually replace entire object I guess
                hdfout.flush()
                logger.debug('Processed chromosome {}'.format(chrom))
        hdfout.put('metadata', metadata, format='table')
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
