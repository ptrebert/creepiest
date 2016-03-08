# coding=utf-8

"""
Module to handle conversion of region files into HDF5 format
"""

import pandas as pd
import multiprocessing as mp
import collections as col
import re as re
import operator as op
import scipy.stats as stats

from crplib.auxiliary.file_ops import text_file_mode
from crplib.metadata.md_region import gen_obj_and_md, MD_SIGNAL_COLDEFS

from crplib.auxiliary.constants import DIV_B_TO_GB


def assemble_worker_args(args):
    """
    :param args:
    :return:
    """
    cols = (0, 1, 2)
    if args.namecol != -1:
        cols += args.namecol,
    if args.scorecol != -1:
        cols += args.scorecol
    arglist = []
    for fp in args.inputfile:
        commons = dict()
        commons['inputfile'] = fp
        commons['keepchroms'] = args.keepchroms
        commons['keeptop'] = args.keeptop
        commons['scorecol'] = args.scorecol
        commons['columns'] = cols
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
        if this[0] == step[0]:
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


def process_regions(kwargs):
    """
    :param kwargs:
    :return:
    """
    mypid = mp.current_process().pid
    fpath = kwargs['inputfile']
    chr_match = re.compile(kwargs['keepchroms'])
    getvals = op.itemgetter(*kwargs['columns'])
    regions = []
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        for line in infile:
            if not line or chr_match.match(line) is None:
                continue
            regions.append(getvals(line.split()))
    if kwargs['scorecol'] != -1 and kwargs['keeptop'] < 100.:
        # convention here: score is always last index by construction
        # and assume ranking where highest score is first one
        scores = [float(reg[-1]) for reg in regions]
        thres = stats.scoreatpercentile(scores, 100 - kwargs['keeptop'])
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
        out.append(reg + name,)
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
                metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
            mapres = pool.map_async(process_regions, arglist, chunksize=1)
            all_regions = col.deque()
            chroms = set()
            for pid, regobj in mapres.get():
                logger.debug('Worker (PID {}) completed chromosome {}'.format(pid))
                # collect all chromosomes in dataset(s)
                [chroms.add(reg[0]) for reg in regobj]
                all_regions.extend(regobj)
            logger.debug('All input files processed, sorting regions...')
            all_regions = sorted(all_regions)
            if len(args.inputfile) == 1:
                pass  # nothing more to do ?
            else:
                all_regions = merge_overlapping_regions(all_regions)
            if args.namecol == -1:
                all_regions = add_names(all_regions)
            for chrom in sorted(chroms):
                grp, valobj, metadata = gen_obj_and_md(metadata, args.grouproot, chrom, args.inputfile,
                                                       [reg for reg in all_regions if reg[0] == chrom])
                hdfout.put(grp, valobj, format='table')
                hdfout.flush()
                logger.debug('Processed chromosome {}'.format(chrom))
        _, _, metadata = gen_obj_and_md(metadata, args.grouproot, 'wg', args.inputfile, all_regions)
        hdfout.put('metadata', metadata, format='table')
    logger.debug('HDF file closed: {}'.format(args.outputfile))
    return 0
