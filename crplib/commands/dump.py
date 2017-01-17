# coding=utf-8

"""
Small module to dump data from HDF format files into text formats (BED-like)
"""

import sys as sys

import numpy as np
import pandas as pd

from crplib.auxiliary.file_ops import text_file_mode
from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups


def set_output(outfile):
    """
    :param outfile:
    :return:
    """
    if outfile == 'stdout':
        return sys.stdout
    else:
        opn, mode = text_file_mode(outfile, read=False)
        return opn(outfile, mode)


def get_index_column_order(mapref, full):
    """
    :param mapref:
    :param full:
    :return:
    """
    if mapref == 'target' and full:
        col_order = ['tchrom', 'tstart', 'tend', 'tstrand', 'blockid',
                     'qchrom', 'qstart', 'qend', 'qstrand']
        sort_cols = ['tstart', 'tend']
    elif mapref == 'target' and not full:
        col_order = ['tchrom', 'tstart', 'tend', 'blockid']
        sort_cols = ['tstart', 'tend']
    elif mapref == 'query' and full:
        col_order = ['qchrom', 'qstart', 'qend', 'tstrand', 'blockid',
                     'tchrom', 'tstart', 'tend', 'qstrand']
        sort_cols = ['qstart', 'qend']
    elif mapref == 'query' and not full:
        col_order = ['qchrom', 'qstart', 'qend', 'blockid']
        sort_cols = ['qstart', 'qend']
    else:
        raise ValueError('Cannot handle combination of map reference {} and {}'.format(mapref, full))
    return col_order, sort_cols


def dump_index(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    with pd.HDFStore(args.inputfile, 'r') as hdf:
        block_groups = get_valid_hdf5_groups(args.inputfile, '/qt/blocks')
        assert block_groups, 'Map index file does not contain standard groups: /qt/blocks/qchrom/tchrom'
        logger.debug('Identified {} blocks in map file'.format(len(block_groups)))
        query_chroms = sorted(set([b.rsplit('/', 2)[1] for b in block_groups]))
        target_chroms = sorted(set([b.rsplit('/', 1)[1] for b in block_groups]))
        if args.mapreference == 'target':
            block_filter = '/qt/blocks/{bottom}/{top}'
            top_chroms = target_chroms
            bottom_chroms = query_chroms
        else:
            block_filter = '/qt/blocks/{top}/{bottom}'
            top_chroms = query_chroms
            bottom_chroms = target_chroms
        iter_done = False
        try:
            out_dest = set_output(args.outputfile)
            for top in top_chroms:
                for bottom in bottom_chroms:
                    selector = block_filter.format(**{'top': top, 'bottom': bottom})
                    if selector not in block_groups:
                        continue
                    logger.debug('Dumping block: {}'.format(selector))
                    _, qchrom, tchrom = selector.rsplit('/', 2)
                    blocks = hdf[selector]
                    blocks['tchrom'] = tchrom
                    blocks['qchrom'] = qchrom
                    blocks['tstrand'] = '+'
                    blocks.replace({'qstrand': {1: '+', -1: '-'}}, inplace=True)
                    blocks['blockid'] = blocks.index
                    order, sort_by = get_index_column_order(args.mapreference, args.fullblocks)
                    blocks = blocks[order]
                    blocks.sort_values(sort_by, axis=0, inplace=True)
                    logger.debug('Writing {} rows...'.format(blocks.shape[0]))
                    blocks.to_csv(out_dest, sep='\t', header=False, index=False)
            iter_done = True
            out_dest.close()
        except ValueError:
            if not iter_done:
                logger.error('Error raised before dumping map file completed')
                raise
    logger.debug('Dumped map index')
    return


def dump_signal(args, logger):
    """
    :return:
    """
    stat = {'mean': np.mean, 'max': np.max, 'min': np.min,
            'median': np.median, 'sum': np.sum, 'product': np.product}
    raise NotImplementedError()


def run_dump_data(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    logger.debug('Dumping data from HDF file {}'.format(args.inputfile))
    if args.mapreference:
        logger.debug('Map reference specified as {} - assuming map index file.'.format(args.mapreference))
        dump_index(args, logger)
    elif args.summstat:
        logger.debug('Summary statistic set to {} - assuming signal file.'.format(args.summstat))
        dump_signal(args, logger)
    else:
        raise NotImplementedError('Dumping region files not implemented')

    return 0
