# coding=utf-8

"""
Small module to dump data from HDF format files into text formats (BED-like)
"""

import sys as sys

import numpy as np
import pandas as pd

from crplib.auxiliary.file_ops import text_file_mode
from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups, determine_crp_filetype


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


def rearrange_columns(dataset, chrom_prefix):
    """
    :param dataset:
    :return:
    """
    idx_count = 9
    known_cols = {'chrom': 0, 'start': 1, 'end': 2, 'name': 3,
                  'score': 4, 'strand': 5, 'signalValue': 6,
                  'pValue': 7, 'qValue': 8}
    col_order = []
    for c in dataset.columns:
        try:
            idx = known_cols[c]
            col_order.append((idx, c))
        except KeyError:
            col_order.append((idx_count, c))
            idx_count += 1
    col_order = sorted(col_order, key=lambda x: x[0])
    col_order = [x[1] for x in col_order]
    dataset = dataset[col_order]
    if chrom_prefix:
        prefix_cols = dataset.columns.tolist()
        prefix_cols[0] = '#chrom'
        dataset.columns = prefix_cols
    return dataset


def dump_regions(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    load_groups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    dump = []
    with pd.HDFStore(args.inputfile, 'r') as hdf:
        for lg in load_groups:
            _, chrom_name = lg.rsplit('/', 1)
            data = hdf[lg]
            data['chrom'] = chrom_name
            assert 'crp_group_index' not in data.columns, 'Column name duplicate: crp_group_index'
            data['crp_group_index'] = data.index.tolist()
            assert 'crp_group_path' not in data.columns, 'Column name duplicate: crp_group_path'
            data['crp_group_path'] = lg
            dump.append(data)
    logger.debug('Concatenating all data groups...')
    dump = pd.concat(dump, axis=0, join='outer', ignore_index=True)
    dump.sort_values(by=['chrom', 'start', 'end'], inplace=True)
    dump.reset_index(drop=True, inplace=True)
    logger.debug('Final dataset size: {}'.format(dump.shape))
    dump = rearrange_columns(dump, args.commentheader)
    col_delim = args.delimiter.strip('"')
    if args.outputfile in ['stdout', '-']:
        logger.debug('Dumping data to stdout')
        out = sys.stdout
    else:
        logger.debug('Dumping data to file')
        out = args.outputfile
    if args.noheader and not args.rtable:
        dump.to_csv(out, sep=col_delim, header=False, index=False)
    elif args.noheader and args.rtable:
        dump.to_csv(out, sep=col_delim, header=False, index=True, index_label=None)
    elif not args.noheader and not args.rtable:
        dump.to_csv(out, sep=col_delim, header=True, index=False, index_label=None)
    elif not args.noheader and args.rtable:
        dump.to_csv(out, sep=col_delim, header=True, index=True, index_label=None)
    else:
        raise RuntimeError('Impossible combination of parameters...')
    return


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
    else:
        logger.debug('Determining datatype of input file...')
        dtype = determine_crp_filetype(args.inputfile)
        if dtype is None:
            # this implies a valid metadata record that does not match
            # any of the hard-coded column definitions - possibly, somebody
            # has been tinkering with this file
            raise RuntimeError('Could not determine datatype despite valid metadata record.\n'
                               'Have you manually edited this file? If not, please contact the developer.')
        elif dtype == 'signal':
            dump_signal(args, logger)
        elif dtype == 'regions':
            dump_regions(args, logger)
        elif dtype == 'features':
            raise NotImplementedError('DTYPE {} not implemented'.format(dtype))
        else:
            raise ValueError('Unknown datatype in inputfile: {}'.format(dtype))
    return 0
