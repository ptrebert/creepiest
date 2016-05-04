# coding=utf-8

"""
Map a signal track from one species/assembly to another
"""

import numpy as np
import pandas as pd
import re as re

from crplib.metadata.md_signal import gen_obj_and_md, MD_SIGNAL_COLDEFS
from crplib.auxiliary.text_parsers import read_chromosome_sizes, get_chain_iterator
from crplib.auxiliary.file_ops import text_file_mode, create_filepath
from crplib.auxiliary.hdf_ops import get_valid_hdf5_groups


def allocate_chrom_arrays(chroms):
    """
    :param chroms:
    :return:
    """
    ret = dict()
    for name, size in chroms:
        ret[name] = np.zeros(size, dtype='float64')
    return ret


def get_alignment_blocks(chainfile, chroms, chromre, logger):
    """
    :param chainfile:
     :type: str
    :param chroms:
     :type: list of tuples of (str, int)
    :param chromre:
     :type: str
    :return:
    """
    opn, mode = text_file_mode(chainfile)
    accept_query = [cname + '$' for cname, size in chroms]
    query_re = re.compile('(' + '|'.join(accept_query) + ')')
    logger.debug('Filtering for query chroms: {}'.format(query_re.pattern))
    blocks = []
    check_target = re.compile(chromre)
    read_targets = set()
    read_queries = set()
    with opn(chainfile, mode) as cf:
        # TODO eventually
        # right now, chain iterator can only filter for one or all targets
        # not for several... could be changed
        for alnblock in get_chain_iterator(cf, qcheck=query_re):
            if check_target.match(alnblock[0]) is not None:
                read_targets.add(alnblock[0])
                read_queries.add(alnblock[4])
                blocks.append(alnblock)
    # returns sorted by target/reference chromosome - start - end...
    logger.debug('Read target chromosomes: {}'.format(read_targets))
    logger.debug('Read query chromosomes: {}'.format(read_queries))
    return sorted(blocks)


def map_signal_data(sigfile, grouproot, alnblocks, qchroms):
    """
    :param sigfile:
    :param grouproot:
    :param alnblocks:
    :param qchroms:
    :return:
    """
    datagroups = get_valid_hdf5_groups(sigfile, grouproot)
    sig_to_map = None
    tchrom = None
    # note here: if alnblocks is empty, returns all-zero query chroms
    for block in alnblocks:
        if block[0] != tchrom:
            chrom_group = [g for g in datagroups if g.endswith(block[0])]
            assert len(chrom_group) == 1, 'No or ambiguous reference chromosome: {}'.format(block[0])
            with pd.HDFStore(sigfile, 'r') as infile:
                sig_to_map = infile[chrom_group[0]]
            tchrom = block[0]
        if block[7] == '-':
            qchroms[block[4]][block[5]:block[6]] = sig_to_map[block[1]:block[2]][::-1]
        else:
            qchroms[block[4]][block[5]:block[6]] = sig_to_map[block[1]:block[2]]
    return qchroms


def run_map_signal(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    args.__dict__['keepchroms'] = args.keepchroms.strip('"')
    logger.debug('Chromosome select pattern: {}'.format(args.keepchroms))
    csizes = read_chromosome_sizes(args.querychroms, args.keepchroms)
    csizes = list(csizes.items())
    step = args.allocate
    _ = create_filepath(args.outputfile, logger)
    logger.debug('Processing {} chromosomes at a time'.format(step))
    for idx in range(0, len(csizes), step):
        try:
            proc_chroms = csizes[idx:idx+step]
        except IndexError:
            proc_chroms = csizes[idx:]
        logger.debug('Processing chromosomes: {}'.format(proc_chroms))
        alnblocks = get_alignment_blocks(args.chainfile, proc_chroms, args.keepchroms, logger)
        if not alnblocks:
            logger.warning('No alignment blocks for query chromosomes: {}'.format(proc_chroms))
        else:
            logger.debug('Read {} alignment blocks from chain file'.format(len(alnblocks)))
        logger.debug('Allocating memory for mapped signal')
        chroms = allocate_chrom_arrays(proc_chroms)
        logger.debug('Mapping signal data')
        chroms = map_signal_data(args.inputfile, args.inputgroup, alnblocks, chroms)
        logger.debug('Mapping complete')
        with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdfout:
            if 'metadata' in hdfout:
                metadata = hdfout['metadata']
            else:
                metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
            for chrom, valobj in chroms.items():
                logger.debug('Saving mapped signal for chromosome {}'.format(chrom))
                grp, valobj, metadata = gen_obj_and_md(metadata, args.outputgroup, chrom, args.inputfile, valobj)
                hdfout.put(grp, valobj, format='fixed')
                hdfout.flush()
            hdfout.put('metadata', metadata, format='table')
            hdfout.flush()
    return 0
