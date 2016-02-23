# coding=utf-8

"""
Module to convert various data formats to HDF5
"""

import os as os
import importlib as imp


def convert_fasta_genome(args, logger):
    """
    :param args:
    :return:
    """
    assert os.path.isfile(args.inputfile), 'Invalid path to input file: {}'.format(args.input)
    assert os.path.isfile(args.chromsizes), 'Invalid path to chromosome sizes file: {}'.format(args.chromsizes)
    # remove double quotes from pattern
    args.__dict__['keepchroms'] = args.keepchroms.strip('"')
    logger.debug('Chromosome select pattern: {}'.format(args.keepchroms))
    mod = imp.import_module('crplib.commands.convert_fasta')
    rv = mod.run_fasta_conversion(args, logger)
    assert os.path.isfile(args.outputfile), 'No output file created - conversion failed? {}'.format(args.outputfile)
    return rv


def convert_bedgraph_signal(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    assert os.path.isfile(args.chromsizes),\
        'Invalid path to chromosome sizes file: {}'.format(args.chromsizes)
    assert all([os.path.isfile(f) for f in args.inputfile]),\
        'Invalid path(s) to input file(s): {}'.format(args.inputfile)
    args.__dict__['keepchroms'] = args.keepchroms.strip('"')
    logger.debug('Chromosome select pattern: {}'.format(args.keepchroms))
    mod = imp.import_module('crplib.commands.convert_bedgraph')
    rv = mod.run_bedgraph_conversion(args, logger)
    assert os.path.isfile(args.outputfile), 'No output file created - conversion failed? {}'.format(args.outputfile)
    return rv


def run_conversion(args):
    """
    :param args: command line parameters
     :type: Namespace object
    :return: exit code
     :rtype: int
    """
    logger = args.module_logger
    rv = 0
    try:
        logger.debug('Running conversion type: {}'.format(args.subparser_name))
        convtype = {'convfna': convert_fasta_genome,
                    'convbg': convert_bedgraph_signal,
                    'convreg': lambda x: x}
        convcall = convtype[args.subparser_name]
        rv = convcall(args, logger)
        logger.debug('Conversion complete')
        return rv
    except Exception as e:
        logger.error('Error: {}'.format(str(e)))
        raise e
