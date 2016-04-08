# coding=utf-8

"""
Module to convert various data formats to HDF5
"""

import os as os
import importlib as imp


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
    assert 0. <= args.clip <= 100., 'Clip value outside of range 0...100: {}'.format(args.clip)
    mod = imp.import_module('crplib.commands.convert_bedgraph')
    rv = mod.run_bedgraph_conversion(args, logger)
    assert os.path.isfile(args.outputfile), 'No output file created - conversion failed? {}'.format(args.outputfile)
    return rv


def convert_genomic_region(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    assert all([os.path.isfile(f) for f in args.inputfile]), \
        'Invalid path(s) to input file(s): {}'.format(args.inputfile)
    assert 1. <= args.keeptop <= 100., \
        'Keeping top N percent value outside of range 1...100: {}'.format(args.keeptop)
    mod = imp.import_module('crplib.commands.convert_region')
    rv = mod.run_region_conversion(args, logger)
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
    try:
        logger.debug('Running conversion type: {}'.format(args.subparser_name))
        convtype = {'signal': convert_bedgraph_signal,
                    'region': convert_genomic_region}
        convcall = convtype[args.task]
        args.__dict__['keepchroms'] = args.keepchroms.strip('"')
        logger.debug('Chromosome select pattern: {}'.format(args.keepchroms))
        rv = convcall(args, logger)
        logger.debug('Conversion complete')
        return rv
    except Exception as e:
        logger.error('Error during conversion: {}'.format(str(e)))
        raise e
