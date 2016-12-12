# coding=utf-8

"""
Module to convert various data formats to HDF5
"""

import os as os
import importlib as imp

from crplib.auxiliary.hdf_ops import check_path_infos


def convert_bedgraph_signal(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    assert os.path.isfile(args.chromsizes),\
        'Invalid path to chromosome sizes file: {}'.format(args.chromsizes)
    assert all([os.path.isfile(f) for f in args.inputfiles]),\
        'Invalid path(s) to input file(s): {}'.format(args.inputfiles)
    assert 0. <= args.clip <= 100., 'Clip value outside of range 0...100: {}'.format(args.clip)
    _, grp, fp = check_path_infos(args.outputfile, args.outputgroup)
    setattr(args, 'outputfile', fp)
    setattr(args, 'outputgroup', grp)
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


def convert_chain_file(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    assert os.path.isfile(args.inputfile[0]), 'Invalid path to input file: {}'.format(args.inputfile)
    mod = imp.import_module('crplib.commands.convert_chain')
    rv = mod.run_chain_conversion(args, logger)
    assert os.path.isfile(args.outputfile), 'No output file created - conversion failed? {}'.format(args.outputfile)
    return rv


def convert_motifdb_file(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    assert os.path.isfile(args.inputfile[0])
    mod = imp.import_module('crplib.commands.convert_motifdb')
    rv = mod.run_motifdb_conversion(args, logger)
    assert os.path.isfile(args.outputfile), 'No output file created - conversion failed? {}'.format(args.outputfile)
    return rv


def convert_map_file(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    assert os.path.isfile(args.inputfiles[0]), 'Invalid path to input file: {}'.format(args.inputfiles)
    mod = imp.import_module('crplib.commands.convert_map')
    rv = mod.run_map_conversion(args, logger)
    assert os.path.isfile(args.outputfile), 'No output file created - conversion failed? {}'.format(args.outputfile)
    return rv


def run_conversion(args):
    """
    :param args: command line parameters
    :type args: Namespace object
    :return: exit code
    :rtype: int
    """
    logger = args.module_logger
    try:
        logger.debug('Running conversion type: {}'.format(args.task))
        convtype = {'signal': convert_bedgraph_signal,
                    'region': convert_genomic_region,
                    'chain': convert_chain_file,
                    'tfmotif': convert_motifdb_file,
                    'map': convert_map_file}
        convcall = convtype[args.task]
        setattr(args, 'selectchroms', args.selectchroms.strip('"'))
        logger.debug('Chromosome select pattern: {}'.format(args.selectchroms))
        rv = convcall(args, logger)
        logger.debug('Conversion complete')
        return rv
    except Exception as e:
        logger.error('Error during conversion: {}'.format(str(e)))
        raise e
