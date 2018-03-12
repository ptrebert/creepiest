# coding=utf-8

"""
This module defines all sub-parsers for the different commands
and contains the appropriate execute functions
"""

import importlib as implib

from crplib.auxiliary.cmdline_parameters import *


def add_sub_parsers(main_parser):
    """
    :param main_parser:
    :return:
    """
    subparsers = main_parser.add_subparsers(dest='subparser_name', title='Subcommands')
    subparsers = _add_tests_command(subparsers)
    subparsers = _add_info_command(subparsers)
    subparsers = _add_dump_command(subparsers)
    subparsers = _add_convert_command(subparsers)
    subparsers = _add_compfeat_command(subparsers)
    subparsers = _add_train_command(subparsers)
    subparsers = _add_mapsig_command(subparsers)
    subparsers = _add_apply_command(subparsers)
    subparsers = _add_correlation_command(subparsers)
    subparsers = _add_match_command(subparsers)
    subparsers = _add_merge_command(subparsers)
    subparsers = _add_eval_command(subparsers)
    subparsers = _add_norm_command(subparsers)
    return main_parser


def _add_tests_command(subparsers):
    """
    :param subparser:
    :return:
    """
    parser_tests = subparsers.add_parser('tests',
                                         help='Run unit tests',
                                         description='This command executes all unit tests, among them is an '
                                                     'import test that attempts to import all modules necessary '
                                                     'to run any CREEPIEST tool. It is highly recommended to run '
                                                     'unit tests prior to any analysis run.')
    parser_tests.set_defaults(execute=_tests_execute)
    return subparsers


def _tests_execute(args):
    """
    :param args:
    :return:
    """
    tests = implib.import_module('crplib.commands.tests')
    module_path = tests.__file__
    retval = tests.run_tests(args, module_path)
    return retval


def _add_info_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_info = subparsers.add_parser('info',
                                        help='Print info about HDF file',
                                        description='... to be updated ...')
    comgroup = parser_info.add_argument_group('Print info')
    comgroup.add_argument(*single_input['args'], **single_input['kwargs'])
    parser_info.set_defaults(execute=_info_execute)
    return subparsers


def _info_execute(args):
    """
    :param args:
    :return:
    """
    info = implib.import_module('crplib.commands.info')
    retval = info.run_print_info(args)
    return retval


def _add_dump_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_dump = subparsers.add_parser('dump',
                                        help='Dump data from HDF file to (compressed) text',
                                        description='...to be updated...')
    comgroup = parser_dump.add_argument_group('Module I/O parameter')
    comgroup.add_argument(*single_input['args'], **single_input['kwargs'])
    comgroup.add_argument(*input_group['args'], **input_group['kwargs'])
    comgroup.add_argument(*single_output['args'], **single_output['kwargs'])

    comgroup = parser_dump.add_argument_group('Module runtime parameter: map index')
    comgroup.add_argument(*map_reference['args'], **map_reference['kwargs'])
    comgroup.add_argument('--full-blocks', '-fb', action='store_true', default=False, dest='fullblocks',
                          help='Dump full blocks, i.e., blocks contain info about both target and query assembly.')

    comgroup = parser_dump.add_argument_group('Module runtime parameter: signal track')
    comgroup.add_argument('--skip-ends', '-ske', type=int, default=0, dest='skipends',
                          help='Skip this number of bases at the beginning/end of each chromosome. Default: 0')
    comgroup.add_argument('--resolution', '-res', type=int, default=1, dest='resolution',
                          help='Summarize the signal value over bins of this size and dump those. Default: 1bp')
    comgroup.add_argument('--statistic', '-stat', type=str, choices=['mean', 'max', 'min', 'sum', 'median', 'product'],
                          dest='summstat', default='mean',
                          help='Use this function to summarize the signal values. Default: mean')

    comgroup = parser_dump.add_argument_group('Module runtime parameter: region file')
    comgroup.add_argument('--delimiter', '-delim', type=str, default="\t", dest='delimiter',
                          help='Column delimiter enclosed in double quotes. Default: "\t"')
    comgroup.add_argument('--comment-header', '-cmh', action='store_true', default=False, dest='commentheader',
                          help='If set, first character of header row will be set to "#". Default: False')
    comgroup.add_argument('--no-header', '-noh', action='store_true', default=False, dest='noheader',
                          help='Output just data and no header line. Default: False')
    comgroup.add_argument('--R-table', '-rtab', action='store_true', default=False, dest='rtable',
                          help='Print the row indices as first column w/o header to make R import straightforward.'
                               ' Default: False')

    parser_dump.set_defaults(execute=_dump_execute)
    return subparsers


def _dump_execute(args):
    """
    :param args:
    :return:
    """
    dump = implib.import_module('crplib.commands.dump')
    retval = dump.run_dump_data(args)
    return retval


def _add_convert_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_convert = subparsers.add_parser('convert',
                                           help='Convert various file formats to HDF5 format'
                                                ' for faster I/O. HDF5 is the default data'
                                                ' format of the CREEPIEST tool.',
                                           description='... to be updated ...')

    comgroup = parser_convert.add_argument_group('General parameters')

    tmp = dict(default_task)
    tmp['kwargs']['choices'] = ['signal', 'region', 'map', 'tfmotif']
    tmp['kwargs']['help'] = 'Specify task to execute: convert signal (bedGraph),' \
                            ' region (BED), map (tsv) or TF motif data.'
    comgroup.add_argument(*tmp['args'], **tmp['kwargs'])
    comgroup.add_argument(*multi_input['args'], **multi_input['kwargs'])
    comgroup.add_argument(*single_hdfout['args'], **single_hdfout['kwargs'])
    comgroup.add_argument(*output_group['args'], **output_group['kwargs'])
    comgroup.add_argument(*select_chroms['args'], **select_chroms['kwargs'])

    comgroup = parser_convert.add_argument_group('Signal conversion')

    comgroup.add_argument('--chrom-sizes', '-cs', type=str, default='', dest='chromsizes',
                          help='Full path to UCSC-style 2 column file with chromosome sizes')
    comgroup.add_argument('--no-qnorm', '-nq', action='store_true', default=False, dest='noqnorm',
                          help='Do not perform quantile normalization before merging several input files.'
                               ' This will decrease the run time. When merging several replicate'
                               ' experiments, performing quantile normalization is recommended. Default: FALSE')
    comgroup.add_argument('--merge-stat', '-ms', type=str, default='mean', dest='mergestat',
                          choices=['mean', 'median', 'max', 'min'],
                          help='Use this statistic to merge several input files: mean, median, min, max. Default: mean')
    comgroup.add_argument('--clip', '-cl', type=float, default=99.95, dest='clip',
                          help='Clip signal values above this percentile to reduce the effect of strong'
                               ' outliers. Default: 99.95')
    comgroup.add_argument('--dec-ranks', '-drk', action='store_true', default=False, dest='decranks',
                          help='Transform non-zero signal to decile ranks (0 <= 10 <= 20 etc.)')

    comgroup = parser_convert.add_argument_group('Region conversion')

    comgroup.add_argument('--name-idx', '-nix', type=int, default=-1, dest='nameidx',
                          help='Specify column index (0-based) with region names. If set to -1,'
                               ' generic names will be assigned based on genomic sort order'
                               ' (if the file does not contain a header; see --use-header). Default: -1')
    comgroup.add_argument('--score-idx', '-six', type=int, default=-1, dest='scoreidx',
                          help='Specify column index (0-based) with score to rank regions. If set to'
                               ' -1 no ranking can be performed. Assumes ranking from high to low score. Default: -1')
    comgroup.add_argument('--use-header', action='store_true', default=False, dest='useheader',
                          help='If set, the BED file to convert must have a full header (starting # will'
                               ' be ignored). The first three fields are always assumed to be "chromosome",'
                               ' "region start" and "region end", no matter what their actual name is.'
                               ' If the header includes the field "name", it will be used automatically,'
                               ' i.e. setting --name-idx appropriately is not required.')
    comgroup.add_argument('--keep-top', '-topk', type=float, default=95, dest='keeptop',
                          help='Specify top N percent of regions to keep after ranking. Requires --score-idx'
                               ' to be set to a valid column index. Default: 95')
    comgroup.add_argument('--filter-size', '-fs', type=int, default=0, dest='filtersize',
                          help='Remove regions smaller than this value. Default: 0')

    comgroup = parser_convert.add_argument_group('Map conversion')

    comgroup.add_argument('--target-assembly', '-tassm', type=str, default='', dest='target',
                          help='Specify the target [from] assembly name, e.g., hg19.'
                               ' If the name is not specified, it will be inferred'
                               ' from the chromosome sizes file (string only consisting'
                               ' of alphanumerical characters from the beginning of the filename).')
    comgroup.add_argument('--target-chrom', '-tchr', type=str, default='', dest='targetchrom',
                          help='2 column text file with chromosome sizes for the target assembly.')
    comgroup.add_argument('--query-assembly', '-qassm', type=str, default='', dest='query',
                          help='Specify the query [to] assembly name, e.g., mm9.'
                               ' If the name is not specified, it will be inferred'
                               ' from the chromosome sizes file (string only consisting'
                               ' of alphanumerical characters from the beginning of the filename).')
    comgroup.add_argument('--query-chrom', '-qchr', type=str, default='', dest='querychrom',
                          help='2 column text file with chromosome sizes for the query assembly.')
    comgroup.add_argument('--index-col', '-idxc', type=int, dest='indexcol',
                          help='Specify the column number (starting at 0) that contains the'
                               ' numerical index.')
    comgroup.add_argument('--min-map', '-min', action='store_true', default=False, dest='minmap',
                          help='Generate minimal map file - this will result in ~15% smaller'
                               ' files but comes at the cost of slightly longer computation times'
                               ' in downstream analysis.')

    comgroup = parser_convert.add_argument_group('TF motif parameter')
    comgroup.add_argument('--motif-db', '-mdb', type=str, default='', dest='motifdb')
    comgroup.add_argument('--db-format', '-dbf', type=str, choices=['meme', 'map', 'list'], dest='dbformat')
    comgroup.add_argument('--query', '-qry', action='store_true', default=False, dest='query',
                          help='Build index for query')
    parser_convert.set_defaults(execute=_convert_execute)
    return subparsers


def _convert_execute(args):
    """
    :param args:
    :return:
    """
    convert = implib.import_module('crplib.commands.convert')
    retval = convert.run_conversion(args)
    return retval


def _add_norm_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_norm = subparsers.add_parser('norm',
                                        help='Perform quantile normalization',
                                        description='...to be updated...')
    parser_norm.add_argument(*multi_input['args'], **multi_input['kwargs'])
    parser_norm.add_argument('--output-dir', '-od', type=str, dest='outdir', default='as_input',
                             help='Specify output directory for all files or "as_input"'
                                  ' to put normalized files next to the originals.'
                                  ' Default: as_input')
    parser_norm.add_argument('--replace', '-rep', type=str, dest='replace', default='.h5',
                             help='Specify what to replace in the old filename to create the'
                                  ' new one. Default: .h5')
    parser_norm.add_argument('--suffix', '-suf', type=str, dest='suffix', default='.norm.h5',
                             help='Specify the new suffix, i.e., the string that should replace'
                                  ' the part of the old filename specified above.'
                                  ' Default: .norm.h5')
    parser_norm.add_argument('--sym-link', '-sym', action='store_true', default=False, dest='symlink',
                             help='If normalization is run in batch mode on a diverse dataset, create a'
                                  ' symbolic link (using the derived output filename) for samples'
                                  ' w/o replicates that would otherwise not give a new output and thus,'
                                  ' potentially, lead to failed runs in the workflow. Default: False')
    parser_norm.add_argument(*output_group['args'], **output_group['kwargs'])
    parser_norm.set_defaults(execute=_norm_execute)
    return subparsers


def _norm_execute(args):
    """
    :param args:
    :return:
    """
    norm = implib.import_module('crplib.commands.norm')
    retval = norm.run_normalization(args)
    return retval


def _add_compfeat_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_compfeat = subparsers.add_parser('compfeat',
                                            help='Compute features for dataset',
                                            description='...to be updated...')

    comgroup = parser_compfeat.add_argument_group('General parameters')

    tmp = dict(default_task)
    tmp['kwargs']['choices'] = ['regress', 'classify', 'matched']
    tmp['kwargs']['help'] = 'Specify task to execute: compute features for regression, classification,' \
                            ' or matched classification (each positive sample is matched to a negative sample)'
    comgroup.add_argument(*tmp['args'], **tmp['kwargs'])
    comgroup.add_argument(*single_input['args'], **single_input['kwargs'])
    comgroup.add_argument(*input_group['args'], **input_group['kwargs'])
    comgroup.add_argument(*single_hdfout['args'], **single_hdfout['kwargs'])
    comgroup.add_argument(*output_group['args'], **output_group['kwargs'])
    comgroup.add_argument(*select_chroms['args'], **select_chroms['kwargs'])
    comgroup.add_argument(*hdf_mapfile['args'], **hdf_mapfile['kwargs'])
    comgroup.add_argument(*map_reference['args'], **map_reference['kwargs'])

    comgroup.add_argument('--features', '-ft', type=str, nargs='+', default=[], dest='features')
    comgroup.add_argument('--kmers', '-km', type=int, nargs='+', default=[], dest='kmers')
    comgroup.add_argument('--seq-file', '-sf', type=str, default='', dest='seqfile')
    comgroup.add_argument('--tf-motifs', '-tfm', type=str, default='', dest='tfmotifs')

    comgroup.add_argument('--signal-file', '-sigf', type=str, default=[], nargs='+', dest='sigfile')
    comgroup.add_argument('--asc-regions', '-asc', type=str, default=[], nargs='+', dest='ascregions')
    comgroup.add_argument('--roi-file', '-roif', type=str, default=[], nargs='+', dest='roifile')
    comgroup.add_argument('--roi-quant', '-roiq', type=str, default=['all'], nargs='+',
                          choices=['all', 'binary', 'counts', 'coverage'], dest='roiquant')

    comgroup = parser_compfeat.add_argument_group('Parameter for signal regression (regsig)')
    comgroup.add_argument('--num-samples', '-smp', type=int, default=20000, dest='numsamples',
                          help='Number of training samples to collect')
    comgroup.add_argument('--resolution', '-res', type=int, default=25, dest='resolution')

    comgroup = parser_compfeat.add_argument_group('Parameter for binary classification (groups)')
    comgroup.add_argument('--pos-ingroup', '-pig', type=str, default='', dest='posingroup')
    comgroup.add_argument('--neg-ingroup', '-nig', type=str, default='', dest='negingroup')
    comgroup.add_argument('--add-seq', '-ads', action='store_true', default=False, dest='addseq')
    comgroup.add_argument('--pos-outgroup', '-pog', type=str, default='', dest='posoutgroup')
    comgroup.add_argument('--neg-outgroup', '-nog', type=str, default='', dest='negoutgroup')

    comgroup.add_argument('--window', '-win', type=int, default=0, dest='window')
    comgroup.add_argument('--stepsize', '-stp', type=int, default=0, dest='stepsize')

    parser_compfeat.set_defaults(execute=_compfeat_execute)
    return subparsers


def _compfeat_execute(args):
    """
    :param args:
    :return:
    """
    compfeat = implib.import_module('crplib.commands.compfeat')
    retval = compfeat.run_compute_features(args)
    return retval


def _add_merge_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_merge = subparsers.add_parser('merge',
                                         help='Merge two datasets or add new data',
                                         description='... to be updated ...')
    comgroup = parser_merge.add_argument_group('General parameters')
    comgroup.add_argument(*multi_input['args'], **multi_input['kwargs'])
    comgroup.add_argument(*single_hdfout['args'], **single_hdfout['kwargs'])
    comgroup.add_argument(*output_group['args'], **output_group['kwargs'])

    comgroup.add_argument('--merge-on', '-mrg', type=str, default=['name'], nargs='+', dest='mergeon')
    comgroup.add_argument('--add-values', '-val', type=str, nargs='*', default=[], dest='valfile')
    comgroup.add_argument('--from-column', '-col', type=str, nargs='*', default=[], dest='valcolumn')
    comgroup.add_argument('--wg-dataset', '-wgd', action='store_true', default=False, dest='wgdataset',
                          help='The dataset to load the values from is is not split by chromosome.'
                               ' Currently, this is only supported with a single'
                               ' column name specified via "--merge-on".')
    comgroup.add_argument('--ignore-missing', '-igm', action='store_true', default=False, dest='ignoremissing',
                          help='When merging, ignore missing data (= drop those entries).')

    comgroup = parser_merge.add_argument_group('Stack datasets vertically')
    comgroup.add_argument('--just-stack', '-jst', action='store_true', default=False, dest='stack')
    comgroup.add_argument('--add-indicator', '-ind', action='store_true', default=False, dest='indicator')

    parser_merge.set_defaults(execute=_merge_execute)
    return subparsers


def _merge_execute(args):
    """
    :param args:
    :return:
    """
    merge = implib.import_module('crplib.commands.merge')
    retval = merge.run_merge_datasets(args)
    return retval


def _add_train_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_train = subparsers.add_parser('train',
                                         help='Train a model on a given training dataset',
                                         description='...to be updated...')
    comgroup = parser_train.add_argument_group('Train model')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile',
                          help='Full path to training data file')
    comgroup.add_argument('--input-group', '-ig', type=str, default='', dest='inputgroup',
                          help='Specify group root in training data file to be loaded.')
    comgroup.add_argument('--model-spec', '-ms', type=str, required=True, dest='modelspec',
                          help='Full path to model specification in JSON format.')
    comgroup.add_argument('--no-tuning', '-nt', action='store_true', default=False, dest='notuning',
                          help='Do not search for better model parameters via cross validation. Default: False')
    comgroup.add_argument('--use-features', '-uft', type=str, nargs='+', default=[], dest='usefeatures')
    mutgrp = parser_train.add_mutually_exclusive_group()
    mutgrp.add_argument('--drop-features', '-drop', type=str, default='', dest='dropfeatures')
    mutgrp.add_argument('--keep-features', '-keep', type=str, default='', dest='keepfeatures')
    comgroup.add_argument('--target-var', '-var', type=str, dest='targetvar',
                          help='Name of the dependant variable (labels/output/target) column in the dataset.')
    comgroup.add_argument('--derive-target', '-drv', type=str, dest='derivetarget',
                          help='Calculate the target variable using this rule/formula. '
                               'Individual columns in the dataset can be accessed '
                               'via "data.<COLUMN_NAME>". The derived variable will be named "target".')
    comgroup.add_argument('--balance', '-bal', action='store_true', default=False, dest='balance',
                          help='Balance classes, i.e., sample from majority class to reach number of'
                               ' samples in minority class (only downsampling supported).')
    comgroup.add_argument('--sub-sample', '-smp', type=int, default=2000, dest='subsample',
                          help='Limit number of samples per class. Default: 2000')
    comgroup.add_argument('--cv-folds', '-cv', type=int, default=10, dest='cvfolds',
                          help='Number of folds in cross validation. Default: 10')
    comgroup.add_argument('--model-output', '-mo', type=str, dest='modelout',
                          help='Specify file path to store the serialized representation of the trained model.')
    comgroup.add_argument('--metadata-output', '-do', type=str, dest='metadataout',
                          help='Specify file path to store the metadata of the trained model. If not specified,'
                               ' the model output path will be used and the file extension replaced with ".json".'
                               ' Default: <empty>')
    comgroup.add_argument('--calc-weights', '-cwt', action='store_true', default=False, dest='calcweights')
    comgroup.add_argument('--load-weights', '-swt', type=str, default='', dest='loadweights')
    comgroup.add_argument('--subset', '-sub', type=str, default='', dest='subset')
    comgroup.add_argument('--crp-metadata', '-cmd', type=str, default=None, dest='crpmetadata',
                          help='The JSON file specified for "subset" is a CREEPIEST metadata'
                               ' file. Evaluate this expression to select the sample names'
                               ' from the file. Default: None')
    parser_train.set_defaults(execute=_train_execute)
    return subparsers


def _train_execute(args):
    """
    :param args:
    :return:
    """
    train = implib.import_module('crplib.commands.train')
    retval = train.run_train_model(args)
    return retval


def _add_mapsig_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_mapsig = subparsers.add_parser('mapsig',
                                          help='Map a signal track from one assembly to another',
                                          description='...to be updated...')
    comgroup = parser_mapsig.add_argument_group('Map signal track')

    comgroup.add_argument(*single_input['args'], **single_input['kwargs'])
    comgroup.add_argument(*input_group['args'], **input_group['kwargs'])
    comgroup.add_argument(*hdf_mapfile['args'], **hdf_mapfile['kwargs'])
    comgroup.add_argument(*select_chroms['args'], **select_chroms['kwargs'])
    comgroup.add_argument('--allocate-chroms', '-ac', type=int, default=2, dest='allocate',
                          help='Number of chromosomes of the query assembly to be allocated'
                               ' in memory simultaneously. Note that each worker process still'
                               ' has to load the actual data of the target assembly to be mapped.'
                               ' If you are tight on memory, set this to 1 and also set the number'
                               ' of worker processes to 1. Default: 2')
    comgroup.add_argument(*single_hdfout['args'], **single_hdfout['kwargs'])
    comgroup.add_argument(*output_group['args'], **output_group['kwargs'])
    parser_mapsig.set_defaults(execute=_mapsig_execute)
    return subparsers


def _mapsig_execute(args):
    """
    :param args:
    :return:
    """
    mapsig = implib.import_module('crplib.commands.map_signal')
    retval = mapsig.run_map_signal(args)
    return retval


def _add_apply_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_apply = subparsers.add_parser('apply',
                                         help='Apply a trained model to a dataset',
                                         description='... to be updated ...')
    parser_apply.add_argument('--task', '-t', type=str, choices=['test', 'est'], dest='task',
                              help='Specify task', required=True)
    comgroup = parser_apply.add_argument_group('General parameters')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile',
                          help='Full path to input file in HDF5 format.')
    comgroup.add_argument('--input-group', '-ig', type=str, default='', dest='inputgroup',
                          help='Group root path for input. Default: <empty>')
    comgroup.add_argument('--output', '-o', type=str, required=True, dest='outputfile',
                          help='Full path to output file in HDF5 format.')
    comgroup.add_argument('--output-group', '-og', type=str, default='', dest='outputgroup',
                          help='Group root path for output. Default: <empty>')
    comgroup.add_argument('--model-file', '-mdf', type=str, required=True, dest='modelfile',
                          help='Specify full path to model file produced with train command.')
    comgroup.add_argument('--model-metadata', '-mdd', type=str, default='', dest='modelmetadata',
                          help='Path to JSON file with model metadata. If left empty, use the same'
                               ' path as for the model file and replace extension with ".json".'
                               ' Default: <empty>')
    comgroup.add_argument('--reduce-classes', '-red', type=list, nargs='*', default=[], dest='reduce')
    comgroup.add_argument('--subset', '-sub', type=str, default='', dest='subset')
    comgroup.add_argument('--crp-metadata', '-cmd', type=str, default=None, dest='crpmetadata',
                          help='The JSON file specified for "subset" is a CREEPIEST metadata'
                               ' file. Evaluate this expression to select the sample names'
                               ' from the file. Default: None')
    comgroup.add_argument('--num-perm', '-nump', type=int, default=0, dest='numperm',
                          help='This runs permutation tests if a value larger than'
                               ' zero is specified. The significance of the default'
                               ' CV score is described as Test 1 in'
                               ' Ojala and Garriga. Journal of ML Research (2010) vol. 11.'
                               ' Note that this can take a substantial amount of time.'
                               ' Default: 0')
    comgroup.add_argument('--cv-perm', '-cvp', type=int, default=10, dest='cvperm',
                          help='Run this number of cross-validation folds in the permutation test.')
    comgroup.add_argument('--num-rand', '-numr', type=int, default=0, dest='numrand',
                          help='This executes a simplified randomization test by just randomizing'
                               ' the output labels and applying the trained model w/o new cross-'
                               ' validation (as it is done for the permuation above). This'
                               ' is considerably faster as no new (temporary) model is trained.'
                               ' Default: 0')
    comgroup.add_argument('--scoring', '-sc', type=str, default='', dest='scoring')
    comgroup.add_argument('--seq-file', '-seq', type=str, dest='seqfile',
                          help='Full path to genomic sequence file in 2bit format.')
    comgroup.add_argument('--target-index', '-idx', type=str, default='', dest='targetindex',
                          help='Full path to target index created with "convert" command.'
                               ' Only required for task "cons". Default: <empty>')
    comgroup.add_argument('--no-smoothing', '-nosm', action='store_true', default=False, dest='nosmooth',
                          help='Do no smooth signal estimate at the end. Default: False')

    comgroup = parser_apply.add_argument_group('Deprecated parameters')
    comgroup.add_argument('--target-var', '-var', type=str, default='', dest='targetvar',
                          help='DEPRECATED: set value will be ignored - information is loaded'
                               ' from model metadata file. Exists for backward compatibility.')
    comgroup.add_argument('--derive-target', '-drv', type=str, default='', dest='derivetarget',
                          help='DEPRECATED: set value will be ignored - information is loaded'
                               ' from model metadata file. Exists for backward compatibility.')
    parser_apply.set_defaults(execute=_apply_execute)
    return subparsers


def _apply_execute(args):
    """
    :param args:
    :return:
    """
    apply = implib.import_module('crplib.commands.apply')
    retval = apply.run_apply_model(args)
    return retval


def _add_correlation_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_corr = subparsers.add_parser('corr',
                                        help='Compute pairwise similarity statistics between signal tracks.',
                                        description='... to be updated ...')
    comgroup = parser_corr.add_argument_group('Compute correlation')
    comgroup.add_argument('--task', '-t', type=str, choices=['cons', 'active', 'full', 'roi'], dest='task',
                          help='Specify task...')
    comgroup.add_argument('--measure', '-ms', type=str, choices=['pearson', 'spearman'], nargs='+',
                          required=True, dest='measure', help='Specify statistic(s) to compute')
    comgroup.add_argument('--limit-to-roi', action='store_true', default=False, dest='roilimit')
    comgroup.add_argument('--skip-size', '-sks', type=int, default=0, dest='skipsize',
                          help='Skip regions below that size. Default: 0')
    comgroup.add_argument('--roi-file', '-roi', type=str, dest='roifile', default='')
    comgroup.add_argument(*hdf_mapfile['args'], **hdf_mapfile['kwargs'])
    comgroup.add_argument('--map-reference', '-mref', type=str, choices=['target', 'query'],
                          dest='mapreference',
                          help='Specify whether to use target or query coordinates from map file.')
    comgroup.add_argument('--input-a', '-ia', type=str, required=True, dest='inputfilea')
    comgroup.add_argument('--input-group-a', '-iga', type=str, default='default', dest='inputgroupa')
    comgroup.add_argument('--input-b', '-ib', type=str, required=True, dest='inputfileb')
    comgroup.add_argument('--input-group-b', '-igb', type=str, default='default', dest='inputgroupb')
    comgroup.add_argument(*single_jsonout['args'], **single_jsonout['kwargs'])
    parser_corr.set_defaults(execute=_correlation_execute)
    return subparsers


def _correlation_execute(args):
    """
    :param args:
    :return:
    """
    corr = implib.import_module('crplib.commands.correlation')
    retval = corr.run_compute_correlation(args)
    return retval


def _add_eval_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_eval = subparsers.add_parser('eval',
                                        help='Evaluate stuff',
                                        description='...to be updated...')
    comgroup = parser_eval.add_argument_group('Evaluate things')
    comgroup.add_argument('--task', '-t', type=str, choices=['overlap'], dest='task',
                          help='Specify task')
    comgroup.add_argument('--input', '-i', type=str, nargs='+', dest='inputfile')
    comgroup.add_argument('--roi-file', '-roi', type=str, default='', dest='roifile')

    comgroup.add_argument('--output', '-o', type=str, dest='outputfile')

    parser_eval.set_defaults(execute=_eval_execute)
    return subparsers


def _eval_execute(args):
    """
    :param args:
    :return:
    """
    eval = implib.import_module('crplib.commands.eval')
    retval = eval.run_evaluation(args)
    return retval


def _add_match_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_match = subparsers.add_parser('match',
                                         help='Search for matching genomic regions based on sequence features',
                                         description='This command works on a set of genomic regions. It computes'
                                                     ' the basic sequence features requested by the user and then'
                                                     ' searches in the genomic complement (the whole genome minus'
                                                     ' the input regions) for similar regions. A common usecase is'
                                                     ' the search for a set of background regions closely matching'
                                                     ' the set of input (foreground) regions.')
    comgroup = parser_match.add_argument_group('Module I/O parameter')
    comgroup.add_argument(*single_input['args'], **single_input['kwargs'])
    comgroup.add_argument(*input_group['args'], **input_group['kwargs'])
    comgroup.add_argument(*single_hdfout['args'], **single_hdfout['kwargs'])
    comgroup.add_argument(*output_group['args'], **output_group['kwargs'])
    comgroup.add_argument(*twobit_genome['args'], **twobit_genome['kwargs'])

    comgroup = parser_match.add_argument_group('Module runtime parameter')
    comgroup.add_argument('--features', '-ft', type=str, nargs='+', default=['len', 'gc', 'rep'], dest='features')
    comgroup.add_argument('--kmers', '-km', type=int, nargs='+', default=[], dest='kmers')
    comgroup.add_argument('--timeout', '-to', type=int, default=10, dest='timeout',
                          help='Timeout in minutes for each chromosome to stop searching'
                               ' even if not all foreground regions were matched: Default: 10 min.')
    comgroup.add_argument('--allow-nomatch', '-anm', action='store_true', default=False, dest='allownomatch',
                          help='Allow that no match is found. The corresponding set of input regions will be'
                               ' discarded. Default: False')
    comgroup.add_argument('--relax-init', '-ri', type=float, default=1.0, dest='relaxinit',
                          help='Initial value for relaxed matching in percentage points. Default: 1.0')
    comgroup.add_argument('--relax-step', '-rs', type=float, default=0.5, dest='relaxstep',
                          help='Increment relaxation after each iteration by this value: Default: 0.5')
    comgroup.add_argument('--relax-limit', '-rl', type=float, default=3.0, dest='relaxlimit',
                          help='Maximum allowed relaxation, reset to initial value for next iteration. Default: 3.0')
    comgroup.add_argument('--save-pairs', '-sp', action='store_true', default=False, dest='savepairs',
                          help='Save only pairs of matched foreground - background regions instead'
                               ' of all foreground plus matched background regions. Default: False')
    parser_match.set_defaults(execute=_match_execute)
    return subparsers


def _match_execute(args):
    """
    :param args:
    :return:
    """
    regmatch = implib.import_module('crplib.commands.match')
    retval = regmatch.run_background_match(args)
    return retval
