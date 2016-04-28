# coding=utf-8

"""
This module defines all sub-parsers for the different commands
and contains the appropriate execute functions
"""

import importlib as implib


def add_sub_parsers(main_parser):
    """
    :param main_parser:
    :return:
    """
    subparsers = main_parser.add_subparsers(dest='subparser_name', title='Subcommands')
    subparsers = _add_tests_command(subparsers)
    subparsers = _add_info_command(subparsers)
    subparsers = _add_convert_command(subparsers)
    subparsers = _add_traindata_command(subparsers)
    subparsers = _add_train_command(subparsers)
    subparsers = _add_mapsig_command(subparsers)
    subparsers = _add_apply_command(subparsers)
    subparsers = _add_correlation_command(subparsers)
    subparsers = _add_match_command(subparsers)
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
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile')
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


def _add_convert_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_convert = subparsers.add_parser('convert',
                                           help='Convert bedGraph or BED files to HDF5 format.',
                                           description='... to be updated ...')
    comgroup = parser_convert.add_argument_group('General parameters')
    comgroup.add_argument('--task', '-tk', type=str, required=True, choices=['signal', 'region', 'chain', 'tfmotif'], dest='task',
                          help='Specify task to execute, convert signal (bedGraph), region (BED-like) or chain files.')
    comgroup.add_argument('--keep-chroms', '-kc', type=str, default='"(chr)?[0-9]+(\s|$)"', dest='keepchroms',
                          help='Regular expression pattern (needs to be double quoted) matching'
                               ' chromosome names to keep in the target. Default: "(chr)?[0-9]+(\s|$)" (i.e. autosomes)')
    comgroup.add_argument('--query-check', '-qc', type=str, default='"(chr)?[0-9]+(\s|$)"', dest='qcheck',
                          help='Regular expression pattern to check for valid chromosome names in the query'
                               ' assembly. Default: "(chr)?[0-9]+(\s|$)" (i.e. autosomes)')
    comgroup.add_argument('--input', '-i', type=str, nargs='+', required=True, dest='inputfile',
                          help='Full path to one or more files to convert.')
    comgroup.add_argument('--output', '-o', type=str, required=True, dest='outputfile',
                          help='Full path to output file to be created')
    comgroup.add_argument('--output-group', '-og', type=str, default='', dest='outputgroup',
                          help='Group prefix to store individual chromosomes. Default: <empty>')
    comgroup = parser_convert.add_argument_group('bedGraph parameters')
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
    comgroup = parser_convert.add_argument_group('BED parameters')
    comgroup.add_argument('--name-idx', '-nix', type=int, default=-1, dest='nameidx',
                          help='Specify column index (0-based) with region names. If set to -1,'
                               ' new names will be assigned based on genomic sort order. Default: -1')
    comgroup.add_argument('--score-idx', '-six', type=int, default=-1, dest='scoreidx',
                          help='Specify column index (0-based) with score to rank regions. If set to'
                               ' -1 no ranking can be performed. Default: -1')
    comgroup.add_argument('--keep-top', '-topk', type=float, default=95, dest='keeptop',
                          help='Specify top N percent of regions to keep after ranking. Requires --score-idx'
                               ' to be set to a valid column index. Default: 95')
    comgroup = parser_convert.add_argument_group('TF motif parameter')
    comgroup.add_argument('--motif-db', '-mdb', type=str, default='', dest='motifdb')
    comgroup.add_argument('--db-format', '-dbf', type=str, choices=['meme', 'map', 'list'], dest='dbformat')

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


def _add_traindata_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_traindata = subparsers.add_parser('traindata',
                                             help='Generate training data',
                                             description='...to be updated...')
    parser_traindata.add_argument('--task', '-tk', type=str, choices=['regsig', 'clsreg'], dest='task')
    comgroup = parser_traindata.add_argument_group('General parameters')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile')
    comgroup.add_argument('--output', '-o', type=str, dest='outputfile')
    comgroup.add_argument('--target-index', '-idx', type=str, default='', dest='targetindex',
                          help='Full path to target index created with "convert" command.'
                               ' Only required for task "cons". Default: <empty>')
    comgroup.add_argument('--keep-chroms', '-c', type=str, default='"(chr)?[0-9]+(\s|$)"', dest='keepchroms',
                          help='Regular expression pattern (needs to be double quoted) matching'
                               ' chromosome names to keep. Default: "(chr)?[0-9]+(\s|$)" (i.e. autosomes)')
    comgroup.add_argument('--features', '-ft', type=str, nargs='+', default=[], dest='features')
    comgroup.add_argument('--kmers', '-km', type=int, nargs='+', default=[], dest='kmers')
    comgroup.add_argument('--seq-file', '-sf', type=str, default='', dest='seqfile')

    comgroup = parser_traindata.add_argument_group('Parameter for signal regression (regsig)')
    comgroup.add_argument('--num-samples', '-smp', type=int, default=20000, dest='numsamples',
                          help='Number of training samples to collect')
    comgroup.add_argument('--resolution', '-res', type=int, default=25, dest='resolution')
    comgroup.add_argument('--input-group', '-ig', type=str, default='', dest='inputgroup')
    comgroup.add_argument('--output-group', '-og', type=str, default='', dest='outputgroup',
                          help='Specify a root path to store the individual chromosomes in the HDF5. Default: <empty>')

    comgroup = parser_traindata.add_argument_group('Parameter for region classification (clsreg)')
    comgroup.add_argument('--pos-ingroup', '-pig', type=str, default='', dest='posingroup')
    comgroup.add_argument('--neg-ingroup', '-nig', type=str, default='', dest='negingroup')
    comgroup.add_argument('--add-seq', '-ads', action='store_true', default=False, dest='addseq')
    comgroup.add_argument('--pos-outgroup', '-pog', type=str, default='', dest='posoutgroup')
    comgroup.add_argument('--neg-outgroup', '-nog', type=str, default='', dest='negoutgroup')

    parser_traindata.set_defaults(execute=_traindata_execute)
    return subparsers


def _traindata_execute(args):
    """
    :param args:
    :return:
    """
    traindata = implib.import_module('crplib.commands.traindata')
    retval = traindata.run_collect_traindata(args)
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
    comgroup.add_argument('--train-data', '-td', type=str, required=True, dest='traindata',
                          help='Full path to training data file')
    comgroup.add_argument('--train-group', '-tg', type=str, default='', dest='traingroup',
                          help='Specify group root in training data file to be loaded.')
    comgroup.add_argument('--model-spec', '-ms', type=str, required=True, dest='modelspec',
                          help='Full path to model specification in JSON format.')
    comgroup.add_argument('--no-tuning', '-nt', action='store_true', default=False, dest='notuning',
                          help='Do not search for better model parameters via cross validation. Default: False')
    comgroup.add_argument('--cv-folds', '-fld', type=int, default=10, dest='cvfolds',
                          help='Number of folds in cross validation. Default: 10')
    comgroup.add_argument('--model-output', '-mo', type=str, dest='modelout',
                          help='Specify file path to store the serialized representation of the trained model.')
    comgroup.add_argument('--metadata-output', '-do', type=str, dest='metadataout',
                          help='Specify file path to store the metadata of the trained model. If not specified,'
                               ' the model output path will be used and the file extension replaced with ".json".'
                               ' Default: <empty>')
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
    comgroup.add_argument('--chain-file', '-ch', type=str, required=True, dest='chainfile',
                          help='Full path to chain file with pairwise alignment information between'
                               ' reference (in liftOver parlance: target) and query.')
    comgroup.add_argument('--query-chroms', '-qch', type=str, required=True, dest='querychroms',
                          help='Full path to file with chromosome sizes of the query assembly.')
    comgroup.add_argument('--keep-chroms', '-c', type=str, default='"(chr)?[0-9]+(\s|$)"', dest='keepchroms',
                          help='Regular expression pattern (needs to be double quoted) matching'
                               ' chromosome names to keep. Default: "(chr)?[0-9]+(\s|$)" (i.e. autosomes)')
    comgroup.add_argument('--allocate-chroms', '-ac', type=int, default=2, dest='allocate',
                          help='Number of chromosomes to be hold in memory simultaneously. Default: 2')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile',
                          help='Full path to input file in HDF5 format.')
    comgroup.add_argument('--input-group', '-ig', type=str, default='', dest='inputgroup',
                          help='Group root in input file for signal to map. Default: <empty>')
    comgroup.add_argument('--output', '-o', type=str, dest='outputfile',
                          help='Full path to output file.')
    comgroup.add_argument('--output-group', '-og', type=str, default='', dest='outputgroup',
                          help='Specify group root path for mapped signal in output file.')
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


def _add_compfeat_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_compfeat = subparsers.add_parser('compfeat',
                                            help='Compute various features for genomic regions',
                                            description='... to be updated ...')
    comgroup = parser_compfeat.add_argument_group('Compute genomic features')

    comgroup.add_argument('--desc-feat', '-df', action='store_true', default=False, dest='descfeat',
                          help='Print a short description of the computable features to stdout and exit.')
    comgroup.add_argument('--features', '-ft', type=str, nargs='+', dest='features',
                          help='State list of features to be computed')
    comgroup.add_argument('--kmers', '-km', type=int, nargs='+', default=[2, 3, 4], dest='kmers',
                          help='Values of k for k-mer frequency. Default: 2, 3, 4')
    comgroup.add_argument('--group-root', '-gr', type=str, dest='grouproot', required=True,
                          help='Specify the group of regions in the HDF5 for the feature computation.'
                               ' All subgroups will be considered automatically.')
    comgroup.add_argument('--add-seq', '-ads', type=str, default='', dest='addseq',
                          help='If a valid path to a 2bit file is specified, the genomic sequence will'
                               ' be added to the regions with the key "seq". This is required for'
                               ' sequence-derived features. Default: <empty>')
    comgroup.add_argument('--tfmotifdb', '-tfdb', type=str, default='', dest='tfmotifdb',
                          help='Path to a folder with TF motif matches in separate files'
                               ' per chromosome (naming: ASSEMBLY_CHROM.ext). Default: <empty>')
    comgroup.add_argument('--signal', '-sg', type=str, default='', dest='signal',
                          help='Path to a signal track in HDF5 format to compute signal'
                               ' features for the genomic regions. Default: <empty>')
    comgroup.add_argument('--pwaln', '-pwa', type=str, default='', dest='pwaln',
                          help='Path to file containing information about pairwise alignable regions'
                               ' between reference and query assembly. These regions are reciprocal best.')
    comgroup.add_argument('--input', '-i', type=str, dest='inputfile',
                          help='Full path to input file.')
    parser_compfeat.set_defaults(execute=_compfeat_execute)
    return subparsers


def _compfeat_execute(args):
    """
    :param args:
    :return:
    """
    compfeat = implib.import_module('crplib.commands.comp_feat')
    retval = compfeat.run_compute_features(args)
    return retval


def _add_apply_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_apply = subparsers.add_parser('apply',
                                         help='Apply a trained model to a dataset',
                                         description='... to be updated ...')
    parser_apply.add_argument('--task', '-t', type=str, choices=['estsig'], dest='task',
                              help='Specify task')
    comgroup = parser_apply.add_argument_group('General parameters')
    comgroup.add_argument('--model-file', '-mdf', type=str, required=True, dest='modelfile',
                          help='Specify full path to model file produced with train command.')
    comgroup.add_argument('--model-metadata', type=str, default='', dest='modelmetadata',
                          help='Path to JSON file with model metadata. If left empty, use the same'
                               ' path as for the model file and replace extension with ".json".'
                               ' Default: <empty>')
    comgroup.add_argument('--seq-file', '-seq', type=str, required=True, dest='seqfile',
                          help='Full path to genomic sequence file in 2bit format.')
    comgroup.add_argument('--chain-file', '-chf', type=str, required=True, dest='chainfile',
                          help='Full path to liftOver chain file with reciprocal best chains'
                               ' between target (from/reference) and query (to) assembly.')
    comgroup.add_argument('--no-smoothing', '-nosm', action='store_true', default=False, dest='nosmooth',
                          help='Do no smooth signal estimate at the end. Default: False')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile',
                          help='Full path to input file in HDF5 format.')
    comgroup.add_argument('--input-group', '-ig', type=str, default='', dest='inputgroup',
                          help='Group root path for input. Default: <empty>')
    comgroup.add_argument('--output', '-o', type=str, required=True, dest='outputfile',
                          help='Full path to output file in HDF5 format.')
    comgroup.add_argument('--output-group', '-og', type=str, default='', dest='outputgroup',
                          help='Group root path for output. Default: <empty>')
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
    comgroup.add_argument('--measure', '-ms', type=str, choices=['pearson', 'spearman', 'r2'], nargs='+',
                          required=True, dest='measure', help='Specify statistic(s) to compute')
    comgroup.add_argument('--target-index', '-idx', type=str, default='', dest='targetindex',
                          help='Full path to target index created with "convert" command.'
                               ' Only required for task "cons". Default: <empty>')
    comgroup.add_argument('--index-group', '-idxg', type=str, default='', dest='indexgroup',
                          help='Group root path for target index. Subpaths to conservation mask'
                               ' and split indices will be appended automatically. Default: <empty>')
    comgroup.add_argument('--input-a', '-ia', type=str, required=True, dest='inputfilea')
    comgroup.add_argument('--input-group-a', '-iga', type=str, default='', dest='inputgroupa')
    comgroup.add_argument('--input-b', '-ib', type=str, required=True, dest='inputfileb')
    comgroup.add_argument('--input-group-b', '-igb', type=str, default='', dest='inputgroupb')

    comgroup.add_argument('--output', '-o', type=str, required=True, dest='outputfile',
                          help='Full path to output file in JSON format.')
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
    comgroup = parser_match.add_argument_group('Parameters for background search')
    comgroup.add_argument('--seq-file', '-seq', type=str, required=True, dest='seqfile',
                          help='Full path to genomic sequence file in 2bit format.')
    comgroup.add_argument('--features', '-ft', type=str, nargs='+', default=['len', 'gc', 'rep'], dest='features')
    comgroup.add_argument('--kmers', '-km', type=int, nargs='+', default=[], dest='kmers')
    comgroup.add_argument('--timeout', '-to', type=int, default=10, dest='timeout',
                          help='Timeout in minutes for each chromosome to stop searching'
                               ' even if not all foreground regions were matched: Default: 10 min.')
    comgroup.add_argument('--relax-init', '-ri', type=float, default=1.0, dest='relaxinit',
                          help='Initial value for relaxed matching in percentage points. Default: 1.0')
    comgroup.add_argument('--relax-step', '-rs', type=float, default=0.5, dest='relaxstep',
                          help='Increment relaxation after each iteration by this value: Default: 0.5')
    comgroup.add_argument('--relax-limit', '-rl', type=float, default=3.0, dest='relaxlimit',
                          help='Maximum allowed relaxation, reset to initial value for next iteration. Default: 3.0')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile',
                          help='Full path to input file in HDF5 format.')
    comgroup.add_argument('--input-group', '-ig', type=str, default='', dest='inputgroup',
                          help='Group root path for input. Default: <empty>')
    comgroup.add_argument('--output', '-o', type=str, required=True, dest='outputfile',
                          help='Full path to output file in HDF5 format.')
    comgroup.add_argument('--output-group', '-og', type=str, default='', dest='outputgroup',
                          help='Group root path for output. Default: <empty>')
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
