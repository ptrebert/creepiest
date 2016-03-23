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
    subparsers = _add_convbg_command(subparsers)
    subparsers = _add_convreg_command(subparsers)
    subparsers = _add_traindata_command(subparsers)
    #subparsers = _add_regmatch_command(subparsers)
    #subparsers = _add_mkmap_command(subparsers)
    #subparsers = _add_sigmap_command(subparsers)
    #subparsers = _add_tfscan_command(subparsers)
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


def _add_convfna_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_convfna = subparsers.add_parser('convfna',
                                           help='Convert assembly in (FASTA nucleotide file) to HDF5 format',
                                           description='...to be updated...')
    comgroup = parser_convfna.add_argument_group('Convert FASTA parameters')
    comgroup.add_argument('--assembly', '-a', type=str, required=True, dest='assembly',
                          help='Specify name of assembly.')
    comgroup.add_argument('--chrom-sizes', '-s', type=str, required=True, dest='chromsizes',
                          help='Full path to UCSC-style 2 column file with chromosome sizes')
    comgroup.add_argument('--keep-chroms', '-c', type=str, default='"(chr)?[0-9]+(\s|$)"', dest='keepchroms',
                          help='Regular expression pattern (needs to be double quoted) matching'
                               ' chromosome names to keep. Default: "(chr)?[0-9]+(\s|$)" (i.e. autosomes)')
    comgroup.add_argument('--no-replace', '-norp', action='store_true', default=False, dest='noreplace',
                          help='Replace all non ACGTN letters in the sequence with N (case sensitive)')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile',
                          help='Full path to input file to be converted')
    comgroup.add_argument('--output', '-o', type=str, required=True, dest='outputfile',
                          help='Full path to output file')
    parser_convfna.set_defaults(execute=_convert_execute)
    return subparsers


def _add_convbg_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_convbg = subparsers.add_parser('convbg',
                                          help='Convert bedGraph signal tracks to HDF5. If several signal'
                                               ' tracks are specified as input, build a single merged track.',
                                          description='...to be updated...')
    comgroup = parser_convbg.add_argument_group('Convert bedGraph parameters')
    comgroup.add_argument('--chrom-sizes', '-s', type=str, required=True, dest='chromsizes',
                          help='Full path to UCSC-style 2 column file with chromosome sizes')
    comgroup.add_argument('--keep-chroms', '-c', type=str, default='"(chr)?[0-9]+(\s|$)"', dest='keepchroms',
                          help='Regular expression pattern (needs to be double quoted) matching'
                               ' chromosome names to keep. Default: "(chr)?[0-9]+(\s|$)" (i.e. autosomes)')
    comgroup.add_argument('--no-qnorm', '-nq', action='store_true', default=False, dest='noqnorm',
                          help='Do not perform quantile normalization before merging several input files.'
                               ' This will substantially decrease the run time. When merging several replicate'
                               ' experiments, performing quantile normalization is recommended. Default: FALSE')
    comgroup.add_argument('--merge-stat', '-ms', type=str, default='mean', choices=['mean', 'median', 'max', 'min'],
                          dest='mergestat',
                          help='Use this statistic to merge several input files: mean, median, min, max. Default: mean')
    comgroup.add_argument('--group-root', '-gr', type=str, default='', dest='grouproot',
                          help='Specify a root path to store the individual chromosomes in the HDF5. Default: <empty>')
    comgroup.add_argument('--clip', '-cl', type=float, default=99.95, dest='clip',
                          help='Clip signal values above this percentile. Default: 99.95')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile', nargs='+',
                          help='Full path to input file(s) to be converted.')
    comgroup.add_argument('--output', '-o', type=str, required=True, dest='outputfile',
                          help='Full path to output file')
    parser_convbg.set_defaults(execute=_convert_execute)
    return subparsers


def _add_convreg_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_convreg = subparsers.add_parser('convreg',
                                           help='Convert genomic regions (e.g. peaks) to HDF5. If serveral files'
                                                ' are specified as input, build a single merged set.',
                                           description='... to be updated ...')
    comgroup = parser_convreg.add_argument_group('Convert region parameters')
    comgroup.add_argument('--keep-chroms', '-c', type=str, default='"(chr)?[0-9]+(\s|$)"', dest='keepchroms',
                          help='Regular expression pattern (needs to be double quoted) matching'
                               ' chromosome names to keep. Default: "(chr)?[0-9]+(\s|$)" (i.e. autosomes)')
    comgroup.add_argument('--name-idx', '-n', type=int, default=-1, dest='nameidx',
                          help='Specify column index (0-based) with region names. If set to -1,'
                               ' new names will be assigned based on genomic sort order. Default: -1')
    comgroup.add_argument('--score-idx', '-s', type=int, default=-1, dest='scoreidx',
                          help='Specify column index (0-based) with score to rank regions. If set to'
                               ' -1 no ranking can be performed. Default: -1')
    comgroup.add_argument('--keep-top', '-k', type=float, default=95, dest='keeptop',
                          help='Specify top N percent of regions to keep after ranking. Requires --score-idx'
                               ' to be set. Default: 95')
    comgroup.add_argument('--group-root', '-gr', type=str, default='', dest='grouproot',
                          help='Specify a root path to store the individual chromosomes in the HDF5. Default: <empty>')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile', nargs='+',
                          help='Full path to input file(s) to be converted.')
    comgroup.add_argument('--output', '-o', type=str, required=True, dest='outputfile',
                          help='Full path to output file.')
    parser_convreg.set_defaults(execute=_convert_execute)
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
    comgroup = parser_traindata.add_argument_group('Generate training data')
    comgroup.add_argument('--task', '-tk', type=str, choices=['regsig'], dest='task')
    comgroup.add_argument('--chrom-sizes', '-s', type=str, required=True, dest='chromsizes',
                          help='Full path to UCSC-style 2 column file with chromosome sizes')
    comgroup.add_argument('--keep-chroms', '-c', type=str, default='"(chr)?[0-9]+(\s|$)"', dest='keepchroms',
                          help='Regular expression pattern (needs to be double quoted) matching'
                               ' chromosome names to keep. Default: "(chr)?[0-9]+(\s|$)" (i.e. autosomes)')
    comgroup.add_argument('--num-samples', '-smp', type=int, default=50000, dest='numsamples',
                          help='Number of training samples to collect')
    comgroup.add_argument('--resolution', '-res', type=int, default=25, dest='resolution')
    comgroup.add_argument('--two-pass', '-tp', action='store_true', default=False, dest='twopass')
    comgroup.add_argument('--features', '-ft', type=str, nargs='+', default=['prm', 'kmers', 'gc'], dest='features')
    comgroup.add_argument('--kmers', '-km', type=int, nargs='+', default=[2, 3], dest='kmers')
    comgroup.add_argument('--chain-file', '-cf', type=str, required=True, dest='chainfile')
    comgroup.add_argument('--seq-file', '-sf', type=str, required=True, dest='seqfile')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile')
    comgroup.add_argument('--input-group', '-ig', type=str, default='', dest='inputgroup')
    comgroup.add_argument('--group-root', '-gr', type=str, default='', dest='grouproot',
                          help='Specify a root path to store the individual chromosomes in the HDF5. Default: <empty>')
    comgroup.add_argument('--output', '-o', type=str, dest='outputfile')
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
                          help='Full path to chain file for mapping information from reference (target) to query')
    comgroup.add_argument('--query-group', '-qgr', type=str, required=True, dest='querygroup',
                          help='Specify group root path for mapped signal')
    comgroup.add_argument('--query-chroms', '-qch', type=str, required=True, dest='querychroms',
                          help='Full path to file with chromosome sizes of the query assembly.')
    comgroup.add_argument('--input', '-i', type=str, required=True, dest='inputfile',
                          help='Full path to input file in HDF5 format.')
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
                          help='Print a short description of the computable features to stdout.')
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


def _add_regmatch_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_regmatch = subparsers.add_parser('regmatch',
                                            help='Search for matching genomic regions based on sequence features',
                                            description='This command works on a set of genomic regions. It computes'
                                                        ' the basic sequence features requested by the user and then'
                                                        ' searches in the genomic complement (the whole genome minus'
                                                        ' the input regions) for similar regions. A common usecase is'
                                                        ' the search for a set of background regions closely matching'
                                                        ' the set of input (foreground) regions.')
    mutex_group = parser_regmatch.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--crpk-source', '-crpkf', dest='crpksource', nargs='+',
                             help='List full paths to individual CRPK files or give a list of paths to folders'
                                  ' containing CRPK files.')
    mutex_group.add_argument('--bed-source', '-bed', dest='bedsource', nargs='+',
                             help='List full path to individual BED files or give a list of paths to folders'
                                  ' containing the input BED files.')
    parser_regmatch.add_argument('--genome', '-gen', dest='genome', type=str, required=True,
                                 help='The genomic sequence as one FASTA file (can be gzipped)')
    parser_regmatch.add_argument('--chrom-sizes', '-chrs', dest='chromsizes', type=str, required=True,
                                 help='Path to file containing chromosome sizes (as provided by UCSC) compatible'
                                      ' to the FASTA sequence file (e.g. concerning chromosome names)')
    parser_regmatch.add_argument('--filter-ambig', '-fltn', dest='filtern', type=float, default=5.0,
                                 help='Filter out sequences with too many ambiguous positions (letter N).'
                                      ' Default is 5 percent.')
    parser_regmatch.add_argument('--features', '-feat', dest='features', nargs='+', required=True,
                                 help='List the features that should be used to select a match for an input region.')
    parser_regmatch.add_argument('--dist-chr', '-dchr', dest='distchrom', action='store_true', default=False,
                                 help='Set if the chromosomal distribution of the input regions should be matched, i.e.'
                                 ' a match is only searched for on the same chromosome. This can increase the runtime'
                                 ' quite considerably.')
    parser_regmatch.add_argument('--relaxation', '-rel', dest='relax', type=float, default=2.,
                                 help='Initial relaxation as percentage points (default: 2.0), i.e. a value of 50'
                                 ' percent will be matched by any value between 48 and 52 whereas an absolute'
                                 ' value of 50 will be matched by any value roughly between 49 and 51.')
    parser_regmatch.add_argument('--increment', '-incr', dest='increment', type=float, default=1.0,
                                 help='Increase the relaxation by this value in each iteration. Default is 1.0')
    parser_regmatch.add_argument('--limit-relax', '-lim', dest='limitrelax', type=float, default=5.0,
                                 help='Limit the relaxation to that value. In the next iteration, it will be set'
                                 ' to its initial value again.')
    parser_regmatch.add_argument('--run-limit', '-rl', dest='runlimit', type=float, default=0.,
                                 help='Limit the runtime to this number of hours. Default runtime is unlimited.')
    parser_regmatch.add_argument('--no-clobber', '-noclb', dest='noclobber', action='store_true', default=False,
                                 help='If set, make a backup copy of the input file instead of replacing it.')
    parser_regmatch.set_defaults(execute=_regmatch_execute)
    return subparsers


def _regmatch_execute(args):
    """
    :param args:
    :return:
    """
    regmatch = implib.import_module('lib.commands.regmatch')
    retval = regmatch.run_region_matching(args)
    return retval


def _add_mkmap_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_mkmap = subparsers.add_parser('mkmap',
                                         help='Create a CREEP for signal mapping based on a pairwise alignment'
                                              ' in UCSC axt format or UCSC chain file.')
    parser_mkmap.add_argument('--server-cfg', '-srv', dest='srvconfig', required=True, type=str,
                              help='Specify the full path to the CREEPIEST annotation server configuration file.')
    parser_mkmap.add_argument('--reference', '-r', dest='reference', type=str, required=True,
                              help='Specify reference assembly, e.g. hg19 or mm10 (also: primary assembly)')
    parser_mkmap.add_argument('--target', '-t', dest='target', type=str, required=True,
                              help='Specify target assembly, e.g. mm10 or hg19 (also: aligning assembly)')
    parser_mkmap.add_argument('--chrom', '-c', dest='chrom', type=str, default='.*',
                              help='Specify a regular expression to restrict the output to chromosomes'
                                   ' matching the expression. This filtering applies to both reference'
                                   ' and target species.')
    grp = parser_mkmap.add_mutually_exclusive_group(required=True)
    grp.add_argument('--axt-input', '-ai', dest='axtinput', type=str,
                     help='Specify full path to input file in axt format')
    grp.add_argument('--chain-input', '-ci', dest='chaininput', type=str,
                     help='Specify full path to chain file (reference to target)')
    parser_mkmap.add_argument('--gap', '-g', dest='gap', type=int, default=25,
                              help='Specify the acceptable gap length between two consecutive alignment blocks.'
                                   ' That length should be smaller than the (average) region size in the signal'
                                   ' tracks to be mapped.')
    parser_mkmap.add_argument('--outdir', '-o', dest='outdir', type=str, required=True,
                              help='Specify path to output directory (file name will be assembled automatically).')
    parser_mkmap.set_defaults(execute=_mkmap_execute)
    return subparsers


def _mkmap_execute(args):
    """
    :param args:
    :return:
    """
    mapbed = implib.import_module('lib.commands.mkmap')
    retval = mapbed.run_mkmap(args)
    return retval


def _add_sigmap_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_sigmap = subparsers.add_parser('sigmap',
                                          help='Take an arbitrary number of signal CREEPs and one mapping CREEP'
                                               ' with compatible reference and map all signal tracks from the'
                                               ' reference to the target assembly.')
    parser_sigmap.add_argument('--map-creep', '-crmp', dest='mapcreep', required=True, type=str,
                               help='Specify full path to map CREEP with compatible resolution and correct'
                                    ' reference assembly.')
    parser_sigmap.add_argument('--sig-creep', '-crsg', dest='sigcreeps', required=True, nargs='+', type=str,
                               help='Specify full path to single or multiple signal CREEPs (as space separated list)'
                                    ' or provide folder path(s) containing an arbitrary number of signal CREEPs.')
    parser_sigmap.add_argument('--self-map', '-self', dest='selfmap', default=False, action='store_true',
                               help='If the map CREEP represents a self alignment of an assembly, this option'
                                    ' must be set to ensure proper mapping of the signal values. This option'
                                    ' assumes that the map CREEP only contains non-trivial alignment regions'
                                    ' for the self alignment.')
    parser_sigmap.set_defaults(execute=_sigmap_execute)
    return subparsers


def _sigmap_execute(args):
    """
    :param args:
    :return:
    """
    sigmap = implib.import_module('lib.commands.sigmap')
    retval = sigmap.run_sigmap(args)
    return retval


def _add_tfscan_command(subparsers):
    """
    :param subparsers:
    :return:
    """
    parser_tfscan = subparsers.add_parser('tfscan',
                                          help='Scan a set of sequences for transcription factor binding sites using'
                                               ' the MEME suite (Fimo). The Fimo executable must be available in PATH.')
    parser_tfscan.add_argument('--genome', '-g', dest='genome', required=True, type=str,
                               help='Identifier of genome assembly, e.g. hg19, used for file/folder naming.')
    parser_tfscan.add_argument('--source', '-s', dest='source', required=True, type=str,
                               help='Folder with sequence files (FASTA format) of genome assembly.')
    parser_tfscan.add_argument('--tf-motifs', '-t', dest='motiffile', required=True, type=str,
                               help='Full path to TF motif file in MEME format containing PWMs and character'
                                    ' frequencies for A, C, G and T.')
    parser_tfscan.add_argument('--match-files', '-mf', dest='matchfile', default='all', type=str,
                               help='Regular expression to match sequence files; defaults to match all files'
                                    ' with .fa/.fasta file extension. Enclose expression with double quotes > " <')
    parser_tfscan.add_argument('--output', '-o', dest='output', type=str, required=True,
                               help='Base output folder, sub-folder <genome> will be created and'
                                    'output files stored; base directory will be created if it does not exist.')
    parser_tfscan.add_argument('--pvalue', '-p', dest='pthresh', type=str, default='1e-5',
                               help='Threshold for p-value to report motif match, can be in scientific'
                                    '>e< notation (default: 1e-5).')
    parser_tfscan.set_defaults(execute=_tfscan_execute)
    return subparsers


def _tfscan_execute(args):
    """
    :param args:
    :return:
    """
    tfscan = implib.import_module('lib.commands.tfscan')
    retval = tfscan.run_tfscan(args)
    return retval

