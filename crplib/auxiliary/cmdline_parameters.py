# coding=utf-8

"""
This is a convenience module to centrally define command line
parameters including help messages for commonly used options
"""

# file parameters

single_input = {'args': ('--input', '-i'),
                'kwargs': {'required': True,
                           'dest': 'inputfile', 'type': str,
                           'help': 'Specify full path to single input file, e.g., '
                                   '/long/path/to/file.ext'}}

multi_input = {'args': ('--input', '-i'),
               'kwargs': {'nargs': '+', 'required': True,
                          'dest': 'inputfiles', 'type': str,
                          'help': 'Specify full path to one or more input files, e.g., '
                                  '/long/path/to/fileA.ext /longer/path/to/fileB.ext'}}

single_hdfout = {'args': ('--output', '-o'),
                 'kwargs': {'required': True,
                            'dest': 'outputfile', 'type': str,
                            'help': 'Specify full path to HDF output file, e.g., '
                                    '/long/path/to/file.h5 - note that the path will '
                                    'be created if it does not exist.'}}

single_jsonout = {'args': ('--output', '-o'),
                  'kwargs': {'required': True,
                             'dest': 'outputfile', 'type': str,
                             'help': 'Specify full path to JSON output file, e.g., '
                                     '/long/path/to/file.json - note that the path will '
                                     'be created if it does not exist.'}}

hdf_indexfile = {'args': ('--index-file', '-idx'),
                 'kwargs': {'dest': 'indexfile', 'type': str,
                            'help': 'Specify full path to HDF index file, e.g., '
                                    '/long/path/to/index_file.h5'}}

twobit_genome = {'args': ('--2bit-genome', '-2bg'),
                 'kwargs': {'dest': '2bitgenome', 'type': str,
                            'help': 'Specify full path to 2bit genome file, e.g., '
                                    '/long/path/to/genome_file.2bit'}}

# HDF group parameters

input_group = {'args': ('--input-group', '-ig'),
               'kwargs': {'dest': 'inputgroup', 'type': str,
                          'default': '',
                          'help': 'Specify root path to data groups in input file, e.g., '
                                  '/group/root/to/data - note that the group information '
                                  'can also be prepended to the file name separated by a colon, e.g., '
                                  '/group/root/to/data:/long/path/to/file.h5. Default group: <empty>'}}

output_group = {'args': ('--output-group', '-og'),
                'kwargs': {'dest': 'outputgroup', 'type': str,
                           'default': '',
                           'help': 'Specify root path to data groups in output file, e.g., '
                                   '/group/root/to/data - note that the group information '
                                   'can also be prepended to the file name separated by a colon, e.g., '
                                   '/group/root/to/data:/long/path/to/file.h5. The data will be stored'
                                   ' splitted by individual chromosomes. Default group: <empty>'}}

# misc parameters

filter_chroms = {'args': ('--select-chroms', '-slc'),
                 'kwargs': {'dest': 'selectchroms', 'type': str,
                            'default': '"(chr)?[0-9][0-9A-Z]?(\s|$)"',
                            'help': 'Specify regular expression to select chromosomes by name. '
                                    'The regular expression needs to be double-quoted. '
                                    'Default: autosomes (i.e.: "(chr)?[0-9][0-9A-Z]?(\s|$)"'}}
