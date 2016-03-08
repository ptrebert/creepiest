# coding=utf-8

import os as os

from ruffus import *


def build_pipeline(args, config, sci_obj):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """

    pipe = Pipeline(name=config.get('Pipeline', 'name'))

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    genome_temp = config.get('Pipeline', 'gentemp')
    chrom_temp = config.get('Pipeline', 'chromtemp')
    outdir = config.get('Pipeline', 'outdir')

    inputfiles = []
    for path in config.get('Pipeline', 'inputfiles').split(':'):
        allentries = os.listdir(path)
        for entry in allentries:
            fpath = os.path.join(path, entry)
            if os.path.isfile(fpath):
                inputfiles.append(fpath)

    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, output=inputfiles)

    cmd = config.get('Pipeline', 'filtchrom')
    regexp = '(?P<ASSEMBLY>\w+)\.chrom\.sizes'
    # step 1: select only autosomes from chrom sizes file and create single dummy files
    filtchrom = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                               name='filtchrom',
                               input=output_from(init),
                               filter=formatter(regexp),
                               output=os.path.join(chrom_temp, '{ASSEMBLY[0]}_chr*'),
                               extras=[chrom_temp, '{ASSEMBLY[0]}_chr*', cmd, jobcall]).mkdir(chrom_temp)

    cmd = config.get('Pipeline', 'conv2bit')
    # step 2: convert 2bit files into gzipped fasta files
    conv2bit = pipe.transform(task_func=sci_obj.get_jobf('inref_out'),
                              name='conv2bit',
                              input=output_from(filtchrom),
                              filter=formatter('(?P<ASSEMBLY>\w+)_(?P<CHROM>chr[0-9]+)\.sg\.txt'),
                              output=os.path.join(genome_temp, '{ASSEMBLY[0]}_{CHROM[0]}.fa'),
                              add_inputs=add_inputs(os.path.join(outdir, '{ASSEMBLY[0]}.2bit')),
                              extras=[cmd, '.2bit', jobcall])

    # step 3: run FIMO on genomes
    cmd = config.get('Pipeline', 'tfscan')
    tfscan = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='tfscan',
                            input=output_from(conv2bit),
                            filter=suffix('.fa'),
                            output='.tsv.gz',
                            extras=[cmd, jobcall]).mkdir(os.path.join(genome_temp, 'tfscan'))

    return pipe