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

    tempdir = config.get('Pipeline', 'tempdir')
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
    # step 1: select only autosomes from chrom sizes file
    filtchrom = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='filtchrom',
                               input=output_from(init),
                               filter=suffix('.chrom.sizes'),
                               output='.selected',
                               output_dir=tempdir,
                               extras=[cmd, jobcall]).mkdir(tempdir)

    cmd = config.get('Pipeline', 'conv2bit')
    # step 2: convert 2bit files into gzipped fasta files
    conv2bit = pipe.collate(task_func=sci_obj.get_jobf('inref_out'),
                            name='conv2bit',
                            input=output_from(init, filtchrom),
                            filter=formatter('(?P<ASSEMBLY>\w+)\.(2bit|selected)'),
                            output=os.path.join(tempdir, '{ASSEMBLY[0]}.fa'),
                            extras=[cmd, '.selected', jobcall])

    # step 3: run FIMO on genomes
    cmd = config.get('Pipeline', 'tfscan')
    base_out = os.path.join(tempdir, 'tfscan')
    tfscan = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                            name='tfscan',
                            input=output_from(conv2bit),
                            filter=formatter('(?P<ASSEMBLY>\w+)\.fa'),
                            output=os.path.join(base_out, '{ASSEMBLY[0]}', 'fimo*'),
                            extras=[os.path.join(base_out, '{ASSEMBLY[0]}'), 'fimo*', cmd, jobcall]).mkdir(os.path.join(tempdir, 'tfscan'))

    return pipe