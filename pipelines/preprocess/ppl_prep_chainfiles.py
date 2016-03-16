# coding=utf-8

import os as os
import fnmatch as fnm

from ruffus import *


def link_input_files(indir, trgdir):
    """
    :param indir:
    :param trgdir:
    :return:
    """
    if not os.path.isdir(trgdir):
        os.makedirs(trgdir, exist_ok=True)
    all_files = os.listdir(indir)
    all_chains = fnm.filter(all_files, '*chain.gz')
    retfiles = []
    for chain in all_chains:
        assemblies = chain.split('.')[0]
        # no need to check, this will throw if >To< is not unique
        target, query = assemblies.split('To')
        # is there a decapitalize?
        target = target[0].lower() + target[1:]
        query = query[0].lower() + query[1:]
        new_name = target + '_to_' + query + '.over.chain.gz'
        new_path = os.path.join(trgdir, new_name)
        if os.path.isfile(new_path):
            retfiles.append(new_path)
        else:
            os.link(os.path.join(indir, chain), new_path)
            retfiles.append(new_path)
    return retfiles


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

    indir = config.get('Pipeline', 'indir')
    tempdir = config.get('Pipeline', 'tempdir')
    outdir = config.get('Pipeline', 'outdir')

    inputfiles = link_input_files(indir, tempdir)

    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, output=inputfiles)

    chain_re = '(?P<TARGET>\w+)_to_(?P<QUERY>\w+)(?P<EXT>\.[\w\.]+)'

    cmd = config.get('Pipeline', 'swap')
    swap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                          name='swap',
                          input=output_from(init),
                          filter=formatter(chain_re),
                          output=os.path.join(tempdir, '{QUERY[0]}_to_{TARGET[0]}.tbest.chain.gz'),
                          extras=[cmd, jobcall]).mkdir(tempdir)

    cmd = config.get('Pipeline', 'prenet')
    prenet = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='prenet',
                            input=output_from(swap),
                            filter=formatter(chain_re),
                            output=os.path.join(tempdir, '{TARGET[0]}_to_{QUERY[0]}.rbest.net.gz'),
                            extras=[cmd, jobcall]).mkdir(tempdir)

    cmd = config.get('Pipeline', 'netchain')
    netchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='netchain',
                              input=output_from(prenet),
                              filter=formatter(chain_re),
                              output=os.path.join(outdir, '{TARGET[0]}_to_{QUERY[0]}.rbest.chain.gz'),
                              extras=[cmd, jobcall]).mkdir(outdir)

    cmd = config.get('Pipeline', 'reswap')
    reswap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='reswap',
                            input=output_from(netchain),
                            filter=formatter(chain_re),
                            output=os.path.join(outdir, '{QUERY[0]}_to_{TARGET[0]}.rbest.chain.gz'),
                            extras=[cmd, jobcall]).mkdir(outdir)



    return pipe