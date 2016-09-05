# coding=utf-8

import os as os
import fnmatch as fnm
import re as re

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


def get_chromosome_names(assm, folder):
    """
    :param assm:
    :param folder:
    :return:
    """
    sizefiles = os.listdir(folder)
    name_match = re.compile('chr[0-9]+$')
    chromnames = []
    for szf in sizefiles:
        if szf == (assm + '.chrom.sizes'):
            with open(os.path.join(folder, szf), 'r') as infile:
                for line in infile:
                    if not line.strip():
                        continue
                    cols = line.strip().split()
                    if name_match.match(cols[0]) is not None:
                        chromnames.append(cols[0])
    assert chromnames, 'No chromosome names loaded for assembly {}'.format(assm)
    return ','.join(chromnames)


def make_filter_command(chainfiles, chainre, chromdir, outdir, cmd, jobcall):
    """
    :param chainfiles:
    :param chainre:
    :param chromdir:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    regexp = re.compile(chainre)
    arglist = []
    for chf in chainfiles:
        mobj = regexp.search(chf)
        if mobj is None:
            raise ValueError('Cannot match target/query in file: {}'.format(chf))
        target = mobj.group('TARGET')
        query = mobj.group('QUERY')
        target_chroms = get_chromosome_names(target, chromdir)
        query_chroms = get_chromosome_names(query, chromdir)
        tmp = cmd.format(**{'targetchroms': target_chroms, 'querychroms': query_chroms})
        outfile = os.path.basename(chf).replace('.over.', '.filt.')
        outpath = os.path.join(outdir, outfile)
        arglist.append([chf, outpath, tmp, jobcall])
    assert arglist, 'No call created for chain filtering'
    return arglist


def make_qfilter_command(indir, chainre, chromdir, outbase, targets, cmd, jobcall):
    """
    :param chainfiles:
    :param chainre:
    :param chromdir:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    regexp = re.compile(chainre)
    chainfiles = os.listdir(indir)
    chainfiles = [os.path.join(indir, f) for f in chainfiles if f.endswith('.chain.gz')]
    arglist = []
    for chf in chainfiles:
        mobj = regexp.search(chf)
        if mobj is None:
            raise ValueError('Cannot match target/query in file: {}'.format(chf))
        target = mobj.group('TARGET')
        if target not in targets:
            continue
        query = mobj.group('QUERY')
        query_chroms = get_chromosome_names(query, chromdir)
        for qc in sorted(query_chroms.split(',')):
            tmp = cmd.format(**{'chrom': qc})
            outfile = os.path.basename(chf).replace('.chain.gz', '.qchain.{}.tsv.gz'.format(qc))
            outdir = os.path.join(outbase, '{}_to_{}'.format(target, query))
            os.makedirs(outdir, exist_ok=True)
            outpath = os.path.join(outdir, outfile)
            arglist.append([chf, outpath, tmp, jobcall])
    assert arglist, 'No call created for chain filtering'
    return arglist


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
    outnet = config.get('Pipeline', 'outnet')
    outchain = config.get('Pipeline', 'outchain')
    outmap = config.get('Pipeline', 'outmap')
    outidx = config.get('Pipeline', 'outidx')
    refdir = config.get('Pipeline', 'refdir')
    targets = config.get('Pipeline', 'targets').split()

    inputfiles = link_input_files(indir, tempdir)

    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, output=inputfiles)

    chain_re = '(?P<TARGET>\w+)_to_(?P<QUERY>\w+)(?P<EXT>\.[\w\.]+)'

    cmd = config.get('Pipeline', 'filter')
    filt_params = make_filter_command(inputfiles, chain_re, refdir, os.path.join(tempdir, 'filter'), cmd, jobcall)
    filt = pipe.files(sci_obj.get_jobf('in_out'),
                      filt_params,
                      name='filt').mkdir(os.path.join(tempdir, 'filter')).follows(init)

    cmd = config.get('Pipeline', 'swap')
    swap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                          name='swap',
                          input=output_from(filt),
                          filter=formatter(chain_re),
                          output=os.path.join(tempdir, 'swap', '{QUERY[0]}_to_{TARGET[0]}.tbest.chain.gz'),
                          extras=[cmd, jobcall]).mkdir(os.path.join(tempdir, 'swap'))

    cmd = config.get('Pipeline', 'qrybnet')
    qrybnet = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='qrybnet',
                             input=output_from(swap),
                             filter=formatter(chain_re),
                             output=os.path.join(outnet, '{TARGET[0]}_to_{QUERY[0]}.rbest.net.gz'),
                             extras=[cmd, jobcall]).mkdir(outnet)

    cmd = config.get('Pipeline', 'qrybchain')
    qrybchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='qrybchain',
                               input=output_from(qrybnet),
                               filter=formatter(chain_re),
                               output=os.path.join(outchain, '{TARGET[0]}_to_{QUERY[0]}.rbest.chain.gz'),
                               extras=[cmd, jobcall]).mkdir(outchain)

    cmd = config.get('Pipeline', 'trgbchain')
    trgbchain = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='trgbchain',
                               input=output_from(qrybchain),
                               filter=formatter(chain_re),
                               output=os.path.join(outchain, '{QUERY[0]}_to_{TARGET[0]}.rbest.chain.gz'),
                               extras=[cmd, jobcall]).mkdir(outchain)

    cmd = config.get('Pipeline', 'trgbnet')
    trgbnet = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='trgbnet',
                             input=output_from(trgbchain),
                             filter=formatter(chain_re),
                             output=os.path.join(outnet, '{TARGET[0]}_to_{QUERY[0]}.rbest.net.gz'),
                             extras=[cmd, jobcall]).mkdir(outnet)

    cmd = config.get('Pipeline', 'qsymm')
    qsymm_params = make_qfilter_command(outchain, chain_re, refdir, os.path.join(tempdir, 'qsymm'), targets, cmd, jobcall)
    qsymm = pipe.files(sci_obj.get_jobf('in_out'),
                       qsymm_params,
                       name='qsymm').mkdir(os.path.join(tempdir, 'qsymm')).follows(trgbchain)

    cmd = config.get('Pipeline', 'mergeqsymm')
    mergeqsymm = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                              name='mergeqsymm',
                              input=output_from(qsymm),
                              filter=formatter('(?P<FILEID>[0-9_A-Za-z\.]+)\.chr[0-9]+\.tsv\.gz$'),
                              output=os.path.join('{path[0]}', '{FILEID[0]}.wg.tsv.gz'),
                              extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'symmext')
    symmext = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='symmext',
                             input=output_from(mergeqsymm),
                             filter=formatter('(?P<TARGET>\w+)_to_(?P<QUERY>\w+)\.rbest\.qchain\.wg\.tsv\.gz'),
                             output=os.path.join('{path[0]}', '{TARGET[0]}_to_{QUERY[0]}.rbest.mapext.tsv.gz'),
                             extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sortmap')
    sortmap = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='sortmap',
                             input=output_from(symmext),
                             filter=suffix('mapext.tsv.gz'),
                             output='mapext.sort.tsv.gz',
                             output_dir=outmap,
                             extras=[cmd, jobcall]).mkdir(outmap)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    map_re = '(?P<TARGET>\w+)_to_(?P<QUERY>\w+)(?P<EXT>\.[\w\.]+)'

    cmd = config.get('Pipeline', 'trgidx')
    trgidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='trgidx',
                            input=output_from(sortmap),
                            filter=formatter(map_re),
                            output=os.path.join(outidx, '{TARGET[0]}_to_{QUERY[0]}.rbest.mapext.trgidx.h5'),
                            extras=[cmd, jobcall]).mkdir(outidx)

    cmd = config.get('Pipeline', 'qryidx')
    qryidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='qryidx',
                            input=output_from(sortmap),
                            filter=formatter(map_re),
                            output=os.path.join(outidx, '{TARGET[0]}_to_{QUERY[0]}.rbest.mapext.qryidx.h5'),
                            extras=[cmd, jobcall]).mkdir(outidx)

    cmd = config.get('Pipeline', 'runall')
    runall = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                        name='runall',
                        input=output_from(trgidx, qryidx, sortmap, trgbnet, qrybnet),
                        output=os.path.join(tempdir, 'runall_prep_chain.chk'),
                        extras=[cmd, jobcall])

    return pipe