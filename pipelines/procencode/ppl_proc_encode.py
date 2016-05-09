# coding=utf-8

import os as os
import fnmatch as fnm

from ruffus import *


def make_input_pairs(datadir, chaindir, targets, queries):
    """
    :param datadir:
    :param chaindir:
    :return:
    """
    datafiles = os.listdir(datadir)
    datafiles = fnm.filter(datafiles, '*srcsig.h5')
    chainfiles = os.listdir(chaindir)
    chainfiles = fnm.filter(chainfiles, '*.chain.gz')
    targets = targets.split()
    queries = queries.split()
    pairs = []
    for df in datafiles:
        df_path = os.path.join(datadir, df)
        assembly = df.split('_')[1]
        if assembly not in targets:
            continue
        for chf in chainfiles:
            if chf.startswith(assembly):
                query = chf.split('.')[0].split('_')[2]
                if query in queries:
                    pairs.append([df_path, os.path.join(chaindir, chf)])
    return pairs


def make_corr_pairs(inpairs, targets, datadir, outdir, cmd, jobcall):
    """
    :param inpairs:
    :param targets:
    :param datadir:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    targets = targets.split()
    mates = []
    for root, dirs, files in os.walk(datadir):
        if files:
            basedir, subdir = os.path.split(root)
            if all([t in subdir for t in targets]):
                for f in files:
                    mates.append(os.path.join(root, f))
    infiles = [ip[0] for ip in inpairs]
    params = []
    # assayed: ENCDS059CRP_mm9_MEL_H3K27ac_L08.srcsig.h5
    # mapped/estimated: ENCDS006CRP_mm9_K562_H3K27ac_L00.trg.hg19.mapsig.h5
    for inf in infiles:
        _, qry, qcell, lib, qlab = inf.split('.')[0].split('_')
        if lib.startswith('IgG') or lib in ['H3K9me1']:
            continue
        matched_files = fnm.filter(mates, '*_' + qry + '_*' + '_' + lib + '_*')
        assert matched_files, 'No matches for input file {}'.format(inf)
        for mf in matched_files:
            _, mfname = os.path.split(mf)
            _, qry, tcell, lib, tlab = mfname.split('.')[0].split('_')
            _, _, trg, datatype, _ = mfname.split('.')
            tmp = cmd.format(**{'query': qry, 'target': trg})
            outname = 'corr_{}_{}_{}_{}_vs_{}_{}.trg.{}.{}.json'.format(qry, lib, qcell, qlab, tcell, tlab, trg, datatype)
            outpath = os.path.join(outdir, outname)
            params.append([[inf, mf], outpath, tmp, jobcall])
    return params


def build_pipeline(args, config, sci_obj):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """

    pipe = Pipeline(name=config.get('Pipeline', 'name'))

    tempdir = config.get('Pipeline', 'tempdir')
    outdir = config.get('Pipeline', 'outdir')
    chaindir = config.get('Pipeline', 'chaindir')
    datadir = config.get('Pipeline', 'datadir')

    inputpairs = make_input_pairs(datadir, chaindir, config.get('Pipeline', 'targets'),
                                  config.get('Pipeline', 'queries'))

    sig_re = '(?P<DSID>ENCDS[0-9]+CRP)_(?P<ASSEMBLY>[a-zA-Z0-9]+)_(?P<CELL>[0-9A-Za-z]+)'\
             '_(?P<DTYPE>[0-9A-Za-z]+)_(?P<LAB>L[0-9]+)\.srcsig\.h5'

    chain_re = '(?P<TARGET>\w+)_to_(?P<QUERY>\w+)\.rbest\.chain\.gz'


    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, output=inputpairs)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'mapsig').replace('\n', ' ')  # deal with multiline value
    # step 1: map signal tracks
    stepdir_map = os.path.join(outdir, 'output', 'mappings', 'signal')
    mapsig = pipe.transform(task_func=sci_obj.get_jobf('inpair_out'),
                            name='mapsig',
                            input=output_from(init),
                            filter=formatter(sig_re, chain_re),
                            output=os.path.join(stepdir_map, '{QUERY[1]}_from_{TARGET[1]}', '{DSID[0]}_{QUERY[1]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.trg.{TARGET[1]}.mapsig.h5'),
                            extras=[cmd, jobcall]).mkdir(stepdir_map)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'tdatsig').replace('\n', ' ')
    stepdir = os.path.join(outdir, 'traindata', 'regression', 'signal')
    tdatsig = pipe.transform(task_func=sci_obj.get_jobf('inpair_out'),
                             name='tdatsig',
                             input=output_from(init),
                             filter=formatter(sig_re, chain_re),
                             output=os.path.join(stepdir, '{ASSEMBLY[0]}_to_{QUERY[1]}', '{DSID[0]}_{ASSEMBLY[0]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.qry.{QUERY[1]}.tdat.sig.h5'),
                             extras=[cmd, jobcall]).mkdir(stepdir)

    tdatsig_re = '(?P<DSID>ENCDS[0-9]+CRP)_(?P<ASSEMBLY>[a-zA-Z0-9]+)_(?P<CELL>[0-9A-Za-z]+)'\
                 '_(?P<DTYPE>[0-9A-Za-z]+)_(?P<LAB>L[0-9]+)\.qry\.(?P<QUERY>[a-zA-Z0-9]+)\.tdat\.sig\.h5'

    stepdir = os.path.join(outdir, 'models', 'regression', 'signal')
    cmd = config.get('Pipeline', 'trainsig').replace('\n', ' ')
    trainsig = pipe.subdivide(task_func=sci_obj.get_jobf('in_outpair'),
                              name='trainsig',
                              input=output_from(tdatsig),
                              filter=formatter(tdatsig_re),
                              output=[os.path.join(stepdir, '{ASSEMBLY[0]}_to_{QUERY[0]}', '{DSID[0]}_{ASSEMBLY[0]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.qry.{QUERY[0]}.rfreg.sig.pck'),
                                      os.path.join(stepdir, '{ASSEMBLY[0]}_to_{QUERY[0]}', '{DSID[0]}_{ASSEMBLY[0]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.qry.{QUERY[0]}.rfreg.sig.json')],
                              extras=[cmd, jobcall]).mkdir(stepdir)

    trainsig_re = '(?P<DSID>ENCDS[0-9]+CRP)_(?P<TARGET>[a-zA-Z0-9]+)_(?P<CELL>[0-9A-Za-z]+)'\
                  '_(?P<DTYPE>[0-9A-Za-z]+)_(?P<LAB>L[0-9]+)\.qry\.(?P<QUERY>[a-zA-Z0-9]+)\.rfreg\.sig\.pck'

    cmd = config.get('Pipeline', 'estsig').replace('\n', ' ')
    stepdir_est = os.path.join(outdir, 'output', 'estimates', 'signal')
    reffile = os.path.join(stepdir_map, '{QUERY[0]}_from_{TARGET[0]}', '{DSID[0]}_{QUERY[0]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.trg.{TARGET[0]}.mapsig.h5')
    estsig = pipe.transform(task_func=sci_obj.get_jobf('in_out_ref'),
                            name='estsig',
                            input=output_from(trainsig),
                            filter=formatter(trainsig_re),
                            output=os.path.join(stepdir_est, '{QUERY[0]}_from_{TARGET[0]}', '{DSID[0]}_{QUERY[0]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.trg.{TARGET[0]}.estsig.h5'),
                            extras=[reffile, cmd, jobcall]).mkdir(stepdir_est).follows(mapsig)

    # make_corr_pairs(inpairs, targets, datadir, outdir, cmd, jobcall):

    cmd = config.get('Pipeline', 'corrmap').replace('\n', ' ')
    stepdir = os.path.join(outdir, 'eval', 'signal', 'correlation', 'mapping')
    corrmap_params = make_corr_pairs(inputpairs, config.get('Pipeline', 'targets'), stepdir_map, stepdir, cmd, jobcall)
    corrmap = pipe.files(sci_obj.get_jobf('inpair_out'),
                         corrmap_params,
                         name='corrmap').mkdir(stepdir).follows(mapsig)

    cmd = config.get('Pipeline', 'correst').replace('\n', ' ')
    stepdir = os.path.join(outdir, 'eval', 'signal', 'correlation', 'estimate')
    correst_params = make_corr_pairs(inputpairs, config.get('Pipeline', 'targets'), stepdir_est, stepdir, cmd, jobcall)
    correst = pipe.files(sci_obj.get_jobf('inpair_out'),
                         correst_params,
                         name='correst').mkdir(stepdir).follows(estsig)

    cmd = config.get('Pipeline', 'runall')
    runall = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                        name='runall',
                        input=output_from(corrmap, correst),
                        output=os.path.join(outdir, 'runall_proc_encode.chk'),
                        extras=[cmd, jobcall])

    return pipe
