# coding=utf-8

import os as os
import csv as csv
import fnmatch as fnm
import operator as op
import collections as col
import json as js

from ruffus import *


def collect_full_paths(rootdir, pattern):
    """
    :param rootdir:
    :param pattern:
    :return:
    """
    all_files = []
    for root, dirs, files in os.walk(rootdir):
        if files:
            filt = fnm.filter(files, pattern)
            for f in filt:
                all_files.append(os.path.join(root, f))
    return all_files


def link_encode_files(metadata, indir, outdir, ffmt):
    """
    :param indir:
    :param outdir:
    :param metadata:
    :return:
    """
    encode_files = collect_full_paths(indir, '*.' + ffmt)
    fp_map = dict((os.path.basename(fp.split('.')[0]), fp) for fp in encode_files)

    linked = []
    getvals = op.itemgetter(*('Experiment accession', 'File accession', 'Assembly',
                              'Biosample term name', 'Experiment target', 'Lab'))
    with open(metadata, 'r', newline='') as mdfile:
        rows = csv.DictReader(mdfile, delimiter='\t')
        for r in rows:
            if r['File format'] != ffmt:
                continue
            if r['File accession'] in fp_map:
                if r['Biosample term name'] == 'liver' and r['Biosample life stage'] == 'em':
                    # manual filtering here since there is not a full epigenome available
                    # for that combination
                    continue
                target_name = '_'.join(list(getvals(r)))
                rep = r['Biological replicate(s)']
                if rep == 'n/a':
                    rep = '0'
                else:
                    try:
                        rep = str(int(rep))
                    except ValueError:
                        if rep == '1, 2':  # some weird info for a mouse ChIP
                            rep = '0'
                        else:
                            print(r)
                            raise
                target_name += '_' + 'b' + rep + '.bw'
                target_path = os.path.join(outdir, target_name)
                if not os.path.islink(target_path):
                    srcpath = fp_map[r['File accession']]
                    os.symlink(srcpath, target_path)
                linked.append(target_path)
    return linked


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
    refdir = config.get('Pipeline', 'refdir')
    datadir = config.get('Pipeline', 'datadir')

    enc_mdfile = config.get('Pipeline', 'encmd')

    inputfiles = []
    tmpf = os.listdir(refdir)
    inputfiles.extend([os.path.join(refdir, tf) for tf in tmpf])
    tmpf = link_encode_files(enc_mdfile, datadir, tempdir, 'bigWig')
    inputfiles.extend(tmpf)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, output=inputfiles, name='init')

    cmd = config.get('Pipeline', 'convbw')
    # convert bigWig signal to bedGraph, select only autosomes
    convbw = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='convbw',
                            input=output_from(init),
                            filter=suffix('.bw'),
                            output='.bg.gz',
                            output_dir=tempdir,
                            extras=[cmd, jobcall]).mkdir(tempdir)

    cmd = config.get('Pipeline', 'convbb')
    convbb = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='convbb',
                            input=output_from(init),
                            filter=suffix('.bigBed'),
                            output='.bed.gz',
                            output_dir=tempdir,
                            extras=[cmd, jobcall]).mkdir(tempdir)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # ENCSR857MYS_ENCFF001ZHE_mm9_ESE14_H3K9me3_RHPSU_b1.bg.gz
    regexp = '(?P<EXPID>[A-Z0-9]+)_(?P<FFID>[0-9A-Z]+)_(?P<ASSEMBLY>[a-zA-Z0-9]+)_' \
             '(?P<CELL>[0-9A-Za-z]+)_(?P<LIB>[0-9A-Za-z]+)_(?P<LAB>[A-Z]+)_(?P<REP>b[0-9])\.bg\.gz'
    cmd = config.get('Pipeline', 'convbg')
    convbg = pipe.collate(task_func=sci_obj.get_jobf('ins_out_ref'),
                          name='convbg',
                          input=output_from(convbw),
                          filter=formatter(regexp),
                          output=os.path.join(outdir, '{EXPID[0]}_{ASSEMBLY[0]}_{CELL[0]}_{LIB[0]}_{LAB[0]}.srcsig.h5'),
                          extras=[os.path.join(refdir, '{ASSEMBLY[0]}.chrom.sizes'), cmd, jobcall]).mkdir(outdir)

    # ====
    # last update changed behavior up to here; everything below is potentially outdated
    # ====

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # ENCDS062CRP_ENCFF001MLR_mm9_CH12_InpControl_L05_R1.broad/narrow.bed.gz
    regexp = '(?P<DSID>ENCDS[0-9]+CRP)_(?P<FFID>ENCFF[0-9A-Z]+)_(?P<ASSEMBLY>[a-zA-Z0-9]+)_' \
             '(?P<CELL>[0-9A-Za-z]+)_(?P<DTYPE>[0-9A-Za-z]+)_(?P<LAB>L[0-9]+)_(?P<REP>R[0-9])\.(broad|narrow)\.bed\.gz'
    cmd = config.get('Pipeline', 'convbed')
    convbed = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                           name='convbed',
                           input=output_from(convbb),
                           filter=formatter(regexp),
                           output=os.path.join(outdir, '{DSID[0]}_{ASSEMBLY[0]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.srcpk.h5'),
                           extras=[cmd, jobcall]).mkdir(outdir)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # ENCDS059CRP_mm9_MEL_H3K27ac_L08.srcpk.h5
    regexp = '(?P<DSID>ENCDS[0-9]+CRP)_(?P<ASSEMBLY>[a-zA-Z0-9]+)_' \
             '(?P<CELL>[0-9A-Za-z]+)_(?P<DTYPE>H3K[0-9A-Za-z]+)_(?P<LAB>L[0-9]+)\.srcpk\.h5'
    cmd = config.get('Pipeline', 'matchbed').replace('\n', ' ')
    matchbed = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='matchbed',
                              input=output_from(convbed),
                              filter=formatter(regexp),
                              output=os.path.join(outdir, '{DSID[0]}_{ASSEMBLY[0]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.fgbg_pk.h5'),
                              extras=[cmd, jobcall])


    cmd = config.get('Pipeline', 'runall')
    runall = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                        name='runall',
                        input=output_from(convbg, convbed),
                        output=os.path.join(tempdir, 'runall_prep_encode.chk'),
                        extras=[cmd, jobcall]).jobs_limit(1)

    return pipe
