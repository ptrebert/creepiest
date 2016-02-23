# coding=utf-8

import os as os
import csv as csv
import json as js

from ruffus import *


def link_encode_files(mdfile, jsfile, indir, trgdir):
    """
    :param mdfile:
    :param indir:
    :param trgdir:
    :return:
    """
    from preprocess.__init__ import ENCODE_BIOSAMPLES_MAP as SAMPLE_MAP
    from preprocess.__init__ import ENCODE_LAB_MAP as LAB_MAP
    # map: all ENCODE files belonging to the same sample (key)
    smp_exp_map = js.load(open(jsfile, 'r'))

    exp_smp_map = dict()
    for smpid, vals in smp_exp_map.items():
        for dpbid, encid in vals:  # list of lists: DeepBlue_ExpID, ENCODE file acc
            if encid in exp_smp_map:
                raise ValueError('ENCODE file acc already in mapping: {} - {}'.format(encid, smpid))
            exp_smp_map[encid] = smpid

    dlfiles = os.listdir(indir)
    accnums = set([f.split('.')[0] for f in dlfiles])
    linked_files = []
    with open(mdfile, 'r') as infile:
        rows = csv.DictReader(infile, delimiter='\t')
        for r in rows:
            if r['File accession'] in accnums and r['File format'] == 'bigWig':
                acc = r['File accession']
                try:
                    smpid = exp_smp_map[acc]
                except KeyError:
                    print(r)
                    raise
                asmbl = r['Assembly']
                labid = LAB_MAP[r['Lab']]
                sample = SAMPLE_MAP[r['Biosample term name']]
                mark = r['Experiment target'].rsplit('-', 1)[0]
                if 'IgG' in mark:
                    mark = 'IgGControl'
                elif mark == 'Control':
                    mark = 'InpControl'
                elif mark.startswith('H3K'):
                    pass
                else:
                    continue  # ignore DNaseI and other for now
                repnum = 'R0' if not r['Biological replicate(s)'] else 'R' + r['Biological replicate(s)']
                fname = '_'.join([acc, smpid, asmbl, sample, mark, labid, repnum]) + '.bigWig'
                tpath = os.path.join(trgdir, fname)
                if os.path.islink(tpath):
                    linked_files.append(tpath)
                    continue
                else:
                    fpath = os.path.join(indir, acc + '.bigWig')
                    os.link(fpath, tpath)
                    linked_files.append(tpath)
    return linked_files


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
    refdir = config.get('Pipeline', 'refdir')
    datadir = config.get('Pipeline', 'datadir')

    enc_mdfile = config.get('Pipeline', 'encmd')
    enc_dsets = config.get('Pipeline', 'encds')

    inputfiles = []
    tmpf = os.listdir(refdir)
    inputfiles.extend([os.path.join(refdir, tf) for tf in tmpf])
    tmpf = link_encode_files(enc_mdfile, enc_dsets, datadir, tempdir)
    inputfiles.extend(tmpf)

    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, output=inputfiles)

    cmd = config.get('Pipeline', 'convbw')
    # step 1: select only autosomes from chrom sizes file
    convbw = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='convbw',
                            input=output_from(init),
                            filter=suffix('.bigWig'),
                            output='.bg.gz',
                            output_dir=tempdir,
                            extras=[cmd, jobcall]).mkdir(tempdir)


    # step 3: convert fasta files into HDF5 files
    #convert_fasta = ruf.transform()

    return pipe