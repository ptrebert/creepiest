# coding=utf-8

import os as os
import csv as csv
import collections as col
import json as js

from ruffus import *


def link_encode_files(mdfile, jsfile, indir, trgdir):
    """ One of those magic functions one needs to deal with metadata...
    :param mdfile:
    :param indir:
    :param trgdir:
    :return:
    """
    from preprocess.__init__ import ENCODE_BIOSAMPLES_MAP as SAMPLE_MAP
    from preprocess.__init__ import ENCODE_LAB_MAP as LAB_MAP
    # map: ENCODE experiment ids (ENCSR...) mapped to custom made
    # ENCODE dataset ids (ENCDS...); see Jupyter Notebook make_ENCODE_ids.ipynb
    map_encsr_encds = js.load(open(jsfile, 'r'))

    dlfiles = os.listdir(indir)
    accnums = set([f.split('.')[0] for f in dlfiles])
    linked_files = []
    dscounter = col.Counter()
    with open(mdfile, 'r') as infile:
        rows = csv.DictReader(infile, delimiter='\t')
        for r in rows:
            if r['File accession'] in accnums and r['File format'] == 'bigWig':
                facc = r['File accession']
                expid = r['Experiment accession']
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
                elif not mark.strip() and r['Assay'] == 'DNase-seq':
                    mark = 'DNaseI'
                else:
                    continue  # ignore others for now
                try:
                    dsid = map_encsr_encds[expid]
                except KeyError:
                    if r['Assay'] == 'ChIP-seq':
                        # this should not happen for histone samples
                        if mark.startswith('H3K'):
                            print('Error?!')
                            print(r)
                        continue
                    dsid = 'ENCDS000CRP'
                dscounter[dsid] += 1
                repnum = 'R0' if not r['Biological replicate(s)'] else 'R' + r['Biological replicate(s)']
                fname = '_'.join([dsid, facc, asmbl, sample, mark, labid, repnum]) + '.bigWig'
                tpath = os.path.join(trgdir, fname)
                if os.path.islink(tpath):
                    linked_files.append(tpath)
                    continue
                else:
                    fpath = os.path.join(indir, facc + '.bigWig')
                    os.link(fpath, tpath)
                    linked_files.append(tpath)
    # remove controls for non-histone ChIPseq etc and datasets with only
    # single file [exception: DNaseI]; apparently, a problem with matching
    # controls to experiments
    filtered_files = []
    for lf in linked_files:
        dsid = os.path.basename(lf).split('_')[0]
        if dscounter[dsid] > 1:
            filtered_files.append(lf)
        else:
            os.unlink(lf)

    return filtered_files


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

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'convbg')
    # ENCDS062CRP_ENCFF001MLR_mm9_CH12_InpControl_L05_R1.bg.gz
    regexp = '(?P<DSID>ENCDS[0-9]+CRP)_(?P<FFID>ENCFF[0-9A-Z]+)_(?P<ASSEMBLY>[a-zA-Z0-9]+)_' \
             '(?P<CELL>[0-9A-Za-z]+)_(?P<DTYPE>[0-9A-Za-z]+)_(?P<LAB>L[0-9]+)_(?P<REP>R[0-9])\.bg\.gz'
    # step 2: convert fasta files into HDF5 files
    convbg = pipe.collate(task_func=sci_obj.get_jobf('ins_out_ref'),
                          name='convbg',
                          input=output_from(convbw),
                          filter=formatter(regexp),
                          output=os.path.join(outdir, '{DSID[0]}_{ASSEMBLY[0]}_{CELL[0]}_{DTYPE[0]}_{LAB[0]}.h5'),
                          extras=[os.path.join(refdir, '{ASSEMBLY[0]}.chrom.sizes'), cmd, jobcall])

    return pipe