# coding=utf-8

import os as os
import csv as csv
import gzip as gz
import re as re
import collections as col
import json as js
import sys as sys
import numpy as np

from ruffus import *


def link_encode_files(mdfile, indir, trgdir):
    """ One of those magic functions one needs to deal with metadata...
    :param mdfile:
    :param indir:
    :param trgdir:
    :return:
    """
    from preprocess.__init__ import ENCODE_BIOSAMPLES_MAP as SAMPLE_MAP
    from preprocess.__init__ import ENCODE_LAB_MAP as LAB_MAP
    from preprocess.__init__ import GENCODE_MAP as GENC_MAP

    dlfiles = os.listdir(indir)
    accnums = set([f.split('.')[0] for f in dlfiles])
    done = set()
    linked_files = []
    with open(mdfile, 'r') as infile:
        rows = csv.DictReader(infile, delimiter='\t')
        for r in rows:
            facc = r['File accession']
            if facc in done:
                continue
            ffmt = r['File format']
            if facc in accnums and ffmt in ['gtf', 'bigWig']:
                try:
                    sample = SAMPLE_MAP[r['Biosample term name']]
                    if ffmt == 'gtf':
                        gencver = GENC_MAP[facc]
                except KeyError:
                    # irrelevant sample or not using predefined GENCODE version
                    continue
                compartment = r['Biosample subcellular fraction term name']
                if compartment == 'nucleus' or compartment == 'cytosol':
                    continue
                if ffmt == 'bigWig':
                    ext = '.bigWig'
                    quant = 'signal'
                elif ffmt == 'gtf':
                    ext = '.gtf.gz'
                    quant = 'genc-' + gencver
                else:
                    raise ValueError(ffmt)
                outtype = r['Output type']
                if outtype not in ['gene quantifications', 'signal']:
                    continue
                assay = r['Assay']
                libtype = r['Library made from']
                treatment = r['Biosample treatments']
                if assay != 'RNA-seq' or libtype != 'polyadenylated mRNA' or treatment:
                    continue
                mark = 'mRNA'
                expid = r['Experiment accession']
                asmbl = r['Assembly']
                labid = LAB_MAP[r['Lab']]
                repnum = 'R0' if not r['Biological replicate(s)'] else 'R' + r['Biological replicate(s)']
                fname = '_'.join([expid, facc, asmbl, sample, mark, labid, repnum]) + '.' + quant + ext
                tpath = os.path.join(trgdir, fname)
                if os.path.islink(tpath):
                    linked_files.append(tpath)
                    continue
                else:
                    fpath = os.path.join(indir, facc + ext)
                    os.link(fpath, tpath)
                    linked_files.append(tpath)
                done.add(facc)
    return linked_files


def _parse_gtf_line(line):
    """
    :param line:
    :return:
    """
    chrom, _, entry, start, end, _, strand, _, attributes = line.split('\t')
    try:
        infos = dict([(t.split()[0], t.split()[1].strip('"')) for t in attributes.split(';') if t and len(t) > 1])
    except IndexError:
        sys.stderr.write('\n{}\n'.format(attributes))
        raise
    return chrom, int(start), int(end), strand, infos


def convert_gtf(inpath, outpath):
    """
    :param inpath:
    :param outpath:
    :return:
    """
    genes = dict()
    chrom_match = re.compile('chr[0-9]+')
    with gz.open(inpath, 'rb') as gtf:
        for line in gtf:
            line = line.decode('ascii').strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if chrom_match.match(parts[0]) is None:
                continue
            if parts[2] == 'gene':
                c, s, e, strand, infos = _parse_gtf_line(line)
                if infos['gene_type'] != 'protein_coding':
                    continue
                this_gene = {'chrom': c, 'start': s, 'end': e, 'symbol': infos['gene_name'],
                             'ensid': infos['gene_id'], 'transcripts': [], 'strand': strand}
                if strand == '+':
                    this_gene['promoter'] = (s - 1500, s + 500)
                    this_gene['body'] = (s + 500, e)
                else:
                    this_gene['promoter'] = (e - 500, e + 1500)
                    this_gene['body'] = (s, e - 500)
                genes[this_gene['ensid']] = this_gene
            elif parts[2] == 'transcript':
                c, s, e, strand, infos = _parse_gtf_line(line)
                this_trans = {'chrom': c, 'start': s, 'end': e, 'strand': strand, 'ensid': infos['transcript_id']}
                gid = infos['gene_id']
                if gid not in genes:
                    continue
                genes[gid]['transcripts'].append(this_trans)
    genes = sorted(genes.values(), key=lambda d: (int(d['chrom'].strip('chr')), d['start'], d['end']))
    with open(outpath, 'w') as outf:
        _ = js.dump(genes, outf, indent=1)
    return outpath


def merge_expquant_files(inputfiles, outputfile, genemodel):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    filter_genes = set()
    with open(genemodel) as ref:
        all_genes = js.load(ref)
        for g in all_genes:
            filter_genes.add(g['ensid'])
            filter_genes.add(g['ensid'].rsplit('.', 1)[0])

    collect_estimates = col.defaultdict(list)
    for fp in inputfiles:
        with gz.open(fp, 'rb') as gtf:
            for line in gtf:
                line = line.decode('ascii').strip()
                if not line or line.startswith('#'):
                    continue
                c, s, e, strand, infos = _parse_gtf_line(line)
                query_id = infos['gene_id'].rsplit('.', 1)[0]
                if query_id in filter_genes:
                    try:
                        collect_estimates[query_id].append(float(infos['FPKM']))
                    except KeyError:
                        collect_estimates[query_id].append(float(infos['RPKM1']) + float(infos['RPKM2']))

    all_exps = []
    with open(outputfile, 'w') as outf:
        _ = outf.write('\t'.join(['#chrom', 'start', 'end', 'ensid', 'symbol', 'expr']) + '\n')
        for gene in all_genes:
            line = '\t'.join([gene['chrom'], str(gene['start']), str(gene['end']), gene['ensid'], gene['symbol'], ''])
            try:
                vals = collect_estimates[gene['ensid'].rsplit('.', 1)[0]]
                avg = np.average(vals)
                if np.isnan(avg):
                    avg = 0.
                all_exps.append(avg)
                line += str(avg) + '\n'
                _ = outf.write(line)
            except KeyError:
                line += '0' + '\n'
                _ = outf.write(line)
    return outputfile


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

    refdir = config.get('Pipeline', 'refdir')
    tempdir = config.get('Pipeline', 'tempdir')
    bedout = config.get('Pipeline', 'bedout')
    hdfout = config.get('Pipeline', 'hdfout')
    datadir = config.get('Pipeline', 'datadir')

    enc_mdfile = config.get('Pipeline', 'encmd')

    inputfiles = [os.path.join(refdir, f) for f in os.listdir(refdir)]
    tmpf = link_encode_files(enc_mdfile, datadir, tempdir)
    inputfiles.extend(tmpf)

    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, name='init', output=inputfiles).mkdir(tempdir)

    convgenc = pipe.transform(task_func=convert_gtf,
                              name='convgenc',
                              input=output_from(init),
                              filter=formatter('gencode\.(?P<GENCVER>v[0-9a-zA-Z]+)\.+'),
                              output=os.path.join(refdir, 'gencode_{GENCVER[0]}_main.json')).jobs_limit(2)

    mergequant = pipe.collate(task_func=merge_expquant_files,
                              name='mergequant',
                              input=output_from(init),
                              # ENCSR000COQ_ENCFF000EWE_hg19_GM12878_mRNA_L13_R0.genc-v7.gtf.gz
                              filter=formatter('(?P<ENCEXP>\w+)_(?P<ENCF>\w+)_(?P<ASS>\w+)_(?P<CELL>\w+)_mRNA_(?P<LAB>\w+)_(?P<REP>\w+)\.genc-(?P<GENCVER>\w+)\.gtf\.gz'),
                              output=os.path.join(bedout, '{ENCEXP[0]}_{ASS[0]}_{CELL[0]}_mRNA_{LAB[0]}.genc-{GENCVER[0]}.bed'),
                              extras=[os.path.join(refdir, 'gencode_{GENCVER[0]}_main.json')]).mkdir(bedout).jobs_limit(1)

    return pipe
