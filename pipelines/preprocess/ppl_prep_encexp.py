# coding=utf-8

import os as os
import csv as csv
import gzip as gz
import re as re
import io as io
import collections as col
import fnmatch as fnm
import json as js
import sys as sys
import numpy as np

from ruffus import *

# This is from GENCODE, based on the description of the transcriptome fasta files:
# Protein-coding transcript sequences
# Transcript biotypes: protein_coding, nonsense_mediated_decay,
# non_stop_decay, IG_*_gene, TR_*_gene, polymorphic_pseudogene
ACCEPT_TRANSCRIPT = {'protein_coding', 'nonsense_mediated_decay', 'non_stop_decay', 'polymorphic_pseudogene'}
# plus the following: IG_*_gene, TR_*_gene


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


def link_fastq_files(metadata, indir, outdir):
    """
    :param indir:
    :param outdir:
    :param metadata:
    :return:
    """
    from preprocess.__init__ import ENCODE_LAB_MAP as LAB_MAP

    all_files = os.listdir(indir)
    bam_files = fnm.filter(all_files, '*.fastq.gz')
    fp_map = dict((fn, os.path.join(indir, fn)) for fn in bam_files)

    linked = []
    with open(metadata, 'r') as mdfile:
        rows = csv.DictReader(mdfile, delimiter='\t')
        for r in rows:
            if not bool(int(r['use'])):
                continue
            fn = os.path.basename(r['fastq'])
            encff = fn.split('.')[0]
            old_path = fp_map[fn]
            new_name = '_'.join([r['experiment'], encff, r['assembly'], r['cell'],
                                 'mRNA', LAB_MAP[r['lab'].strip('"')], 'R' + str(r['rep'])])
            new_name += '.se.fastq.gz'
            new_path = os.path.join(outdir, new_name)
            if not os.path.islink(new_path):
                os.link(old_path, new_path)
            linked.append(new_path)
    return linked


def convert_gtf(inpath, outpath):
    """
    :param inpath:
    :param outpath:
    :return:
    """
    genes = dict()
    chrom_match = re.compile('chr[0-9XY]+')
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
                biotype = infos['gene_type']
                if biotype in ACCEPT_TRANSCRIPT:
                    pass
                elif biotype.startswith('IG_') and biotype.endswith('_gene'):
                    pass
                elif biotype.startswith('TR_') and biotype.endswith('_gene'):
                    pass
                else:
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
    assert genes, 'No genes selected from GTF file'
    genes = sorted(genes.values(), key=lambda d: (d['chrom'], d['start'], d['end']))
    with open(outpath, 'w') as outf:
        _ = js.dump(genes, outf, indent=1)
    return outpath


def create_gtf_white_black_list(inputfile, outputfiles):
    """ Create a white- and blacklist from an input GTF file following
    as much as possible the definitions used by GENCODE for gene model FASTA files
    see: http://www.gencodegenes.org/releases/
    :param inputfile:
    :param outputfiles:
    :return:
    """
    white = io.StringIO()
    black = io.StringIO()
    chrom_re = re.compile("chr[0-9XY]+$")
    trans_re = re.compile("\stranscript_type\s\"(?P<TRANSTYPE>\w+)\";\s")
    with gz.open(inputfile, 'r') as gtf:

        for line in gtf:
            line = line.decode('utf-8')
            if line.startswith('#'):
                white.write(line)
                black.write(line)
                continue
            elif not line:
                continue
            else:
                cols = line.strip().split('\t')
                if cols[2] == 'transcript':
                    if chrom_re.match(cols[0]) is None:
                        black.write(line)
                        continue
                    else:
                        mobj = trans_re.search(cols[8])
                        assert mobj is not None, 'Could not identify transcript type in line: {}'.format(line)
                        ttype = mobj.group('TRANSTYPE').strip()
                        if ttype in ACCEPT_TRANSCRIPT:
                            white.write(line)
                        elif ttype.startswith('IG_') and ttype.endswith('_gene'):
                            white.write(line)
                        elif ttype.startswith('TR_') and ttype.endswith('_gene'):
                            white.write(line)
                        else:
                            black.write(line)
    with open(outputfiles[0], 'w') as out:
        _ = out.write(white.getvalue())
    with open(outputfiles[1], 'w') as out:
        _ = out.write(black.getvalue())
    return outputfiles


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


def dump_gene_regions(inputfile, outputfile, regtype):
    """
    :param inputfile:
    :param outputfile:
    :param regtype:
    :return:
    """
    with open(inputfile, 'r') as infile:
        genes = js.load(infile)
    outbuf = io.StringIO()
    # NB: by construction, genes are sorted
    # by genomic coordinate
    for gene in genes:
        start, end = gene[regtype]
        _ = outbuf.write('\t'.join([gene['chrom'], str(start), str(end), gene['ensid']]) + '\n')
    with open(outputfile, 'w') as outfile:
        _ = outfile.write(outbuf.getvalue())
    return outputfile


def make_musbed_params(quantdir, fastqdir, outdir, gencodefile):
    """
    :param quantdir:
    :param fastqdir:
    :param outdir:
    :param gencodefile:
    :return:
    """
    gencver = gencodefile.split('_')[1]
    all_files = os.listdir(fastqdir)
    all_fastq = fnm.filter(all_files, '*.fastq.gz')
    arglist = []
    for root, dirs, files in os.walk(quantdir):
        if files:
            for f in files:
                if f == 'quant.sf':
                    _, encexp = os.path.split(root)
                    expfiles = fnm.filter(all_fastq, '*' + encexp + '*')
                    assert len(expfiles) == 2, 'Missing files for ENCODE experiment: {}'.format(encexp)
                    expfile = expfiles[0]
                    components = expfile.split('.')[0].split('_')
                    new_name = '_'.join([encexp, components[2], components[3], 'mRNA', components[5]]) + '.genc-' + gencver + '.bed'
                    new_path = os.path.join(outdir, new_name)
                    arglist.append([os.path.join(root, f), new_path, gencodefile])
    return arglist


def salmon_to_bed_genes(inputfile, outputfile, gencode):
    """
    :param inputfile:
    :param outputfile:
    :param gencode:
    :return:
    """
    gene_expression = col.defaultdict(float)
    with open(inputfile, 'r') as inf:
        rows = csv.DictReader(inf, delimiter='\t')
        for r in rows:
            name_parts = r['Name'].split('|')
            gene_symbol = name_parts[5]
            if gene_symbol.startswith('mt-'):
                # ignore chrM genes
                continue
            tpm = float(r['TPM'])
            gene_expression[gene_symbol] += tpm

    annotation = js.load(open(gencode, 'r'))
    annotation = dict((g['symbol'], (g['chrom'], g['start'], g['end'], g['ensid'])) for g in annotation)
    outbuf = []
    for symbol, tpm in gene_expression.items():
        chrom, start, end, ensid = annotation[symbol]
        gene_region = chrom, start, end, ensid, symbol, str(tpm)
        outbuf.append(gene_region)

    outbuf = sorted(outbuf, key=lambda x: (x[0], x[1], x[2]))
    with open(outputfile, 'w') as outf:
        _ = outf.write('\t'.join(['#chrom', 'start', 'end', 'ensid', 'symbol', 'expr']) + '\n')
        outbuf = ['\t'.join(map(str, gene)) for gene in outbuf]
        _ = outf.write('\n'.join(outbuf) + '\n')
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
    tmpquant = config.get('Pipeline', 'tmpquant')
    bedout = config.get('Pipeline', 'bedout')
    hdfout = config.get('Pipeline', 'hdfout')
    datadir = config.get('Pipeline', 'datadir')

    enc_mdfile = config.get('Pipeline', 'encmd')
    bam_mdfile = config.get('Pipeline', 'bammd')

    inputfiles = [os.path.join(refdir, f) for f in os.listdir(refdir)]
    tmpf = link_encode_files(enc_mdfile, datadir, tempdir)
    inputfiles.extend(tmpf)
    tmpf = link_fastq_files(bam_mdfile, datadir, tempdir)
    inputfiles.extend(tmpf)

    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, name='init', output=inputfiles).mkdir(tempdir)

    convgenc = pipe.transform(task_func=convert_gtf,
                              name='convgenc',
                              input=output_from(init),
                              filter=formatter('gencode\.(?P<GENCVER>v[0-9a-zA-Z]+)\.annotation\.(?P<ASSEMBLY>\w+)\.+'),
                              output=os.path.join(refdir, 'gencode_{GENCVER[0]}_{ASSEMBLY[0]}_main.json')).jobs_limit(2)

    mergequant = pipe.collate(task_func=merge_expquant_files,
                              name='mergequant',
                              input=output_from(init),
                              # ENCSR000COQ_ENCFF000EWE_hg19_GM12878_mRNA_L13_R0.genc-v7.gtf.gz
                              filter=formatter('(?P<ENCEXP>\w+)_(?P<ENCF>\w+)_(?P<ASS>\w+)_(?P<CELL>\w+)_mRNA_(?P<LAB>\w+)_(?P<REP>\w+)\.genc-(?P<GENCVER>\w+)\.gtf\.gz'),
                              output=os.path.join(bedout, '{ENCEXP[0]}_{ASS[0]}_{CELL[0]}_mRNA_{LAB[0]}.genc-{GENCVER[0]}.bed'),
                              extras=[os.path.join(refdir, 'gencode_{GENCVER[0]}_{ASS[0]}_main.json')]).mkdir(bedout).jobs_limit(1)

    promdump = pipe.transform(task_func=dump_gene_regions,
                              name='promdump',
                              input=output_from(convgenc),
                              filter=suffix('.json'),
                              output='.2kbprom.bed',
                              output_dir=refdir,
                              extras=['promoter']).jobs_limit(2)

    bodydump = pipe.transform(task_func=dump_gene_regions,
                              name='bodydump',
                              input=output_from(convgenc),
                              filter=suffix('.json'),
                              output='.body.bed',
                              output_dir=refdir,
                              extras=['body']).jobs_limit(2)

    cmd = config.get('Pipeline', 'convdump')
    convdump = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='convdump',
                              input=output_from(promdump, bodydump),
                              filter=formatter('gencode_(?P<GENCVER>\w+)_(?P<ASSEMBLY>\w+)_main\.(?P<REGTYPE>\w+)\.bed'),
                              output=os.path.join(refdir, 'gencode_{GENCVER[0]}_{ASSEMBLY[0]}_main.{REGTYPE[0]}.h5'),
                              extras=[cmd, jobcall])

    makelists = pipe.subdivide(task_func=create_gtf_white_black_list,
                               name='makelists',
                               input=output_from(init),
                               filter=formatter('gencode\.(?P<GENCVER>v[0-9a-zA-Z]+)\.annotation\.(?P<ASSEMBLY>\w+)\.+'),
                               output=[os.path.join(refdir, 'gencode_{GENCVER[0]}_{ASSEMBLY[0]}_pctranscripts_white.gtf'),
                                       os.path.join(refdir, 'gencode_{GENCVER[0]}_{ASSEMBLY[0]}_pctranscripts_black.gtf')]).jobs_limit(2)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'salmonidx')
    salmonidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                               name='salmonidx',
                               input=output_from(init),
                               filter=formatter('(?P<TRANSCRIPTOME>[\w\.]+)\.fa\.gz$'),
                               output=os.path.join(refdir, '{TRANSCRIPTOME[0]}.idx', 'hash.bin'),
                               extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'musquant').replace('\n', ' ')
    musquant = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                            name='musquant',
                            input=output_from(init),
                            filter=formatter('(?P<ENCEXP>\w+)_(?P<ENCF>\w+)_(?P<ASS>\w+)_(?P<CELL>\w+)_mRNA_(?P<LAB>\w+)_(?P<REP>\w+)\.se\.fastq'),
                            output=os.path.join(tmpquant, '{ENCEXP[0]}', 'quant.sf'),
                            extras=[cmd, jobcall]).mkdir(tmpquant).follows(salmonidx)

    musbed_params = make_musbed_params(tmpquant, tempdir, bedout, os.path.join(refdir, 'gencode_vM1_mm9_main.json'))

    musbed = pipe.files(salmon_to_bed_genes, musbed_params, name='musbed').jobs_limit(2).follows(musquant).follows(convgenc)

    return pipe
