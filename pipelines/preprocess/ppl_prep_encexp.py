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


def link_fastq_files(metadata, indir, outdir):
    """
    :param indir:
    :param outdir:
    :param metadata:
    :return:
    """
    from preprocess.__init__ import ENCODE_LAB_MAP as LAB_MAP

    fastq_files = collect_full_paths(indir, '*.fastq.gz')
    fp_map = dict((os.path.basename(fp), fp) for fp in fastq_files)

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


def link_sra_files(metadata, indir, outdir):
    """
    :param metadata:
    :param indir:
    :param outdir:
    :return:
    """
    # this is a misnomer, should simply have consecutive IDs for all labs
    from preprocess.__init__ import ENCODE_LAB_MAP as LAB_MAP

    sra_files = collect_full_paths(indir, '*.sra')
    fp_map = dict((os.path.basename(fp), fp) for fp in sra_files)
    linked = []
    with open(metadata, 'r') as mdfile:
        rows = csv.DictReader(mdfile, delimiter='\t')
        for r in rows:
            if not bool(int(r['use'])):
                continue
            fn = os.path.basename(r['sra'])
            sraff = fn.split('.')[0]
            old_path = fp_map[fn]
            new_name = '_'.join([sraff, r['assembly'], r['cell'],
                                 'mRNA', LAB_MAP[r['lab'].strip('"')], 'R' + str(r['rep'])])
            new_name += '.pe.sra'
            new_path = os.path.join(outdir, new_name)
            if not os.path.islink(new_path):
                os.link(old_path, new_path)
            linked.append(new_path)
    return linked


def make_pe_fastq_quant(indir, species, baseout, cmd, jobcall):
    """
    :param indir:
    :param species:
    :param baseout:
    :return:
    """
    readfiles_1 = sorted(collect_full_paths(indir, '*' + species + '*1.fastq.gz'))
    readfiles_2 = sorted(collect_full_paths(indir, '*' + species + '*2.fastq.gz'))

    ctype = None
    arglist = []
    reads1 = []
    reads2 = []
    for rf1 in readfiles_1:
        fn = os.path.basename(rf1)
        cell = fn.split('_')[2]
        rep = fn.split('.')[0].split('_')[-1]
        if cell != ctype and ctype is not None:
            fmt = {'reads1': ' '.join(reads1), 'reads2': ' '.join(reads2),
                   'species': species, 'cell': ctype}
            tmp = cmd.format(**fmt)
            outfile = os.path.join(baseout, species + '_' + ctype, 'quant.sf')
            arglist.append([reads1 + reads2, outfile, tmp, jobcall])
            reads1 = []
            reads2 = []
        ctype = cell
        reads1.append(rf1)
        mates = fnm.filter(readfiles_2, '*' + species + '_' + ctype + '*' + rep + '*')
        assert len(mates) == 1, 'Unmated read: {}'.format(rf1)
        reads2.append(mates[0])
    if reads1:
        fmt = {'reads1': ' '.join(reads1), 'reads2': ' '.join(reads2),
                   'species': species, 'cell': ctype}
        tmp = cmd.format(**fmt)
        outfile = os.path.join(baseout, species + '_' + ctype, 'quant.sf')
        arglist.append([reads1 + reads2, outfile, tmp, jobcall])
    return arglist


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
    _ = outbuf.write('\t'.join(['#chrom', 'start', 'end', 'ensemblID', 'symbol']) + '\n')
    # NB: by construction, genes are sorted
    # by genomic coordinate
    for gene in genes:
        start, end = gene[regtype]
        _ = outbuf.write('\t'.join([gene['chrom'], str(start), str(end), gene['ensemblID'], gene['symbol']]) + '\n')
    with open(outputfile, 'w') as outfile:
        _ = outfile.write(outbuf.getvalue())
    return outputfile


def make_convbed_params(quantdir, fastqdir, outdir, gencodefile):
    """
    :param quantdir:
    :param fastqdir:
    :param outdir:
    :param gencodefile:
    :return:
    """
    gencver = gencodefile.split('_')[1]
    assembly = gencodefile.split('.')[0].split('_')[2]
    all_files = os.listdir(fastqdir)
    all_fastq = fnm.filter(all_files, '*.fastq.gz')
    arglist = []
    for root, dirs, files in os.walk(quantdir):
        if files:
            for f in files:
                if f == 'quant.sf':
                    _, encexp = os.path.split(root)
                    expfiles = fnm.filter(all_fastq, '*' + encexp + '*')
                    assert len(expfiles) >= 2, 'Missing files for ENCODE experiment: {} ({})'.format(encexp, expfiles)
                    if encexp.startswith('ENCSR'):
                        expfile = expfiles[0]
                        components = expfile.split('.')[0].split('_')
                        if components[2] != assembly:
                            continue
                        new_name = '_'.join([encexp, components[2], components[3], 'mRNA', components[5]]) + '.genc-' + gencver + '.bed'
                        new_path = os.path.join(outdir, new_name)
                        arglist.append([os.path.join(root, f), new_path, gencodefile])
                    elif encexp.startswith(assembly):
                        components = expfiles[0].split('_')
                        cell, lab = components[2], components[4]
                        all_sra_ids = set([(expf.split('_')[0]) for expf in expfiles])
                        all_sra_ids = list(map(list, all_sra_ids))
                        joined_id = ''
                        for chars in zip(*all_sra_ids):
                            if len(set(chars)) == 1:
                                joined_id += chars[0]
                            elif 'u' in joined_id:
                                joined_id += ''.join(sorted(set(chars)))
                            else:
                                joined_id += 'u' + ''.join(sorted(set(chars)))
                        new_name = '_'.join([joined_id, assembly, cell, 'mRNA', lab]) + '.genc-' + gencver + '.bed'
                        new_path = os.path.join(outdir, new_name)
                        arglist.append([os.path.join(root, f), new_path, gencodefile])
                    else:
                        # assembly does not match
                        pass
    return arglist


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


def _read_fasta_header(line):
    """
    :param line:
    :return:
    """
    gene_id = ''
    line = line.decode('ascii').strip()
    if line.startswith('>'):
        gene_id = line.split('|')[1]
    return gene_id


def filter_gencode(annotation, outpath, transcriptome):
    """
    :param annotation:
    :param outpath:
    :param transcriptome:
    :return:
    """
    with gz.open(transcriptome, 'rb') as fasta:
        genes_in_transcriptome = set([_read_fasta_header(line) for line in fasta])
    gencode_genes = dict()
    #print(os.path.basename(annotation))
    skipped = 0
    with gz.open(annotation, 'rb') as gtf:
        for line in gtf:
            line = line.decode('ascii').strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if parts[2] == 'gene':
                c, s, e, strand, infos = _parse_gtf_line(line)
                if infos['gene_id'] not in genes_in_transcriptome:
                    continue
                if e - s < 250:
                    print('Skipping - size: {} / {}'.format(infos['gene_name'], e - s))
                    skipped += 1
                    continue
                this_gene = {'chrom': c, 'start': s, 'end': e, 'symbol': infos['gene_name'],
                             'ensemblID': infos['gene_id'], 'transcripts': [], 'strand': strand}
                if strand == '+':
                    this_gene['promoter'] = s - 1250, s + 250
                    if e - (s + 250) < 250:
                        this_gene['body'] = s, e
                    elif e - (s + 250) < 500:
                        this_gene['body'] = s + 125, e
                    else:
                        this_gene['body'] = s + 250, e
                else:
                    this_gene['promoter'] = e - 250, e + 1250
                    if (e - 250) - s < 250:
                        this_gene['body'] = s, e
                    elif (e - 250) - s < 500:
                        this_gene['body'] = s, e - 125
                    else:
                        this_gene['body'] = (s, e)
                assert this_gene['body'][0] < this_gene['body'][1], 'Malformed gene: {}'.format(this_gene)
                gencode_genes[this_gene['ensemblID']] = this_gene
            elif parts[2] == 'transcript':
                c, s, e, strand, infos = _parse_gtf_line(line)
                this_trans = {'chrom': c, 'start': s, 'end': e, 'strand': strand, 'ensemblID': infos['transcript_id']}
                gid = infos['gene_id']
                if gid not in gencode_genes:
                    continue
                gencode_genes[gid]['transcripts'].append(this_trans)
    #print('Total genes skipped: {}'.format(skipped))
    assert gencode_genes, 'No genes selected from GTF file'
    genes = sorted(gencode_genes.values(), key=lambda d: (d['chrom'], d['start'], d['end']))
    with open(outpath, 'w') as outf:
        _ = js.dump(genes, outf, indent=1)
    return outpath


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
            gene_id = name_parts[1]
            gene_expression[gene_id] += float(r['TPM'])

    annotation = js.load(open(gencode, 'r'))
    annotation = dict((g['ensemblID'], (g['chrom'], g['start'], g['end'], g['symbol'])) for g in annotation)
    outbuf = []
    for geneid, tpm in gene_expression.items():
        try:
            chrom, start, end, symbol = annotation[geneid]
            if tpm < 1:
                level = '0'
            elif tpm < 10:
                level = '1'
            else:
                level = '2'
            gene_region = chrom, start, end, geneid, symbol, str(tpm), level
            outbuf.append(gene_region)
        except KeyError:
            pass  # filtered small genes

    outbuf = sorted(outbuf, key=lambda x: (x[0], x[1], x[2]))
    with open(outputfile, 'w') as outf:
        _ = outf.write('\t'.join(['#chrom', 'start', 'end', 'ensemblID', 'symbol', 'tpm', 'y_depvar']) + '\n')
        outbuf = ['\t'.join(map(str, gene)) for gene in outbuf]
        _ = outf.write('\n'.join(outbuf) + '\n')
    return outputfile


def convert_hcop_table(inputfile, outputfile):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    fn = os.path.basename(inputfile)
    from_species = fn.split('_')[0]
    to_species = fn.split('_')[1]
    outbuffer = io.StringIO()
    _ = outbuffer.write('\t'.join([from_species, to_species, 'support', 'assert_ids']) + '\n')
    with gz.open(inputfile, 'rb') as inf:
        # skip header
        _ = inf.readline()
        for line in inf:
            line = line.decode('ascii')
            if not line:
                continue
            _, from_ens, from_assert, _, to_ens, to_assert = line.strip().split()
            if not (from_ens.startswith('ENS') and to_ens.startswith('ENS')):
                continue
            support = str(min(len(from_assert.split(',')), len(to_assert.split(','))))
            _ = outbuffer.write('\t'.join([from_ens, to_ens, support, from_assert + '@' + to_assert]) + '\n')
    with open(outputfile, 'w') as outf:
        _ = outf.write(outbuffer.getvalue())
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
    orthdir = config.get('Pipeline', 'orthdir')
    tempdir = config.get('Pipeline', 'tempdir')
    tmpquant = config.get('Pipeline', 'tmpquant')
    hdfout = config.get('Pipeline', 'hdfout')
    datadir = config.get('Pipeline', 'datadir')

    fastq_mdfile = config.get('Pipeline', 'fastqmd')
    sra_mdfile = config.get('Pipeline', 'sramd')

    inputfiles = [os.path.join(refdir, f) for f in os.listdir(refdir)]
    inputfiles = [inpf for inpf in inputfiles if os.path.isfile(inpf)]
    tmpf = link_fastq_files(fastq_mdfile, datadir, tempdir)
    inputfiles.extend(tmpf)
    tmpf = link_sra_files(sra_mdfile, datadir, tempdir)
    inputfiles.extend(tmpf)

    # step 0: initiate pipeline with input files
    init = pipe.originate(task_func=lambda x: x, name='init', output=inputfiles).mkdir(tempdir)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'musidx')
    musidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='musidx',
                            input=output_from(init),
                            filter=formatter('(?P<TRANSCRIPTOME>gencode\.vM1\.\w+)\.fa\.gz$'),
                            output=os.path.join(refdir, '{TRANSCRIPTOME[0]}_k21.idx', 'hash.bin'),
                            extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'hsaidx')
    hsaidx = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='hsaidx',
                            input=output_from(init),
                            filter=formatter('(?P<TRANSCRIPTOME>gencode\.v19\.\w+)\.fa\.gz$'),
                            output=os.path.join(refdir, '{TRANSCRIPTOME[0]}_k31.idx', 'hash.bin'),
                            extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'musquant').replace('\n', ' ')
    musquant = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                            name='musquant',
                            input=output_from(init),
                            filter=formatter('(?P<ENCEXP>\w+)_(?P<ENCF>\w+)_mm9_(?P<CELL>\w+)_mRNA_(?P<LAB>\w+)_(?P<REP>\w+)\.se\.fastq'),
                            output=os.path.join(tmpquant, '{ENCEXP[0]}', 'quant.sf'),
                            extras=[cmd, jobcall]).mkdir(tmpquant).follows(musidx)

    cmd = config.get('Pipeline', 'hsaquant').replace('\n', ' ')
    hsaquant = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                            name='hsaquant',
                            input=output_from(init),
                            filter=formatter('(?P<ENCEXP>\w+)_(?P<ENCF>\w+)_hg19_(?P<CELL>\w+)_mRNA_(?P<LAB>\w+)_(?P<REP>\w+)\.se\.fastq'),
                            output=os.path.join(tmpquant, '{ENCEXP[0]}', 'quant.sf'),
                            extras=[cmd, jobcall]).mkdir(tmpquant).follows(musidx)

    cmd = config.get('Pipeline', 'sradump')
    sradump = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                             name='sradump',
                             input=output_from(init),
                             filter=formatter('(?P<SRARUN>SRR[0-9]+)_.+\.sra'),
                             output=os.path.join(tempdir, '{SRARUN[0]}*.fastq.gz'),
                             extras=[tempdir, '{SRARUN[0]}*.fastq.gz', cmd, jobcall])

    cmd = config.get('Pipeline', 'musquantpe').replace('\n', ' ')
    muspe_params = make_pe_fastq_quant(tempdir, 'mm9', tmpquant, cmd, jobcall)
    musquantpe = pipe.files(sci_obj.get_jobf('ins_out'),
                            muspe_params,
                            name='musquantpe')


    genc_re = 'gencode\.(?P<GENCVER>v[M0-9]+)\.\w+\.(?P<ASSEMBLY>\w+)\.gtf\.gz'
    gencdump = pipe.transform(task_func=filter_gencode,
                              name='gencdump',
                              input=output_from(init),
                              filter=formatter(genc_re),
                              output=os.path.join(refdir, 'gencode_{GENCVER[0]}_{ASSEMBLY[0]}.pc_transcripts.json'),
                              extras=[os.path.join(refdir, 'gencode.{GENCVER[0]}.pc_transcripts.fa.gz')])

    mus_bedparams = make_convbed_params(tmpquant, tempdir, tmpquant,
                                        os.path.join(refdir, 'gencode_vM1_mm9.pc_transcripts.json'))
    musbed = pipe.files(salmon_to_bed_genes, mus_bedparams, name='musbed').jobs_limit(2)

    hsa_bedparams = make_convbed_params(tmpquant, tempdir, tmpquant,
                                        os.path.join(refdir, 'gencode_v19_hg19.pc_transcripts.json'))
    hsabed = pipe.files(salmon_to_bed_genes, hsa_bedparams, name='hsabed').jobs_limit(2)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    promdump = pipe.transform(task_func=dump_gene_regions,
                              name='promdump',
                              input=output_from(gencdump),
                              filter=suffix('.json'),
                              output='.promoter.bed',
                              output_dir=refdir,
                              extras=['promoter'])

    bodydump = pipe.transform(task_func=dump_gene_regions,
                              name='bodydump',
                              input=output_from(gencdump),
                              filter=suffix('.json'),
                              output='.body.bed',
                              output_dir=refdir,
                              extras=['body'])

    # gencode_v19_hg19.pc_transcripts.body.bed
    cmd = config.get('Pipeline', 'convreg')
    convregions = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                 name='convregions',
                                 input=output_from(promdump, bodydump),
                                 filter=formatter('gencode_(?P<GENCVER>\w+)_(?P<ASSEMBLY>\w+)\.pc_transcripts\.(?P<REGTYPE>\w+)\.bed'),
                                 output=os.path.join(refdir, 'gencode_{GENCVER[0]}_{ASSEMBLY[0]}.coding.{REGTYPE[0]}.h5'),
                                 extras=[cmd, jobcall])

    # ENCSR000CGU_mm9_ESB4_mRNA_L08.genc-vM1.bed
    cmd = config.get('Pipeline', 'convexpr')
    convexpr = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='convexpr',
                              input=output_from(musbed, hsabed),
                              filter=formatter('(?P<SAMPLE>\w+)\.genc-(?P<GENCVER>\w+)\.bed'),
                              output=os.path.join(hdfout, '{SAMPLE[0]}.genc-{GENCVER[0]}.h5'),
                              extras=[cmd, jobcall]).mkdir(hdfout)


    orth_files = os.listdir(orthdir)
    orth_files = [os.path.join(orthdir, f) for f in orth_files]

    initorth = pipe.originate(lambda x: x, orth_files, name='initorth')

    convorth = pipe.transform(task_func=convert_hcop_table,
                              name='convorth',
                              input=output_from(initorth),
                              filter=suffix('_six_column.txt.gz'),
                              output='_orthologs.tsv',
                              output_dir=orthdir).jobs_limit(2)

    cmd = config.get('Pipeline', 'runall')
    runall = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                          name='runall',
                          input=output_from(musidx, hsaidx,
                                            musquant, hsaquant,
                                            gencdump, sradump,
                                            musbed, hsabed,
                                            promdump, bodydump,
                                            convregions, convexpr,
                                            convorth,
                                            musquantpe),
                          filter=formatter(),
                          output=os.path.join(tempdir, 'runall_encexp.chk'),
                          extras=[cmd, jobcall])





    return pipe
