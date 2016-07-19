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
import re as re
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
    fastq_files = collect_full_paths(indir, '*.fastq.gz')
    fp_map = dict((os.path.basename(fp), fp) for fp in fastq_files)
    fastq_md = dict()

    linked = []
    pc = 1
    pairs = dict()
    with open(metadata, 'r', newline='') as mdfile:
        rows = csv.DictReader(mdfile, delimiter='\t')
        for r in rows:
            if 'use' in r:
                if not bool(int(r['use'])):
                    continue
            if r['File format'] != 'fastq':
                continue
            fastq_md[r['File accession']] = r
            if r['Run type'] == 'pe':
                if r['File accession'] in pairs or r['Paired with'] in pairs:
                    continue
                pairs[r['File accession']] = str(pc)
                pairs[r['Paired with']] = str(pc)
                pc += 1

    for fn, fp in fp_map.items():
        assert os.path.isfile(fp), 'File missing: {}'.format(fp)
        acc = fn.split('.')[0]
        md = fastq_md[acc]
        target_name = '_'.join([md['Experiment accession'], md['File accession'], md['Biosample organism'],
                                md['Biosample term name'], 'mRNA', md['Lab']])
        if md['Run type'] == 'pe':
            target_name += '.' + md['Run type'] + '-' + md['Read length']
            pnum = pairs[acc]
            target_name += '.b' + md['Biological replicate(s)'] + 'p' + pnum + 'r_' + md['Paired end']
        else:
            target_name += '.' + md['Run type'] + '-' + md['Read length']
            target_name += '.b' + md['Biological replicate(s)']
        target_name += '.fastq.gz'
        target_path = os.path.join(outdir, target_name)
        if not os.path.islink(target_path):
            os.symlink(fp, target_path)
        linked.append(target_path)
    return linked


def link_sra_files(metadata, indir, outdir):
    """
    :param metadata:
    :param indir:
    :param outdir:
    :return:
    """
    sra_files = collect_full_paths(indir, '*.sra')
    fp_map = dict((os.path.basename(fp), fp) for fp in sra_files)
    linked = []
    pc = 0
    with open(metadata, 'r', newline='') as mdfile:
        rows = csv.DictReader(mdfile, delimiter='\t')
        for r in rows:
            if not bool(int(r['use'])):
                continue
            pc += 1
            fn = os.path.basename(r['sra'])
            expnum = r['sra'].split('/')[-3]
            assert expnum.startswith('SRX'), 'Wrong entry selected: {}'.format(expnum)
            expnum = int(expnum.strip('SRX'))
            if r['rep'] == '1':
                expnum += 1
            elif r['rep'] == '3':
                expnum -= 1
            else:
                pass
            expid = 'SRX' + str(expnum)
            acc = fn.split('.')[0]
            old_path = fp_map[fn]
            new_name = '_'.join([expid, acc, r['species'], r['cell'], 'mRNA', r['lab']])
            new_name += '.' + r['runtype'] + '-' + r['readlength'] + '.b' + r['rep'] + 'p' + str(pc) + 'r'
            new_name += '.sra'
            new_path = os.path.join(outdir, new_name)
            if not os.path.islink(new_path):
                os.symlink(old_path, new_path)
            linked.append(new_path)
    return linked


def normalize_btau_annotation(inpath, outpath):
    """
    :param inpath:
    :param outpath:
    :return:
    """
    entries = []
    chrom = re.compile('^chr[0-9XY]$')
    with gz.open(inpath, 'rb') as gff:
        for line in gff:
            line = line.decode('ascii').strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if cols[2] not in ['gene', 'transcript']:
                continue
            if '_protein_coding' not in cols[1]:
                continue
            if chrom.match(cols[0]) is None:
                continue
            entry = {'chrom': cols[0], 'source': 'Ensembl', 'biotype': cols[2],
                     'start': cols[3], 'end': cols[4], 'score': cols[5],
                     'strand': cols[6], 'frame': cols[7], 'attributes': cols[8]}
            entries.append(entry)
    assert entries, 'No genes/transcripts selected from file {}'.format(inpath)
    entries = sorted(entries, key=lambda x: (x['chrom'], int(x['start']), x['biotype'], int(x['end'])))
    fieldnames = ['chrom', 'source', 'biotype', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    outbuffer = io.StringIO()
    writer = csv.DictWriter(outbuffer, fieldnames=fieldnames, delimiter='\t')
    writer.writerows(entries)
    with gz.open(outpath, 'wb') as out:
        _ = out.write(outbuffer.getvalue().encode('utf-8'))
    return outpath


def _parse_gtf_line(line):
    """
    :param line:
    :return:
    """
    chrom, _, entry, start, end, _, strand, _, attributes = line.split('\t')
    if 'Dbxref' in attributes:
        infos = dict()
        try:
            for entry in attributes.split(';'):
                k, v = entry.split('=')
                if k == 'Dbxref':
                    k, v = v.split(':')
                infos[k.strip()] = v.strip()
        except Exception:
            sys.stderr.write('\n{}\n'.format(attributes))
            raise
    else:
        try:
            infos = dict([(t.split()[0], t.split()[1].strip('"')) for t in attributes.split(';') if t and len(t) > 1])
        except IndexError:
            sys.stderr.write('\n{}\n'.format(attributes))
            raise
    return chrom, int(start), int(end), strand, infos


def make_genemap(inputfile, outputpair):
    """ Pure convenience function to build gene maps (transcript to gene)
    for Salmon to directly perform quantification on the gene level.
    Requires to clean up the transcriptome files s.t. the FASTA header
    only contains the transcript name
    :param inputfile:
    :param outputpair:
    :return:
    """
    norm_fasta = io.StringIO()
    gene_map = io.StringIO()
    with gz.open(inputfile, 'rb') as fa:
        seqbuf = io.StringIO()
        c = 0
        for line in fa:
            line = line.decode('ascii').strip()
            if line.startswith('>'):
                if c > 0:
                    for mobj in re.finditer('(\w{1,80})', seqbuf.getvalue()):
                        _ = norm_fasta.write(mobj.group(0) + '\n')
                    c = 0
                    seqbuf = io.StringIO()
                parts = line.split('|')
                trans, gene = parts[0], parts[1]
                _ = norm_fasta.write(trans.strip() + '\n')
                trans = trans.strip('> ')
                _ = gene_map.write(trans + '\t' + gene.strip() + '\n')
            else:
                c += seqbuf.write(line)
    out_fasta, out_map = outputpair
    with gz.open(out_fasta, 'wb') as fa:
        _ = fa.write(norm_fasta.getvalue().encode('ascii'))

    with open(out_map, 'w') as out:
        _ = out.write(gene_map.getvalue())

    return out_fasta, out_map


def extract_exp_mate_info(filepaths):
    """
    :param filepaths:
    :return:
    """
    res = []
    for fp in filepaths:
        path, fn = os.path.split(fp)
        expid, _, _, _, _, libinfo, _ = fn.split('_')
        mate_id = libinfo.split('.')[-1]
        res.append((expid, mate_id, fp))
    return sorted(res)


def make_pe_quant(indir, baseout, species, index, cmd, jobcall):
    """
    :param indir:
    :param baseout:
    :param species:
    :param index:
    :param cmd:
    :param jobcall:
    :return:
    """
    readfiles_1 = sorted(collect_full_paths(indir, '*' + species + '*_1.fastq.gz'))
    readfiles_2 = sorted(collect_full_paths(indir, '*' + species + '*_2.fastq.gz'))

    readfiles_1 = extract_exp_mate_info(readfiles_1)
    readfiles_2 = extract_exp_mate_info(readfiles_2)

    curr_exp = readfiles_1[0][0]
    genemap = index.rsplit('_', 1)[0] + '_map.tsv'
    arglist = []
    reads1 = []
    reads2 = []
    for r1, r2 in zip(readfiles_1, readfiles_2):
        assert r1[0] == r2[0], 'Experiment mismatch: {} - {}'.format(r1, r2)
        assert r1[1] == r2[1], 'Mate mismatch: {} - {}'.format(r1, r2)
        if r1[0] == curr_exp:
            reads1.append(r1[2])
            reads2.append(r2[2])
        else:
            tmp = cmd.format(**{'reads1': ' '.join(reads1), 'reads2': ' '.join(reads2),
                                'expid': curr_exp, 'index': index, 'genemap': genemap})
            outpath = os.path.join(baseout, curr_exp, 'quant.genes.sf')
            arglist.append([reads1 + reads2, outpath, tmp, jobcall])
            curr_exp = r1[0]
            reads1 = [r1[2]]
            reads2 = [r2[2]]
    tmp = cmd.format(**{'reads1': ' '.join(reads1), 'reads2': ' '.join(reads2),
                        'expid': curr_exp, 'index': index, 'genemap': genemap})
    outpath = os.path.join(baseout, curr_exp, 'quant.genes.sf')
    arglist.append([reads1 + reads2, outpath, tmp, jobcall])
    return arglist


def salmon_to_bed(inputfile, outputfile, genemodel, datadir, expid):
    """
    :param inputfile:
    :param outputfile:
    :param genemodel:
    :return:
    """
    ra_genemodel = dict([(d['id'], d) for d in js.load(open(genemodel, 'r'))])
    with open(inputfile, 'r', newline='') as quantfile:
        rows = csv.DictReader(quantfile, delimiter='\t')
        for r in rows:
            this_gene = ra_genemodel[r['Name']]
            tpm = float(r['TPM'])
            lvl = 0
            if 1 <= tpm < 10:
                lvl = 1
            if 10 <= tpm:
                lvl = 2
            this_gene['tpm'] = tpm
            this_gene['level'] = lvl

    fp, fn = os.path.split(genemodel)
    parts = fn.split('_')
    species, assembly, auth, ver = parts[0], parts[1], parts[2], parts[3].split('.')[0]

    assoc_files = os.listdir(datadir)
    assoc_files = fnm.filter(assoc_files, expid + '*' + '.fastq.gz')
    assert assoc_files, 'No files for experiment ID: {}'.format(expid)
    common_name = ''
    for af in assoc_files:
        parts = af.split('_')
        new_name = '_'.join([parts[0], assembly, parts[3], parts[4], parts[5].split('.')[0]])
        if common_name:
            assert new_name == common_name, 'Information mismatch: {} \n {}'.format(common_name, assoc_files)
        common_name = new_name

    outpath = os.path.join(datadir, 'tmp', common_name + '.' + auth + '-' + ver + '.bed')
    genes = sorted(ra_genemodel.values(), key=lambda d: (d['chrom'], d['start'], d['end'], d['id']))
    with open(outpath, 'w') as out:
        writer = csv.DictWriter(out, fieldnames=['chrom', 'start', 'end', 'id', 'tpm', 'strand', 'symbol', 'level'],
                                extrasaction='ignore', delimiter='\t')
        writer.writeheader()
        writer.writerows(genes)
    return outputfile


def dump_gene_regions(inputfile, outputfile, regtype):
    """
    :param inputfile:
    :param outputfile:
    :param regtype:
    :return:
    """
    regtypes = {'core': {'+': {'refpoint': 'start', 'from': -250, 'to': 250},
                         '-': {'refpoint': 'end', 'from': -250, 'to': 250}},
                'uprr': {'+': {'refpoint': 'start', 'from': -500, 'to': -5500},
                         '-': {'refpoint': 'end', 'from': 500, 'to': 5500}}}
    fieldnames = ['chrom', 'start', 'end', 'id', 'symbol']
    with open(inputfile, 'r') as infile:
        genes = js.load(infile)
    outbuf = []
    for gene in genes:
        select = {k: gene[k] for k in fieldnames}
        if regtype == 'body':
            pass
        else:
            adapt = regtypes[regtype][gene['strand']]
            select['start'] = gene[adapt['refpoint']] + adapt['from']
            select['end'] = gene[adapt['refpoint']] + adapt['to']
        outbuf.append(select)
    with open(outputfile, 'w') as out:
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(outbuf)
    return outputfile


def _read_fasta_header(line):
    """ This is specific to the transcriptome annotation released by GENCODE
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
                assert s < e, 'Malformed gene: {}'.format(line.strip())
                this_gene = {'chrom': c, 'start': s, 'end': e, 'symbol': infos['gene_name'],
                             'id': infos['gene_id'], 'transcripts': [], 'strand': strand, 'authority': 'Ensembl'}
                assert this_gene['id'].startswith('ENS'), 'Unknown gene id: {}'.format(this_gene)
                gencode_genes[this_gene['id']] = this_gene
            elif parts[2] == 'transcript':
                c, s, e, strand, infos = _parse_gtf_line(line)
                this_trans = {'chrom': c, 'start': s, 'end': e, 'strand': strand,
                              'id': infos['transcript_id'], 'authority': 'Ensembl'}
                assert this_trans['id'].startswith('ENS'), 'Unknown transcript id: {}'.format(this_trans)
                gid = infos['gene_id']
                if gid not in gencode_genes:
                    continue
                gencode_genes[gid]['transcripts'].append(this_trans)
    assert gencode_genes, 'No genes selected from GTF file'
    genes = sorted(gencode_genes.values(), key=lambda d: (d['chrom'], d['start'], d['end']))
    with open(outpath, 'w') as outf:
        _ = js.dump(genes, outf, indent=1)
    return outpath


# the following functions may be outdated


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

    encode_mdfile = config.get('Pipeline', 'encmd')
    sra_mdfile = config.get('Pipeline', 'sramd')

    inputfiles = [os.path.join(refdir, f) for f in os.listdir(refdir)]
    inputfiles = [inpf for inpf in inputfiles if os.path.isfile(inpf)]
    tmpf = link_fastq_files(encode_mdfile, datadir, tempdir)
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

    cmd = config.get('Pipeline', 'sradump')
    sradump = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                             name='sradump',
                             input=output_from(init),
                             filter=formatter('(?P<SRARUN>SRR[0-9]+)_.+\.sra'),
                             output=os.path.join(tempdir, '*{SRARUN[0]}*.fastq.gz'),
                             extras=[tempdir, '*{SRARUN[0]}*.fastq.gz', cmd, jobcall])

    # gencode.v19.annotation.hg19.gtf.gz
    # gencode.v19.pc_transcripts.fa.gz
    hsagencode = pipe.transform(task_func=filter_gencode,
                                name='hsagencode',
                                input=os.path.join(refdir, 'gencode.v19.annotation.gtf.gz'),
                                filter=formatter('.+'),
                                output=os.path.join(refdir, 'hsa_hg19_gencode_v19.pc_transcripts.json'),
                                extras=[os.path.join(refdir, 'gencode.v19.pc_transcripts.fa.gz')])

    # gencode.v19.annotation.hg19.gtf.gz
    # gencode.v19.pc_transcripts.fa.gz
    mmugencode = pipe.transform(task_func=filter_gencode,
                                name='mmugencode',
                                input=os.path.join(refdir, 'gencode.vM1.annotation.gtf.gz'),
                                filter=formatter('.+'),
                                output=os.path.join(refdir, 'mmu_mm9_gencode_vM1.pc_transcripts.json'),
                                extras=[os.path.join(refdir, 'gencode.vM1.pc_transcripts.fa.gz')])

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'mrgens').replace('\n', ' ')
    mrgens = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='mrgens',
                            input=output_from(init),
                            filter=formatter('ensGene_UCSC_(?P<SPECIES>[a-z]+)_(?P<ASSEMBLY>[A-Za-z0-9]+)\.txt\.gz'),
                            output=os.path.join(refdir, '{SPECIES[0]}_{ASSEMBLY[0]}_ensembl_vN.pc_transcripts.gff.gz'),
                            extras=[cmd, jobcall])

    normbtau = pipe.transform(task_func=normalize_btau_annotation,
                              name='normbtau',
                              input=os.path.join(refdir, 'Ensembl75_liftOver_Btau_4.6.1_genes.gff3.gz'),
                              filter=formatter('.+'),
                              output=os.path.join(refdir, 'bta_bosTau7_ensembl_v75.pc_transcripts.gff.gz'))

    cmd = config.get('Pipeline', 'mktrans').replace('\n', ' ')
    mktrans = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='mktrans',
                             input=output_from(normbtau, mrgens),
                             filter=formatter('(?P<SPECIES>[a-z]+)_(?P<ASSEMBLY>[0-9A-Za-z]+)_(?P<VERSION>[0-9A-Za-z_]+)\.pc_transcripts\.gff\.gz'),
                             output=os.path.join(refdir, '{SPECIES[0]}_{ASSEMBLY[0]}_{VERSION[0]}.pc_transcripts.fa.gz'),
                             extras=[cmd, jobcall])

    mkhsamap = pipe.subdivide(task_func=make_genemap,
                              name='mkhsamap',
                              input=os.path.join(refdir, 'gencode.v19.pc_transcripts.fa.gz'),
                              filter=formatter('.+'),
                              output=[os.path.join(refdir, 'hsa_hg19_gencode_v19.pctr_norm.fa.gz'),
                                      os.path.join(refdir, 'hsa_hg19_gencode_v19.pctr_map.tsv')])

    mkmmumap = pipe.subdivide(task_func=make_genemap,
                              name='mkmmumap',
                              input=os.path.join(refdir, 'gencode.vM1.pc_transcripts.fa.gz'),
                              filter=formatter('.+'),
                              output=[os.path.join(refdir, 'mmu_mm9_gencode_vM1.pctr_norm.fa.gz'),
                                      os.path.join(refdir, 'mmu_mm9_gencode_vM1.pctr_map.tsv')])

    mkgenemap = pipe.subdivide(task_func=make_genemap,
                               name='mkgenemap',
                               input=output_from(mktrans),
                               filter=formatter('(?P<TRANSCRIPTOME>[\w\.]+).pc_transcripts.fa.gz'),
                               output=[os.path.join(refdir, '{TRANSCRIPTOME[0]}.pctr_norm.fa.gz'),
                                       os.path.join(refdir, '{TRANSCRIPTOME[0]}.pctr_map.tsv')]).jobs_limit(2)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # building all Salmon indices
    cmd = config.get('Pipeline', 'hsaidx31')
    hsaidx31 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='hsaidx31',
                              input=output_from(mkhsamap),
                              filter=formatter('(?P<TRANSCRIPTOME>hsa_hg19_gencode_v19\.pctr_norm)\.fa\.gz$'),
                              output=os.path.join(refdir, '{TRANSCRIPTOME[0]}.k31.idx', 'hash.bin'),
                              extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'mmuidx13')
    mmuidx13 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='mmuidx13',
                              input=output_from(mkmmumap),
                              filter=formatter('(?P<TRANSCRIPTOME>mmu_mm9_gencode_vM1\.pctr_norm)\.fa\.gz$'),
                              output=os.path.join(refdir, '{TRANSCRIPTOME[0]}.k13.idx', 'hash.bin'),
                              extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'genidx19')
    mmuidx19 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='mmuidx19',
                              input=output_from(mkmmumap),
                              filter=formatter('(?P<TRANSCRIPTOME>mmu_mm9_gencode_vM1\.pctr_norm)\.fa\.gz$'),
                              output=os.path.join(refdir, '{TRANSCRIPTOME[0]}.k19.idx', 'hash.bin'),
                              extras=[cmd, jobcall])

    genidx19 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                              name='genidx19',
                              input=output_from(mkgenemap),
                              filter=formatter('(?P<TRANSCRIPTOME>[\w\.]+)\.fa\.gz$'),
                              output=os.path.join(refdir, '{TRANSCRIPTOME[0]}.k19.idx', 'hash.bin'),
                              extras=[cmd, jobcall])

    # quantification
    cmd = config.get('Pipeline', 'qmmuse').replace('\n', ' ')
    qmmuse = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                          name='qmmuse',
                          input=output_from(init),
                          filter=formatter('(?P<EXPID>\w+)_(?P<RUNID>\w+)_mmu_(?P<CELL>\w+)_mRNA_(?P<LAB>\w+)\.se.+'),
                          output=os.path.join(tmpquant, '{EXPID[0]}', 'quant.genes.sf'),
                          extras=[cmd, jobcall]).mkdir(tmpquant)

    cmd = config.get('Pipeline', 'qallpe').replace('\n', ' ')
    qhsape_params = make_pe_quant(tempdir, tmpquant, 'hsa',
                                  os.path.join(refdir, 'hsa_hg19_gencode_v19.pctr_norm.k31.idx'), cmd, jobcall)
    qhsape = pipe.files(sci_obj.get_jobf('ins_out'),
                        qhsape_params,
                        name='qhsape').follows(hsaidx31)

    qmmupe_params = make_pe_quant(tempdir, tmpquant, 'mmu',
                                  os.path.join(refdir, 'mmu_mm9_gencode_vM1.pctr_norm.k19.idx'), cmd, jobcall)
    qmmupe = pipe.files(sci_obj.get_jobf('ins_out'),
                        qmmupe_params,
                        name='qmmupe').follows(mmuidx19)

    qsscpe_params = make_pe_quant(tempdir, tmpquant, 'ssc',
                                  os.path.join(refdir, 'ssc_susScr2_ensembl_vN.pctr_norm.k19.idx'), cmd, jobcall)
    qsscpe = pipe.files(sci_obj.get_jobf('ins_out'),
                        qsscpe_params,
                        name='qsscpe').follows(genidx19)

    qbtape_params = make_pe_quant(tempdir, tmpquant, 'bta',
                                   os.path.join(refdir, 'bta_bosTau7_ensembl_v75.pctr_norm.k19.idx'), cmd, jobcall)
    qbtape = pipe.files(sci_obj.get_jobf('ins_out'),
                        qbtape_params,
                        name='qbtape').follows(genidx19)

    # Salmon segfaults here
    #qcfape_params = make_pe_quant(tempdir, tmpquant, 'cfa',
    #                              os.path.join(refdir, 'cfa_canFam3_ensembl_vN.pctr_norm.k19.idx'), cmd, jobcall)
    #qcfape = pipe.files(sci_obj.get_jobf('ins_out'),
    #                    qcfape_params,
    #                    name='qcfape').follows(genidx19)

    # conversion to BED
    cvbedhsa = pipe.subdivide(task_func=salmon_to_bed,
                              name='cvbedhsa',
                              input=output_from(qhsape),
                              filter=formatter('.+/(?P<EXPID>[A-Z0-9]+)/quant.genes.sf'),
                              output=os.path.join(tmpquant, '{EXPID[0]}_*'),
                              extras=[os.path.join(refdir, 'hsa_hg19_gencode_v19.pc_transcripts.json'),
                                      tempdir, '{EXPID[0]}']).follows(qhsape).jobs_limit(2)

    cvbedmmu = pipe.subdivide(task_func=salmon_to_bed,
                              name='cvbedmmu',
                              input=output_from(qmmupe, qmmuse),
                              filter=formatter('.+/(?P<EXPID>[A-Z0-9]+)/quant.genes.sf'),
                              output=os.path.join(tmpquant, '{EXPID[0]}_*'),
                              extras=[os.path.join(refdir, 'mmu_mm9_gencode_vM1.pc_transcripts.json'),
                                      tempdir, '{EXPID[0]}']).follows(qmmupe).jobs_limit(2)

    cvbedbta = pipe.subdivide(task_func=salmon_to_bed,
                              name='cvbedbta',
                              input=output_from(qbtape),
                              filter=formatter('.+/(?P<EXPID>[A-Z0-9]+)/quant.genes.sf'),
                              output=os.path.join(tmpquant, '{EXPID[0]}_*'),
                              extras=[os.path.join(refdir, 'bta_bosTau7_ensembl_v75.pc_transcripts.json'),
                                      tempdir, '{EXPID[0]}']).follows(qbtape).jobs_limit(2)

    cvbedssc = pipe.subdivide(task_func=salmon_to_bed,
                              name='cvbedssc',
                              input=output_from(qsscpe),
                              filter=formatter('.+/(?P<EXPID>[A-Z0-9]+)/quant.genes.sf'),
                              output=os.path.join(tmpquant, '{EXPID[0]}_*'),
                              extras=[os.path.join(refdir, 'ssc_susScr2_ensembl_vN.pc_transcripts.json'),
                                      tempdir, '{EXPID[0]}']).follows(qsscpe).jobs_limit(2)


    # dump various gene regions to BED for conversion

    all_models = collect_full_paths(refdir, '*pc_transcripts.json')

    dumpbody = pipe.transform(task_func=dump_gene_regions,
                              name='dumpbody',
                              input=all_models,
                              filter=suffix('.pc_transcripts.json'),
                              output='.body.bed',
                              extras=['body']).follows(mkgenemap)

    dumpcore = pipe.transform(task_func=dump_gene_regions,
                              name='dumpcore',
                              input=all_models,
                              filter=suffix('.pc_transcripts.json'),
                              output='.core.bed',
                              extras=['core']).follows(mkgenemap)

    dumpuprr = pipe.transform(task_func=dump_gene_regions,
                              name='dumpuprr',
                              input=all_models,
                              filter=suffix('.pc_transcripts.json'),
                              output='.uprr.bed',
                              extras=['uprr']).follows(mkgenemap)


    cmd = config.get('Pipeline', 'runall')
    runall = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                        name='runall',
                        input=output_from(init, sradump,
                                          hsagencode, mmugencode,
                                          mrgens, normbtau, mktrans,
                                          mkhsamap, mkmmumap, mkgenemap,
                                          hsaidx31, mmuidx13, mmuidx19, genidx19,
                                          qmmuse, qhsape, qmmupe, qsscpe, qbtape,
                                          cvbedhsa, cvbedmmu, cvbedbta, cvbedssc,
                                          dumpbody, dumpcore, dumpuprr),
                        output=os.path.join(tempdir, 'runall_encexp.chk'),
                        extras=[cmd, jobcall])

    orth_files = os.listdir(orthdir)
    orth_files = [os.path.join(orthdir, f) for f in orth_files]

    initorth = pipe.originate(lambda x: x, orth_files, name='initorth')

    convorth = pipe.transform(task_func=convert_hcop_table,
                              name='convorth',
                              input=output_from(initorth),
                              filter=suffix('_six_column.txt.gz'),
                              output='_orthologs.tsv',
                              output_dir=orthdir).jobs_limit(2)







    return pipe
