#!/usr/bin/env python3

import sys as sys
import io as io
import json as json
import gzip as gz
import argparse as argp
import re as re
import twobitreader as tbr


def determine_open_mode(fp):
    if fp.endswith('.gz'):
        return gz.open, 'rt'
    else:
        return open, 'r'


def get_seqrc(sequence):
    """
    :param sequence:
    :return:
    """
    rcmap = {'A': 'T', 'a': 't',
             'C': 'G', 'c': 'g',
             'G': 'C', 'g': 'c',
             'T': 'A', 't': 'a',
             'N': 'N', 'n': 'n'}
    return ''.join(map(lambda x: rcmap[x], reversed(list(sequence))))


def parse_arguments():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--assembly', type=str, dest='assembly', required=True)
    parser.add_argument('--seq-file', type=str, dest='seqfile', required=True)
    parser.add_argument('--annotation', type=str, dest='annotation', required=True)
    parser.add_argument('--out-model', type=str, dest='outmodel')
    parser.add_argument('--out-fasta', type=str, dest='outfasta')
    parser.add_argument('--coord-base', type=int, default=0, dest='coordbase')
    args = parser.parse_args()
    return args


def parse_annotation_line(line):
    """
    :param line:
    :return:
    """
    chrom, source, biotype, start, end, _, strand, _, attributes = line.strip().split('\t')
    info = dict()
    for entry in attributes.split(';'):
        k, v = entry.split('=')
        if k.lower() == 'dbxref':
            try:
                k, v = v.split(':')
            except ValueError:
                try:
                    k, v = v.split(' ')
                except ValueError:
                    pass
        info[k.strip()] = v.strip()
    return source, chrom, int(start), int(end), strand, biotype, info


def parse_gff(fp, source='.+'):
    """
    :param fp:
    :return:
    """
    opn, mode = determine_open_mode(fp)
    genemodel = dict()
    chrom_match = re.compile('chr[0-9XYZW]+(\s|$)')
    source_match = re.compile(source)
    num_trans = 0
    with opn(fp, mode, encoding='ascii') as infile:
        for line in infile:
            if not line or line.startswith('#'):
                continue
            src, c, s, e, strand, biotype, attr = parse_annotation_line(line)
            if source_match.match(src) is None:
                continue
            if biotype not in ['gene', 'mRNA', 'transcript']:
                continue
            if chrom_match.match(c) is None:
                continue
            assert s < e, 'Malformed entry: {}'.format(line.strip())
            if biotype == 'gene':
                this_gene = {'chrom': c, 'start': s, 'end': e, 'strand': strand, 'transcripts': []}
                if 'symbol_ncbi' in attr:
                    this_gene['symbol'] = attr['symbol_ncbi']
                    this_gene['authority'] = 'NCBI'
                    this_gene['id'] = attr['NCBI_gene']
                elif 'ID' in attr:
                    if 'Gene_name' in attr:
                        this_gene['symbol'] = attr['Gene_name']
                    else:
                        this_gene['symbol'] = 'no_symbol'
                    this_gene['authority'] = 'Ensembl'
                    this_gene['id'] = attr['ID']
                    assert this_gene['id'].startswith('ENS'), 'This is not an Ensembl identifier: {}'.format(attr['ID'])
                else:
                    raise ValueError('Cannot identify gene name/id: {}'.format(line))
                genemodel[this_gene['id']] = this_gene
            else:  # biotype is mRNA/transcript
                if 'GeneID' in attr:
                    gene_id = attr['GeneID']
                elif 'Parent' in attr:
                    gene_id = attr['Parent']
                else:
                    raise ValueError('Cannot identify gene for transcript: {}'.format(attr))
                this_trans = {'chrom': c, 'start': s, 'end': e, 'id': attr['Name']}
                trans_id = this_trans['id']
                if trans_id.startswith('ENS'):
                    this_trans['authority'] = 'Ensembl'
                    if 'Transcript_name' in attr:
                        this_trans['name'] = attr['Transcript_name']
                    else:
                        this_trans['name'] = 'no_name'
                else:
                    this_trans['authority'] = 'RefSeq'
                genemodel[gene_id]['transcripts'].append(this_trans)
                num_trans += 1
    genemodel = sorted(genemodel.values(), key=lambda x: (x['chrom'], x['start'], x['end']))
    return genemodel


def check_need_rc(seq, is_minus):
    """
    :param seq:
    :param is_minus:
    :return:
    """
    if is_minus:
        return get_seqrc(seq)
    return seq


def build_transcriptome(genemodel, coordbase, seqfile, outpath):
    """
    :param genemodel:
    :param seqfile:
    :param outpath:
    :return:
    """
    outbuffer = io.StringIO()
    current_chrom = None
    chrom_seq = ''
    wg_seq = tbr.TwoBitFile(seqfile)
    for gene in genemodel:
        if gene['chrom'] != current_chrom:
            current_chrom = gene['chrom']
            chrom_seq = wg_seq[current_chrom]
        minus = gene['strand'] == '-'
        for trans in gene['transcripts']:
            _ = outbuffer.write('>{}|{}'.format(trans['id'], gene['id']) + '\n')
            # NB: GTF/GFF coordinates 1-based, sequence 0-based
            trans_seq = check_need_rc(chrom_seq[trans['start'] - coordbase:trans['end']], minus)
            tlen = trans['end'] - trans['start']
            assert len(trans_seq) == tlen, 'Length mismatch: {} - {} vs {}'.format(trans, tlen, len(trans_seq))
            for idx in range(0, len(trans_seq), 80):
                _ = outbuffer.write(trans_seq[idx:idx+80] + '\n')

    with gz.open(outpath, 'wt') as out:
        _ = out.write(outbuffer.getvalue())

    return 0


if __name__ == '__main__':
    args = parse_arguments()
    source_select = {}
    genemodel = parse_gff(args.annotation, source_select.get(args.assembly, '.+'))
    _ = build_transcriptome(genemodel, args.coordbase, args.seqfile, args.outfasta)

    if not args.outmodel:
        outmodel = args.outfasta.replace('.fa.gz', '.json')
    else:
        outmodel = args.outmodel

    with open(outmodel, 'w') as outf:
        _ = json.dump(genemodel, outf, indent=1, ensure_ascii=True, check_circular=True)

    sys.exit(0)
