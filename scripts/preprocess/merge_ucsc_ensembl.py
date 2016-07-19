#!/usr/bin/env python3

import sys as sys
import os as os
import csv as csv
import argparse as argp
import io as io
import gzip as gz


def parse_arguments():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--locations', type=str, dest='locations')
    parser.add_argument('--to-biotype', type=str, dest='tobiotype')
    parser.add_argument('--to-name', type=str, dest='toname')
    parser.add_argument('--out-gff', type=str, dest='outgff')
    args = parser.parse_args()
    return args


def read_mapping(fp, header, select_field='', select_value=''):
    """
    :param fp:
    :param header:
    :return:
    """
    content = []
    with gz.open(fp, 'rt') as mapping:
        rows = csv.DictReader(mapping, fieldnames=header, delimiter='\t')
        for r in rows:
            if select_field:
                if r[select_field] != select_value:
                    continue
            content.append(r)
    return content


def read_ucsc_dump(fp):
    """
    :param fp:
    :return:
    """
    header = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
              'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    content = []
    with gz.open(fp, 'rt') as dump:
        rows = csv.DictReader(dump, fieldnames=header, delimiter='\t')
        for r in rows:
            # filter for ...?
            content.append(r)
    return content


if __name__ == '__main__':
    args = parse_arguments()
    base_info = read_ucsc_dump(args.locations)
    base_info = sorted(base_info, key=lambda x: (x['chrom'], int(x['txStart']), int(x['txEnd'])))
    biotype_map = read_mapping(args.tobiotype, ['transcript_name', 'biotype'], 'biotype', 'protein_coding')
    all_pc_transcripts = set(d['transcript_name'] for d in biotype_map)
    #print('Protein coding transcripts: {}'.format(len(all_pc_transcripts)))
    name_map = read_mapping(args.toname, ['transcript_name', 'gene_name'])
    name_map = dict([(d['transcript_name'], d['gene_name']) for d in name_map if d['transcript_name'] in all_pc_transcripts])
    #print('Named entities: {}'.format(len(name_map)))

    genes = dict()
    for entry in base_info:
        if entry['name'] not in all_pc_transcripts:
            continue
        gene_id = entry['name2']
        if entry['name'] in name_map:
            symbol = name_map[entry['name']]
        else:
            symbol = 'no_symbol'
        if gene_id in genes:
            genes[gene_id]['start'].append(int(entry['txStart']))
            genes[gene_id]['end'].append(int(entry['txEnd']))
            genes[gene_id]['transcripts'].append(entry['name'])
            genes[gene_id]['symbol'].append(symbol)
        else:
            new_gene = {'start': [int(entry['txStart'])], 'end': [int(entry['txEnd'])],
                        'symbol': [symbol], 'transcripts': [entry['name']], 'chrom': entry['chrom'],
                        'strand': entry['strand'], 'id': gene_id}
            genes[gene_id] = new_gene

    outbuffer = io.StringIO()
    base_info = dict([(d['name'], d) for d in base_info])
    genes = sorted(genes.values(), key=lambda x: (x['chrom'], min(x['start']), max(x['end'])))
    skipped_dups = 0
    for g in genes:
        s, e = str(min(g['start'])), str(max(g['end']))
        line = '\t'.join([g['chrom'], 'Ensembl', 'gene', s, e, '.', g['strand'], '.', ''])
        assert len(set(g['symbol'])) == 1, 'Mixed symbol for gene: {}'.format(g)
        symbol = set(g['symbol']).pop()
        attributes = 'ID={};Gene_name={}'.format(g['id'], symbol)
        line += attributes + '\n'
        _ = outbuffer.write(line)
        dups = set()
        for t in sorted(g['transcripts'], reverse=True):
            # start with highest/most recent transcript number
            tinfo = base_info[t]
            pos = tinfo['chrom'], tinfo['txStart'], tinfo['txEnd']
            if pos in dups:
                skipped_dups += 1
                continue
            dups.add(pos)
            line = '\t'.join([tinfo['chrom'], 'Ensembl', 'transcript', tinfo['txStart'],
                              tinfo['txEnd'], '.', tinfo['strand'], '.', ''])
            attributes = 'ID={};Parent={};Name={};Transcript_name={}'.format(tinfo['name'], g['id'], tinfo['name'], name_map.get(tinfo['name'], 'no_name'))
            line += attributes + '\n'
            _ = outbuffer.write(line)

    #print('Skipped positional duplicates: {}'.format(skipped_dups))
    with gz.open(args.outgff, 'wt') as outf:
        _ = outf.write(outbuffer.getvalue())

    sys.exit(0)


