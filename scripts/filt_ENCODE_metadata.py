#!/usr/bin/env python3
# coding=utf-8

"""
Stand-alone script to parse the raw ENCODE metadata and create
reasonable sample annotation files
"""

import sys as sys
import operator as op

import collections as col

from pitlib.readers.rowdata import rows_to_buffer


ENCODE_GET = ['File accession', 'Output type', 'Biosample term name', 'Biosample term id',
              'Experiment target', 'Biosample sex', 'Assembly', 'File download URL', 'md5sum',
              'Paired with', 'Technical replicate', 'Biological replicate(s)']

ENCODE_FILETYPES = ['bigWig', 'bed broadPeak', 'bed narrowPeak']

ENCODE_OUTPUTS = ['peaks', 'signal', 'raw signal']

ENCODE_SAMPLES = {'hg19': ['GM12878', 'H1-hESC', 'K562', 'HepG2', 'HUVEC'],
                  'mm9': ['CH12.LX', 'ES-Bruce4', 'ES-E14', 'MEL cell line']}

ENCODE_TARGETS = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3', 'CTCF', 'EP300', 'POLR2A', 'Control']


if __name__ == '__main__':
    mdfile = sys.argv[1]
    outfile = sys.argv[2]
    rows = rows_to_buffer(mdfile)
    itemgetter = op.itemgetter(*tuple(ENCODE_GET))
    assembly = None
    select_samples = None
    header = '\t'.join(ENCODE_GET)
    buffer = []
    count_signals = col.defaultdict(col.Counter)
    count_peaks = col.defaultdict(col.Counter)
    for row in rows:
        if row['File format'] in ENCODE_FILETYPES and \
           row['Output type'] in ENCODE_OUTPUTS:
            if assembly is None:
                assembly = row['Assembly']
                select_samples = ENCODE_SAMPLES[assembly]
            if row['Biosample term name'] in select_samples and \
               any([row['Experiment target'].startswith(t + '-') for t in ENCODE_TARGETS]):
                if row['Output type'] == 'peaks':
                    count_peaks[row['Biosample term name']][row['Experiment target']] += 1
                else:
                    count_signals[row['Biosample term name']][row['Experiment target']] += 1
                infos = itemgetter(row)
                infos = '\t'.join(['NA' if not i else i for i in infos])
                buffer.append(infos)
    with open(outfile, 'w') as outf:
        _ = outf.write(header + '\n')
        _ = outf.write('\n'.join(buffer))
    sys.exit(0)
