#!/usr/bin/env python3
# coding=utf-8

"""
Isolated script to parse ENCODE metadata file, select relevant information
and create a file listing for batch download of data
"""

import os as os
import sys as sys
import csv as csv
import traceback as trb


# since this is a stand-alone script, simple config right here

BIOSAMPLES = ['HepG2', 'K562', 'GM12878', 'H1-hESC', 'MEL cell line', 'CH12.LX', 'ES-Bruce4', 'ES-E14']
FORMATS = ['bigWig', 'gtf', 'bigBed narrowPeak', 'bigBed broadPeak']
OUTPUTS = ['peaks', 'signal', 'raw signal', 'base overlap signal',
           'exon quantifications', 'transcript quantifications', 'gene quantifications']

DL_FOLDER = '/TL/deep/fhgfs/projects/pebert/thesis/biodata/dlfolder/encode'

MDFILE = '/home/pebert/work/code/mpggit/creepiest/datasrc/encode/20160211_ENCODE_metadata.tsv'
LISTFILE = '/TL/deep/fhgfs/projects/pebert/thesis/biodata/dlfolder/listing_encode.txt'

# end of config


if __name__ == '__main__':
    try:
        existing_files = os.listdir(DL_FOLDER)
        existing_acc = set([f.split('.')[0] for f in existing_files])

        dlfiles = []
        with open(MDFILE, 'r') as infile:
            rowdata = csv.DictReader(infile, delimiter='\t')
            for row in rowdata:
                # ignore all data treated with some... treatment...
                if not row['Biosample treatments'].strip():
                    acc = row['File accession']
                    if acc not in existing_acc:
                        smp, form, out = row['Biosample term name'], row['File format'], row['Output type']
                        if smp in BIOSAMPLES and form in FORMATS and out in OUTPUTS:
                            dlfiles.append(row['File download URL'])

        if dlfiles:
            print('Writing file listing...')
            with open(LISTFILE, 'w') as outfile:
                _ = outfile.write('\n'.join(dlfiles) + '\n')
        else:
            print('Nothing to do')

    except Exception as err:
        trb.print_exc()
        sys.exit(1)
    else:
        sys.exit(0)
