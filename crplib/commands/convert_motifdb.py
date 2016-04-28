# coding=utf-8

"""
Module to convert TF motif databases
and associated binned motif counts to HDF5
"""

import pandas as pd
import numpy as np
import operator as op
import multiprocessing as mp

from crplib.auxiliary.text_parsers import text_file_mode, get_meme_iterator, get_binned_motifs_iterator
from crplib.mlfeat.featdef import FEAT_TFMOTIF_PREFIX
from crplib.metadata.md_helpers import normalize_group_path
from crplib.auxiliary.constants import DIV_B_TO_MB


def build_motifmap(fpath, dbfmt):
    """
    :param fpath:
    :param dbfmt:
    :return:
    """
    opn, mode = text_file_mode(fpath)
    if dbfmt == 'list':
        with opn(fpath, mode) as infile:
            motifs = infile.read().strip().split()
            motifmap = dict((m, m) for m in motifs)
    elif dbfmt == 'map':
        motifmap = dict()
        with opn(fpath, mode) as infile:
            for line in infile:
                if not line:
                    continue
                k, v = line.strip().split()
                motifmap[k] = v
    elif dbfmt == 'meme':
        motifmap = dict()
        with opn(fpath, mode) as infile:
            memit = get_meme_iterator(infile)
            for tfid, tfname in memit:
                motifmap[tfid] = tfname
    else:
        # NB: this would raise if another choice is added to the subparser w/o
        # changing the code here
        raise ValueError('Unknown format for motif DB: {}'.format(dbfmt))
    return motifmap


def process_merged_file(params):
    """
    :param params:
    :return:
    """
    motifmap = params['motifmap']
    getcounts = op.itemgetter(*tuple(sorted(motifmap.keys())))
    opn, mode = text_file_mode(params['inputfile'])
    with opn(params['inputfile'], mode) as infile:
        tfmit = get_binned_motifs_iterator(infile)
        chrom_counts = []
        chrom_indices = []
        curr_chrom = None
        for chrom, index, counts in tfmit:
            if chrom != curr_chrom and curr_chrom is not None:
                yield curr_chrom, chrom_indices, chrom_counts
                curr_chrom = chrom
                chrom_indices = [index]
                chrom_counts = [np.array(getcounts(counts), dtype=np.int32)]
                continue
            chrom_indices.append(index)
            chrom_counts.append(np.array(getcounts(counts), dtype=np.int32))
        yield curr_chrom, np.array(chrom_indices, dtype=np.int64), np.array(chrom_counts)
    return


def process_split_files(params):
    """
    :param params:
    :return:
    """
    motifmap = params['motifmap']
    getcounts = op.itemgetter(*tuple(sorted(motifmap.keys())))
    opn, mode = text_file_mode(params['inputfile'])
    with opn(params['inputfile'], mode) as infile:
        tfmit = get_binned_motifs_iterator(infile)
        chrom_counts = []
        chrom_indices = []
        mychrom = None
        for chrom, index, counts in tfmit:
            mychrom = chrom
            chrom_indices.append(index)
            chrom_counts.append(np.array(getcounts(counts), dtype=np.int32))
    return mychrom, np.array(chrom_indices, dtype=np.int64), np.array(chrom_counts)


def build_dataframe_colnames(motifmap, fmt):
    """
    :param motifmap:
    :return:
    """
    motif_order = sorted(motifmap.keys())
    colnames = []
    if fmt == 'list':
        for motif in motif_order:
            colnames.append(FEAT_TFMOTIF_PREFIX + motif)
    else:
        for motif in motif_order:
            colnames.append(FEAT_TFMOTIF_PREFIX + motif + '-' + motifmap[motif])
    return colnames


def run_motifdb_conversion(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    logger.debug('Building name/ID map')
    motifmap = build_motifmap(args.motifdb, args.dbformat)
    logger.debug('Identified {} motifs in DB file'.format(len(motifmap)))
    colnames = build_dataframe_colnames(motifmap, args.dbformat)
    arglist = [{'inputfile': inpf, 'motifmap': motifmap} for inpf in args.inputfile]
    # TODO unify the approaches...
    logger.debug('Start processing...')
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdfout:
        chroms_seen = set()
        with mp.Pool(args.workers) as pool:
            if len(arglist) == 1:
                resit = process_merged_file(arglist[0])
            else:
                resit = pool.imap_unordered(process_split_files, arglist, chunksize=1)
            for chrom, indices, counts in resit:
                assert chrom not in chroms_seen, 'Encountered {} twice in dataset'.format(chrom)
                chroms_seen.add(chrom)
                logger.debug('Processed chromosome {}'.format(chrom))
                df = pd.DataFrame(counts, index=indices, columns=colnames, dtype=np.int32)
                group = normalize_group_path(args.outputgroup, chrom)
                hdfout.put(group, df, format='fixed')
                hdfout.flush()
                logger.debug('Data saved')
    return 0
