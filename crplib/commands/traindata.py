# coding=utf-8

"""
Module to collect training data
Data are stored in separate files in the form of Pandas dataframes
"""

import os as os
import random as rand
import numpy as np
import re as re
import multiprocessing as mp
import pandas as pd

from crplib.metadata.md_traindata import gen_obj_and_md, MD_TRAINDATA_COLDEFS
from crplib.auxiliary.hdf_ops import load_masked_sigtrack, get_valid_hdf5_groups, \
    get_trgindex_groups, get_valid_chrom_group
from crplib.auxiliary.text_parsers import read_chromosome_sizes
from crplib.auxiliary.file_ops import create_filepath
from crplib.auxiliary.seq_parsers import get_twobit_seq, add_seq_regions
from crplib.mlfeat.featdef import feat_mapsig, feat_tf_motifs,\
    get_online_version, check_online_available
from crplib.auxiliary.constants import CHROMOSOME_BOUNDARY


def sample_signal_traindata(params):
    """
    :param params:
    :return:
    """
    chrom = params['chrom']
    res = params['resolution']
    # TODO
    # ad-hoc values...
    stepsize = res * 2
    step_bw = res * 4
    step_fw = res * 6
    # memory-wise large objects
    index_groups = get_trgindex_groups(params['targetindex'], '')
    with pd.HDFStore(params['targetindex'], 'r') as idx:
        mask = idx[index_groups[chrom]['mask']]
    chromseq = get_twobit_seq(params['seqfile'], params['chrom'])
    signal = load_masked_sigtrack(params['inputfile'], '', params['group'],
                                  params['chrom'], '', mask=mask)
    samples = []
    # make function available in local namespace
    mapfeat = feat_mapsig
    lolim = params['lolim']
    hilim = len(chromseq) - lolim  # which is chromosome boundary constant
    for n in range(params['numsamples']):
        pos = rand.randint(lolim, hilim)
        for move in range(pos - step_bw, pos + step_fw, stepsize):
            seq_reg = {'sample_n': n, 'start': move, 'end': move + res}
            y_dep = np.average(signal[move:move + res].data)
            seq_reg['y_depvar'] = y_dep
            seq_reg['seq'] = chromseq[move:move + res]
            seq_reg.update(mapfeat(signal[move:move + res]))
            samples.append(seq_reg)
    comp_seqfeat = get_online_version(params['features'], params['kmers'])
    samples = list(map(comp_seqfeat, samples))
    return chrom, samples


def get_region_traindata(params):
    """
    :param params:
    :return:
    """
    chrom = params['chrom']
    with pd.HDFStore(params['inputfile'], 'r') as hdfin:
        neg_samples = hdfin[params['negingroup']]
        pos_samples = hdfin[params['posingroup']].iloc[neg_samples.index, ]
        neg_samples = neg_samples.assign(index=neg_samples.index)
        neg_samples = neg_samples.assign(y_depvar=np.zeros(neg_samples.shape[0]))
        cleanup = [c for c in neg_samples.columns if c.startswith('ft')]
        neg_samples.drop(cleanup, axis='columns', inplace=True)
        pos_samples = pos_samples.assign(index=pos_samples.index)
        pos_samples = pos_samples.assign(y_depvar=np.ones(pos_samples.shape[0]))
        cleanup = [c for c in pos_samples.columns if c.startswith('ft')]
        pos_samples.drop(cleanup, axis='columns', inplace=True)
        pos_samples = pos_samples.to_dict('record')
        neg_samples = neg_samples.to_dict('record')
    if params['addseq']:
        assert os.path.isfile(params['seqfile']), 'Invalid path to sequence file: {}'.format(params['seqfile'])
        pos_samples = add_seq_regions(pos_samples, params['seqfile'], chrom)
        neg_samples = add_seq_regions(neg_samples, params['seqfile'], chrom)
    if len(check_online_available(params['features'])) > 0:
        compfeat = get_online_version(params['features'], params['kmers'])
        pos_samples = list(map(compfeat, pos_samples))
        neg_samples = list(map(compfeat, neg_samples))
    if 'tfm' in params['features']:
        assert os.path.isfile(params['tfmotifs']), 'Invalid path to TF motifs file specified: {}'.format(params['tfmotifs'])
        with pd.HDFStore(params['tfmotifs'], 'r') as tffile:
            chrgroup = get_valid_chrom_group(params['tfmotifs'], chrom)
            tfdata = tffile[chrgroup]
            pos_samples = feat_tf_motifs(pos_samples, tfdata)
            neg_samples = feat_tf_motifs(neg_samples, tfdata)
    if 'msig' in params['features']:
        assert params['signallabel'], 'Need a signal label for feature msig'
        label = params['signallabel']
        index_groups = get_trgindex_groups(params['targetindex'], '')
        with pd.HDFStore(params['targetindex'], 'r') as idx:
            mask = idx[index_groups[chrom]['mask']]
        if not params['signalgroup']:
            signal = load_masked_sigtrack(params['signalfile'], '', '', chrom, mask=mask)
        else:
            signal = load_masked_sigtrack(params['signalfile'], '', params['signalgroup'], chrom, mask=mask)
        mapfeat = feat_mapsig
        for pos, neg in zip(pos_samples, neg_samples):
            pos.update(mapfeat(signal[pos['start']:pos['end']], label))
            neg.update(mapfeat(signal[neg['start']:neg['end']], label))
    pos_samples = rebuild_dataframe(pos_samples)
    neg_samples = rebuild_dataframe(neg_samples)
    posgrp = params['posoutgroup']
    neggrp = params['negoutgroup']
    return chrom, pos_samples, posgrp, neg_samples, neggrp


def rebuild_dataframe(samples):
    """
    :param samples:
    :return:
    """
    assert samples, 'Received empty list of samples'
    samples = sorted(samples, key=lambda d: d['index'])
    index = [d['index'] for d in samples]
    df = pd.DataFrame.from_dict(samples, orient='columns')
    df.index = index
    df.drop('index', axis='columns', inplace=True)
    return df


def assemble_regsig_args(chromlim, args):
    """
    :param chroms:
    :param chromlim:
    :param args:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfile
    commons['targetindex'] = args.targetindex
    commons['seqfile'] = args.seqfile
    commons['group'] = args.inputgroup
    commons['resolution'] = args.resolution
    commons['features'] = args.features
    commons['kmers'] = args.kmers

    index_groups = get_trgindex_groups(args.targetindex, '')

    num_chrom = len(index_groups.keys())
    smp_per_chrom = args.numsamples // num_chrom
    rest = args.numsamples
    dist_rest = smp_per_chrom * num_chrom < args.numsamples
    for name in index_groups.keys():
        tmp = dict(commons)
        tmp['chrom'] = name
        if dist_rest:
            # last chromosome loses a few sample points...
            tmp['numsamples'] = min(rest, smp_per_chrom + 1)
            rest -= min(rest, smp_per_chrom + 1)
        else:
            tmp['numsamples'] = smp_per_chrom
        tmp['lolim'] = chromlim
        arglist.append(tmp)
    return arglist


def collect_regsig_trainsamples(args, chromlim, logger):
    """
    :param args:
    :param csizes:
    :param chromlim:
    :param logger:
    :return:
    """
    arglist = assemble_regsig_args(chromlim, args)
    logger.debug('Collecting training data')
    with pd.HDFStore(args.outputfile, 'w', complevel=9, complib='blosc', encoding='utf-8') as hdfout:
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(sample_signal_traindata, arglist)
            metadata = pd.DataFrame(columns=MD_TRAINDATA_COLDEFS)
            for chrom, samples in resit:
                logger.debug('Processed chromosome {}'.format(chrom))
                grp, dataobj, metadata = gen_obj_and_md(metadata, args.outputgroup, chrom, args, samples)
                hdfout.put(grp, dataobj, format='fixed')
                hdfout.flush()
            hdfout.put('metadata', metadata, format='table')
    return 0


def assemble_clsreg_args(args, logger):
    """
    :param args:
    :return:
    """
    commons = dict()
    commons['inputfile'] = args.inputfile
    commons['addseq'] = args.addseq
    commons['seqfile'] = args.seqfile
    commons['features'] = args.features
    commons['kmers'] = args.kmers
    commons['tfmotifs'] = args.tfmotifs
    commons['signalfile'] = args.signalfile
    commons['signalgroup'] = args.signalgroup
    commons['signallabel'] = args.signallabel
    commons['targetindex'] = args.targetindex
    check = re.compile(args.keepchroms)
    posgroups = get_valid_hdf5_groups(args.inputfile, args.posingroup)
    neggroups = get_valid_hdf5_groups(args.inputfile, args.negingroup)
    arglist = []
    for grp in posgroups:
        prefix, chrom = os.path.split(grp)
        if check.match(chrom) is None:
            logger.debug('Skipping chromosome {}'.format(chrom))
            continue
        neggrp = list(filter(lambda x: x.endswith(chrom), neggroups))
        assert len(neggrp) == 1, 'Cannot find negative group for positive {}'.format(grp)
        neggrp = neggrp[0]
        tmp = dict(commons)
        tmp['posingroup'] = grp
        tmp['negingroup'] = neggrp
        tmp['posoutgroup'] = os.path.join(args.posoutgroup, chrom)
        tmp['negoutgroup'] = os.path.join(args.negoutgroup, chrom)
        tmp['chrom'] = chrom
        arglist.append(tmp)
    return arglist


def collect_clsreg_trainsamples(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    arglist = assemble_clsreg_args(args, logger)
    logger.debug('Argument list of size {} to process'.format(len(arglist)))
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9, encoding='utf-8') as hdfout:
        metadata = pd.DataFrame(columns=MD_TRAINDATA_COLDEFS)
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(get_region_traindata, arglist, chunksize=1)
            for chrom, pos_samples, posgrp, neg_samples, neggrp in resit:
                grp, dataobj, metadata = gen_obj_and_md(metadata, posgrp, chrom, args, pos_samples)
                hdfout.put(grp, dataobj, format='fixed')
                grp, dataobj, metadata = gen_obj_and_md(metadata, neggrp, chrom, args, neg_samples)
                hdfout.put(grp, dataobj, format='fixed')
                hdfout.flush()
                logger.debug('Processed chromosome {}'.format(chrom))
        hdfout.put('metadata', metadata, format='table')
        hdfout.flush()
    logger.debug('Collecting training data done')
    return 0


def run_collect_traindata(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    args.__dict__['keepchroms'] = args.keepchroms.strip('"')
    _ = create_filepath(args.outputfile, logger)
    logger.debug('Chromosome select pattern: {}'.format(args.keepchroms))
    if args.task == 'regsig':
        logger.debug('Collecting training data for task {}'.format(args.task))
        # "magic number" following common limits, e.g., in ChromImpute
        chromlim = CHROMOSOME_BOUNDARY
        _ = collect_regsig_trainsamples(args, chromlim, logger)
    elif args.task == 'clsreg':
        logger.debug('Collecting training data for task {}'.format(args.task))
        assert args.posingroup and args.negingroup, 'Need to specify HDF groups for positive and negative class'
        _ = collect_clsreg_trainsamples(args, logger)
    else:
        raise NotImplementedError('Task unknown: {}'.format(args.task))
    return 0
