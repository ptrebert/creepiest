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

from crplib.metadata.md_featdata import gen_obj_and_md, MD_FEATDATA_COLDEFS
from crplib.auxiliary.hdf_ops import load_masked_sigtrack, get_valid_hdf5_groups, \
    get_mapindex_groups, get_valid_chrom_group, get_default_group
from crplib.auxiliary.file_ops import create_filepath, check_array_serializable, load_mmap_array
from crplib.auxiliary.seq_parsers import get_twobit_seq, add_seq_regions
from crplib.mlfeat.featdef import feat_mapsig, feat_tf_motifs,\
    get_online_version, check_online_available, feat_roi, verify_sample_integrity
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
    index_groups = get_mapindex_groups(params['mapfile'], params['mapreference'])
    with pd.HDFStore(params['mapfile'], 'r') as idx:
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


def add_signal_features(samples, params):
    """
    :param samples:
    :param params:
    :return:
    """
    siglabels = params['siglabels']
    sigfiles = params['sigfiles']
    siggroups = params['siggroups']
    chrom = params['chrom']
    index_groups = get_mapindex_groups(params['mapfile'], params['mapreference'])
    with pd.HDFStore(params['mapfile'], 'r') as idx:
        mask = idx[index_groups[chrom]['mask']]
    mapfeat = feat_mapsig
    for sgf in sigfiles:
        this_group = siggroups[sgf]
        this_label = siglabels[sgf]
        signal = load_masked_sigtrack(sgf, '', this_group, chrom, mask=mask)
        for smp in samples:
            smp.update(mapfeat(signal[smp['start']:smp['end']], this_label))
    return samples


def add_roi_features(samples, params):
    """
    :param samples:
    :param params:
    :return:
    """
    roilabels = params['roilabels']
    roifiles = params['roifiles']
    roigroups = params['roigroups']
    chrom = params['chrom']
    for rf in roifiles:
        with pd.HDFStore(rf, 'r') as hdf:
            grp = roigroups[rf]
            if not grp.endswith(chrom):
                grp = os.path.join(grp, chrom)
            rois = hdf[grp]
            samples = feat_roi(samples, rois, roilabels[rf], params['roiquant'])
    return samples


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
    verify_sample_integrity(pos_samples, True)
    verify_sample_integrity(neg_samples, True)
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
        pos_samples = add_signal_features(pos_samples, params)
        neg_samples = add_signal_features(neg_samples, params)
    if 'roi' in params['features']:
        pos_samples = add_roi_features(pos_samples, params)
        neg_samples = add_roi_features(neg_samples, params)
    pos_samples = rebuild_dataframe(pos_samples)
    neg_samples = rebuild_dataframe(neg_samples)
    posgrp = params['posoutgroup']
    neggrp = params['negoutgroup']
    return chrom, pos_samples, posgrp, neg_samples, neggrp


def prep_scan_regions(params):
    """
    :param params:
    :return:
    """
    chrom = params['chrom']
    with pd.HDFStore(params['inputfile'], 'r') as hdfin:
        samples = hdfin[params['inputgroup']]
        cleanup = [c for c in samples.columns if c.startswith('ft')]
        samples.drop(cleanup, axis='columns', inplace=True)
        if params['window'] > 0:
            assert 'name' in samples.columns,\
                'Need name column to perform sliding window scan - file {} / chrom {}'.format(params['inputfile'], chrom)
            assert samples.name.unique().size == samples.name.size,\
                'Region names have to be unique per chromosome - file {} / chrom {}'.format(params['inputfile'], chrom)
        samples = samples.to_dict('record')
    if params['addseq']:
        assert os.path.isfile(params['seqfile']), 'Invalid path to sequence file: {}'.format(params['seqfile'])
        samples = add_seq_regions(samples, params['seqfile'], chrom)
    if params['window'] > 0:
        assert params['stepsize'] > 0,\
            'Step size invalid, cannot proceed: {} (window: {})'.format(params['step'], params['window'])
        samples = split_regions(samples, params['window'], params['stepsize'])
    verify_sample_integrity(samples, True)
    if len(check_online_available(params['features'])) > 0:
        compfeat = get_online_version(params['features'], params['kmers'])
        samples = list(map(compfeat, samples))
    if 'tfm' in params['features']:
        assert os.path.isfile(params['tfmotifs']), 'Invalid path to TF motifs file specified: {}'.format(params['tfmotifs'])
        with pd.HDFStore(params['tfmotifs'], 'r') as tffile:
            chrgroup = get_valid_chrom_group(params['tfmotifs'], chrom)
            tfdata = tffile[chrgroup]
            samples = feat_tf_motifs(samples, tfdata)
    if 'msig' in params['features']:
        samples = add_signal_features(samples, params)
    if 'roi' in params['features']:
        samples = add_roi_features(samples, params)
    samples = rebuild_dataframe(samples)
    # DEBUG free some bytes in the dataframe
    samples.drop(['seq'], axis='columns', inplace=True)
    samples = check_array_serializable(samples)
    outgroup = params['outputgroup']
    return chrom, samples, outgroup


def split_regions(samples, window, stepsize):
    """
    :param samples:
    :param window:
    :param stepsize:
    :return:
    """
    subsamples = []
    # this is just to avoid running over the end
    # if "gapped" windows are required at some point,
    # need to check subsamples
    assert window >= stepsize, 'Step sizes larger than window not supported: {} < {}'.format(window, stepsize)
    for smp in samples:
        added = False
        e = smp['end'] // window * window
        seq = smp['seq']
        offset = smp['start']
        # Note to self here:
        # if e <= smp['start'], won't fire
        # added is False, add entire sample
        for idx in range(smp['start'], e, stepsize):
            if idx + window > smp['end']:
                break
            subsamples.append({'source': smp['name'], 'start': idx, 'end': idx + window,
                               'seq': seq[idx - offset:idx + window - offset]})
            added = True
        if not added:
            # it could be that for very small regions the window size is too large,
            # breaking the loop immediately and thus skipping the sample
            subsamples.append({'source': smp['name'], 'start': smp['start'], 'end': smp['end'],
                               'seq': seq})
    return subsamples


def rebuild_dataframe(samples):
    """
    :param samples:
    :return:
    """
    assert samples, 'Received empty list of samples'
    if 'index' in samples[0]:
        samples = sorted(samples, key=lambda d: d['index'])
        index = [d['index'] for d in samples]
        df = pd.DataFrame.from_dict(samples, orient='columns')
        df.index = index
        df.drop('index', axis='columns', inplace=True)
    else:
        samples = sorted(samples, key=lambda d: (d['start'], d['end']))
        df = pd.DataFrame.from_dict(samples, orient='columns')
    return df


def build_featurefile_info(files):
    """
    :param files:
    :return:
    """
    labels = dict()
    groups = dict()
    featfiles = []
    for f in files:
        label, grp, fp = f.split(':')
        assert os.path.isfile(fp), 'Path to file {} invalid'.format(fp)
        if not grp or grp in ['default', 'auto']:
            grp = get_default_group(fp)
        label = label.strip('_')
        labels[fp] = label
        groups[fp] = grp
        featfiles.append(fp)
    return featfiles, labels, groups


def assemble_regsig_args(chromlim, args):
    """
    :param chromlim:
    :param args:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfile
    commons['mapfile'] = args.mapfile
    commons['mapreference'] = args.mapreference
    commons['seqfile'] = args.seqfile
    commons['group'] = args.inputgroup
    commons['resolution'] = args.resolution
    commons['features'] = args.features
    commons['kmers'] = args.kmers

    index_groups = get_mapindex_groups(args.mapfile,args.mapreference)

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


def prepare_regsig_samples(args, chromlim, logger):
    """
    :param args:
    :param chromlim:
    :param logger:
    :return:
    """
    arglist = assemble_regsig_args(chromlim, args)
    logger.debug('Collecting training data')
    _ = create_filepath(args.outputfile, logger)
    with pd.HDFStore(args.outputfile, args.filemode, complevel=9, complib='blosc', encoding='utf-8') as hdfout:
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(sample_signal_traindata, arglist)
            metadata = pd.DataFrame(columns=MD_FEATDATA_COLDEFS)
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
    sigfiles, siglabels, siggroups = build_featurefile_info(args.sigfile)
    commons['siglabels'] = siglabels
    commons['sigfiles'] = sigfiles
    commons['siggroups'] = siggroups
    roifiles, roilabels, roigroups = build_featurefile_info(args.roifile)
    commons['roilabels'] = roilabels
    commons['roifiles'] = roifiles
    commons['roigroups'] = roigroups
    commons['roiquant'] = args.roiquant
    commons['mapfile'] = args.mapfile
    commons['mapreference'] = args.mapreference
    check = re.compile(args.selectchroms)
    posgroups = get_valid_hdf5_groups(args.inputfile, args.posingroup)
    neggroups = get_valid_hdf5_groups(args.inputfile, args.negingroup)
    arglist = []
    for grp in posgroups:
        prefix, chrom = os.path.split(grp)
        if check.match(chrom) is None:
            logger.debug('Skipping group {}'.format(chrom))
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


def assemble_scnreg_args(args, logger):
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
    commons['window'] = args.window
    commons['stepsize'] = args.stepsize
    sigfiles, siglabels, siggroups = build_featurefile_info(args.sigfile)
    commons['siglabels'] = siglabels
    commons['sigfiles'] = sigfiles
    commons['siggroups'] = siggroups
    commons['mapfile'] = args.mapfile
    commons['mapreference'] = args.mapreference
    roifiles, roilabels, roigroups = build_featurefile_info(args.roifile)
    commons['roilabels'] = roilabels
    commons['roifiles'] = roifiles
    commons['roigroups'] = roigroups
    commons['roiquant'] = args.roiquant
    check = re.compile(args.selectchroms)
    ingroups = get_valid_hdf5_groups(args.inputfile, args.inputgroup)
    arglist = []
    for grp in ingroups:
        prefix, chrom = os.path.split(grp)
        if check.match(chrom) is None:
            logger.debug('Skipping group {}'.format(chrom))
            continue
        tmp = dict(commons)
        tmp['inputgroup'] = grp
        tmp['outputgroup'] = os.path.join(args.outputgroup, chrom)
        tmp['chrom'] = chrom
        arglist.append(tmp)
    return arglist


def prepare_clsreg_samples(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    arglist = assemble_clsreg_args(args, logger)
    logger.debug('Argument list of size {} to process'.format(len(arglist)))
    _ = create_filepath(args.outputfile, logger)
    with pd.HDFStore(args.outputfile, args.filemode, complib='blosc', complevel=9, encoding='utf-8') as hdfout:
        metadata = pd.DataFrame(columns=MD_FEATDATA_COLDEFS)
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


def prepare_scnreg_samples(args, logger):
    """
    :param args:
    :param logger:
    :return:
    """
    arglist = assemble_scnreg_args(args, logger)
    logger.debug('Argument list of size {} to process'.format(len(arglist)))
    with pd.HDFStore(args.outputfile, args.filemode, complib='blosc', complevel=9, encoding='utf-8') as hdfout:
        metadata = pd.DataFrame(columns=MD_FEATDATA_COLDEFS)
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(prep_scan_regions, arglist, chunksize=1)
            for chrom, samples, group in resit:
                logger.debug('Received data for chromosome {}'.format(chrom))
                if isinstance(samples, tuple):
                    logger.debug('Large dataset for chromosome {}, reading from file buffer'.format(chrom))
                    tmpfn = samples[0]
                    samples, rm_tmp = load_mmap_array(tmpfn, 'pandas')
                    if not rm_tmp:
                        logger.warning('Could not remove temp buffer {}'.format(tmpfn))
                grp, dataobj, metadata = gen_obj_and_md(metadata, group, chrom, args, samples)
                hdfout.put(grp, dataobj, format='fixed')
                hdfout.flush()
                logger.debug('Data flushed to file')
        hdfout.put('metadata', metadata, format='table')
        hdfout.flush()
    logger.debug('Collecting training data done')
    return 0


def run_compute_features(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    setattr(args, 'selectchroms', args.selectchroms.strip('"'))
    _ = create_filepath(args.outputfile, logger)
    logger.debug('Chromosome select pattern: {}'.format(args.selectchroms))
    if args.task == 'regress':
        logger.debug('Computing features/sampling data for task {}'.format(args.task))
        # "magic number" following common limits, e.g., in ChromImpute
        chromlim = CHROMOSOME_BOUNDARY
        _ = prepare_regsig_samples(args, chromlim, logger)
    elif args.task == 'groups':
        logger.debug('Computing features for task {}'.format(args.task))
        assert args.posingroup and args.negingroup, 'Need to specify HDF groups for positive and negative class'
        _ = prepare_clsreg_samples(args, logger)
    elif args.task == 'classify':
        logger.debug('Computing region features for task: {}'.format(args.task))
        _ = prepare_scnreg_samples(args, logger)
    else:
        raise NotImplementedError('Task unknown: {}'.format(args.task))
    return 0
