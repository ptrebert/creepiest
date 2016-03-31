# coding=utf-8

"""
Module to collect training data
Data are stored in separate files in the form of Pandas dataframes
"""

import random as rand
import numpy as np
import multiprocessing as mp
import pandas as pd

from crplib.metadata.md_traindata import gen_obj_and_md, MD_TRAINDATA_COLDEFS
from crplib.auxiliary.file_ops import load_masked_sigtrack
from crplib.auxiliary.text_parsers import read_chromosome_sizes
from crplib.auxiliary.seq_parsers import get_twobit_seq
from crplib.mlfeat.featdef import feat_mapsig, get_online_version


def sample_signal_traindata(params):
    """
    :param params:
    :return:
    """
    mypid = mp.current_process().pid
    lolim = params['lolim']
    hilim = params['hilim']

    res = params['resolution']
    # TODO
    # ad-hoc values... currently, no strategy
    # justifying the selection...
    stepsize = res * 2
    step_bw = res * 4
    step_fw = res * 6
    # memory-wise large objects
    chromseq = get_twobit_seq(params['seqfile'], params['chrom'])
    signal = load_masked_sigtrack(params['inputfile'], params['chainfile'], params['group'],
                                  params['chrom'], params['size'])
    samples = []
    # make function available in local namespace
    mapfeat = feat_mapsig
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
    return mypid, params['chrom'], samples


def assemble_worker_args(chroms, chromlim, args):
    """
    :param chroms:
    :param chromlim:
    :param args:
    :return:
    """
    arglist = []
    commons = dict()
    commons['inputfile'] = args.inputfile
    commons['chainfile'] = args.chainfile
    commons['seqfile'] = args.seqfile
    commons['group'] = args.inputgroup
    commons['resolution'] = args.resolution
    commons['features'] = args.features
    commons['kmers'] = tuple(args.kmers)

    num_chrom = len(chroms)
    smp_per_chrom = args.numsamples // num_chrom
    rest = args.numsamples
    dist_rest = smp_per_chrom * num_chrom < args.numsamples
    for name, size in chroms.items():
        tmp = dict(commons)
        tmp['chrom'] = name
        tmp['size'] = size
        if dist_rest:
            # last chromosome loses a few sample points...
            tmp['numsamples'] = min(rest, smp_per_chrom + 1)
            rest -= min(rest, smp_per_chrom + 1)
        else:
            tmp['numsamples'] = smp_per_chrom
        tmp['lolim'] = chromlim
        tmp['hilim'] = size - chromlim
        arglist.append(tmp)
    return arglist


def collect_sigres_trainsamples(args, csizes, chromlim, logger):
    """
    :param args:
    :param csizes:
    :param chromlim:
    :param logger:
    :return:
    """

    arglist = assemble_worker_args(csizes, chromlim, args)

    with pd.HDFStore(args.outputfile, 'w', complevel=9, complib='blosc') as hdfout:
        with mp.Pool(args.workers) as pool:
            mapres = pool.map_async(sample_signal_traindata, arglist)
            metadata = pd.DataFrame(columns=MD_TRAINDATA_COLDEFS)
            for pid, chrom, samples in mapres.get():
                logger.debug('Process {} finished chromosome {}'.format(pid, chrom))
                grp, dataobj, metadata = gen_obj_and_md(metadata, args.grouproot, chrom, 'seq',
                                                        args.inputfile, args.chainfile, samples)
                hdfout.put(grp, dataobj, format='fixed')
                hdfout.flush()
            hdfout.put('metadata', metadata, format='table')
    return 0


def run_collect_traindata(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    args.__dict__['keepchroms'] = args.keepchroms.strip('"')
    logger.debug('Chromosome select pattern: {}'.format(args.keepchroms))
    if args.task == 'regsig':
        logger.debug('Collecting training data for task {}'.format(args.task))
        csizes = read_chromosome_sizes(args.chromsizes, args.keepchroms)
        # "magic number" following common limits, e.g., in ChromImpute
        chromlim = 10000
        _ = collect_sigres_trainsamples(args, csizes, chromlim, logger)
    else:
        raise NotImplementedError('Task unknown: {}'.format(args.task))
    return 0
