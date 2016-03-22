# coding=utf-8

"""
Module to collect training data
Data are stored in separate files in the form of Pandas dataframes
Allows to train new/different learning algorithms
"""

import random as rand
import numpy as np
import multiprocessing as mp
import pandas as pd

from crplib.metadata.md_traindata import gen_obj_and_md, MD_TRAINDATA_COLDEFS
from crplib.auxiliary.file_ops import text_file_mode
from crplib.auxiliary.text_parsers import read_chromosome_sizes
from crplib.auxiliary.text_parsers import get_chain_iterator
from crplib.auxiliary.seq_parsers import get_twobit_seq
from crplib.mlfeat.featdef import feat_dist, feat_mapsig, get_online_version


def build_conservation_mask(chainfile, chrom, csize):
    """ Build a mask that indicates
    1: is not conserved (= is masked)
    0: is conserved (= is not masked)
    :param chainfile:
    :param chrom:
    :param csize:
    :return:
    """
    mask = np.ones(csize, dtype=np.bool)
    opn, mode = text_file_mode(chainfile)
    with opn(chainfile, mode) as cf:
        chainit = get_chain_iterator(cf, select=chrom)
        for aln in chainit:
            mask[aln[1]:aln[2] + 1] = 0
    return mask


def load_masked_sigtrack(hdfsig, chainfile, group, chrom, csize):
    """
    :return:
    """
    mask = build_conservation_mask(chainfile, chrom, csize)
    with pd.HDFStore(hdfsig, 'r', complib='blosc', complevel=9) as hdf:
        signal = np.ma.array(hdf[group + '/' + chrom].values, mask=mask)
    return signal


def sample_twopass(params):
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
    local_win = res * 5
    far_win = res * 100
    # memory-wise large objects
    chromseq = get_twobit_seq(params['seqfile'], params['chrom'])
    signal = load_masked_sigtrack(params['inputfile'], params['chainfile'], params['group'],
                                  params['chrom'], params['size'])
    pass_one = []
    pass_two = []
    dstfeat = feat_dist
    mapfeat = feat_mapsig
    for n in range(params['numsamples']):
        pos = rand.randint(lolim, hilim)
        dist_reg = {'sample_n': n, 'start': pos, 'end': pos + res}
        dist_reg.update(dstfeat(signal[pos - local_win:pos + local_win], 'local'))
        dist_reg.update(dstfeat(signal[pos - far_win:pos], 'upst'))
        dist_reg.update(dstfeat(signal[pos:pos + far_win], 'dwnst'))
        # for dist estimator/2nd pass, get full information
        y_dep = np.average(signal[pos:pos + res].data)
        dist_reg['y_depvar'] = y_dep
        pass_two.append(dist_reg)
        for move in range(pos - step_bw, pos + step_fw, stepsize):
            seq_reg = {'sample_n': n, 'start': move, 'end': move + res}
            y_dep = np.average(signal[move:move + res].data)
            seq_reg['y_depvar'] = y_dep
            seq_reg['seq'] = chromseq[move:move + res]
            seq_reg.update(feat_mapsig(signal[move:move + res]))
            pass_one.append(seq_reg)
    seqfeat = get_online_version(params['features'], params['kmers'])
    pass_one = list(map(seqfeat, pass_one))
    return mypid, params['chrom'], pass_one, pass_two


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
            if args.twopass:
                mapres = pool.map_async(sample_twopass, arglist)
                metadata = pd.DataFrame(columns=MD_TRAINDATA_COLDEFS)
                for pid, chrom, pass_one, pass_two in mapres.get():
                    logger.debug('Process {} finished chromosome {}'.format(pid, chrom))
                    grp, dataobj, metadata = gen_obj_and_md(metadata, args.grouproot, chrom, 'seq',
                                                            args.inputfile, args.chainfile, pass_one)
                    hdfout.put(grp, dataobj, format='fixed')
                    hdfout.flush()
                    grp, dataobj, metadata = gen_obj_and_md(metadata, args.grouproot, chrom, 'dist',
                                                            args.inputfile, args.chainfile, pass_two)
                    hdfout.put(grp, dataobj, format='fixed')
                    hdfout.flush()
                hdfout.put('metadata', metadata, format='table')
            else:
                raise NotImplementedError('One-pass sampling not available')
    return 0


def run_collect_traindata(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    args.__dict__['keepchroms'] = args.keepchroms.strip('"')
    if args.task == 'regsig':
        logger.debug('Collecting training data for task {}'.format(args.task))
        csizes = read_chromosome_sizes(args.chromsizes, args.keepchroms)
        if args.twopass:
            chromlim = max(args.resolution * 100, 10000)
            _ = collect_sigres_trainsamples(args, csizes, chromlim, logger)
        else:
            chromlim = 10000
    else:
        raise NotImplementedError('Task unknown: {}'.format(args.task))
    return 0
