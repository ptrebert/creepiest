# coding=utf-8

"""
Configuration module to hold feature definitions and
isolated functions to compute these features

"""

import numpy as np
import scipy.stats as stats
import re as re
import copy as cp
import itertools as itt
import functools as funct
import collections as col

FEAT_LENGTH = 'ftlen_abs_length'
FEAT_OECPG = 'ftgc_rat_oeCpG'
FEAT_GC = 'ftgc_pct_GC'
FEAT_CPG = 'ftgc_pct_CpG'
FEAT_REPCONT = 'ftrep_pct_repcon'

# the k avoids singular occurrence of string NA
# can cause problems when feat name is stripped
# to just the kmer: ftkm_pct_NA --> NA
FEAT_KMERFREQ_PREFIX = 'ftkm_pct_k'
FEAT_COREPROM_PREFIX = 'ftprm_pct_'

# all features for signal regression
FEAT_DIST_MEAN = 'ftdst_abs_mean_'
FEAT_DIST_VAR = 'ftdst_abs_var_'
FEAT_DIST_MIN = 'ftdst_abs_min_'
FEAT_DIST_MAX = 'ftdst_abs_max_'
FEAT_DIST_SKEW = 'ftdst_abs_skew_'
FEAT_DIST_KURT = 'ftdst_abs_kurt_'

FEAT_MAPSIG_PREFIX = 'ftmap_'


def _get_feat_fun_map():
    feat_fun_map = {'len': feat_region_length,
                    'prm': feat_coreprom_motifs,
                    'gc': feat_gccpg_content,
                    'cpg': feat_gccpg_content,
                    'oecpg': feat_gccpg_content,
                    'rep': feat_repetitive_content,
                    'kmers': feat_kmer_frequency,
                    'tfm': None}
    return feat_fun_map


def _make_kmer_dict(k, alphabet='ACGTN'):
    """
    :param k:
    :param alphabet:
    :return:
    """
    kmerit = itt.product(alphabet, repeat=k)
    kmers = dict()
    while 1:
        try:
            kmers[FEAT_KMERFREQ_PREFIX + ("".join(next(kmerit)))] = 0
        except StopIteration:
            break
    return kmers


def get_online_version(features, kmers=None):
    """ Return a closure encompassing all feature
    functions - intended to use is with
    multiprocessing.Pool.map() or similar
    """
    if 'kmers' in features:
        assert kmers is not None, 'No values for k-mer frequency specified'
    funmap = _get_feat_fun_map()
    exec_functions = set()
    for ft in features:
        if ft == 'kmers':
            for k in kmers:
                kd = _make_kmer_dict(k)
                part = funct.partial(feat_single_kmer, *(k, kd))
                exec_functions.add(part)
        else:
            exec_functions.add(funmap[ft])

    def compute_features(d):
        for f in exec_functions:
            d = f(d)
        return d
    return compute_features


def feat_region_length(region):
    """ overengineering, anyone?
    :param region:
    :return:
    """
    reglen = region['end'] - region['start']
    assert reglen > 0, 'Malformed genomic region (length): {}'.format(region)
    region[FEAT_LENGTH] = reglen
    return region


def feat_repetitive_content(region):
    """
    :param region:
    :return:
    """
    try:
        _ = region[FEAT_REPCONT]
        return region
    except KeyError:
        tmpseq = region['seq']
        seqlen = float(len(tmpseq))
        repmasked = tmpseq.count('a') + tmpseq.count('c') + tmpseq.count('g') + tmpseq.count('t')
        region[FEAT_REPCONT] = round((repmasked / seqlen) * 100., 2)
        return region


def feat_single_kmer(k, kmers, region):
    """
    :param k:
    :param kmers:
    :param region:
    :return:
    """
    # TODO
    # the deep copy is probably necessary for the intended
    # use case, see get_online_version
    mykmers = cp.deepcopy(kmers)
    tmpseq = region['seq'].upper()  # kmers are not case-sensitive
    seqlen = len(tmpseq)
    total_klen = float(seqlen - k + 1)
    wordfreqs = col.defaultdict(int)
    for i in range(0, k):
        # curly braces in literal part of format string
        # need to be escaped with curly braces
        words = re.findall('.{{{}}}'.format(k), tmpseq[i:])
        for w in words:
            wordfreqs[FEAT_KMERFREQ_PREFIX + w] += 1
    for key, val in wordfreqs.items():
        val = round((val / total_klen) * 100., 2)
        mykmers[key] = val
    region.update(mykmers)
    return region


def feat_kmer_frequency(kmers, region):
    """
    :param region:
    :param kmers:
    :return:
    """
    tmpseq = region['seq'].upper()
    seqlen = len(tmpseq)
    for k in kmers:
        # normalization factor
        # number of possible substrings of length k in the seq
        # that is: n - k + 1 (k len substr in seq of len n)
        total_klen = float(seqlen - k + 1)
        kmerit = itt.product('ACGTN', repeat=k)
        kmerdict = dict()
        while 1:
            try:
                kmerdict[FEAT_KMERFREQ_PREFIX + ''.join(next(kmerit))] = 0
            except StopIteration:
                break
        wordfreqs = col.defaultdict(int)
        for i in range(0, k):
            # curly braces in literal part of format string
            # need to be escaped with curly braces
            words = re.findall('.{{{}}}'.format(k), tmpseq[i:])
            for w in words:
                wordfreqs[FEAT_KMERFREQ_PREFIX + w] += 1
        for key, val in wordfreqs.items():
            val = round((val / total_klen) * 100., 2)
            kmerdict[key] = val
        region.update(kmerdict)
    return region


def feat_gccpg_content(region):
    """
    :param region:
    :return:
    """
    try:
        _ = region[FEAT_OECPG]
        _ = region[FEAT_GC]
        _ = region[FEAT_CPG]
        return region
    except KeyError:
        seq = region['seq'].lower()
        seqlen = float(len(seq))  # for division
        total_G = seq.count('g')
        total_C = seq.count('c')
        total_CpG = seq.count('cg')
        region[FEAT_GC] = round(((total_G + total_C) / seqlen) * 100., 2)
        region[FEAT_CPG] = round((total_CpG / (seqlen / 2.)) * 100., 2)
        # this definition of the obs-exp ratio is taken from UCSC
        region[FEAT_OECPG] = round((total_CpG / (max(1, total_G) * max(1, total_C))) * seqlen, 2)
        return region


def feat_coreprom_motifs(region):
    """
    :param region:
     :type: dict
    :return:
    """

    core_motifs = []

    # DOI:10.1016/j.gene.2006.09.029
    # Yang et al., 2007 Gene
    # Refer specifically to TATA-less promoters

    core_motifs.append(('elkM3', '[GC]CGGAAG[CT]'))  # not sure if that is a reasonable one
    core_motifs.append(('sp1M6', 'GGGCGG[AG]'))  # not sure if that is a reasonable one
    core_motifs.append(('novelM22', 'TGCGCA[ACGTN][GT]'))

    # DOI:10.1016/j.ydbio.2009.08.009
    # Juven-Gershon, Kadonaga 2010, Dev. Bio
    # DOI:10.1101/gad.1026202
    # Butler, Kadonaga 2002, Genes & Dev.

    core_motifs.append(('tataM3', 'TATA[AT]AA[AG]'))  # TATA box
    core_motifs.append(('inrM4', '[TC][TC]A[ACGTN][AT][TC][TC]'))  # initiator (Inr)
    core_motifs.append(('breMx', '[GC][GC][AG]CGCC'))  # TFIIB recognition element (BRE)
    core_motifs.append(('dpeM9', '[AG]G[AT][TC][GAC](T)?'))  # downstream core promoter element (DPE)
    core_motifs.append(('mteM10', 'C[GC]A[AG]C[GC][GC]AACG[GC]'))  # motif ten (MTE)

    # DOI:10.1093/nar/gkv1032
    # Marbach-Bar et al., 2016 NAR
    core_motifs.append(('dtieMx', 'G[CGT][CGT][AG][AGT][ACGTN][ACT]GG'))  # Downstream Transcription Initiation Element (DTIE)
    tmpseq = region['seq'].upper()
    try:
        reglen = region[FEAT_LENGTH]
        assert len(tmpseq) == reglen, 'Malformed genomic region (length): {}'.format(region)
    except KeyError:
        try:
            reglen = region['end'] - region['start']
            assert len(tmpseq) == reglen, 'Malformed genomic region (length): {}'.format(region)
        except KeyError:
            reglen = len(tmpseq)
    for name, motifre in core_motifs:
        bpcov = sum(len(m) for m in re.findall(motifre, tmpseq))
        region[FEAT_COREPROM_PREFIX + name] = round((bpcov / reglen) * 100, 2)
    return region


def feat_dist(sample, suffix):
    """
    :param sample:
    :param suffix:
    :return:
    """
    desc = stats.describe(sample)
    kurt = stats.kurtosis(sample, fisher=False)
    ret = dict()
    ret[FEAT_DIST_MEAN + suffix] = desc.mean
    ret[FEAT_DIST_VAR + suffix] = desc.variance
    ret[FEAT_DIST_SKEW + suffix] = desc.skewness
    ret[FEAT_DIST_MIN + suffix] = desc.minmax[0]
    ret[FEAT_DIST_MAX + suffix] = desc.minmax[1]
    ret[FEAT_DIST_KURT + suffix] = kurt
    return ret


def feat_mapsig(sample):
    """
    :param sample:
     :type: numpy.MaskedArray
    :return:
    """
    conserved = np.ma.count(sample)
    reglen = sample.size
    ret = dict()
    ret[FEAT_MAPSIG_PREFIX + 'pct_cons'] = round((conserved / reglen) * 100, 2)
    if conserved == 0:
        cons_mean = 0
        cons_max = 0
        cons_min = 0
    else:
        cons_mean = np.ma.mean(sample)
        cons_max = np.ma.max(sample)
        cons_min = np.ma.min(sample)
    ret[FEAT_MAPSIG_PREFIX + 'abs_mean'] = cons_mean
    ret[FEAT_MAPSIG_PREFIX + 'abs_max'] = cons_max
    ret[FEAT_MAPSIG_PREFIX + 'abs_min'] = cons_min
    return ret
