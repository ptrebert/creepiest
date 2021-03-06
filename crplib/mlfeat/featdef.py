# coding=utf-8

"""
Configuration module to hold feature definitions and
isolated functions to compute these features

"""

import sys as sys
import re as re
import copy as cp
import itertools as itt
import functools as funct
import collections as col

import numpy as np
import scipy.stats as stats
import pandas as pd

FEAT_LENGTH = 'ftlen_abs_length'
FEAT_RELLENGTH = 'ftlen_pct_length'
FEAT_OECPG = 'ftoecpg_rat_oeCpG'
FEAT_GC = 'ftgc_pct_GC'
FEAT_CPG = 'ftcpg_pct_CpG'
FEAT_REPCONT = 'ftrep_pct_repcon'

# the k avoids singular occurrence of string NA
# can cause problems when feat name is stripped
# to just the kmer: ftkm_pct_NA --> NA
FEAT_KMERFREQ_PREFIX = 'ftkmf_pct_k'

FEAT_COREPROM_PREFIX = 'ftprm_'

FEAT_TFMOTIF_PREFIX = 'fttfm_pct_'

# all features for signal regression
# distribution features are (atm) no longer supported
FEAT_DIST_MEAN = 'ftdst_abs_mean_'
FEAT_DIST_VAR = 'ftdst_abs_var_'
FEAT_DIST_MIN = 'ftdst_abs_min_'
FEAT_DIST_MAX = 'ftdst_abs_max_'
FEAT_DIST_SKEW = 'ftdst_abs_skew_'
FEAT_DIST_KURT = 'ftdst_abs_kurt_'

FEAT_MAPSIG_PREFIX = 'ftmsig_'
FEAT_ROI_PREFIX = 'ftroi_'
FEAT_ASCREG_PREFIX = 'ftasc_'

FEAT_DNASE_MEDPROB = 'ftdnm_abs_medprob'
FEAT_DNASE_MAXPROB = 'ftdnm_abs_maxprob'
FEAT_DNASE_MEDAD = 'ftdnm_abs_medad'

# priors are not computed, but manually added to the training data
FEAT_DERIVED = 'ftdrv'

FEAT_PREFIX_MAP = {'len': [FEAT_LENGTH, FEAT_RELLENGTH], 'prm': [FEAT_COREPROM_PREFIX],
                   'gc': [FEAT_GC], 'cpg': [FEAT_CPG], 'oecpg': [FEAT_OECPG],
                   'rep': [FEAT_REPCONT], 'kmf': [FEAT_KMERFREQ_PREFIX], 'tfm': [FEAT_TFMOTIF_PREFIX],
                   'msig': [FEAT_MAPSIG_PREFIX], 'roi': [FEAT_ROI_PREFIX],
                   'dnm ': [FEAT_DNASE_MEDPROB, FEAT_DNASE_MAXPROB, FEAT_DNASE_MEDAD],
                   'drv': [FEAT_DERIVED], 'asc': [FEAT_ASCREG_PREFIX]}

FEAT_CLASS_MAP = dict()
for ftclass, prefixes in FEAT_PREFIX_MAP.items():
    for pref in prefixes:
        FEAT_CLASS_MAP[pref] = ftclass
        FEAT_CLASS_MAP[pref.strip('_')] = ftclass
        if pref.startswith('ftkmf_pct'):
            FEAT_CLASS_MAP[pref.strip('_k')] = ftclass


def _format_malformed_region(reg):
    """ Remove unnecessary information from
    region upon error
    :param reg:
    :return:
    """
    infos = []
    for k, v in reg.items():
        if k == FEAT_LENGTH or k == FEAT_RELLENGTH:
            infos.append((k, v))
            continue
        if k.startswith('ft'):
            continue
        if k == 'seq':
            infos.append((k, 'seq len: {}'.format(len(v))))
            continue
        infos.append((k, v))
    return sorted(infos)


def _get_feat_fun_map():
    feat_fun_map = {'len': feat_region_length,
                    'prm': feat_coreprom_motifs,
                    'gc': feat_gc_content,
                    'cpg': feat_cpg_content,
                    'oecpg': feat_oecpg_content,
                    'rep': feat_repetitive_content,
                    'kmf': feat_kmer_frequency,
                    'tfm': feat_tf_motifs,
                    'msig': feat_mapsig,
                    'roi': feat_roi,
                    'dnm': feat_dnase_motif,
                    'asc': feat_ascreg}
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


def verify_sample_integrity(samples, check_seq=False):
    """
    :param samples:
    :param check_seq:
    :return:
    """
    for smp in samples:
        try:
            reglen = smp['end'] - smp['start']
            assert reglen > 0, 'Malformed region (length): {}'.format(_format_malformed_region(smp))
            if check_seq:
                seqlen = len(smp['seq'])
                assert reglen == seqlen, 'Malformed region (seq. length): {}'.format(_format_malformed_region(smp))
        except KeyError as ke:
            raise AssertionError('Was expecting key {} while checking region integrity'.format(str(ke)))
    return


def get_prefix_list(features):
    """
    :param features:
    :return:
    """
    relevant_prefixes = []
    for k, v in FEAT_PREFIX_MAP.items():
        if k in features:
            relevant_prefixes.extend(v)
    return relevant_prefixes


def get_classes_from_names(featnames):
    """
    :param featnames:
    :return:
    """
    q1 = col.deque(featnames)
    q1.append('SENTINEL')
    q2 = col.deque()
    splits = 0
    ftclasses = set()
    while 1:
        for ft in iter(q1.popleft, 'SENTINEL'):
            prefix = ft.rsplit('_', splits)[0]
            if prefix in FEAT_CLASS_MAP:
                ftclasses.add(FEAT_CLASS_MAP[prefix])
            else:
                if splits > ft.count('_'):
                    raise AssertionError('Cannot find feature class for feature name: {} (splits: {})'.format(ft, splits))
                q2.append(ft)
        splits += 1
        q2.append('SENTINEL')
        if len(q2) <= 1:
            break
        q1, q2 = q2, q1
    return sorted(ftclasses)


def check_online_available(reqfeat):
    """
    :param reqfeat:
    :return:
    """
    avfeat = filter(lambda f: f in ['len', 'prm', 'gc', 'cpg', 'oecpg', 'rep', 'kmf', 'dnm'], reqfeat)
    return list(avfeat)


def get_online_version(features, kmers=None, yardstick=0):
    """ Return a closure encompassing all feature
    functions - intended to use is with
    multiprocessing.Pool.map() or similar
    """
    if 'kmf' in features:
        assert kmers is not None, 'No values for k-mer frequency specified'
    funmap = _get_feat_fun_map()
    exec_functions = set()
    avfeat = check_online_available(features)
    for ft in avfeat:
        if ft == 'kmf':
            for k in kmers:
                kd = _make_kmer_dict(k)
                part = funct.partial(feat_single_kmer, *(k, kd))
                exec_functions.add(part)
        elif ft == 'len':
            part = funct.partial(feat_region_length, *(yardstick, ))
            exec_functions.add(part)
        else:
            exec_functions.add(funmap[ft])

    def compute_features(d):
        for f in exec_functions:
            d = f(d)
        return d
    return compute_features


def feat_region_length(yardstick, region):
    """
    :param yardstick:
    :param region:
    :return:
    """
    reglen = region['end'] - region['start']
    assert reglen > 0, 'Malformed genomic region (len): {}'.format(_format_malformed_region(region))
    if yardstick > 0:
        region[FEAT_RELLENGTH] = reglen / yardstick * 100
    else:
        region[FEAT_LENGTH] = reglen
    return region


def feat_repetitive_content(region):
    """
    :param region:
    :return:
    """
    tmpseq = region['seq']
    seqlen = len(tmpseq)
    assert seqlen > 0, 'Malformed genomic region (len): {}'.format(_format_malformed_region(region))
    repmasked = tmpseq.count('a') + tmpseq.count('c') + tmpseq.count('g') + tmpseq.count('t')
    region[FEAT_REPCONT] = (repmasked / seqlen) * 100.
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
        val = (val / total_klen) * 100.
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
            val = (val / total_klen) * 100.
            kmerdict[key] = val
        region.update(kmerdict)
    return region


def feat_gc_content(region):
    """
    :param region:
    :return:
    """
    seq = region['seq'].lower()
    seqlen = len(seq)
    assert seqlen > 0, 'Malformed genomic region (len): {}'.format(_format_malformed_region(region))
    total_G = seq.count('g')
    total_C = seq.count('c')
    region[FEAT_GC] = ((total_G + total_C) / seqlen) * 100.
    return region


def feat_cpg_content(region):
    """
    :param region:
    :return:
    """
    seq = region['seq'].lower()
    seqlen = len(seq)
    assert seqlen > 0, 'Malformed genomic region (len): {}'.format(_format_malformed_region(region))
    total_CpG = seq.count('cg')
    region[FEAT_CPG] = (total_CpG / (seqlen / 2.)) * 100.
    return region


def feat_oecpg_content(region):
    """
    :param region:
    :return:
    """
    seq = region['seq'].lower()
    seqlen = len(seq)  # for division
    assert seqlen > 0, 'Malformed genomic region (len): {}'.format(_format_malformed_region(region))
    total_G = seq.count('g')
    total_C = seq.count('c')
    total_CpG = seq.count('cg')
    # this definition of the obs-exp ratio is taken from UCSC
    region[FEAT_OECPG] = (total_CpG / (max(1, total_G) * max(1, total_C))) * seqlen
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
    reglen = len(tmpseq)
    for name, motifre in core_motifs:
        occurrences = re.findall(motifre, tmpseq)
        bpcov = sum(len(m) for m in occurrences)
        region[FEAT_COREPROM_PREFIX + 'pct_' + name] = (bpcov / reglen) * 100
        region[FEAT_COREPROM_PREFIX + 'abs_' + name] = len(occurrences)
    return region


@funct.lru_cache(maxsize=4096, typed=False)
def _pwm_motif_prob(match):
    """
    :param match:
    :param pwm:
    :return:
    """
    # PWM for preferred DNaseI cleavage sites
    # taken from
    # Herrera and Chaires, J. Mol. Biol. (1994) 236, 405-411
    pwm = pd.DataFrame([[44, 4, 28, 19, 31, 24],
                        [26, 39, 46, 42, 2, 20],
                        [7, 41, 15, 19, 33, 35],
                        [22, 17, 11, 20, 33, 20]],
                       index=list('ATCG'), columns=[0, 1, 2, 3, 4, 5], dtype=np.float64)
    pwm /= 100.
    p = (pwm.lookup(list(match), [0, 1, 2, 3, 4, 5])).prod()
    return p


def feat_dnase_motif(region):
    """ Since the computations for this function are a bit more involved,
    this is the first candidate to remove from the set of online feature
    functions
    :param region:
    :return:
    """
    tmpseq = region['seq'].upper()
    region[FEAT_DNASE_MEDPROB] = 0
    region[FEAT_DNASE_MAXPROB] = 0
    region[FEAT_DNASE_MEDAD] = 0
    #motif_prob = funct.partial(_pwm_motif_prob, *(pwm, np.arange(6)))
    motif_prob = _pwm_motif_prob
    match_probs = []
    for idx in range(0, 6):
        match_probs.extend([motif_prob(m) for m in re.findall('[ACGT]{6}', tmpseq[idx:])])
        # Note to self:
        # benchmarking indicates that findall is faster than finditer here
        # I assume because of mobj.group() call for each object
        #match_probs.extend([motif_prob(mobj.group(0)) for mobj in re.finditer('[ACGT]{6}', tmpseq[idx:])])
    if match_probs:
        match_probs = np.array(match_probs, dtype=np.float64)
        medprob = np.median(match_probs)
        region[FEAT_DNASE_MEDPROB] = medprob
        region[FEAT_DNASE_MAXPROB] = np.max(match_probs)
        region[FEAT_DNASE_MEDAD] = np.median(np.absolute(match_probs - medprob))
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


def feat_mapsig(sample, infix=''):
    """
    :param sample:
     :type: numpy.MaskedArray
    :return:
    """
    if infix and not infix.endswith('_'):
        infix += '_'
    # numpy.ma.MaskedArray.count
    # => count the non-masked elements
    # Usage: mask out positions in signal track
    # that are not conserved
    # is masked = is not conserved
    # is not masked = is conserved
    conserved = np.ma.count(sample)
    reglen = sample.size
    ret = dict()
    ret[FEAT_MAPSIG_PREFIX + infix + 'pct_cons'] = (conserved / reglen) * 100
    if conserved == 0:
        cons_mean = 0
        cons_median = 0
        cons_max = 0
        cons_min = 0
    else:
        cons_mean = np.ma.mean(sample)
        cons_median = np.ma.median(sample)
        cons_max = np.ma.max(sample)
        cons_min = np.ma.min(sample)
    ret[FEAT_MAPSIG_PREFIX + infix + 'abs_mean'] = cons_mean
    ret[FEAT_MAPSIG_PREFIX + infix + 'abs_median'] = cons_median
    ret[FEAT_MAPSIG_PREFIX + infix + 'abs_max'] = cons_max
    ret[FEAT_MAPSIG_PREFIX + infix + 'abs_min'] = cons_min
    return ret


def feat_ascreg(sample, ascregions, asclabel, signal, siglabel):
    """
    :param sample:
    :param ascregions:
    :param asclabel:
    :param signal:
    :param siglabel:
    :return:
    """
    dist = []
    sig = []
    d_min = sys.maxsize
    nxt_idx = -1
    sample_mid = int(sample['start'] + (sample['end'] - sample['start']) // 2)
    for idx, (s, e) in enumerate(zip(ascregions['starts'], ascregions['ends'])):
        asc_mid = int(s + (e - s) // 2)
        asc_sig = signal[s:e]
        c = np.ma.count(asc_sig)
        if c == 0:
            sig.append((0, 0, 0, 0))
        else:
            sig.append((np.ma.min(asc_sig), np.ma.max(asc_sig),
                        np.ma.median(asc_sig), np.ma.mean(asc_sig)))
        d = abs(asc_mid - sample_mid)
        dist.append(d)
        if d < d_min:
            d_min = d
            nxt_idx = idx
    assert nxt_idx > -1, 'Invalid index for closest associated region'
    dist = np.array(dist, dtype=np.int32)
    ret = dict()
    # features irrespective of signal
    label_prefix = FEAT_ASCREG_PREFIX + asclabel + '_'
    ret[label_prefix + 'abs_num'] = dist.size
    ret[label_prefix + 'abs_dist_mean'] = dist.mean()
    ret[label_prefix + 'abs_dist_min'] = dist.min()
    ret[label_prefix + 'abs_dist_max'] = dist.max()
    # signal features
    label_prefix = FEAT_ASCREG_PREFIX + asclabel + '_' + siglabel + '_'
    mins = np.array([x[0] for x in sig], dtype=np.float32)
    maxs = np.array([x[1] for x in sig], dtype=np.float32)
    means = np.array([x[2] for x in sig], dtype=np.float32)
    medians = np.array([x[3] for x in sig], dtype=np.float32)
    ret[label_prefix + 'abs_min'] = mins.mean()
    ret[label_prefix + 'abs_max'] = maxs.mean()
    ret[label_prefix + 'abs_median'] = medians.mean()
    ret[label_prefix + 'abs_mean'] = means.mean()
    # same again weighted
    ret[label_prefix + 'abs_wt-min'] = np.average(mins, weights=ascregions['weights'])
    ret[label_prefix + 'abs_wt-max'] = np.average(maxs, weights=ascregions['weights'])
    ret[label_prefix + 'abs_wt-median'] = np.average(medians, weights=ascregions['weights'])
    ret[label_prefix + 'abs_wt-mean'] = np.average(means, weights=ascregions['weights'])
    # same again for closest only
    ret[label_prefix + 'abs_nx-min'] = mins[nxt_idx]
    ret[label_prefix + 'abs_nx-max'] = maxs[nxt_idx]
    ret[label_prefix + 'abs_nx-median'] = medians[nxt_idx]
    ret[label_prefix + 'abs_nx-mean'] = means[nxt_idx]
    return ret


def feat_ascreg_default(asclabel, siglabel):
    """
    :param asclabel:
    :param siglabel:
    :return:
    """
    ret = dict()
    # features irrespective of signal
    label_prefix = FEAT_ASCREG_PREFIX + asclabel + '_'
    ret[label_prefix + 'abs_num'] = 0
    ret[label_prefix + 'abs_dist_mean'] = 0
    ret[label_prefix + 'abs_dist_min'] = 0
    ret[label_prefix + 'abs_dist_max'] = 0
    # signal features
    label_prefix = FEAT_ASCREG_PREFIX + asclabel + '_' + siglabel + '_'
    ret[label_prefix + 'abs_min'] = 0
    ret[label_prefix + 'abs_max'] = 0
    ret[label_prefix + 'abs_median'] = 0
    ret[label_prefix + 'abs_mean'] = 0
    # same again weighted
    ret[label_prefix + 'abs_wt-min'] = 0
    ret[label_prefix + 'abs_wt-max'] = 0
    ret[label_prefix + 'abs_wt-median'] = 0
    ret[label_prefix + 'abs_wt-mean'] = 0
    # same again for closest only
    ret[label_prefix + 'abs_nx-min'] = 0
    ret[label_prefix + 'abs_nx-max'] = 0
    ret[label_prefix + 'abs_nx-median'] = 0
    ret[label_prefix + 'abs_nx-mean'] = 0
    return ret


def _feat_roi_default(infix, quantify):
    """
    :param infix:
    :param quantify:
    :return:
    """
    feats = {}
    if len(quantify & {'binary', 'all'}) > 0:
        feats.update({FEAT_ROI_PREFIX + infix + 'bin_obs': 0})
    if len(quantify & {'density', 'all'}) > 0:
        feats.update({FEAT_ROI_PREFIX + infix + 'abs_num': 0.0})
    if len(quantify & {'coverage', 'all'}) > 0:
        feats.update({FEAT_ROI_PREFIX + infix + 'pct_cov': 0.0})
    return feats


def _feat_roi_bpcov(row, start, end):
    return max(row.start, start) - min(row.end, end)


def feat_roi(regions, rois, infix='', quantify={'all'}):
    """
    :param regions:
    :param rois:
    :param infix:
    :param quantify:
    :return:
    """
    quantify = set(quantify)
    if infix and not infix.endswith('_'):
        infix += '_'
    default = _feat_roi_default(infix, quantify)
    calc_bpcov = _feat_roi_bpcov
    [reg.update(default) for reg in regions]
    for reg in regions:
        ovl = rois.query('start < {} and end > {}'.format(reg['end'], reg['start']))
        if ovl.empty:
            continue
        s, e = reg['start'], reg['end']
        reglen = e - s
        bpcov = ovl.apply(calc_bpcov, axis=1, args=(s, e))
        tot_cov = bpcov.sum()
        if 'all' in quantify:
            reg[FEAT_ROI_PREFIX + infix + 'bin_obs'] = 1
            reg[FEAT_ROI_PREFIX + infix + 'pct_cov'] = tot_cov / reglen * 100.
            reg[FEAT_ROI_PREFIX + infix + 'abs_num'] = ovl.shape[0]
        # TODO here:
        # write smaller functions that allow to avoid this clunky construct
        # by setting defaults as necessary (as soon as performance becomes an issue)
        else:
            try:
                reg[FEAT_ROI_PREFIX + infix + 'bin_obs'] = 1
            except KeyError:
                pass
            try:
                reg[FEAT_ROI_PREFIX + infix + 'abs_num'] = ovl.shape[0]
            except KeyError:
                pass
            try:
                reg[FEAT_ROI_PREFIX + infix + 'pct_cov'] = tot_cov / reglen * 100.
            except KeyError:
                pass
    return regions


def _weighted_motifs(row, start, end, binsize):
    """
    :param row:
    :param start:
    :param end:
    :param binsize:
    :return:
    """
    # Important: this works with the TF motifs
    # because these don't have a name column by construction,
    # i.e. row.name gives the index (start coordinate by construction).
    # Iff there is a name column in the row that needs to be accessed,
    # it has to be done as row['name']
    ridx = row.name
    if start >= ridx and end <= ridx + binsize:
        w = (end - start) / binsize
        return row * w
    elif start < ridx and end > ridx + binsize:
        w = (ridx + binsize - ridx) / binsize
        return row * w
    elif start < ridx and end <= ridx + binsize:
        # left overlap
        w = (end - ridx) / binsize
        return row * w
    elif start >= ridx and end > ridx + binsize:
        # right overlap
        w = (ridx + binsize - start) / binsize
        return row * w
    else:
        # first time this throws, need a UnitTest...
        raise AssertionError('Unable to compute motif weights for values: {}, {} - {} {}'.format(row, start, end, binsize))


def feat_tf_motifs(regions, motifcounts):
    """
    :param regions:
     :type: list of dict
    :param motifcounts: Pandas DataFrame
    :return:
    """
    motifnames = list(motifcounts.columns)
    base_counts = col.Counter(motifnames)  # all zero counts
    mc = motifcounts
    binsize = mc.index[1] - mc.index[0]
    assert binsize > 0, 'Unable to determine bin size for indices: {} and {}'.format(mc.index[1], mc.index[0])
    wm = _weighted_motifs
    for reg in regions:
        s, e = reg['start'], reg['end']
        ovl = mc.loc[(mc.index < e) & (mc.index + binsize > s), ]
        if ovl.empty:
            reg.update(base_counts)
            continue
        ovl = ovl.apply(wm, axis=1, args=(s, e, binsize))
        ovl = ovl.sum(axis=0) / (e - s) * 100
        reg.update(ovl.to_dict())

    return regions

