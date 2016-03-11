# coding=utf-8

"""
Configuration module to hold feature definitions and
isolated functions to compute these features
"""

import re as re
import itertools as itt
import collections as col

FEAT_LENGTH = 'feat_scn_abs_length'
FEAT_OECPG = 'feat_scn_rat_oeCpG'
FEAT_GC = 'feat_scn_pct_GC'
FEAT_CPG = 'feat_scn_pct_CpG'
FEAT_REPCONT = 'feat_scn_pct_repcont'

FEAT_KMER_VALUES = 2, 3, 4

# the k avoids singular occurrence of string NA
# can cause problems when feat name is stripped
# to just the kmer: feat_scn_pct_NA --> NA
FEAT_KMERFREQ_PREFIX = 'feat_scn_pct_k'
FEAT_COREPROM_PREFIX = 'feat_scn_pct_'


def feat_region_length(region):
    """ over engineering, anyone?
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


def feat_kmer_frequency(region, kmers=FEAT_KMER_VALUES):
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

    core_motifs = core_motifs.append(('elkM3', '[GC]CGGAAG[CT]'))  # not sure if that is a reasonable one
    core_motifs = core_motifs.append(('sp1M6', 'GGGCGG[AG]'))  # not sure if that is a reasonable one
    core_motifs = core_motifs.append(('novelM22', 'TGCGCA[ACGTN][GT]'))

    # DOI:10.1016/j.ydbio.2009.08.009
    # Juven-Gershon, Kadonaga 2010, Dev. Bio
    # 10.1101/gad.1026202
    # Butler, Kadonaga 2002, Genes & Dev.

    core_motifs = core_motifs.append(('tataM3', 'TATA[AT]AA[AG]'))  # TATA box
    core_motifs = core_motifs.append(('inrM4', '[TC][TC]A[ACGTN][AT][TC][TC]'))  # initiator (Inr)
    core_motifs = core_motifs.append(('breMx', '[GC][GC][AG]CGCC'))  # TFIIB recognition element (BRE)
    core_motifs = core_motifs.append(('dpeM9', '[AG]G[AT][TC][GAC](T)?'))  # downstream core promoter element (DPE)
    core_motifs = core_motifs.append(('mteM10', 'C[GC]A[AG]C[GC][GC]AACG[GC]'))  # motif ten (MTE)

    # DOI:10.1093/nar/gkv1032
    # Marbach-Bar et al., 2016 NAR
    core_motifs = core_motifs.append(('dtieMx', 'G[CGT][CGT][AG][AGT][ACGTN][ACT]GG'))  # Downstream Transcription Initiation Element (DTIE)

    try:
        reglen = region[FEAT_LENGTH]
    except KeyError:
        reglen = region['end'] - region['start']
    tmpseq = region['seq'].upper()
    assert len(tmpseq) == reglen, 'Malformed genomic region (length): {}'.format(region)
    for name, motifre in core_motifs:
        bpcov = sum(len(m) for m in re.findall(motifre, tmpseq))
        region[FEAT_COREPROM_PREFIX + name] = round((bpcov / reglen) * 100, 2)

    return region


