# coding=utf-8

"""
Some helper functions to parse text-based files
with more or less no well-defined format
"""

import re as re
import collections as col

from crplib.auxiliary.file_ops import text_file_mode


def read_chromosome_sizes(fpath, keep='.+'):
    """
    :param fpath:
    :param keep:
    :return:
    """
    chroms = dict()
    keeper = re.compile(keep)
    opn, mode = text_file_mode(fpath)
    with opn(fpath, mode=mode, encoding='ascii') as infile:
        for line in infile:
            if not line.strip():
                continue
            cols = line.strip().split()
            cname = cols[0].strip()
            m = keeper.match(cname)
            if m is not None:
                csize = int(cols[1])
                chroms[cname] = csize
    assert chroms, 'No chromosomes from file {} selected with pattern {}'.format(fpath, keep)
    return chroms


def _read_chain_header(line):
    _, _, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, _ = line.strip().split()
    return tName, int(tSize), tStrand, int(tStart), int(tEnd), qName, int(qSize), qStrand, int(qStart), int(qEnd)


def _read_aln_dataline(line):
    size, dt, dq = line.strip().split()
    return int(size), int(dt), int(dq)


def _read_block_end(line):
    size = line.strip().split()[0]
    return int(size)


def get_chain_iterator(fobj, tselect='all', qcheck=None):
    """ Returns an iterator over chain files as used by
    UCSC liftOver tool. The design assumes a simple parallelization, i.e.
    many processes can read the same chain file, each one filtering
    out 1 target chain and many query chains. This always assumes that there
    are some blocks in the chain file, otherwise raises AssertionError.

    :param fobj: object supporting line iteration/read
     :type: file-like object opened in text mode
    :param tselect: exact name match for allowed target chromosome
     :type: str
    :param qcheck: compiled regex object to check for allowed query chromosomes
     :type: None or re.regex object
    :return:
    """
    tchrom = None
    qchrom = None
    tstrand = '+'
    qstrand = None
    trun = -1
    qrun = -1
    exp_tend = -1
    exp_qend = -1
    read_head = _read_chain_header
    read_aln = _read_aln_dataline
    read_end = _read_block_end
    skip = False
    bc = 0
    for line in fobj:
        if line.strip():
            if line.startswith('chain'):
                assert trun == exp_tend or exp_tend == -1, 'Missed expected block end: {} vs {}'.format(trun, exp_tend)
                assert qrun == exp_qend or exp_qend == -1, 'Missed expected block end: {} vs {}'.format(qrun, exp_qend)
                parts = read_head(line)
                assert parts[2] == '+', 'Reverse target chain is unexpected: {}'.format(line)
                tchrom = parts[0]
                if not (tchrom == tselect or tselect == 'all'):
                    skip = True
                    continue
                qchrom = parts[5]
                if qcheck is not None:
                    mobj = qcheck.match(qchrom)
                    if mobj is None:
                        skip = True
                        continue
                skip = False
                trun = parts[3]
                exp_tend = parts[4]
                qstrand = parts[7]
                if qstrand == '+':
                    qrun = parts[8]
                    exp_qend = parts[9]
                elif qstrand == '-':
                    qrun = parts[6] - parts[9]
                    exp_qend = parts[6] - parts[8]
                else:
                    raise ValueError('Unknown qStrand specified: {}'.format(line))
            else:
                if skip:
                    continue
                try:
                    size, dt, dq = read_aln(line)
                    yield tchrom, trun, trun + size, tstrand, qchrom, qrun, qrun + size, qstrand
                    bc += 1
                    trun += size + dt
                    qrun += size + dq
                except ValueError:
                    size = read_end(line)
                    yield tchrom, trun, trun + size, tstrand, qchrom, qrun, qrun + size, qstrand
                    bc += 1
                    trun += size
                    qrun += size
    assert bc > 0, 'No aln. blocks read from chain file for target select {}'.format(tselect)
    return


def chromsize_from_chain(chainfile, chrom):
    """
    :param chainfile:
    :param chrom:
    :return:
    """
    read_head = _read_chain_header
    opn, mode = text_file_mode(chainfile)
    chrom_size = 0
    with opn(chainfile, mode=mode, encoding='ascii') as chf:
        for line in chf:
            if line.strip() and line.startswith('chain'):
                parts = read_head(line)
                if parts[0] == chrom:
                    chrom_size = parts[1]
                    break
    assert chrom_size > 0, 'No entry in chain file {} for chromosome: {}'.format(chainfile, chrom)
    return chrom_size


def get_meme_iterator(fobj):
    """
    :param fobj:
    :return:
    """

    header = fobj.readline()
    assert header.startswith('MEME'), 'Supposed MEME motif DB does not start with field "MEME"'
    for line in fobj:
        if line.startswith('MOTIF'):
            _, tfid, tfname = line.strip().split()
            yield tfid, tfname
    return


def _read_ovl_line(line):
    """ chr4	0	250	chr4	1	12	MA0041.1	11
    :param line:
    :return:
    """
    cols = line.strip().split()
    return cols[0], int(cols[1]), cols[6], int(cols[7])


def get_binned_motifs_iterator(fobj):
    """
    :param fobj:
    :return:
    """
    rdl = _read_ovl_line
    init = rdl(fobj.readline())
    curr_chrom = init[0]
    curr_index = init[1]
    fobj.seek(0)
    collect = col.Counter()
    for line in fobj:
        if not line:
            continue
        chrom, index, tfid, ovl = rdl(line)
        if index != curr_index:
            yield curr_chrom, curr_index, collect
            if curr_chrom == chrom:
                assert curr_index < index, \
                    'Binned motif file seems to be unsorted: step from {} to {}'.format(curr_index, index)
            curr_chrom = chrom
            curr_index = index
            collect = col.Counter()
        collect[tfid] += ovl
    yield curr_chrom, curr_index, collect
    return
