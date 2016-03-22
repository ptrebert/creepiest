# coding=utf-8

"""
Some helper functions to parse text-based files
with more or less no well-defined format
"""

import re as re

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


def get_chain_iterator(fobj, select='all'):
    """ Returns an iterator over chain files as used by
    UCSC liftOver tool
    :param fobj:
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
    for line in fobj:
        if line.strip():
            if line.startswith('chain'):
                assert trun == exp_tend or exp_tend == -1, 'Missed expected block end: {} vs {}'.format(trun, exp_tend)
                assert qrun == exp_qend or exp_qend == -1, 'Missed expected block end: {} vs {}'.format(qrun, exp_qend)
                parts = read_head(line)
                assert parts[2] == '+', 'Reverse target chain is unexpected: {}'.format(line)
                tchrom = parts[0]
                if not (tchrom == select or select == 'all'):
                    skip = True
                    continue
                skip = False
                trun = parts[3]
                exp_tend = parts[4]
                qchrom = parts[5]
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
                    trun += size + dt
                    qrun += size + dq
                except ValueError:
                    size = read_end(line)
                    yield tchrom, trun, trun + size, tstrand, qchrom, qrun, qrun + size, qstrand
                    trun += size
                    qrun += size
    return
