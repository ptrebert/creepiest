#!/usr/bin/env python3

import os as os
import sys as sys
import argparse as argp
import traceback as trb
import re as re
import gzip as gz
import operator as op
import functools as fnt
import json as js
import numpy as np
import numpy.ma as msk
import io as io
import pickle as pck

from crplib.auxiliary.text_parsers import get_chain_iterator, chromsize_from_chain


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(add_help=True)
    parser.add_argument('--task', '-tk', type=str, choices=['chaintobed', 'qfilter', 'symmext',
                                                            'maptobedgraph', 'swap'], dest='task')

    comgroup = parser.add_argument_group('File arguments')
    comgroup.add_argument('--chain-file', '-chf', type=str, dest='chainfile',
                          help='A standard UCSC chain file')
    comgroup.add_argument('--map-file', '-mpf', type=str, dest='mapfile',
                          help='A one-to-one mapping file created from a chain file')
    comgroup.add_argument('--output-file', '-otf', type=str, dest='output')

    comgroup = parser.add_argument_group('Task: Dump chain to BED [chaintobed]')
    comgroup.add_argument('--no-buffer', '-nob', action='store_true', default=False, dest='nobuffer')

    comgroup = parser.add_argument_group('Task: Filter query for symmetrical coverage [qfilter]')
    comgroup.add_argument('--chrom', '-chr', type=str, dest='chrom')

    comgroup = parser.add_argument_group('Task: Check for mergeable blocks in a mapping [symmext]')
    comgroup.add_argument('--query-sizes', '-qcs', type=str, dest='qsizes')
    comgroup.add_argument('--gap-size', '-gap', type=int, default=200, dest='gapsize')
    comgroup.add_argument('--cov-ratio', '-cov', type=float, default=0.7, dest='covratio')
    comgroup.add_argument('--min-cov', '-mcv', type=int, default=100, dest='mincov')
    comgroup.add_argument('--check-ext', '-cex', type=int, nargs='+', default=[100, 50, 25], dest='checkext')

    comgroup = parser.add_argument_group('Task: Dump mapping to bedGraph [maptobedgraph]')
    comgroup.add_argument('--query', '-qry', action='store_true', default=False,
                          help='Dump (unsorted) query coordinates; otherwise output is target coordinates')

    args = parser.parse_args()
    assert max(args.checkext) < args.gapsize,\
        'Maximal extension has to be shorter than gap size: {} vs {}'.format(max(args.checkext), args.gapsize)
    setattr(args, 'checkext', sorted(args.checkext, reverse=True))

    return args


def chain_query_overlap_filter(args):
    """
    :param args:
    :return:
    """
    chromsize = chromsize_from_chain(args.chainfile, args.chrom, False)
    # tchrom, trun, trun + size, tstrand
    #   0       1        2          3
    # qchrom, qrun, qrun + size, qstrand
    #   4       5        6          7
    # chain_id, chain_score
    #   8          9
    get_target = op.itemgetter(*(0, 1, 2, 3))
    get_query = op.itemgetter(*(4, 5, 6, 7))
    get_blockid = op.itemgetter(*(8, ))
    get_score = op.itemgetter(*(9, ))

    iva = msk.masked_array(np.arange(chromsize, dtype=np.int32))
    iva.mask = False
    contained = 0
    discard_short = 0
    discard_slice = 0
    partial_recov = 0
    symm_blocks = 0
    recov_len = 0
    total = 0
    block_buffer = io.StringIO()
    with gz.open(args.chainfile, 'rt') as infile:
        chainit = get_chain_iterator(infile, qselect=re.compile(args.chrom + '$'), min_size=args.mincov, min_score=1000)
        last_score = sys.maxsize
        for block in chainit:
            assert get_score(block) <= last_score, 'Chains not sorted in descending order for score: {}'.format(block)
            last_score = get_score(block)
            total += 1
            qs, qe = get_query(block)[1], get_query(block)[2]
            if iva[qs:qe].count() == 0:  # all values masked
                contained += 1
            elif msk.count_masked(iva[qs:qe]) == 0:
                iva[qs:qe].mask = 1
                block_buffer.write('\t'.join(map(str, get_target(block))))
                block_buffer.write('\t' + get_blockid(block) + '\t')
                block_buffer.write('\t'.join(map(str, get_query(block))) + '\n')
                symm_blocks += 1
            elif qe - qs <= 25:
                discard_short += 1
            else:
                unmask_slices = msk.clump_unmasked(iva[qs:qe])
                max_slice = sorted(unmask_slices, key=lambda x: x.stop - x.start, reverse=True)[0]
                if max_slice.stop - max_slice.start < 25:
                    discard_slice += 1
                    continue
                else:
                    partial_recov += 1
                    ads = abs(qs - (qs + max_slice.start))
                    ade = abs(qe - (qs + max_slice.stop))
                    rlen = (qe - ade) - (qs + ads)
                    recov_len += rlen
                    iva[qs + ads:qe - ade].mask = 1
                    treg = list(get_target(block))
                    treg[1] += ads
                    treg[2] -= ade
                    qreg = list(get_query(block))
                    qreg[1] += ads
                    qreg[2] -= ade
                    assert qreg[2] - qreg[1] > 1, 'Adapting block edges resulted in block of length <= 1: {}'.format(block)
                    block_buffer.write('\t'.join(map(str, treg)))
                    block_buffer.write('\t' + get_blockid(block) + '\t')
                    block_buffer.write('\t'.join(map(str, qreg)) + '\n')
                    continue
    chrom_cov = int(iva.mask.sum())
    metadata = {'chainfile': args.chainfile, 'chromosome': args.chrom,
                'chrom_size': chromsize, 'chrom_pct_cov': round(chrom_cov / chromsize, 3),
                'total_blocks': total, 'symm_blocks': symm_blocks,
                'symm_recov_cov': chrom_cov, 'recov_blocks': partial_recov,
                'recov_cov': int(recov_len), 'discarded_short': discard_short,
                'discarded_slice': discard_slice}
    stat_out = args.output.replace('.tsv.gz', '.json')
    with open(stat_out, 'w') as outf:
        _ = js.dump(metadata, outf, indent=1, sort_keys=True)

    with gz.open(args.output, 'wt') as outf:
        _ = outf.write(block_buffer.getvalue())

    return 0


def build_query_chrom_masks(chromsizes):
    """
    :param chromsizes:
    :return:
    """
    chrommasks = dict()
    selector = re.compile('chr[0-9]+$')
    with open(chromsizes, 'r') as infile:
        for line in infile.readlines():
            c, s = line.strip().split()[:2]
            if selector.match(c) is not None:
                chrommasks[c] = msk.masked_array(np.arange(int(s)), dtype=np.int32)
                chrommasks[c].mask = False
    assert chrommasks, 'No chromosomes were masked'
    return chrommasks


def get_mapfile_iterator(fobj):
    """
    :param fobj:
    :return:
    """
    pml = parse_map_line
    while 1:
        lines = fobj.readlines(32768)
        if not lines:
            break
        for line in lines:
            yield pml(line)
    return


def mask_query_blocks(mapfile, qmasks):
    """
    :param mapfile:
    :param qmasks:
    :return:
    """
    qc, qs, qe = 5, 6, 7
    with gz.open(mapfile, 'rt') as infile:
        lineiter = get_mapfile_iterator(infile)
        for line in lineiter:
            assert msk.count_masked(qmasks[line[qc]][line[qs]:line[qe]]) == 0,\
                'Overlapping query blocks at position: {}'.format(line)
            qmasks[line[qc]][line[qs]:line[qe]].mask = 1
    return qmasks


def parse_map_line(line):
    """
    :param line:
    :return:
    """
    tchrom, tstart, tend, tstrand, blockid, qchrom, qstart, qend, qstrand = line.strip().split()
    return tchrom, int(tstart), int(tend), tstrand, blockid.split('.')[0], qchrom, int(qstart), int(qend), qstrand


def check_block_merge(blocks, chrom):
    """ A test run of this function on hg19 / mm9 indicates that a straightforward
    merge of consecutive blocks is not possible (due to varying gap sizes)

    Hence, this function is deprecated / useless and left here just as a reminder

    :param blocks:
    :param chrom:
    :return:
    """
    if len(blocks) < 2:
        return blocks, False
    outblocks = []
    merged = False
    for idx in range(len(blocks)):
        subset = blocks[idx:]
        if len(subset) == 1:
            outblocks.append(subset[0])
            break
        elif (subset[-1][2] - subset[0][1]) == (subset[-1][7] - subset[0][6]) and \
                sum(s[7] - s[6] for s in subset) == chrom[subset[0][6]:subset[0][7]].mask.sum():
            # same length in target and query and no intervening blocks in query
            outblocks.append((subset[0][0], subset[0][1], subset[-1][2], subset[0][3],
                              subset[0][4],
                              subset[0][5], subset[0][6], subset[-1][7], subset[0][8]))
            chrom[subset[0][6]:subset[-1][7]].mask = 1
            merged = True
            break
        else:
            outblocks.append(blocks[idx])
    assert outblocks, 'No output generated for input of length {}'.format(len(blocks))
    return outblocks, merged


def check_gap_extendable(covratio, mincov, checkext, blocks, chrom):
    """
    :param blocks:
    :param chrom:
    :return:
    """

    ts = blocks[0][1]
    qs = blocks[0][6]
    te = blocks[-1][2]
    qe = blocks[-1][7]
    max_len = max(te - ts, qe - qs)
    tcov = sum([b[2] - b[1] for b in blocks])
    qcov = sum([b[7] - b[6] for b in blocks])
    assert tcov == qcov, 'Coverage not equal ({} vs {}) for blocks: {}'.format(tcov, qcov, blocks)
    if tcov < mincov:
        return blocks
    outblocks = []
    for ext in checkext:
        new_len = max_len // ext * ext + ext
        if chrom[qs:qs+new_len].mask.sum() == qcov:
            # no overlapping blocks from other chains in this gap, extend
            if qcov / (qs+new_len - qs) < covratio:
                continue
            outblocks = [(blocks[0][0], ts, ts+new_len, blocks[0][3],
                          blocks[0][4],
                          blocks[0][5], qs, qs+new_len, blocks[0][8])]
            chrom[qs:qs+new_len].mask = 1
            break
    if not outblocks:
        if chrom[qs:qs+max_len].mask.sum() == qcov:
            if not qcov / (qs+max_len - qs) < covratio:
                outblocks = [(blocks[0][0], ts, ts+max_len, blocks[0][3],
                              blocks[0][4],
                              blocks[0][5], qs, qs+max_len, blocks[0][8])]
                chrom[qs:qs+max_len].mask = 1
            else:
                outblocks = blocks
        else:
            outblocks = blocks
    return outblocks


def check_gap_extend(gapsize, covratio, mincov, checkext, blocks, chrom):
    """
    :param blocks:
    :param chrom:
    :return:
    """
    check_pos = []
    for idx, block in enumerate(blocks):
        try:
            if blocks[idx+1][1] - block[2] >= gapsize and \
                    blocks[idx+1][6] - block[7] >= gapsize:
                check_pos.append(idx)
        except IndexError:
            continue
    if not check_pos:
        return blocks
    outblocks = []
    start_pos = 0
    check_ext = fnt.partial(check_gap_extendable, *(covratio, mincov, checkext))
    for inc_end in check_pos:
        s, e = start_pos, inc_end + 1
        subset = blocks[s:e]
        if len(subset) == 1:
            # isolated singleton
            outblocks.extend(subset)
            start_pos = inc_end + 1
            continue
        assert subset, 'Invalid subset selection of blocks - indices {} to {}'.format(s, e)
        outblocks.extend(check_ext(subset, chrom))
        start_pos = inc_end + 1
    if start_pos < len(blocks):
        outblocks.extend(blocks[start_pos:])
    assert outblocks, 'No output created for input blocks: {}'.format(blocks)
    assert len(outblocks) == len(set(outblocks)), 'Duplicate blocks created: {}'.format(outblocks)
    return outblocks


def symm_extend_blocks(args):
    """
    :param args:
    :return:
    """
    mf_ext = get_file_extension(args, 'mapfile', True)
    maskfile = args.mapfile.replace(mf_ext, '.pck')
    if not os.path.isfile(maskfile):
        cm = build_query_chrom_masks(args.qsizes)
        qmasks = mask_query_blocks(args.mapfile, cm)
        pck.dump(qmasks, open(maskfile, 'wb'))
    else:
        maskfile = args.mapfile.replace(mf_ext, '.pck')
        qmasks = pck.load(open(maskfile, 'rb'))
    raw_cov = summarize_coverage(qmasks)
    out_ext = get_file_extension(args, 'output', True)
    raw_cov_stats = collect_block_statistics(qmasks, args.output.replace(out_ext, '.rawcov.pck'))
    pml = parse_map_line
    ext = 0
    maps = 0
    # empty output file in case there is
    # an incomplete file from a previous run
    with gz.open(args.output, 'wt') as outf:
        _ = outf.write('')
    outbuffer = []
    all_mrg_cov = 0
    check_gap = fnt.partial(check_gap_extend, *(args.gapsize, args.covratio, args.mincov, args.checkext))
    with gz.open(args.mapfile, 'rt') as infile:
        chainbuffer = []
        chain = ''
        chrom = ''
        for line in infile:
            l = pml(line)
            maps += 1
            if l[4] != chain:
                try:
                    chaincov = sum([b[2] - b[1] for b in chainbuffer])
                    outblocks = check_gap(chainbuffer, qmasks[chrom])
                    mrgcov = sum([b[2] - b[1] for b in outblocks])
                    all_mrg_cov += mrgcov
                    assert mrgcov >= chaincov, 'Coverage shrink - lost blocks: {}\n\n{}'.format(chainbuffer, outblocks)
                    outbuffer.extend(outblocks)
                    if len(outblocks) != len(chainbuffer):
                        ext += 1
                except KeyError:
                    pass
                chainbuffer = [l]
                chrom = l[5]
                chain = l[4]
                if len(outbuffer) > 1000000:
                    with gz.open(args.output, 'at') as outf:
                        for b in outbuffer:
                            _ = outf.write('\t'.join(map(str, b)) + '\n')
                    outbuffer = []
            else:
                chainbuffer.append(l)
        chaincov = sum([b[2] - b[1] for b in chainbuffer])
        outblocks = check_gap(chainbuffer, qmasks[chrom])
        mrgcov = sum([b[2] - b[1] for b in outblocks])
        all_mrg_cov += mrgcov
        assert mrgcov >= chaincov, 'Coverage shrink - lost blocks: {}\n\n{}'.format(chainbuffer, outblocks)
        outbuffer.extend(outblocks)
        if len(outblocks) != len(chainbuffer):
            ext += 1
    with gz.open(args.output, 'at') as outf:
        for b in outbuffer:
            _ = outf.write('\t'.join(map(str, b)) + '\n')
        outbuffer = []
    mrg_cov = summarize_coverage(qmasks)
    mrg_cov_stats = collect_block_statistics(qmasks, args.output.replace(out_ext, '.mrgcov.pck'))
    assert mrg_cov == all_mrg_cov, 'Coverage mismatch between mask and regions: {} - {}'.format(mrg_cov, all_mrg_cov)
    md = dict()
    md['num_maps'] = maps
    md['num_ext_blocks'] = ext
    md['raw_cov_bp'] = int(raw_cov)
    md['mrg_cov_bp'] = int(mrg_cov)
    md['mrg_ratio'] = round(mrg_cov / raw_cov, 5)
    md['input_mapfile'] = args.mapfile
    md['output_mapfile'] = args.output
    md['raw_cov_stats'] = raw_cov_stats
    md['mrg_cov_stats'] = mrg_cov_stats
    md['param_min_cov'] = args.mincov
    md['param_cov_ratio'] = args.covratio
    md['param_gap_size'] = args.gapsize
    md['param_check_ext'] = args.checkext
    with open(args.output.replace('.tsv.gz', '.json'), 'w') as metadata:
        js.dump(md, metadata, indent=1)
    assert not outbuffer, 'Out buffer not empty'
    return 0


def collect_block_statistics(qmasks, outfile):
    """
    :param qmasks:
    :return:
    """
    gapdist = []
    sizedist = []
    for chrom, masked in qmasks.items():
        sizes = msk.clump_masked(masked)
        sizedist.append([s.stop - s.start for s in sizes])
        tmp = msk.masked_array(masked.data, msk.logical_not(masked.mask))
        leftend, rightend = msk.flatnotmasked_edges(tmp)
        # ignore the regions to the chromosome boundaries
        gaps = msk.clump_masked(tmp[leftend:rightend+1])
        gapdist.append([g.stop - g.start for g in gaps])
    gapdist = np.concatenate(gapdist)
    sizedist = np.concatenate(sizedist)
    stats = {'gaps': gapdist, 'blocks': sizedist}
    with open(outfile, 'wb') as outf:
        pck.dump(stats, outf)
    return outfile


def summarize_coverage(chroms):
    """
    :param chroms:
    :return:
    """
    cov = 0
    for k, v in chroms.items():
        cov += v.mask.sum()
    return cov


def _block_to_bed_line(block):
    """ tchrom, trun, trun + size, tstrand, qchrom, qrun, qrun + size, qstrand, chain_id + '.' + str(bc), chain_score
    :param block:
    :return:
    """
    return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(*block)


def chain_to_bed_unbuffered(args):
    """
    :param args:
    :return:
    """
    chainit = get_chain_iterator(gz.open(args.chainfile, 'rt'), min_score=1000, min_size=args.mincov)
    bbl = _block_to_bed_line
    for block in chainit:
        sys.stdout.write(bbl(block))
    return 0


def chain_to_bed(args):
    """
    :param args:
    :return:
    """
    outbuffer = io.StringIO()
    chars = 0
    chainit = get_chain_iterator(gz.open(args.chainfile, 'rt'), min_score=1000, min_size=args.mincov)
    bbl = _block_to_bed_line
    for block in chainit:
        chars += outbuffer.write(bbl(block))
        if chars > 1000000:
            sys.stdout.write(outbuffer.getvalue())
            outbuffer = io.StringIO()
    sys.stdout.write(outbuffer.getvalue())
    return 0


def _target_map_to_bedgraph(line):
    """
    :param line:
    :return:
    """
    cols = line.strip().split()
    return '{}\t{}\t{}\t1\n'.format(cols[0], cols[1], cols[2])


def _query_map_to_bedgraph(line):
    """ Note to self: even though the strand can be minus,
    all coordinates are given in plus orientation by construction
    :param line:
    :return:
    """
    cols = line.strip().split()
    return '{}\t{}\t{}\t1\n'.format(cols[5], cols[6], cols[7])


def map_to_bedgraph(mapfile, query):
    """
    :param mapfile:
    :param query:
    :return:
    """
    if not query:
        readline = _target_map_to_bedgraph
    else:
        readline = _query_map_to_bedgraph
    with gz.open(mapfile, 'rt') as infile:
        for line in infile:
            sys.stdout.write(readline(line))
    return 0


def swap_map(mapfile):
    """
    :param mapfile:
    :return:
    """
    inverse = op.itemgetter(*(5, 6, 7, 3, 4, 0, 1, 2, 8))
    with gz.open(mapfile, 'rt') as infile:
        for line in infile:
            if not line:
                continue
            cols = line.strip().split()
            swl = '\t'.join(inverse(cols)) + '\n'
            sys.stdout.write(swl)
    return


def get_file_extension(args, which, compressed):
    """
    :param args:
    :param which:
    :param compressed:
    :return:
    """
    fpath = getattr(args, which)
    if compressed:
        if not fpath.endswith('.gz') or fpath.endswith('.gzip'):
            fpath += '.gz'
        fileext = '.' + '.'.join(fpath.split('.')[-2:])
    else:
        fileext = '.' + fpath.split('.')[-1]
    setattr(args, which, fpath)
    return fileext


if __name__ == '__main__':
    try:
        args = parse_command_line()
        if args.task == 'chaintobed':
            if args.nobuffer:
                _ = chain_to_bed_unbuffered(args)
            else:
                _ = chain_to_bed(args)
        elif args.task == 'qfilter':
            if not (args.output.endswith('.gz') or args.output.endswith('.gzip')):
                setattr(args, 'output', args.output + '.gz')
            _ = chain_query_overlap_filter(args)
        elif args.task == 'symmext':
            _ = symm_extend_blocks(args)
        elif args.task == 'maptobedgraph':
            _ = map_to_bedgraph(args.mapfile, args.query)
        elif args.task == 'swap':
            _ = swap_map(args.mapfile)
        else:
            raise ValueError('Unknown task: {}'.format(args.task))
    except Exception:
        trb.print_exc()
        sys.exit(1)
    else:
        sys.exit(0)
