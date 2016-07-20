# coding=utf-8

import os as os
import collections as col
import json as json
import multiprocessing as mp
import pandas as pd

from string import ascii_uppercase as asciiupper

from crplib.auxiliary.hdf_ops import get_default_group, get_chrom_list, load_data_group
from crplib.auxiliary.file_ops import create_filepath


def collect_all_chroms(filepaths):
    """
    :param filepaths:
    :return:
    """
    chroms = set()
    for fp in filepaths:
        chroms |= set(get_chrom_list(fp))
    return chroms


def assemble_overlap_params(args):
    """
    :param args:
    :return:
    """
    commons = dict(vars(args))
    del commons['module_logger']
    del commons['execute']
    inputfiles = []
    inputlabels = dict()
    inputgroups = dict()
    deflabels = list(asciiupper)
    defidx = 0
    for inpf in sorted(commons['inputfile']):
        if os.path.isfile(inpf):
            try:
                label = deflabels[defidx]
                defidx += 1
            except IndexError:
                raise AssertionError('List of default labels exhausted; please specify less than 27 unlabeled'
                                     ' files for region overlap evaluation.')
            group = get_default_group(inpf)
            fp = inpf
        else:
            assert inpf.count(':') >= 2, 'Cannot distinguish between label and group for file {} - two ":"' \
                                         ' preceding the filepath are needed'.format(inpf)
            label, group, fp = inpf.split(':', 2)
            assert os.path.isfile(fp), 'Invalid path to file specified: {}'.format(fp)
            if not label.strip():
                try:
                    label = deflabels[defidx]
                    defidx += 1
                except IndexError:
                    raise AssertionError('List of default labels exhausted; please specify less than 27 unlabeled'
                                         ' files for region overlap evaluation.')
            if not group.strip():
                group = get_default_group(fp)
        inputfiles.append(fp)
        inputgroups[fp] = group.strip()
        inputlabels[fp] = label.strip()
    commons['inputfile'] = inputfiles
    commons['inputgroup'] = inputgroups
    commons['inputlabel'] = inputlabels
    if commons['roifile']:
        try:
            grp, fp = commons['roifile'].split(':')
        except ValueError:
            fp = commons['roifile']
            assert os.path.isfile(fp), 'Invalid path to ROI file specified: {}'.format(fp)
            grp = get_default_group(fp)
        commons['roifile'] = fp
        commons['roigroup'] = grp
    all_chroms = collect_all_chroms(inputfiles)
    arglist = []
    for c in all_chroms:
        tmp = dict(commons)
        tmp['chrom'] = c
        arglist.append(tmp)
    return arglist


def reg_ovl(row, start, end):
    """
    :param row:
    :param start:
    :param end:
    :return:
    """
    return min(row.end, end) - max(row.start, start)


def compute_overlap(data_a, label_a, data_b, label_b):
    """
    :param data_a:
    :param data_b:
    :return:
    """

    start_min_b = data_b.start.min()
    end_max_b = data_b.end.max()
    hits = col.defaultdict(set)
    cov_bp = 0
    bpovl = reg_ovl
    for row in data_a.itertuples():
        if row.end <= start_min_b:
            continue
        elif row.start > end_max_b:
            break
        else:
            # geq and leq are necessary for cases where two regions exactly coincide
            ovl = data_b.query('start >= {} and end <= {}'.format(row.start, row.end))
            if ovl.empty:
                continue
            hits[label_a].add(row.name)
            hits[label_b] |= set(ovl.name.tolist())
            ovl_bp = (ovl.apply(bpovl, axis=1, args=(row.start, row.end))).sum()
            cov_bp += ovl_bp
    return len(hits[label_a]), len(hits[label_b]), cov_bp


def get_overlap_stats(data_a, label_a, data_b, label_b, roidata):
    """
    :param data_a:
    :param label_a:
    :param data_b:
    :param label_b:
    :param roidata:
    :return:
    """
    res = {'roi_hits': 0, 'roi_miss': 0,
           label_a + '_num': 0, label_a + '_cov': 0,
           label_b + '_num': 0, label_b + '_cov': 0}
    # BED convention: half open intervals, ( ]
    ovl_qry = 'start > {} and end <= {}'
    if roidata is not None:
        roi_hits = 0
        tmp_a = pd.DataFrame(columns=data_a.columns)
        tmp_b = pd.DataFrame(columns=data_b.columns)
        for roi in roidata.itertuples():
            roihit = False
            tmp = data_a.query(ovl_qry.format(roi.start, roi.end))
            if not tmp.empty:
                roihit = True
            tmp_a = pd.concat([tmp_a, tmp], ignore_index=True, axis=0)
            tmp = data_b.query(ovl_qry.format(roi.start, roi.end))
            if not tmp.empty:
                roihit = True
            tmp_b = pd.concat([tmp_b, tmp], ignore_index=True, axis=0)
            if roihit:
                roi_hits += 1
        tmp_a.drop_duplicates(inplace=True)
        tmp_b.drop_duplicates(inplace=True)
        res['roi_hits'] = roi_hits
        res['roi_miss'] = roidata.shape[0] - roi_hits
    else:
        tmp_a = data_a
        tmp_b = data_b
    res[label_a + '_miss'] = tmp_a.shape[0]
    res[label_b + '_miss'] = tmp_b.shape[0]
    if not tmp_a.empty:
        res[label_a + '_num'] = tmp_a.shape[0]
        res[label_a + '_cov'] = int((tmp_a.end - tmp_a.start).sum())
    if not tmp_b.empty:
        res[label_b + '_num'] = tmp_b.shape[0]
        res[label_b + '_cov'] = int((tmp_b.end - tmp_b.start).sum())
    if not (tmp_a.empty or tmp_b.empty):
        hits_a, hits_b, bpcov = compute_overlap(tmp_a, label_a, tmp_b, label_b)
        res[label_a + '_hits'] = hits_a
        res[label_a + '_miss'] = tmp_a.shape[0] - hits_a
        res[label_b + '_hits'] = hits_b
        res[label_b + '_miss'] = tmp_b.shape[0] - hits_b
        res['union_cov_bp'] = int(res[label_a + '_cov'] + res[label_b + '_cov'] - bpcov)
        res['ovl_cov_bp'] = int(bpcov)
        res['jaccard'] = bpcov / res['union_cov_bp']
    return res


def compute_pairwise_overlap(params):
    """
    :param params:
    :return:
    """
    chrom = params['chrom']
    infiles, ingroups, inlabels = params['inputfile'], params['inputgroup'], params['inputlabel']
    roidata = None
    res = {'chrom': chrom, 'roi_info': {}, 'stats': {}}
    if params['roifile']:
        res['roi_info']['roi_file'] = os.path.basename(params['roifile'])
        res['roi_info']['roi_group'] = os.path.basename(params['roigroup'])
        try:
            roidata = load_data_group(params['roifile'], params['roigroup'], chrom)
        except KeyError:
            res['roi_info']['status'] = 'no_roi_regions'
            return res

    done = set()
    for file_a in infiles:
        label_a = inlabels[file_a]
        data_a = load_data_group(file_a, ingroups[file_a], chrom, True)
        for file_b in infiles:
            label_b = inlabels[file_b]
            if label_a == label_b:
                continue
            if (label_a, label_b) in done:
                continue
            data_b = load_data_group(file_b, ingroups[file_b], chrom, True)
            stats = get_overlap_stats(data_a, label_a, data_b, label_b, roidata)
            res['stats'][label_a + '_vs_' + label_b] = stats
        done.add((label_a, label_b))
        done.add((label_b, label_a))
    return res


def run_overlap_evaluation(args, logger):
    """
    :param args:
    :return:
    """
    arglist = assemble_overlap_params(args)
    logger.debug('Assembled list of size {} for processing'.format(len(arglist)))

    results = []
    wg_counts = col.defaultdict(col.Counter)
    with mp.Pool(args.workers) as pool:
        resit = pool.imap_unordered(compute_pairwise_overlap, arglist, chunksize=1)
        for res in resit:
            logger.debug('Received results for {}'.format(res['chrom']))
            results.append(res)
            for comp, stats in res['stats'].items():
                wg_counts[comp].update(stats)
    logger.debug('All statistics collected')
    norm_dict = dict()
    for comp, stats in wg_counts.items():
        stats['jaccard'] = stats['ovl_cov_bp'] / stats['union_cov_bp']
        norm_dict[comp] = dict(stats)
    results.append({'chrom': 'wg', 'roi_info': {}, 'stats': norm_dict})
    results = sorted(results, key=lambda x: x['chrom'])
    create_filepath(args.outputfile, logger)
    logger.debug('Dumping output json')
    with open(args.outputfile, 'w', encoding='ascii') as outf:
        json.dump(results, outf, indent=1, ensure_ascii=True)
    return 0


def run_evaluation(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    logger.debug('Starting evaluation')
    if args.task == 'overlap':
        assert len(args.inputfile) >= 2, 'Need to specify at least two files for region overlap evaluation.' \
                                         ' You specified {}, human.'.format(len(args.inputfile))
        logger.debug('Assessing overlap between region sets')
        run_overlap_evaluation(args, logger)
    else:
        raise ValueError('Unknown task specified: {}'.format(args.task))
    return 0
