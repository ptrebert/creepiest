# coding=utf-8

import os as os
import collections as col
import numpy as np
import pandas as pd

from crplib.auxiliary.hdf_ops import get_default_group,\
    get_chrom_list
from crplib.numalg.normalization import nonzero_qnorm


def check_input_files(filepaths, old_ext, new_ext, outdir, overwrite, logger):
    """
    :param filepaths:
    :param old_ext:
    :param new_ext:
    :param overwrite:
    :param outdir:
    :param logger:
    :return:
    """
    annotation = col.defaultdict(list)
    for fp in filepaths:
        fn = os.path.basename(fp)
        try:
            root = get_default_group(fp)
        except AssertionError as ae:
            logger.error('File {} contains several groups per chromosome - '
                         'multi-sample normalization supports only one group '
                         'per chromosome per file.'.format(fp))
            raise ae
        else:
            chroms_in_file = get_chrom_list(fp, verify=False)
            new_fn = fn.replace(old_ext, new_ext)
            if outdir == 'as_input':
                outpath = os.path.dirname(fp)
            else:
                outpath = outdir
            new_filepath = os.path.join(outpath, new_fn)
            if fn == new_fn:
                if (fp == new_filepath) and overwrite:
                    logger.error('Path to new file {} is identical to old one {} and filemode is'
                                 ' set to overwrite [w] - this will result in data loss,'
                                 ' cannot proceed.'.format(new_filepath, fp))
                    raise AssertionError('Cannot overwrite original file with new file')
            for chrom in chroms_in_file:
                load_group = os.path.join(root, chrom)
                infos = (fn, fp, load_group, chrom, new_filepath)
                annotation[chrom].append(infos)
    return annotation


def load_data(fileinfos):
    """
    :param fileinfos:
    :return:
    """
    _, path, group, _, _ = fileinfos[0]
    with pd.HDFStore(path, 'r') as hdf:
        data = hdf[group]
        data_matrix = np.zeros((len(fileinfos), data.size), dtype=data.dtype)
        data_matrix[0, :] = data.values
    for idx, (_, path, group, _, _) in enumerate(fileinfos[1:], start=1):
        with pd.HDFStore(path, 'r') as hdf:
            data_matrix[idx, :] = hdf[group].values
    return data_matrix





def run_normalization(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    overwrite = args.filemode == 'w'
    chrom_annot = check_input_files(args.inputfile, args.replace, args.suffix,
                                    args.outdir, overwrite, logger)
    logger.debug('Input files annotated, found {} chromosomes'.format(len(list(chrom_annot.keys()))))
    filemode = args.filemode
    for chrom, fileinfos in chrom_annot.items():
        logger.debug('Processing chromosome {}'.format(chrom))
        ordered = sorted(fileinfos, key=lambda x: x[0])
        if len(ordered) == 1:
            logger.warning('Chromosome {} only present in sample {} - skipping'.format(chrom, ordered[0][0]))
            continue
        logger.debug('Loading and normalizing data...')
        norm_data = nonzero_qnorm(load_data(ordered))
        dump_data(norm_data, ordered, filemode)
        logger.debug('Finished processing chromosome {}'.format(chrom))
        filemode = filemode.replace('w', 'a')
    logger.debug('Normalization complete')
    return 0
