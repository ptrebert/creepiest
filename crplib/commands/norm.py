# coding=utf-8

import os as os
import collections as col
import numpy as np
import pandas as pd

from crplib.auxiliary.hdf_ops import get_default_group,\
    get_chrom_list
from crplib.auxiliary.file_ops import create_filepath
from crplib.numalg.normalization import nonzero_qnorm
from crplib.metadata.md_signal import MD_SIGNAL_COLDEFS, gen_obj_and_md
from crplib.auxiliary.constants import DIV_B_TO_GB


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
            create_filepath(outpath, logger)
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


def load_data(fileinfos, logger):
    """
    :param fileinfos:
    :return:
    """
    _, path, group, _, _ = fileinfos[0]
    with pd.HDFStore(path, 'r') as hdf:
        data = hdf[group]
        data_matrix = np.zeros((len(fileinfos), data.size), dtype=data.dtype)
        data_matrix[0, :] = data.values
    logger.debug('First dataset loaded, memory allocated')
    for idx, (_, path, group, _, _) in enumerate(fileinfos[1:], start=1):
        with pd.HDFStore(path, 'r') as hdf:
            data_matrix[idx, :] = hdf[group].values
    logger.debug('Data loading complete')
    logger.debug('Matrix dimension: {}'.format(data_matrix.shape))
    logger.debug('Approx. size in memory: {} GiB'.format(round(data_matrix.nbytes / DIV_B_TO_GB, 2)))
    return data_matrix


def dump_data(fileinfos, datamatrix, outgroup, filemode, logger):
    """
    :param fileinfos:
    :param datamatrix:
    :param outgroup:
    :return:
    """

    for idx, (oldname, oldpath, oldgroup, chrom, newpath) in enumerate(fileinfos):
        logger.debug('Writing {} to file {}'.format(chrom, newpath))
        with pd.HDFStore(newpath, filemode, complevel=9, complib='blosc') as hdf:
            if 'metadata' in hdf:
                metadata = hdf['metadata']
            else:
                metadata = pd.DataFrame(columns=MD_SIGNAL_COLDEFS)
            grp, valobj, metadata = gen_obj_and_md(metadata, outgroup, chrom, oldpath, datamatrix[idx, :])
            hdf.put(grp, valobj, format='fixed')
            hdf.put('/metadata', metadata, format='table')
            hdf.flush()
    logger.debug('All data saved')
    return


def create_sym_link(args, overwrite, logger):
    """
    :param args:
    :param overwrite:
    :param logger:
    :return:
    """
    fp, fn = os.path.split(args.inputfiles[0])
    new_fn = fn.replace(args.replace, args.suffix)
    if args.outdir == 'as_input':
        new_fp = fp
    else:
        new_fp = args.outdir
    old_inpath = args.inputfiles[0]
    new_outpath = os.path.join(new_fp, new_fn)
    if old_inpath == new_outpath:
        # cannot create link to file itself, leave a message and be done with it
        logger.debug('New file path is identical to old file path - cannot create sym link: '
                     '{} - {}'.format(old_inpath, new_outpath))
    else:
        if os.path.islink(new_outpath):
            if overwrite:
                try:
                    os.unlink(new_outpath)
                    os.link(old_inpath, new_outpath)
                except (IOError, OSError) as err:
                    logger.warning('Could not update link - original error message: {}'.format(err))
            else:
                # running in "append" mode, so leave link untouched
                pass
        elif os.path.isfile(new_outpath):
            logger.debug('Sym link path points to regular file - did you manually copy the input file, human?'
                         ' Here it is: {}'.format(new_outpath))
        else:
            try:
                os.link(old_inpath, new_outpath)
            except (IOError, OSError) as err:
                logger.error('Could not create sym link - original error message: {}'.format(err))
                raise err
    return


def run_normalization(args):
    """
    :param args:
    :return:
    """
    logger = args.module_logger
    overwrite = args.filemode == 'w'
    if len(args.inputfiles) == 1:
        chrom_annot = col.defaultdict(list)
        logger.debug('Single input file - no normalization possible, check sym link option...')
        if args.symlink:
            create_sym_link(args, overwrite, logger)
    else:
        chrom_annot = check_input_files(args.inputfiles, args.replace, args.suffix,
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
        norm_data = nonzero_qnorm(load_data(ordered, logger), args.workers)
        dump_data(ordered, norm_data, args.outputgroup, filemode, logger)
        logger.debug('Finished processing chromosome {}'.format(chrom))
        filemode = filemode.replace('w', 'a')
    logger.debug('Normalization complete')
    return 0
