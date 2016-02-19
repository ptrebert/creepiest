# coding=utf-8

"""
Module to convert FASTA nucleotide sequences to HDF5 format
"""

import os as os
import io as io
import re as re
import functools as fnt
import pandas as pd

from crplib.auxiliary.text_parsers import read_chromosome_sizes
from crplib.auxiliary.file_ops import text_file_mode
from crplib.metadata.md_genomes import gen_obj_and_md, MD_GENOMES_COLDEFS


def _replace_n(mobj):
    """
    :param mobj: match object for non ACGTN chars in sequence
     :type: re.match
    :return:
     :rtype: string
    """
    if mobj.group().islower():
        return 'n'
    return 'N'


def _line_proc_fun(norep):
    """
    :param norep:
     :type: bool
    :return: identity function or function to replace non ACGTN in argument
     :rtype: function, str -> str
    """
    if norep:
        f = lambda s: s
    else:
        f = fnt.partial(re.sub, *('[^ACGTNacgtn]', _replace_n))
    return f


def run_fasta_conversion(args, logger):
    """
    :param args:
    :return:
    """
    csizes = read_chromosome_sizes(args.chromsizes, args.keepchroms)
    opn, mode = text_file_mode(args.inputfile)
    proc = _line_proc_fun(args.noreplace)
    logger.warning('Warning: FASTA to HDF5 conversion is deprecated')
    metadata = pd.DataFrame(columns=MD_GENOMES_COLDEFS)
    with pd.HDFStore(args.outputfile, 'a', complevel=9, complib='blosc') as hdfout:
        with opn(args.inputfile, mode=mode, encoding='ascii') as infile:
            seqbuf = io.StringIO()
            chrom = None
            csize = 0
            skip = False
            for line in infile:
                if line.startswith('>'):
                    slen = len(seqbuf.getvalue())
                    if slen > 0:
                        assert csize == slen, 'Length mismatch: {} - {} - {}'.format(chrom, slen, csize)
                        logger.debug('Handling chromosome: {}'.format(chrom))
                        grp, seqobj, metadata = gen_obj_and_md(metadata, args.assembly, chrom,
                                                               os.path.basename(args.inputfile), seqbuf)
                        hdfout.put(grp, seqobj, format='fixed')
                        seqbuf = io.StringIO()
                    new_chrom = line.strip('>').split()[0].strip()
                    try:
                        csize = csizes[new_chrom]
                        chrom = new_chrom
                        skip = False
                        continue
                    except KeyError:
                        chrom = None
                        csize = 0
                        skip = True
                        continue
                elif skip or not line.strip():
                    continue
                else:
                    seqbuf.write(proc(line.strip()))
            if csize > 0:
                slen = len(seqbuf.getvalue())
                logger.debug('Handling chromosome: {}'.format(chrom))
                assert csize == slen, 'Length mismatch: {} - {} - {}'.format(chrom, slen, csize)
                grp, seqobj, metadata = gen_obj_and_md(metadata, args.assembly, chrom,
                                                       os.path.basename(args.inputfile), seqbuf)
                hdfout.put(grp, seqobj, format='fixed')
        hdfout.put('metadata', metadata, format='table')
    return 0
