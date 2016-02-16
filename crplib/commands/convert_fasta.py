# coding=utf-8

"""
Module to convert FASTA nucleotide sequences to HDF5 format
"""

import io as io
import re as re
import functools as fnt
import pandas as pd

from crplib.auxiliary.text_parsers import read_chromosome_sizes
from crplib.auxiliary.file_ops import text_file_mode
from crplib.metadata.md_genomes import gen_obj_and_md, MD_GENOMES_COLDEFS


def _replace_n(mobj):
    """
    :param mobj:
    :return:
    """
    if mobj.group().islower():
        return 'n'
    return 'N'


def _line_proc_fun(norep):
    """
    :param norep:
    :return:
    """
    if norep:
        f = lambda s: s
    else:
        pat = re.compile('[^ACGTNacgtn]')
        f = fnt.partial(pat.sub, repl=_replace_n)
    return f


def run_fasta_conversion(args):
    """
    :param args:
    :return:
    """
    csizes = read_chromosome_sizes(args.chromsizes, args.keepchroms)
    opn, mode = text_file_mode(args.input)
    proc = _line_proc_fun(args.noreplace)

    metadata = pd.DataFrame(columns=MD_GENOMES_COLDEFS)
    with pd.HDFStore(args.output) as hdfout:
        with opn(args.input, mode=mode, encoding='utf-8') as infile:
            seqbuf = io.StringIO()
            chrom = None
            for line in infile:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    try:
                        cs = csizes[chrom]
                        seql = len(seqbuf.getvalue())
                        assert cs == seql, \
                            'Length mismatch for chrom {}: annotated {} - fasta seq. {}'.format(chrom, cs, seql)
                        grp, seqobj, metadata = gen_obj_and_md(metadata, args.assembly, chrom, seqbuf)
                        hdfout.put(grp, seqobj)
                        seqbuf = io.StringIO()
                        chrom = line.strip().strip('>').split()[0]
                    except KeyError:
                        seqbuf = io.StringIO()
                        chrom = line.strip().strip('>').split()[0]
                else:
                    seqbuf.write(proc(line.strip()))
        hdfout.put('metadata', metadata)
    return 0
