# coding=utf-8

import datetime as dt
import pandas as pd


MD_GENOMES_COLDEFS = ['group', 'name', 'ctime', 'mtime', 'size_mb', 'length', 'num_upper', 'num_lower']


def gen_obj_and_md(mdframe, assembly, chrom, seqbuf):
    """
    :param mdframe:
    :param assembly:
    :param chrom:
    :param seqbuf:
    :return:
    """
    grp = assembly + '/' + chrom
    slen = len(seqbuf.getvalue())
    num_low = sum(1 for c in seqbuf.getvalue() if c.islower())
    num_upp = slen - num_low
    ctime = dt.datetime.now()
    seqobj = pd.Series(data=list(seqbuf.getvalue()))
    size_mem = (seqobj.values.nbytes + seqobj.index.nbytes) / (1024 * 1024)
    mtime = dt.datetime.now()
    entries = [grp, chrom, ctime, mtime, int(size_mem), slen, num_upp, num_low]
    tmp = pd.DataFrame([entries, ], columns=MD_GENOMES_COLDEFS)
    mdframe = mdframe.append(tmp, ignore_index=True)
    return grp, seqobj, mdframe
