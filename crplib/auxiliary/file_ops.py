# coding=utf-8

"""
Convenience module for file operations
"""

import os as os
import gzip as gz
import bz2 as bz


def text_file_mode(fpath):
    """ Naive determination of file type and appropriate selection
    of opening function and mode
    :param fpath:
    :return:
    """
    assert os.path.isfile(fpath), 'Invalid path to file: {}'.format(fpath)
    ext = fpath.split('.')[-1].lower()
    # TODO should read magic numbers instead...
    if ext in ['gz', 'gzip']:
        f, m = gz.open, 'rt'  # Python 3.4+
    elif ext in ['bz', 'bz2', 'bzip', 'bzip2']:
        f, m = bz.open, 'rt'
    else:
        f, m = open, 'r'
    return f, m

