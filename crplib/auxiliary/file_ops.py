# coding=utf-8

"""
Convenience module for file operations
"""

import os as os
import gzip as gz
import bz2 as bz
import tempfile as tempf
import hashlib as hsl
import numpy as np
import pandas as pd

from crplib.auxiliary.constants import LIMIT_SERIALIZATION


def text_file_mode(fpath, read=True):
    """
    Naive determination of file type and respective open function and mode

    :param fpath:
    :param read:
    :return:
    """
    if read:
        assert os.path.isfile(fpath), 'Invalid path to file: {}'.format(fpath)
    ext = fpath.split('.')[-1].lower()
    # TODO should read magic numbers instead...
    if ext in ['gz', 'gzip']:
        f, m = gz.open, 'rt'  # Python 3.4+
    elif ext in ['bz', 'bz2', 'bzip', 'bzip2']:
        f, m = bz.open, 'rt'
    else:
        f, m = open, 'r'
    if not read:
        m = m.replace('r', 'w')
    return f, m


def shm_file_object(fpath):
    """
    :param fpath:
    :return:
    """
    assert os.path.isfile(fpath)
    ext = fpath.split('.')[-1].lower()
    if ext in ['gz', 'gzip']:
        f, t = gz.GzipFile, 'gzip'
    elif ext in ['bz', 'bz2', 'bzip', 'bzip2']:
        f, t = bz.BZ2File, 'bzip'
    else:
        # assume plain binary
        f, t = None, 'raw'
    return f, t


def create_filepath(fpath, logger=None):
    """
    :param fpath:
    :return:
    """
    dirs, filename = os.path.split(fpath)
    if not dirs:  # empty could mean simply CWD
        return fpath
    if not os.path.isdir(dirs):
        try:
            os.makedirs(dirs, exist_ok=True)
        except (IOError, OSError) as excp:
            if logger is not None:
                logger.error('Could not create path to file: {}'.format(fpath))
            raise excp
    return fpath


def check_array_serializable(data):
    """
    :param data:
    :return:
    """
    # to make this more or less generic,
    # no duck-typing this time...
    if hasattr(data, 'nbytes'):
        objsize = data.nbytes
        handle = 'numpy'
    # then presume it is a Pandas DataFrame or Series
    elif hasattr(data, 'values'):
        objsize = data.values.nbytes + data.index.nbytes
        handle = 'pandas'
    else:
        raise TypeError('Cannot determine object size - missing attribute nbytes for type: {}'.format(type(data)))
    if objsize < LIMIT_SERIALIZATION:
        return data
    if handle == 'numpy':
        file_buffer = tempf.NamedTemporaryFile('wb', delete=False, suffix='.tmp', prefix='np_mmap_')
        fp = np.memmap(file_buffer, dtype=data.dtype, mode='w+', shape=data.shape)
        fp[:] = data[:]
        # triggers flushing to disk
        del fp
        # return info needed to read data from file buffer
        recov_info = (file_buffer.name, data.dtype, data.shape)
    elif handle == 'pandas':
        file_buffer = tempf.NamedTemporaryFile('wb', delete=False, suffix='.tmp', prefix='pd_hdf_')
        with pd.HDFStore(file_buffer.name, 'w') as hdf:
            hdf.put('tmp', data, format='fixed')
            hdf.flush()
        recov_info = (file_buffer.name, '', '')
    else:
        raise TypeError('Undefined usecase for file buffer')
    return recov_info


def load_mmap_array(tmpfn, handle, dtype=None, shape=None):
    """
    :param tmpfn:
    :param dtype:
    :param shape:
    :return:
    """
    if handle == 'numpy':
        assert dtype is not None and shape is not None, 'Need datatype and shape for numpy mmap file buffer'
        fp = np.memmap(tmpfn, dtype=dtype, shape=shape, mode='r')
        data = np.zeros(dtype=dtype, shape=shape)
        data[:] = fp[:]
        del fp
        assert (data > 0).any(), 'Read all zero data from buffer {}'.format(tmpfn)
    elif handle == 'pandas':
        with pd.HDFStore(tmpfn, 'r') as hdf:
            data = hdf['/tmp']
    else:
        raise TypeError('Undefined usecase for recovering file buffer')
    rm_success = False
    try:
        os.remove(tmpfn)
        rm_success = True
    except (IOError, OSError, FileNotFoundError):
        pass
    return data, rm_success


def get_checksum(fpath, hashtype='md5', blocksize=4096):
    """
    :param fpath:
    :param hashtype:
    :param blocksize:
    :return:
    """
    assert os.path.isfile(fpath), 'Cannot compute checksum for file {} - invalid path: {}'.format(os.path.basename(fpath), fpath)
    select_hash = {'md5': hsl.md5}
    hasher = select_hash[hashtype]()
    with open(fpath, 'rb') as f:
        buf = f.read(blocksize)
        while buf:
            hasher.update(buf)
            buf = f.read(blocksize)
    return hashtype, hasher.hexdigest()
