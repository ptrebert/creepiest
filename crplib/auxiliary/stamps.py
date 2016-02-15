# coding=utf-8

"""
This module offers some convenience wrappers
around random and datetime functionality to
create date, time and random character stamps
for filenames, object identifiers etc.
"""

import datetime as dt
import random as rand
import string as string


def get_date_stamp(fmt='%Y%m%d'):
    """
    :param fmt:
     :type: str
    :return:
     :rtype: string
    """
    dtobj = dt.datetime.now()
    stamp = dtobj.strftime(fmt)
    return stamp


def get_time_stamp(fmt='%H%M%S'):
    """ No POSIX timestamp
    :param fmt: string format expression
     :type: str
    :return:
     :rtype: string
    """
    dtobj = dt.datetime.now()
    stamp = dtobj.strftime(fmt)
    return stamp


def get_char_stamp(upper=False, length=4):
    """
    :param upper:
     :type: bool
    :param length:
    :return:
     :rtype: string
    """
    if upper:
        charseq = string.ascii_uppercase
    else:
        charseq = string.ascii_lowercase
    stamp = ''.join([rand.choice(charseq) for _ in range(length)])
    return stamp


def get_full_filestamp(sep='_', clen=4):
    """
    :param clen:
    :return:
    """
    dstamp = get_date_stamp()
    tstamp = get_time_stamp()
    chars = get_char_stamp(length=clen)
    full_stamp = sep.join([dstamp, tstamp, chars])
    return full_stamp
