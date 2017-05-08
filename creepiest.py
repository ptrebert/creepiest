#!/usr/bin/env python3
# coding=utf-8

"""
Main executable for CREEPIEST framework
"""

import sys as sys
import os as os
import io as io
import logging as logging
import traceback as trb
import argparse as argp
import configparser as confp
import select as select
import datetime as dt
import random as rand
import cProfile as cprf
from argparse import ArgumentError

from crplib.__init__ import __version__
from crplib.__init__ import __description__ as prog_desc
from crplib.__init__ import __default_log_long__ as log_long
from crplib.__init__ import __default_log_short__ as log_short

from crplib.commands.sub_parsers import add_sub_parsers
from crplib.auxiliary.email_notify import send_notification
from crplib.auxiliary.stamps import get_full_filestamp


def check_stdin_config():
    """ See if we get a configuration string via stdin
    Note: this will definitely NOT work under windows
    :return: a configuration obtained via stdin
     :rtype: list
    """
    stdincfg = ''
    r_ok, _, _ = select.select([sys.stdin], [], [], 0.5)
    if r_ok:
        stdincfg = sys.stdin.read()
        # had a default decode here for bytes data, not sure why, just make it a check
        if isinstance(stdincfg, bytes):
            stdincfg = stdincfg.decode("utf-8")
    return stdincfg


def build_main_parser(stdincfg):
    """ Create the parent parser and check if
    there is a config file to be read
    :param stdincfg:
    :return:
    """
    conf_parser = argp.ArgumentParser(prog='CREEPIEST', add_help=False, allow_abbrev=False)
    fbconf = conf_parser.add_argument_group('File-based configuration')
    fbconf.add_argument('--config-file', '-cfg', dest='configfile', metavar='FILEPATH',
                        help='Specify full path to configuration file')
    if stdincfg:
        args, remain_args = conf_parser.parse_known_args(stdincfg)
    else:
        args, remain_args = conf_parser.parse_known_args()
    main_parser = argp.ArgumentParser(add_help=True, parents=[conf_parser], description=prog_desc,
                                      formatter_class=argp.RawDescriptionHelpFormatter, allow_abbrev=False)
    generics = main_parser.add_argument_group('CREEPIEST runtime parameter')
    generics.add_argument('--version', '-v', action='version', version=__version__,
                          help='Print version information and exit')
    generics.add_argument('--verbose', '-vb', action='store_true', default=False, dest='verbose',
                          help='Print status messages to stderr. Equivalent to --debug.'
                               ' Default: FALSE')
    generics.add_argument('--debug', '-dbg', action='store_true', default=False, dest='debug',
                          help='Print status messages to stderr. Equivalent to --verbose.'
                               ' Default: FALSE')
    generics.add_argument('--outfile-mode', '-ofm', type=str, choices=['a', 'w', 'replace', 'append'], default='w',
                          dest='filemode', help='Specify if output overwrites (replace [w]) or is'
                                                ' appended (append [a]) to existing files. Default: replace')
    generics.add_argument('--notify', '-ntf', type=str, default='', dest='notify', metavar='EMAIL',
                          help='Send notification upon run completion to this email address.'
                               ' This assumes that localhost is configured as SMTP server'
                               ' or knows where to find an SMTP server. Default: <empty>')
    generics.add_argument('--send-log', '-snd', action='store_true', dest='sendlog', default=False,
                          help='When sending the email notification, include the complete log.'
                               ' Warning: this can be a lot of text! Default: FALSE')
    generics.add_argument('--log-format', '-lfm', type=str, default='long',
                          dest='log_format', choices=['long', 'short'],
                          help='Specify log format string as long or short (level of detail). Default format: long')
    generics.add_argument('--profile', '-prf', action='store_true', dest='profiling', default=False,
                          help='Profile this run of the CREEPIEST tool. This only works if the C'
                               ' profiling libraries are available to limit the overhead. Creates a'
                               ' cProfile output file in the working directory. Default: FALSE')
    generics.add_argument('--no-dump', '-nod', action='store_true', default=False, dest='nodump',
                          help='If set to TRUE, do not dump a summary of the current configuration'
                               ' and run parameters. Default: FALSE')
    generics.add_argument('--dump-dir', '-dmp', type=str, default=os.getcwd(), dest='dumpdir',
                          help='Specify folder to dump configuration summary. Default: working dir')
    generics.add_argument('--workers', '-wrk', dest='workers', default=1, type=int, metavar='NUM',
                          help='Number of worker processes allowed to start. Note that for I/O intensive'
                               ' task, e.g. data conversion, this number should only be set to a large'
                               ' value if the storage hardware can handle parallel writes. Default: 1')
    main_parser.set_defaults(execute=lambda arg: main_parser.print_help())
    if args.configfile:
        cfgp = confp.ConfigParser()
        cfgp.read([args.configfile])
        for s in cfgp.sections():
            main_parser.set_defaults(**dict(cfgp.items(s)))
        main_parser.set_defaults(**cfgp.defaults())
    return main_parser, remain_args


def dump_config(args, ignored_args):
    """
    :param args:
     :type: Namespace
    :param ignored_args:
     :type: list
    :return: path to dumped config file
     :rtype: string
    """
    current_env = os.environ.copy()
    dump = confp.ConfigParser(allow_no_value=True)
    dump.optionxform = lambda opt: opt
    envk = 'ENVIRONMENT'
    dump.add_section(envk)
    for key in ['USER', 'PATH', 'PYTHONPATH', 'PWD', 'SHELL', 'HOSTNAME']:
        try:
            dump[envk][key] = current_env[key]
        except KeyError:
            dump[envk][key] = 'N/A'
    dtobj = dt.datetime.now()
    dump[envk]['TIMESTAMP'] = dtobj.strftime('%Y-%m-%d - %H:%M:%S')
    dump[envk]['CREEPIEST_VERSION'] = __version__

    cfgk = 'CONFIGURATION'
    dump.add_section(cfgk)
    params = vars(args)
    for key in sorted(params.keys()):
        dump[cfgk][key] = str(params[key])
    dump[cfgk]['ignored_args'] = ' '.join(ignored_args)

    fname = get_full_filestamp() + '_' + args.subparser_name + '_crp.ini'
    fpath = os.path.join(args.dumpdir, fname)
    try:
        with open(fpath, 'w') as outf:
            dump.write(outf)
    except FileNotFoundError:
        os.makedirs(args.dumpdir, exist_ok=True)
        with open(fpath, 'w') as outf:
            dump.write(outf)
    return fpath


def init_logging_system(args, logbuf):
    """
    :param args: command line arguments
     :type: Namespace object
    :param logbuf: buffer for log messages
     :type: io.StringIO or None
    :return: root logger
     :rtype: logger object
    """
    logger = logging.getLogger(name=None)
    log_formats = {'long': log_long, 'short': log_short}
    formatter = logging.Formatter(fmt=log_formats[args.log_format])
    shdl = logging.StreamHandler(stream=sys.stderr)
    shdl.setFormatter(formatter)
    logger.addHandler(shdl)
    if logbuf:
        bufhdl = logging.StreamHandler(stream=logbuf)
        bufhdl.setFormatter(formatter)
        logger.addHandler(bufhdl)
    if args.verbose or args.debug:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.WARNING
    logger.setLevel(loglevel)
    return logger


def normalize_file_mode(selected):
    """
    :param selected:
    :return:
    """
    if selected in ['a', 'w']:
        pass
    elif selected == 'replace':
        selected = 'w'
    elif selected == 'append':
        selected = 'a'
    else:
        raise AssertionError('This should be impossible - but:\n'
                             'In the face of ambiguity, refuse the temptation to guess.')
    return selected


def run():
    """
    :return: exit code
     :rtype: int
    """
    logger = None
    logbuf = None
    args = None
    # this uses OS randomness sources by default if available
    # and is currently not used for critical analysis, no need to
    # manually set seed for reproducibility
    rand.seed()
    retcode = 0
    conf_dump = ''
    dtobj = dt.datetime.now()
    starttime = dtobj.strftime('%Y%m%d - %H%M%S')
    try:
        recv_stdin = check_stdin_config()
        mainprs, remain_args = build_main_parser(recv_stdin)
        mainprs = add_sub_parsers(mainprs)
        args, remain_args = mainprs.parse_known_args(remain_args)
        if args.sendlog:
            logbuf = io.StringIO()
        logger = init_logging_system(args, logbuf)
        logger.debug('Logging system initialized')
        if remain_args:
            raise ValueError('Unknown parameters detected: {}'.format(remain_args))
        if not args.nodump and args.subparser_name not in ['tests', 'info']:
            conf_dump = dump_config(args, remain_args)
        if conf_dump:
            logger.debug('Configuration dumped to: {}'.format(conf_dump))
        logger.debug('Executing command: {}'.format(args.subparser_name))
        setattr(args, 'module_logger', logging.getLogger(args.subparser_name))
        setattr(args, 'filemode', normalize_file_mode(args.filemode))
        if args.profiling:
            retcode = cprf.runctx('args.execute(args)', {}, {'args': args}, 'crp_' + args.subparser_name + '.prf')
        else:
            retcode = args.execute(args)
    except SystemExit as se:
        retcode = se.args[0]
    except Exception as e:
        trbbuf = io.StringIO()
        trb.print_exc(file=trbbuf)
        if logger is not None:
            logger.error('Error: {}'.format(e))
            logger.error('Exception: {}'.format(trbbuf.getvalue()))
        else:
            # this is just a fallback in case the logger init already fails
            sys.stderr.write('\nError: {}\nTraceback: {}\n'.format(e, trbbuf.getvalue()))
        if os.path.isfile(conf_dump):
            try:
                os.remove(conf_dump)
            except (OSError, IOError, FileNotFoundError):
                pass
        retcode = e.args[0]
    finally:
        dtobj = dt.datetime.now()
        endtime = dtobj.strftime('%Y%m%d - %H%M%S')
        if args and args.notify:
            send_notification(args.notify, args.subparser_name, retcode, starttime, endtime, logbuf)
        logging.shutdown()
        return retcode


if __name__ == '__main__':
    exc = run()
    sys.exit(exc)
