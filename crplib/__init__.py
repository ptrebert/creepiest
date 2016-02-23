# coding=utf-8

__version__ = '0.4.0'
__default_log_long__ = '%(asctime)s [%(levelname)s-%(name)s]: %(module)s::%(funcName)s | %(message)s | %(lineno)d'
__default_log_short__ = '%(asctime)s: %(module)s::%(funcName)s | %(message)s'
__author__ = 'Peter Ebert'
__maintainer__ = 'Peter Ebert'
__author_email__ = 'pebert@mpi-inf.mpg.de'

__description__ = """        The CREEPIEST tool
[CRoss-spEcies EPIgenome ESTimation]
------------------------------------

The CREEPIEST tool report errors always on stderr.
If the verbose switch is set, status messages will print to stderr.

For basic usage information, run:
./creepiest <subcommand> --help
"""