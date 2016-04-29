# coding=utf-8

DIV_B_TO_GB = 1024 * 1024 * 1024

DIV_B_TO_MB = 1024 * 1024

# this constant is used to check if an object
# is too large to be serialized and sent directly
# to the calling process/parent (from the worker)
# see "bug" here: https://bugs.python.org/issue17560
# Seems to exist mainly for backward compatibility!?
LIMIT_SERIALIZATION = 2**31 - 1  # NB to myself: is max. size in bytes, i.e. check obj.nbytes - not obj.size!

# ignore this many bases at the beginning of
# each chromosome; this value is also used in, e.g., ChromImpute
CHROMOSOME_BOUNDARY = 10000

TRGIDX_MASK = 'cons/mask'
TRGIDX_SPLITS = 'cons/splits'
TRGIDX_SELECT = 'cons/select'

FEAT_FP_PREC = 5
