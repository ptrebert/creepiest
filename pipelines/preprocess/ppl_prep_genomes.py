# coding=utf-8

from ruffus import *

def build_pipeline(args, config, sci_obj):


    pipe = Pipeline(name='PREP_GENOMES')

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()


    return pipe