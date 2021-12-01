import json
from os.path import join

import pkg_resources

from uDance.PoolAstralWorker import PoolAstralWorker
import multiprocessing as mp


def refine(options):

    astral_exec = pkg_resources.resource_filename('uDance', "tools/astral.jar")
    partition_worker = PoolAstralWorker()
    partition_worker.set_class_attributes(options, astral_exec)
    partition_worker.worker(options.cluster_dir)




    #insert back removed duplicates
    #besttree = join(options.output_fp,"jobsizes.txt")