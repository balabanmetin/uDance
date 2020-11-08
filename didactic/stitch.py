import json
from os.path import join

import pkg_resources

from didactic.PoolAstralWorker import PoolAstralWorker
import multiprocessing as mp


def stitch(options):

    outmap = join(options.output_fp, "outgroup_map.json")
    with open(outmap) as o:
        j = json.load(o)
    partition_dirs = [x for x in j.keys() if int(x) >= 0]
    astral_exec = pkg_resources.resource_filename('didactic', "tools/astral.jar")
    partition_worker = PoolAstralWorker()
    partition_worker.set_class_attributes(options, astral_exec)
    pool = mp.Pool(options.num_thread)
    results = pool.map(partition_worker.worker, partition_dirs)
    pool.close()
    pool.join()

    #insert back removed duplicates
    #besttree = join(options.output_fp,"jobsizes.txt")