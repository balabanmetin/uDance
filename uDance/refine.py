import pkg_resources

from uDance.PoolAstralWorker import PoolAstralWorker


def refine(options):

    astral_exec = pkg_resources.resource_filename('uDance', "tools/ASTRAL/astral.5.17.2.jar")
    astral_mp_exec = pkg_resources.resource_filename('uDance', "tools/ASTRAL/astralmp.5.17.2.jar")
    partition_worker = PoolAstralWorker()
    partition_worker.set_class_attributes(options, astral_exec, astral_mp_exec)
    partition_worker.worker(options.partition_dir)




    #insert back removed duplicates
    #besttree = join(options.output_fp,"jobsizes.txt")
