import pkg_resources

from uDance.PoolAstralWorker import PoolAstralWorker


def refine(options):

    astral_libdir = pkg_resources.resource_filename('uDance', 'tools/ASTRAL/lib/')
    astral_mp_exec = pkg_resources.resource_filename('uDance', 'tools/ASTRAL/astralmp.5.17.2.jar')
    partition_worker = PoolAstralWorker()
    partition_worker.set_class_attributes(options, astral_mp_exec, astral_libdir)
    partition_worker.worker(options.partition_dir)

    # insert back removed duplicates
    # besttree = join(options.output_fp,"jobsizes.txt")
