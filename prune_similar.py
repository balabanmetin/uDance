import json
from multiprocessing import cpu_count
from optparse import OptionParser
from os.path import join

from uDance.subsample_partition import subsample_partition


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-o", "--output", dest="output_fp",
                      help="path for the output directory where partitions are located",
                      metavar="DIRECTORY")
    parser.add_option("-T", "--threads", type=int, dest="num_thread", default=0,
                      help="number of cores used in placement. "
                           "0 to use all cores in the running machine", metavar="NUMBER")
    parser.add_option("-S", "--size", type=int, dest="minimum_size", metavar='NUMBER', default=9000,
                      help="partition size requirement for pruning.")

    (options, args) = parser.parse_args()
    if options.num_thread == 0:
        options.num_thread = cpu_count()

    outmap = join(options.output_fp, "outgroup_map.json")
    with open(outmap) as o:
        j = json.load(o)
    partition_dirs = [x for x in j.keys() if int(x) >= 0]
    dupmapstrs = []
    for i in partition_dirs:
        partition_output_dir = join(options.output_fp, str(i))
        with open(join(partition_output_dir, "species.txt")) as f:
            species = list(map(lambda x: x.strip(), f.readlines()))
            numspecies = len(species)
            if numspecies < options.minimum_size:
                continue
        print(numspecies)
        res = subsample_partition(partition_output_dir, options.minimum_size)
        if res:
            dupmapstrs.append(res)

    if len(dupmapstrs) > 0:
        with open(join(options.output_fp, "dedupe_map.txt"), "a") as f:
            for st in dupmapstrs:
                f.write(st)



