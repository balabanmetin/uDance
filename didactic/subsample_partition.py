from os.path import join
from glob import glob
import numpy as np
from pathlib import Path
from didactic.fasta2dic import readfq
from scipy.sparse.csgraph import connected_components
import time


def subsample_partition(partition_output_dir, cutoff):
    with open(join(partition_output_dir, "species.txt")) as f:
        species = list(map(lambda x: x.strip(), f.readlines()))
        numspecies = len(species)
    ind_to_name = dict(enumerate(species))
    name_to_ind = {v: k for k, v in ind_to_name.items()}
    counts_ij = np.zeros((numspecies, numspecies), dtype=int)
    counts_i = np.zeros((numspecies, ), dtype=int)
    genes = glob(join(partition_output_dir, "*", ""))
    print("number of genes %d. " % len(genes))

    start = time.time()
    for g in genes:
        print(g)
        dupmap_path = Path(join(g,"dupmap.txt"))
        if not dupmap_path.is_file():
            continue
        with open(join(g, "aln.fa")) as af:
            glabels = [name_to_ind[name] for name, seq, _ in readfq(af)]
            for i in glabels:
                counts_i[i] += 1

        with open(dupmap_path) as f:
            for line in f.readlines():
                things = [name_to_ind[j] for j in line.strip().split("\t")]
                for i in things[1:]:
                    counts_i[i] += 1

                if len(things) <= np.sqrt(numspecies):
                    for i in things:
                        for j in things:
                            if i == j:
                                continue
                            counts_ij[i][j] += 1
                else:
                    dupcounts_i = np.zeros((numspecies,), dtype=int)
                    for i in things:
                        dupcounts_i[i] = 1
                    dupcounts_ij = np.logical_and(dupcounts_i[..., np.newaxis], dupcounts_i[np.newaxis, ...]).astype(int)
                    np.fill_diagonal(dupcounts_ij, 0)
                    counts_ij += dupcounts_ij



    print("counting neigh %.3f." % (time.time() - start))
    # print(counts_ij[52][53], counts_i[52], counts_i[53])
    start = time.time()
    #print (counts_i)
    x = np.minimum(counts_i[...,np.newaxis],counts_i[np.newaxis,...]) - counts_ij
    #print(x)
    y = (x <= cutoff)
    print("redo %.3f." % (time.time() - start))
    start = time.time()
    n, components = connected_components(y)
    print("components %.3f." % (time.time() - start))
    print(n)
    return