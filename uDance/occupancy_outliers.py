
# $1 alignments directory
# $2 clusters file
# $3 prot or nuc

import sys
from uDance.count_occupancy import count_occupancy
from uDance.tc_parser import tc_parser
from kmeans1d import cluster

def occupancy_outliers(alignments_dir, clusters_file, protein):
    occupancy, num_genes = count_occupancy(alignments_dir, protein)
    clusters = tc_parser(clusters_file)
    deletedlist = []
    for n,clus in clusters:
        while True:
            occups = [occupancy[tag] for tag in clus]
            if len(occups) <= 1:
                continue
            cl, cent = cluster(occups, k=2)
            cardi_zero = len(cl) - sum(cl)
            cardi_one = sum(cl)
            if cardi_zero > 0 and cardi_zero < cardi_one and (cent[1]-cent[0])/cent[1] >= 0.5:
                minoccup = min(occups)
                deletes = [n for n,i in enumerate(occups) if i == minoccup]
                deletedlist += [i for n,i in enumerate(clus) if n in deletes]
                clus = [i for n,i in enumerate(clus) if n not in deletes]
            else:
                break
    if deletedlist:
        for i in deletedlist:
            print(i)
