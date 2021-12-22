
from uDance.fasta2dic import fasta2dic
from os import listdir
from os.path import isfile, join, splitext


def count_occupancy(alndir, protein):
    only_files = [f for f in listdir(alndir) if isfile(join(alndir, f)) and not f.startswith(".")]

    num_genes = 0
    occupancy = {}
    for aln in only_files:
        aln_input_file = join(alndir, aln)
        basename = splitext(aln)[0]
        fa_dict = fasta2dic(aln_input_file, protein, False)
        num_genes += 1
        for k in fa_dict.keys():
            if k in occupancy:
                occupancy[k] += 1
            else:
                occupancy[k] = 1
    return occupancy, num_genes
