from os.path import join
from glob import glob
import numpy as np
from pathlib import Path
from uDance.fasta2dic import readfq
from scipy.sparse.csgraph import connected_components
import time
import treeswift as ts


def subsample_partition(partition_output_dir, limit):
    with open(join(partition_output_dir, 'species.txt')) as f:
        species = set(map(lambda x: x.strip(), f.readlines()))
    t = ts.read_tree_newick(join(partition_output_dir, 'astral_constraint.nwk'))
    tlabs = t.labels(internal=False)
    for i in tlabs:
        if i in species:
            species.remove(i)

    numspecies = len(species)
    ind_to_name = dict(enumerate(species))
    name_to_ind = {v: k for k, v in ind_to_name.items()}
    counts_ij = np.zeros((numspecies, numspecies), dtype=np.int16)  # counts in how many genes i and j are are identical
    counts_i = np.zeros((numspecies,), dtype=np.int16)
    genes = glob(join(partition_output_dir, '*', ''))
    print('number of genes %d. ' % len(genes))

    start = time.time()
    for g in genes:
        with open(join(g, 'aln.fa')) as af:
            glabels = [name_to_ind[name] for name, seq, _ in readfq(af) if name in name_to_ind]
            for i in glabels:
                counts_i[i] += 1

        dupmap_path = Path(join(g, 'dupmap.txt'))
        if not dupmap_path.is_file():  # every sequence in the alignment is unique
            continue
        with open(dupmap_path) as f:
            for line in f.readlines():
                things = [name_to_ind[j] for j in line.strip().split('\t')[1:] if j in name_to_ind]
                for i in things:
                    counts_i[i] += 1

                # here, if and else does the same thing. The only purpose of the if..else.. is to improve performance
                if len(things) <= np.sqrt(numspecies):
                    for i in things:
                        for j in things:
                            if i == j:
                                continue
                            counts_ij[i][j] += 1
                else:
                    dupcounts_i = np.zeros((numspecies,), dtype=np.int16)
                    for i in things:
                        dupcounts_i[i] = 1
                    dupcounts_ij = np.logical_and(dupcounts_i[..., np.newaxis], dupcounts_i[np.newaxis, ...]).astype(
                        np.int16
                    )
                    # np.fill_diagonal(dupcounts_ij, 0)
                    counts_ij += dupcounts_ij

    print('counting neigh %.3f.' % (time.time() - start))
    # print(counts_ij[52][53], counts_i[52], counts_i[53])
    start = time.time()
    # print (counts_ij)
    # print((counts_ij/np.minimum(counts_i[...,np.newaxis],counts_i[np.newaxis,...])).max())
    x = counts_ij / np.minimum(counts_i[..., np.newaxis], counts_i[np.newaxis, ...])
    x.dump(join(partition_output_dir, 'adj_mat.pkl'), protocol=4)
    # print(x)

    #    clow = 1
    #    chigh = 100
    #
    #    while chigh > clow or chigh == 2:
    #        cutoff = (chigh + clow + 1)//2
    #        print(clow, cutoff, chigh)
    #        y = (x >= (cutoff/100))
    #        n, components = connected_components(y)
    #        if limit < n:
    #            chigh = cutoff-1
    #        elif limit > n:
    #            clow = cutoff
    #        else:
    #            clow = cutoff
    #            chigh = cutoff # making sure clow=chigh always holds at the exit. we will print that.
    #            break

    cutoff = 100
    curr = numspecies
    while cutoff >= 2:
        cutoff -= 1
        y = x >= (cutoff / 100)
        n, components = connected_components(y)
        curr = n
        if curr < limit:
            break

    print(
        'Partition '
        + partition_output_dir
        + ' is pruned to %d taxa at the automatic cutoff %.2f.' % (curr, cutoff / 100)
    )

    print('components %.3f.' % (time.time() - start))

    organized_components = dict()
    for i, comp_id in enumerate(components):
        if comp_id in organized_components:
            organized_components[comp_id].append(i)
        else:
            organized_components[comp_id] = [i]

    pruned_species = []
    dupmapstr = ''
    for v in organized_components.values():
        components_member_names = list(map(lambda x: ind_to_name[x], v))
        pruned_species += components_member_names[1:]
        if len(components_member_names) > 1:
            dupmapstr += str(len(components_member_names)) + '\t' + ', '.join(components_member_names) + '\n'
    if dupmapstr:
        with open(join(partition_output_dir, 'pruning_dupmap.txt'), 'w') as f:
            f.write(dupmapstr)

    if len(pruned_species) == 0:
        return ''
    pruned_species = set(pruned_species)

    # remove pruned from species.txt
    with open(join(partition_output_dir, 'species.txt')) as f:
        prevspecies = set(map(lambda x: x.strip(), f.readlines()))
    newspecieslst = list(prevspecies.difference(pruned_species))
    with open(join(partition_output_dir, 'species.txt'), 'w') as f:
        f.write('\n'.join(newspecieslst) + '\n')

    for g in genes:
        aln_dict = dict()
        with open(join(g, 'aln.fa')) as af:
            for name, seq, _ in readfq(af):
                aln_dict[name] = seq

        # expands dups
        dupmap_path = Path(join(g, 'dupmap.txt'))
        if dupmap_path.is_file():
            with open(dupmap_path) as f:
                for line in f.readlines():
                    things = line.strip().split('\t')
                    for i in things[1:]:
                        aln_dict[i] = aln_dict[things[0]]
        for i in pruned_species:
            if i in aln_dict:
                del aln_dict[i]

        # deduplicate the alignment
        seq_keyed_dict = {}
        for name, seq in aln_dict.items():
            if seq in seq_keyed_dict:
                seq_keyed_dict[seq].append(name)
            else:
                seq_keyed_dict[seq] = [name]

        if len(seq_keyed_dict) >= 4:
            # write trimmed MSA fasta
            res = []
            duplist = []
            for k, v in seq_keyed_dict.items():
                res.append('>' + v[0])
                res.append(k)
                if len(v) > 1:
                    duplist.append('\t'.join(v))

            aln_output_path = join(g, 'aln.fa')
            with open(aln_output_path, 'w', buffering=100000000) as f:
                f.write('\n'.join(res))
                f.write('\n')
            if duplist:
                dupmap_output_path = join(g, 'dupmap.txt')
                with open(dupmap_output_path, 'w', buffering=100000000) as f:
                    f.write('\n'.join(duplist))
                    f.write('\n')
    return dupmapstr
