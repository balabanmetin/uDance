import json
import sys
import multiprocessing as mp
from collections import deque

from compute_bipartition_alignment import compute_bipartition_alignment
from fasta2dic import fasta2dic
from newick_extended import read_tree_newick
import treeswift as ts
from pathlib import Path
from os import listdir
from os.path import isfile, join, splitext, expanduser, abspath
import time
import numpy as np


# inputs: a placement tree
# max number of things in each cluster
from optparse import OptionParser
from treecluster_sum import min_tree_coloring_sum


def aggregate_placements(index_to_node_map, placements):
    for placement in placements:
        for i, seqname in enumerate(placement['n']):
            index = placement["p"][i][0]
            index_to_node_map[index].placements += [seqname]


def set_closest_three_directions(tree):
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.closest_left = (0, node)
            node.closest_right = (0, node)
        else:
            left, right = node.children
            lb, ln = left.closest_child
            node.closest_left = (lb + left.edge_length, ln)
            rb, rn = right.closest_child
            node.closest_right = (rb + right.edge_length, rn)
        if node.closest_left[0] < node.closest_right[0]:
            node.closest_child = node.closest_left
        else:
            node.closest_child = node.closest_right

    for node in tree.traverse_preorder():
        if node == tree.root:
            node.closest_top = (float("inf"), None)
        else:
            sibling = [c for c in node.parent.children if c != node][0]
            sb, sn = sibling.closest_child
            pb, pn = node.parent.closest_top
            if pb < sibling.edge_length + sb:
                node.closest_top = (pb + node.edge_length, pn)
            else:
                node.closest_top = (sb + sibling.edge_length + node.edge_length, sn)


def undo_resolve_polytomies(tree):
    for e in tree.traverse_postorder():
        if e == tree.root:
            continue
        if e.outgroup == False and e.resolved_randomly:  # edge length check just for safety
            par = e.parent
            par.remove_child(e)
            for c in e.children:
                par.add_child(c)


def partition_worker(i, j):
    partition_output_dir = join(options.output_fp, str(i))
    Path(partition_output_dir).mkdir(parents=True, exist_ok=True)
    undo_resolve_polytomies(j)
    newick_path = join(partition_output_dir, "astral_constraint.nwk")
    j.write_tree_newick(newick_path)

    # find all outgroups
    outgroups = [n.label for n in j.traverse_postorder() if n.outgroup]
    if len(outgroups) >= 4:
        constraint = j.extract_tree_with(outgroups, suppress_unifurcations=True)
        constraint.is_rooted = False
        bipartition_path = join(partition_output_dir, "bipartition.fasta")
        with open(bipartition_path, "w") as f:
            f.write(compute_bipartition_alignment(constraint.__str__()))
        raxml_constraint_path = join(partition_output_dir, "raxml_constraint.nwk")
        constraint.write_tree_newick(raxml_constraint_path)

    species_list_path = join(partition_output_dir, "species.txt")
    species_list = []
    with open(species_list_path, "w") as f:
        for n in j.traverse_postorder():
            if n.is_leaf():
                species_list += [n.label]
                f.write(n.label + "\n")
            if hasattr(n, 'placements'):
                for p in n.placements:
                    species_list += [p]
                    f.write(p + "\n")
    return (i, species_list)


#  VECTORIZED (yay!)
#  TODO raise error if directory exists
def alignment_worker(i, j):
    species = species_dict[i]
    partition_aln = {key: fa_dict[key] for key in species if key in fa_dict}

    if len(partition_aln) < 4:
        return None
    aln_length = len(next(iter(partition_aln.values())))
    not_all_gap = np.array([False]*aln_length)
    for s in partition_aln.values():
        not_all_gap = np.logical_or(not_all_gap,  (s != b'-'))
    for k, v in partition_aln.items():
        partition_aln[k] = v[not_all_gap]

    trimmed_aln_length = len(next(iter(partition_aln.values())))
    if trimmed_aln_length >= options.overlap_length and len(partition_aln) >= 4:
        # write trimmed MSA fasta
        partition_output_dir = join(options.output_fp, str(i))
        aln_outdir = join(partition_output_dir, basename)
        Path(aln_outdir).mkdir(parents=True, exist_ok=True)
        aln_output_path = join(aln_outdir, "aln.fa")
        with open(aln_output_path, "w", buffering=100000000) as f:
            for k, v in partition_aln.items():
                f.write(">" + k + "\n")
                f.write(v.tostring().decode("UTF-8") + "\n")

        constraint_outgroup_tree = join(partition_output_dir, "raxml_constraint.nwk")
        if isfile(constraint_outgroup_tree):
            t = ts.read_tree_newick(constraint_outgroup_tree)
            induced_constraints_tree = t.extract_tree_with(list(partition_aln.keys()), suppress_unifurcations=True)
            induced_constraints_tree.is_rooted = False
            numlabels = induced_constraints_tree.num_nodes(internal=False)
            if numlabels >= 4:
                # write fasttree and raxml constraint
                bipartition_path = join(aln_outdir, "bipartition.fasta")
                with open(bipartition_path, "w") as f:
                    f.write(compute_bipartition_alignment(induced_constraints_tree.__str__()))
                induced_raxml_constraint_path = join(aln_outdir, "raxml_constraint.nwk")
                induced_constraints_tree.write_tree_newick(induced_raxml_constraint_path)

        script = join(aln_outdir, "run.sh")
        with open(script, "w") as f:
            f.write("#! /usr/bin/env bash\n\n")
            f.write("export OMP_NUM_THREADS=1\n\n")
            bipartition_path = join(aln_outdir, "bipartition.fasta")
            fasttree_log = join(aln_outdir, "fasttree.log")
            fasttree_err = join(aln_outdir, "fasttree.err")
            fasttree_nwk = join(aln_outdir, "fasttree.nwk")
            if isfile(bipartition_path):
                f.write("FastTree -constraints %s -log %s < %s > %s 2> %s \n"
                        % (bipartition_path, fasttree_log, aln_output_path, fasttree_nwk, fasttree_err))
            else:
                f.write("FastTree -log %s < %s > %s 2> %s \n"
                        % (fasttree_log, aln_output_path, fasttree_nwk, fasttree_err))

            fasttree_resolved_nwk = join(aln_outdir, "fasttree_resolved.nwk")
            f.write("python3 -c \"import sys, treeswift; "
                    "t=treeswift.read_tree_newick(input()); "
                    "t.resolve_polytomies(); print(t)\" < %s > %s \n" % (fasttree_nwk, fasttree_resolved_nwk))

            induced_raxml_constraint_path = join(aln_outdir, "raxml_constraint.nwk")
            raxml_err = join(aln_outdir, "raxml.err")
            raxml_out = join(aln_outdir, "raxml.out")
            if isfile(induced_raxml_constraint_path):
                f.write("raxml-ng --tree %s --tree-constraint %s "
                        "--msa %s --model LG+G --prefix RUN --seed 12345 "
                        "--threads 1 > %s 2> %s \n"
                        % (fasttree_resolved_nwk, induced_raxml_constraint_path, aln_output_path, raxml_out, raxml_err))
            else:
                f.write("raxml-ng --tree %s "
                        "--msa %s --model LG+G --prefix RUN --seed 12345 "
                        "--threads 1 > %s 2> %s \n"
                        % (fasttree_resolved_nwk, aln_output_path, raxml_out, raxml_err))
            return script
    return None


if __name__ == "__main__":
    mp.set_start_method('fork')
    parser = OptionParser()
    parser.add_option("-j", "--jplace", dest="jplace_fp",
                      help="path to the jplace placement file", metavar="FILE")
    parser.add_option("-f", "--threshold", dest="threshold", default="600",
                      help="maximum number of elements in each cluster")
    parser.add_option("-o", "--output", dest="output_fp",
                      help="path for the output directory where files will be placed",
                      metavar="DIRECTORY")
    parser.add_option("-s", "--alignment-dir", dest="alignment_dir_fp",
                      help="path for input directory which contains "
                           "extended reference alignment files (FASTA), "
                           "containing reference and query sequences.",
                      metavar="FILE")
    parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true', default=False,
                      help="input sequences are protein sequences")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in placement. "
                           "0 to use all cores in the running machine", metavar="NUMBER")
    parser.add_option("-l", "--overlap", dest="overlap_length", default="50",
                      help="minimum alignment overlap length needed to use the subalignment"
                           "in subtree refinement", metavar="NUMBER")

    (options, args) = parser.parse_args()

    options.num_thread = int(options.num_thread)
    options.overlap_length = int(options.overlap_length)
    if not options.num_thread:
        options.num_thread = mp.cpu_count()

    options.output_fp = abspath(expanduser(options.output_fp))

    with open(options.jplace_fp) as f:
        jp = json.load(f)
    tstree = read_tree_newick(jp["tree"])

    index_to_node_map = {}
    for e in tstree.traverse_postorder():
        e.placements = []
        if e != tstree.root:
            index_to_node_map[e.edge_index] = e
    aggregate_placements(index_to_node_map, jp["placements"])

    min_tree_coloring_sum(tstree, float(options.threshold))

    set_closest_three_directions(tstree)

    colors = {}
    for n in tstree.traverse_postorder():
        if n.color not in colors:
            colors[n.color] = [n]
        else:
            colors[n.color] += [n]

    color_spanning_tree = ts.Tree()
    r = color_spanning_tree.root
    r.color = tstree.root.color
    r.label = str(r.color)

    color_to_node_map = {r.color: r}
    for n in tstree.traverse_preorder():
        if n.is_leaf():
            continue
        for c in n.children:
            if c.color != n.color:
                par = color_to_node_map[n.color]
                child = ts.Node()
                par.add_child(child)
                color_to_node_map[c.color] = child
                child.joint = n
                child.color = c.color
                child.label = str(child.color)

    copyts = tstree.extract_tree_with(labels=tstree.labels())
    for n in copyts.traverse_preorder():
        n.outgroup = False

    traversal = list(zip(tstree.traverse_preorder(), copyts.traverse_preorder()))
    tree_catalog = {}

    for n, ncopy in traversal:
        ncopy.resolved_randomly = n.resolved_randomly
        ncopy.placements = n.placements
        if n.is_leaf():
            continue
        cl, cr = n.children
        clcopy, crcopy = ncopy.children

        if cl.color != n.color:
            # replace it with representative
            ncopy.remove_child(clcopy)
            clcopy_repr = ts.Node()
            clcopy_repr.outgroup = True
            lb, ln = n.closest_left
            clcopy_repr.label = ln.label
            clcopy_repr.edge_length = lb
            ncopy.add_child(clcopy_repr)

            # make child a new tree
            newTree = ts.Tree()
            newTree.is_rooted = False
            newTree.root.outgroup = False
            newTree.root.add_child(clcopy)
            pb, pn = n.closest_top
            if pn:
                top_repr = ts.Node()
                top_repr.outgroup = True
                top_repr.label = pn.label
                top_repr.edge_length = pb
                newTree.root.add_child(top_repr)

            crcopy_repr = ts.Node()
            crcopy_repr.outgroup = True
            rb, rn = n.closest_right
            crcopy_repr.label = rn.label
            crcopy_repr.edge_length = rb

            if cl.color == cr.color:
                ncopy.remove_child(crcopy)
                ncopy.add_child(crcopy_repr)
                newTree.root.add_child(crcopy)

            else: # cr.color != n.color and cr.color != cl.color:
                newTree.root.add_child(crcopy_repr)

            tree_catalog[cl.color] = newTree

        if cr.color != n.color and cr.color != cl.color:
            ncopy.remove_child(crcopy)
            crcopy_repr = ts.Node()
            crcopy_repr.outgroup = True
            rb, rn = n.closest_right
            crcopy_repr.label = rn.label
            crcopy_repr.edge_length = rb
            ncopy.add_child(crcopy_repr)

            # make child a new tree
            newTree = ts.Tree()
            newTree.is_rooted = False
            newTree.root.outgroup = False
            newTree.root.add_child(crcopy)
            pb, pn = n.closest_top
            if pn:
                top_repr = ts.Node()
                top_repr.outgroup = True
                top_repr.label = pn.label
                top_repr.edge_length = pb
                newTree.root.add_child(top_repr)

            clcopy_repr = ts.Node()
            clcopy_repr.outgroup = True
            lb, ln = n.closest_left
            clcopy_repr.label = ln.label
            clcopy_repr.edge_length = lb
            newTree.root.add_child(clcopy_repr)
            tree_catalog[cr.color] = newTree



    # stitching algorithm:
    # preorder traversal color_to_node_map
    # for each node n, find the joint j in tstree.
    # find LCA of j.left_closest_child and j.right_closest_child in n
    # replace it with n.children

    Path(options.output_fp).mkdir(parents=True, exist_ok=True)

    pool = mp.Pool(options.num_thread)
    species_dict = dict(pool.starmap(partition_worker, tree_catalog.items()))
    pool.close()
    pool.join()

    only_files = [f for f in listdir(options.alignment_dir_fp) if isfile(join(options.alignment_dir_fp, f))]
    main_script = open(join(options.output_fp, "main_script.sh"), "w")
    main_script.write("#! /usr/bin/env bash\n\n")

    for aln in only_files:
        aln_input_file = join(options.alignment_dir_fp, aln)
        basename = splitext(aln)[0]
        try:
            fa_dict = fasta2dic(aln_input_file, options.protein_seqs, False)
            pool = mp.Pool(options.num_thread)
            scripts = pool.starmap(alignment_worker, tree_catalog.items())
            pool.close()
            pool.join()
            valid_scripts = [s for s in scripts if s is not None]
            main_script.write("\n".join(valid_scripts))
            main_script.write("\n")
        except e:
            print("Alignment %s is not a valid fasta alignment" % aln_input_file, file=sys.stderr)

    main_script.close()

    # TODO a bipartition for each alignment
    for i, j in tree_catalog.items():
        count = 0
        for n in j.traverse_postorder():
            try:
                if n.is_leaf():
                    count +=1
                elif n.outgroup is False and hasattr(n, 'placements'):
                    count += len(n.placements)
            except e:
                pass
        print (i, count)



    print("hello")