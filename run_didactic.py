import json
import sys
import multiprocessing as mp

from didactic.PoolAlignmentWorker import PoolAlignmentWorker
from didactic.fasta2dic import fasta2dic
from didactic.newick_extended import read_tree_newick
import treeswift as ts
from pathlib import Path
from os import listdir
from os.path import isfile, join, splitext
from didactic.options import options_config
from didactic.treecluster_sum import min_tree_coloring_sum
from didactic.PoolPartitionWorker import PoolPartitionWorker

# inputs: a placement tree
# max number of things in each cluster


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


if __name__ == "__main__":
    mp.set_start_method('fork')

    options, args = options_config()

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
    partition_worker = PoolPartitionWorker()
    partition_worker.set_class_attributes(options)

    pool = mp.Pool(options.num_thread)
    species_dict = dict(pool.starmap(partition_worker.worker, tree_catalog.items()))
    pool.close()
    pool.join()

    only_files = [f for f in listdir(options.alignment_dir_fp) if isfile(join(options.alignment_dir_fp, f))]
    main_script = open(join(options.output_fp, "main_script.sh"), "w")

    for aln in only_files:
        aln_input_file = join(options.alignment_dir_fp, aln)
        basename = splitext(aln)[0]
        try:
            fa_dict = fasta2dic(aln_input_file, options.protein_seqs, False)
            alignment_worker = PoolAlignmentWorker()
            alignment_worker.set_class_attributes(options, species_dict, fa_dict, basename)
            pool = mp.Pool(options.num_thread)
            scripts = pool.starmap(alignment_worker.worker, tree_catalog.items())
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