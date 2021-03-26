import json
from os.path import join
import treeswift as ts


def stitch(options):

    outmap_file = join(options.output_fp, "outgroup_map.json")
    with open(outmap_file) as o:
        outmap = json.load(o)

    cg_file = join(options.output_fp, "color_spanning_tree.nwk")
    cg = ts.read_tree_newick(cg_file)

    removed = set()
    def _stitch(node):
        if node.label == "-1":
            mytree = ts.Tree()
            mytree.is_rooted = True
            for c in node.children:
                mytree.root.add_child(_stitch(c).root)
            return mytree

        astral_tree_par = ts.read_tree_newick(join(options.output_fp, node.label, "astral_output.nwk"))
        astral_tree_cons = ts.read_tree_newick(join(options.output_fp, node.label, "astral_constraint.nwk"))
        outmap_par = outmap[node.label]

        uptree = ts.read_tree_newick(outmap_par["up"])
        uptree_labels = set(uptree.labels(internal=False))
        astral_tree_cons_labels = set(astral_tree_cons.labels(internal=False))

        for i in astral_tree_par.traverse_postorder(internal=False):
            if i.label not in uptree_labels and i.label in astral_tree_cons_labels:
                notuptree_species = i   # this has to be a backbone species
                break
        astral_tree_par.is_rooted = True
        astral_tree_par.reroot(notuptree_species.parent)
        astral_tree_par.suppress_unifurcations()
        mrca = astral_tree_par.mrca(list(uptree_labels))
        astral_tree_par.reroot(mrca)
        astral_tree_par.suppress_unifurcations()
        if outmap_par['ownsup']:
            deletelist = []
            for c in astral_tree_par.root.children:
                for llab in c.traverse_postorder(internal=False):
                    if llab.label in uptree_labels:
                        deletelist += [c]
                        break

            for i in deletelist:
                for j in i.traverse_postorder(internal=False):
                    if j.label not in uptree_labels:
                        removed.add(j.label + "\t" + node.label)
                astral_tree_par.root.remove_child(i)
            if len(astral_tree_par.root.children) != 1:
                raise ValueError('Astral tree is not binary.')
            astral_tree_par.root = astral_tree_par.root.children[0]  # get rid of the degree 2 node
        else:
            non_uptree = astral_tree_cons_labels.difference(uptree_labels)
            non_uptree_mrca = astral_tree_par.mrca(list(non_uptree))
            to_be_deleted = astral_tree_par
            non_uptree_mrca.parent.remove_child(non_uptree_mrca)
            for j in to_be_deleted.traverse_postorder(internal=False):
                if j.label not in uptree_labels:
                    removed.add(j.label + "\t" + node.label)
            astral_tree_par = ts.Tree()
            astral_tree_par.is_rooted = True
            astral_tree_par.root = non_uptree_mrca

        for c in node.children:
            ctree = _stitch(c)
            c_rep_tree = ts.read_tree_newick(outmap_par["children"][c.label])
            c_rep_tree_labels = set(c_rep_tree.labels(internal=False))
            c_rep_tree_mrca = astral_tree_par.mrca(list(c_rep_tree_labels))
            for j in c_rep_tree_mrca.traverse_postorder(internal=False):
                if j.label not in c_rep_tree_labels:
                    removed.add(j.label + "\t" + node.label)
            c_rep_tree_mrca_parent = c_rep_tree_mrca.parent
            c_rep_tree_mrca_parent.remove_child(c_rep_tree_mrca)
            c_rep_tree_mrca_parent.add_child(ctree.root)

        return astral_tree_par
    stitched_tree = _stitch(cg.root)
    final_tree = join(options.output_fp, "didactic.nwk")
    stitched_tree.write_tree_newick(final_tree)
    unplaced = join(options.output_fp, "unplaced.csv")
    with open(unplaced, "w") as f:
        f.write("\n".join(removed) + "\n")

    return
    # for n in cg.traverse_postorder():
    #     print(outmap[n.label])
