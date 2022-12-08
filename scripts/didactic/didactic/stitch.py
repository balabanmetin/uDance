import json
from os.path import join
import treeswift as ts

def deroot(tree):
    if len(tree.root.children) > 2:
        tree.is_rooted = False
    elif len(tree.root.children) == 2:
        left, right = tree.root.children
        if left.is_leaf():
            ln = left.edge_length + right.edge_length
            tree.root.remove_child(left)
            right.add_child(left)
            left.edge_length = ln
            tree.root = right
            right.edge_length = None
            tree.is_rooted = False
        else:
            ln = left.edge_length + right.edge_length
            tree.root.remove_child(right)
            left.add_child(right)
            right.edge_length = ln
            tree.root = left
            left.edge_length = None
            tree.is_rooted = False
    return



def safe_midpoint_reroot(tree, node):
    pendant_edge_length = node.edge_length
    node.edge_length = 1
    tree.root.edge_length = None # this prevents a leaf with label "ROOT" from appearing after reroot
    tree.reroot(node, 0.5)
    tree.suppress_unifurcations()
    assert len(tree.root.children) == 2
    for rc in tree.root.children:
        rc.edge_length = pendant_edge_length / 2


def stitch(options):
    outmap_file = join(options.output_fp, "outgroup_map.json")
    with open(outmap_file) as o:
        outmap = json.load(o)

    cg_file = join(options.output_fp, "color_spanning_tree.nwk")
    cg = ts.read_tree_newick(cg_file)

    if len(outmap["-1"]["children"]) == 1:
        rm_root = cg.root
        cg.reroot(cg.root.children[0])
        rm_root.parent.remove_child(rm_root)
        outmap.pop("-1")

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
        astral_tree_cons_labels = set(astral_tree_cons.labels(internal=False))

        outmap_par = outmap[node.label]

        uptreestr = outmap_par["up"]
        non_uptree = astral_tree_cons_labels # default
        if uptreestr: # there is an uptree
            uptree = ts.read_tree_newick(outmap_par["up"])

            uptree_labels = set(uptree.labels(internal=False))
            non_uptree = astral_tree_cons_labels.difference(uptree_labels) # override

            astral_tree_par.root.edge_length = None
            for i in astral_tree_par.traverse_postorder(internal=False):
                if i.label in non_uptree:
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
                non_uptree_mrca = astral_tree_par.mrca(list(non_uptree))
                to_be_deleted = astral_tree_par
                non_uptree_mrca.parent.remove_child(non_uptree_mrca)
                for j in to_be_deleted.traverse_postorder(internal=False):
                    if j.label not in uptree_labels:
                        removed.add(j.label + "\t" + node.label)
                astral_tree_par = ts.Tree()
                astral_tree_par.is_rooted = True
                astral_tree_par.root = non_uptree_mrca
        # if there's no uptree, find a backbone species and root at the middle of it's edge.
        # if there is no backbone species either, there must be two children subsets.
        # root at a representative of first child. find mrca of one of the other children. root at there.
        else:
            childreps = set()
            for chd in node.children:
                c_rep_tree = ts.read_tree_newick(outmap_par["children"][chd.label])
                childreps |= set(c_rep_tree.labels(internal=False))
            notreps = astral_tree_cons_labels.difference(childreps)

            constree_norep_species = None
            for i in astral_tree_par.traverse_postorder(internal=False):
                if i.label in astral_tree_cons_labels and i.label in notreps:
                    constree_norep_species = i   # this has to be a backbone species that is not a children repres.
                    break
            if constree_norep_species:
                safe_midpoint_reroot(astral_tree_par, constree_norep_species)
            else: # no backbone exists. only representatives.
                assert len(node.children) >= 2
                first_c = node.children[0]
                first_c_rep_tree = ts.read_tree_newick(outmap_par["children"][first_c.label])
                first_c_rep_tree_labels = set(first_c_rep_tree.labels(internal=False))
                for i in astral_tree_par.traverse_postorder(internal=False):
                    if i.label in first_c_rep_tree_labels:
                        first_c_rep_tree_rep = i  # this has to be a backbone species that is not a children repres.
                        break
                safe_midpoint_reroot(astral_tree_par, first_c_rep_tree_rep)
                second_c = node.children[1]
                second_c_rep_tree = ts.read_tree_newick(outmap_par["children"][second_c.label])
                second_c_rep_tree_labels = set(second_c_rep_tree.labels(internal=False))
                second_c_rep_tree_mrca = astral_tree_par.mrca(list(second_c_rep_tree_labels))
                safe_midpoint_reroot(astral_tree_par, second_c_rep_tree_mrca)

        for c in node.children:
            ownsup_child = outmap[c.label]["ownsup"]
            ctree = _stitch(c)
            c_rep_tree = ts.read_tree_newick(outmap_par["children"][c.label])
            c_rep_tree_labels = set(c_rep_tree.labels(internal=False))
            try:
                c_rep_tree_mrca = astral_tree_par.mrca(list(c_rep_tree_labels))
            except:
                breakpoint()
            c_rep_tree_mrca_parent = c_rep_tree_mrca.parent
            for j in c_rep_tree_mrca.traverse_postorder(internal=False):
                if j.label not in c_rep_tree_labels:
                    removed.add(j.label + "\t" + node.label)
            # if child does not own its "up" edge, only misplaced queries are under the mrca of child representatives.
            # in that case, the edge length of the mrca node comes from the current subtree, not from the child.
            if not ownsup_child:
                mrca_len = c_rep_tree_mrca.edge_length
                c_rep_tree_mrca_parent.remove_child(c_rep_tree_mrca)
                ctree.root.edge_length = mrca_len
                c_rep_tree_mrca_parent.add_child(ctree.root)
            #  if child owns its "up" edge, there are potentially misplaced queries above the mrca.
            #  starting from the mrca, we traverse on the path from mrca to the root until we find an internal node
            #  with least one descendant from backbone
            else:
                current = c_rep_tree_mrca
                parent = c_rep_tree_mrca_parent
                while parent != astral_tree_par.root:
                    has_non_uptree_species = False
                    for cp in parent.children:
                        if cp == current:
                            continue
                        for i in cp.traverse_postorder():
                            if i.label in non_uptree:
                                has_non_uptree_species = True
                    if has_non_uptree_species:
                        break
                    current = parent
                    parent = current.parent
                for j in current.traverse_postorder(internal=False):
                    if j.label not in c_rep_tree_labels:
                        removed.add(j.label + "\t" + node.label)
                parent.remove_child(current)
                parent.add_child(ctree.root)

        return astral_tree_par
    stitched_tree = _stitch(cg.root)
    deroot(stitched_tree)
    final_tree = join(options.output_fp, "didactic.nwk")
    stitched_tree.write_tree_newick(final_tree)
    unplaced = join(options.output_fp, "unplaced.csv")
    with open(unplaced, "w") as f:
        f.write("\n".join(removed) + "\n")

    return
    # for n in cg.traverse_postorder():
    #     print(outmap[n.label])
