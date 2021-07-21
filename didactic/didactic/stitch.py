import json
from os.path import join
import treeswift as ts


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
        if uptreestr: # there is an uptree
            uptree = ts.read_tree_newick(outmap_par["up"])

            uptree_labels = set(uptree.labels(internal=False))
            non_uptree = astral_tree_cons_labels.difference(uptree_labels)

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
        else: # find a query. root at its MRCA. find mrca of one of the childs (representatives). root at there.
            non_uptree = astral_tree_cons_labels
            for i in astral_tree_par.traverse_postorder(internal=False):
                if i.label not in astral_tree_cons_labels:
                    notconstree_species = i   # this has to be a query species
                    break
            astral_tree_par.is_rooted = True
            astral_tree_par.reroot(notconstree_species.parent)
            astral_tree_par.suppress_unifurcations()

            first_c = node.children[0]
            first_c_rep_tree = ts.read_tree_newick(outmap_par["children"][first_c.label])
            first_c_rep_tree_labels = set(first_c_rep_tree.labels(internal=False))
            first_c_rep_tree_mrca = astral_tree_par.mrca(list(first_c_rep_tree_labels))
            astral_tree_par.reroot(first_c_rep_tree_mrca.parent)
            astral_tree_par.suppress_unifurcations()

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
    final_tree = join(options.output_fp, "didactic.nwk")
    stitched_tree.write_tree_newick(final_tree)
    unplaced = join(options.output_fp, "unplaced.csv")
    with open(unplaced, "w") as f:
        f.write("\n".join(removed) + "\n")

    return
    # for n in cg.traverse_postorder():
    #     print(outmap[n.label])