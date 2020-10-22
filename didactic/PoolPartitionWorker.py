from os.path import join
from pathlib import Path

from didactic.compute_bipartition_alignment import compute_bipartition_alignment


class PoolPartitionWorker:
    options = None

    @classmethod
    def set_class_attributes(cls, options):
        cls.options = options

    @staticmethod
    def _undo_resolve_polytomies(tree):
        for e in tree.traverse_postorder():
            if e == tree.root:
                continue
            if e.outgroup is False and e.resolved_randomly:  # edge length check just for safety
                par = e.parent
                par.remove_child(e)
                for c in e.children:
                    par.add_child(c)

    @classmethod
    def worker(cls, i, j):
        partition_output_dir = join(cls.options.output_fp, str(i))
        Path(partition_output_dir).mkdir(parents=True, exist_ok=True)
        cls._undo_resolve_polytomies(j)
        newick_path = join(partition_output_dir, "astral_constraint.nwk")
        j.write_tree_newick(newick_path)
        with open(newick_path, "a") as a_file:
            a_file.write("\n")

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
            with open(raxml_constraint_path, "a") as a_file:
                a_file.write("\n")

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