from os.path import join
import os
from pathlib import Path
import shutil
import treeswift as ts

from uDance.compute_bipartition_alignment import compute_bipartition_alignment


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
        try:
            os.remove(join(partition_output_dir, 'skip_partition'))
        except OSError:
            pass
        cls._undo_resolve_polytomies(j)
        newick_path = join(partition_output_dir, 'astral_constraint.nwk')
        j.write_tree_newick(newick_path)
        with open(newick_path, 'a') as a_file:
            a_file.write('\n')

        # find all outgroups
        with open(join(cls.options.output_fp, 'all_outgroups.txt')) as file:
            lines = file.readlines()
            outgroups_all = set([line.rstrip() for line in lines])
        outgroups_in_partition = [n.label for n in j.traverse_postorder() if n.label in outgroups_all]
        if len(outgroups_in_partition) >= 4:
            constraint = j.extract_tree_with(outgroups_in_partition, suppress_unifurcations=True)
            constraint.is_rooted = False
            bipartition_path = join(partition_output_dir, 'bipartition.fasta')
            with open(bipartition_path, 'w') as f:
                f.write(compute_bipartition_alignment(constraint.__str__()))
            raxml_constraint_path = join(partition_output_dir, 'raxml_constraint.nwk')
            constraint.write_tree_newick(raxml_constraint_path)
            with open(raxml_constraint_path, 'a') as a_file:
                a_file.write('\n')

        species_list_path = join(partition_output_dir, 'species.txt')
        edgeindices_list_path = join(partition_output_dir, 'edgeindices.txt')
        species_list = []
        pcount = 0
        skip = False
        with open(species_list_path, 'w') as f:
            with open(edgeindices_list_path, 'w') as f2:
                for n in j.traverse_postorder():
                    if not n.is_root() and hasattr(n, 'edge_index'):
                        f2.write(str(n.edge_index) + '\n')
                    if n.is_leaf():
                        species_list += [n.label]
                        f.write(n.label + '\n')
                    if hasattr(n, 'placements'):
                        for p in n.placements:
                            pcount += 1
                            species_list += [p]
                            f.write(p + '\n')
        if pcount <= cls.options.min_placements:
            print(
                'Number of placements on the partition %s is less than or equal to %d. '
                'Returning the backbone tree.' % (str(i), cls.options.min_placements)
            )
            skip = True
            Path(join(partition_output_dir, 'skip_partition')).touch()
            # jlabs = list(j.labels(internal=False))
            # backbone_tree_fp = join(cls.options.output_fp, "../backbone.nwk")
            # backbone_tree = ts.read_tree_newick(backbone_tree_fp)
            # extracted_tree = backbone_tree.extract_tree_with(jlabs)
            # input_tree_path = join(partition_output_dir, "astral_input.trees")
            # extracted_tree.write_tree_newick(input_tree_path)
            # for i in ["incremental", "updates"]:
            #     newick_path = join(partition_output_dir, "astral_output.%s.nwk" % i)
            #     newickbl_path = join(partition_output_dir, "astral_output.%s.nwk.bl" % i)
            #     extracted_tree.write_tree_newick(newick_path)
            #     extracted_tree.write_tree_newick(newickbl_path)
            #     log_path = join(partition_output_dir, "astral.%s.log" % i)
            #     with open(log_path, "w") as f:
            #         f.write("Final quartet score is 1\n")

        return species_list_path, skip
