import os
import stat
import shutil
from os.path import join, isfile, abspath, expanduser
from pathlib import Path
import treeswift as ts
import numpy as np
from didactic.compute_bipartition_alignment import compute_bipartition_alignment


class PoolAlignmentWorker:
    options = None
    species_dict = None
    fa_dict = None
    basename = None

    @classmethod
    def set_class_attributes(cls, options, species_dict, fa_dict, basename):
        cls.options = options
        cls.species_dict = species_dict
        cls.fa_dict = fa_dict
        cls.basename = basename

    #  TODO raise error if directory exists
    @classmethod
    def worker(cls, i, _):
        species = cls.species_dict[i]
        partition_aln = {key: cls.fa_dict[key] for key in species if key in cls.fa_dict}

        if len(partition_aln) < 4:
            return None
        aln_length = len(next(iter(partition_aln.values())))
        not_all_gap = np.array([False] * aln_length)
        for s in partition_aln.values():
            not_all_gap = np.logical_or(not_all_gap, (s != b'-'))
        for k, v in partition_aln.items():
            partition_aln[k] = v[not_all_gap]
        trimmed_aln_length = len(next(iter(partition_aln.values())))

        #deduplicate the alignment
        seq_keyed_dict = {}
        for name, sba in partition_aln.items():
            seq = sba.tostring().decode("UTF-8")
            if seq in seq_keyed_dict:
                seq_keyed_dict[seq].append(name)
            else:
                seq_keyed_dict[seq] = [name]

        if trimmed_aln_length >= cls.options.overlap_length and len(seq_keyed_dict) >= 4:
            # write trimmed MSA fasta
            res = []
            duplist = []
            for k, v in seq_keyed_dict.items():
                res.append(">" + v[0])
                res.append(k)
                if len(v) > 1:
                    duplist.append("\t".join(v))

            partition_output_dir = join(cls.options.output_fp, str(i))
            aln_outdir = join(partition_output_dir, cls.basename)
            Path(aln_outdir).mkdir(parents=True, exist_ok=True)
            aln_output_path = join(aln_outdir, "aln.fa")
            with open(aln_output_path, "w", buffering=100000000) as f:
                f.write("\n".join(res))
                f.write("\n")
            if duplist:
                dupmap_output_path = join(aln_outdir, "dupmap.txt")
                with open(dupmap_output_path, "w", buffering=100000000) as f:
                    f.write("\n".join(duplist))
                    f.write("\n")

            constraint_outgroup_tree = join(partition_output_dir, "raxml_constraint.nwk")
            if isfile(constraint_outgroup_tree) and cls.options.constrain_outgroups:
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
                    with open(induced_constraints_tree, "a") as a_file:
                        a_file.write("\n")

            script = join(aln_outdir, "run.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\n\n")
                f.write("export OMP_NUM_THREADS=1\n\n")
                bipartition_path = join(aln_outdir, "bipartition.fasta")
                fasttree_log = join(aln_outdir, "fasttree.log")
                fasttree_err = join(aln_outdir, "fasttree.err")
                fasttree_nwk = join(aln_outdir, "fasttree.nwk")
                if isfile(bipartition_path) and cls.options.constrain_outgroups and not cls.options.use_iqtree:
                    f.write("FastTreeMP -constraints %s -log %s < %s > %s 2> %s \n"
                            % (bipartition_path, fasttree_log, aln_output_path, fasttree_nwk, fasttree_err))
                elif not cls.options.use_iqtree:
                    f.write("FastTreeMP -log %s < %s > %s 2> %s \n"
                            % (fasttree_log, aln_output_path, fasttree_nwk, fasttree_err))

                fasttree_resolved_nwk = join(aln_outdir, "fasttree_resolved.nwk")
                if not cls.options.use_iqtree:
                    f.write("python3 -c \"import sys, treeswift; "
                            "t=treeswift.read_tree_newick(input()); "
                            "t.resolve_polytomies(); print(t)\" < %s > %s \n" % (fasttree_nwk, fasttree_resolved_nwk))

                induced_raxml_constraint_path = join(aln_outdir, "raxml_constraint.nwk")
                raxml_err = join(aln_outdir, "raxml.err")
                raxml_out = join(aln_outdir, "raxml.out")
                raxml_run = join(aln_outdir, "RUN")
                raxml_constraint_path = join(partition_output_dir, "raxml_constraint.nwk")
                astral_constraint_path = join(partition_output_dir, "astral_constraint.nwk")

                iqtree_err = join(aln_outdir, "iqtree.err")
                iqtree_out = join(aln_outdir, "iqtree.out")
                iqtree_run = join(aln_outdir, "RUN")
                iqtree_constraint_path = join(partition_output_dir, "iqtree_constraint.nwk")
                if isfile(induced_raxml_constraint_path) and cls.options.constrain_outgroups:
                    f.write("raxml-ng --force perf_threads --tree %s --tree-constraint %s "
                            "--msa %s --model LG+G --prefix %s --seed 12345 "
                            "--threads 4 > %s 2> %s \n"
                            % (fasttree_resolved_nwk, induced_raxml_constraint_path, aln_output_path,
                               raxml_run, raxml_out, raxml_err))

                else:
                    if cls.options.use_iqtree:
                        shutil.copy(astral_constraint_path, iqtree_constraint_path)
                        if cls.options.protein_seqs:
                            # f.write("iqtree -t %s -s %s --prefix %s "
                            #         "--seed 12345 -T 4 --model LG+G > %s 2> %s \n"
                            #         % (fasttree_resolved_nwk, aln_output_path,
                            #            iqtree_run, iqtree_out, iqtree_err))
                            if os.path.isfile(iqtree_constraint_path):
                                f.write("iqtree2 -s %s --prefix %s "
                                        "--seed 12345 -T AUTO --model LG+G -g %s > %s 2> %s \n"
                                        % (aln_output_path, iqtree_run,
                                           iqtree_constraint_path,
                                           iqtree_out, iqtree_err))
                            else:
                                f.write("iqtree2 -s %s --prefix %s "
                                        "--seed 12345 -T AUTO --model LG+G -g %s > %s 2> %s \n"
                                        % (aln_output_path, iqtree_run,
                                           iqtree_out, iqtree_err))
                        else:
                            if os.path.isfile(iqtree_constraint_path):
                                f.write("iqtree2 -s %s --prefix %s "
                                        "--seed 12345 -T AUTO --model GTR+F+G4 -g %s > %s 2> %s \n"
                                        % (aln_output_path, iqtree_run,
                                           iqtree_constraint_path,
                                           iqtree_out, iqtree_err))
                            else:
                                f.write("iqtree2 -s %s --prefix %s "
                                        "--seed 12345 -T AUTO --model GTR+F+G4 > %s 2> %s \n"
                                        % (aln_output_path, iqtree_run,
                                           iqtree_out, iqtree_err))
                    else:
                        f.write("raxml-ng --force perf_threads --tree %s "
                                "--msa %s --model LG+G --prefix %s --seed 12345 "
                                "--threads 4 > %s 2> %s \n"
                                % (fasttree_resolved_nwk, aln_output_path,
                                   raxml_run, raxml_out, raxml_err))
            st = os.stat(script)
            os.chmod(script, st.st_mode | stat.S_IEXEC)
            return trimmed_aln_length*len(partition_aln), script
        return None
