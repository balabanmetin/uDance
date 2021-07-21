from os.path import join
from glob import glob
from pathlib import Path
from sys import stderr
import shutil
from subprocess import Popen, PIPE

from didactic.expand_dedupe_newick import expand_dedupe_newick


class PoolAstralWorker:
    options = None
    astral_exec = None

    @classmethod
    def set_class_attributes(cls, options, astral_exec):
        cls.options = options
        cls.astral_exec = astral_exec

    @classmethod
    def worker(cls, i):
        partition_output_dir = join(cls.options.output_fp, str(i))
        genes = glob(join(partition_output_dir, "*", ""))
        for gene in genes:
            if cls.options.use_iqtree:
                best = Path(join(gene, 'RUN.treefile'))
            else:
                best = Path(join(gene, 'RUN.raxml.bestTree'))
            bestCollapsed = Path(join(gene, 'RUN.raxml.bestTreeCollapsed'))
            if bestCollapsed.is_file():
                raxtree = bestCollapsed
            elif best.is_file():
                raxtree = best
            else:
                stderr.write("%s/RUN.raxml.bestTree does not exist. RAxML job is corrupted. \n" % gene)
                continue
            treestr = open(raxtree).readline()
            dupmap_file = Path(join(gene, "dupmap.txt"))
            if dupmap_file.is_file():
                dmp = list(map(lambda x: x.strip().split("\t"), open(dupmap_file).readlines()))
                expanded_tree_str = expand_dedupe_newick(treestr, dmp)
            else:
                expanded_tree_str = treestr
            with open(join(gene, "raxml.expanded.nwk"), "w") as out:
                out.write(expanded_tree_str)
        expanded_trees = glob(join(partition_output_dir, "*", "raxml.expanded.nwk"))
        astral_input_file = join(partition_output_dir, "astral_input.trees")

        with open(astral_input_file, 'wb') as wfd:
            for f in expanded_trees:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

        astral_output_file = join(partition_output_dir, "astral_output.nwk")
        astral_log_file = join(partition_output_dir, "astral.log")
        astral_const_file = join(partition_output_dir, "astral_constraint.nwk")
        s = f'cp {astral_input_file} {astral_output_file}\n'
        # s = ["java", "-jar", cls.astral_exec, "-i", astral_input_file,
        #       "-o", astral_output_file, "-j", astral_const_file, "2>", astral_log_file]
        return s

        # with open(astral_log_file, "w") as lg:
        #     with Popen(s, stdout=PIPE, stdin=PIPE, stderr=lg) as p:
        #         astral_stdout = p.stdout.read().decode('utf-8')
        #         #print(astral_stdout)
