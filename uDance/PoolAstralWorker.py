from os.path import join
from glob import glob
from pathlib import Path
from sys import stderr
import shutil
from subprocess import Popen, PIPE
import treeswift as ts
from statistics import median
from kmeans1d import cluster

from uDance.expand_dedupe_newick import expand_dedupe_newick


class PoolAstralWorker:
    options = None
    astral_exec = None

    @classmethod
    def set_class_attributes(cls, options, astral_exec):
        cls.options = options
        cls.astral_exec = astral_exec

    @classmethod
    def worker(cls, partition_output_dir):
        genes = glob(join(partition_output_dir, "*", ""))
        median_map = dict()
        for gene in genes:
            if cls.options.method == 'iqtree':
                best = Path(join(gene, 'RUN.treefile'))
            elif cls.options.method == 'raxml-ng':
                best = Path(join(gene, 'RUN.raxml.bestTree'))
            elif cls.options.method == 'raxml-8':
                best = Path(join(gene, 'bestTree.nwk'))
            bestCollapsed = Path(join(gene, 'RUN.raxml.bestTreeCollapsed'))
            if bestCollapsed.is_file():
                raxtree = bestCollapsed
            elif best.is_file():
                raxtree = best
            else:
                stderr.write("%s/bestTree.nwk does not exist. RAxML job is corrupted. \n" % gene)
                continue
            with open(raxtree) as f:
                treestr = f.readline()
            tf = ts.read_tree_newick(treestr)
            lpps = [float(i.label) for i in tf.traverse_postorder(leaves=False) if i.label]
            if len(lpps) > 0:
                median_map[gene] = median(lpps)
            dupmap_file = Path(join(gene, "dupmap.txt"))
            if dupmap_file.is_file():
                dmp = list(map(lambda x: x.strip().split("\t"), open(dupmap_file).readlines()))
                expanded_tree_str = expand_dedupe_newick(treestr, dmp)
            else:
                expanded_tree_str = treestr
            with open(join(gene, "raxml.expanded.nwk"), "w") as out:
                out.write(expanded_tree_str)
        # remove outlier genes. outlier is defined as having lower median local posterior probability than majority
        # we use 1d k-means (k=2) for outlier detection.
        clusters, centroids = cluster(list(median_map.values()), k=2)
        if 0.8 < sum(clusters)/len(clusters) < 1:
            min_median = min([v for i, v in enumerate(median_map.values()) if clusters[i] == 0])
            numdiscard = len(clusters) - sum(clusters)
            print("In cluster %s, %d gene tree(s) with lower than "
                  "%.2f median lpp are discarded." % (partition_output_dir, numdiscard, min_median), file=stderr)
            expanded_trees = [join(gene, "raxml.expanded.nwk") for i, gene in enumerate(median_map.keys()) if clusters[i] == 1]
        else:
            expanded_trees = glob(join(partition_output_dir, "*", "raxml.expanded.nwk"))
        astral_input_file = join(partition_output_dir, "astral_input.trees")

        with open(astral_input_file, 'wb') as wfd:
            for f in expanded_trees:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

        astral_output_file, astral_log_file, astral_const_file = [dict(), dict(), dict()]
        astral_const_file["incremental"] = join(partition_output_dir, "astral_constraint.nwk")
        astral_const_file["updates"] = join(partition_output_dir, "raxml_constraint.nwk")

        for mtd in ["incremental", "updates"]:
            astral_output_file[mtd] = Path(join(partition_output_dir, "astral_output.%s.nwk" % mtd))
            if mtd == "updates" and not Path(astral_const_file["updates"]).is_file() and not Path(astral_const_file["incremental"]).is_file():
                shutil.copyfile(astral_output_file["incremental"], astral_output_file["updates"])
                break
            astral_log_file[mtd] = join(partition_output_dir, "astral.%s.log" % mtd)
            s = ["java", "-Xmx%sG" % cls.options.memory, "-jar", cls.astral_exec, "-i", astral_input_file,
                 "-o", astral_output_file[mtd]]
            if Path(astral_const_file[mtd]).is_file():
                s += ["-j", astral_const_file[mtd]]
            with open(astral_log_file[mtd], "w") as lg:
                with Popen(s, stdout=PIPE, stdin=PIPE, stderr=lg) as p:
                    astral_stdout = p.stdout.read().decode('utf-8')
                    # print(astral_stdout)
        # if cls.options.use_gpu:
        #     gpu_opt = ""
        # else:
        #     gpu_opt = "-C"

        # s = f'cp {astral_const_file} {astral_output_file}\n'
        # s = ["cp", astral_const_file, astral_output_file]
