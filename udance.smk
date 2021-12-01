
import os

configfile: "config.yaml"

#include: "workflows/decompose.smk"

outdir=config["output_dir"]

localrules: all, clean

rule all:
    input: t="%s/udance.nwk" % outdir

rule clean:
    params: outdir
    shell:
        """
            rm -r {params} data2/count.txt
        """

checkpoint decompose:
    input: j=config["input_jplace"], ind=config["input_dir"]
    output: cst=os.path.join(outdir,"color_spanning_tree.nwk")
    params:
        size=config["decompose_config"]["cluster_size"],
        method=config["decompose_config"]["infer"],
        outd=outdir,
    threads: config["decompose_config"]["cores"]
    shell:
        """
            python run_udance.py decompose -s {input.ind} -o {params.outd} -f {params.size} -j {input.j} -m {params.method} -T {threads}
        """

# phy inf
rule genetreeinfer:
    input:
        "%s/{cluster}/{gene}/run.sh" % outdir
    output:
        "%s/{cluster}/{gene}/RAxML_bestTree.file" % outdir
    shell:
        "bash {input}"


def aggregate_refine_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "%s/{j}/run.sh" % wildcards.cluster))
    return ["%s/%s/%s/RAxML_bestTree.file" % (outdir, wildcards.cluster,j) for j in wc.j]


rule refine:
    input: aggregate_refine_input
    output: "%s/{cluster}/astral_output.nwk" % outdir
    params: o=outdir, method=config["decompose_config"]["infer"]
    shell:
        """
            python run_udance.py refine -c {params.o}/{wildcards.cluster} -m {params.method} -M 1
        """

def aggregate_stitch_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{i}/species.txt"))
    return ["%s/%s/astral_output.nwk" % (outdir, i) for i in wc.i]


rule stitch:
    input: aggregate_stitch_input
    output: t="%s/udance.nwk" % outdir
    params: o=outdir
    shell:
        """
           python run_udance.py stitch -o {params.o}
        """

