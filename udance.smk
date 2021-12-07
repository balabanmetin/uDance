
import os

# configfile: "config.yaml"

#include: "workflows/decompose.smk"

wdr = config["workdir"]
outdir = os.path.join(wdr, "output")
alndir = os.path.join(wdr, "alignments")
bbone = os.path.join(wdr, "backbone.nwk")

localrules: all, clean

rule all:
    input: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])

rule clean:
    params: outdir
    shell:
        """
            rm -r {params} data2/count.txt
        """

rule placement:
    input: b = bbone, ind = alndir
    output: j=os.path.join(outdir,"placement.jplace"), d=directory(os.path.join(outdir,"placement"))
    params: o=outdir,
            f=config["apples_config"]["filter"],
            m=config["apples_config"]["method"], b=config["apples_config"]["base"]
    log: out=os.path.join(outdir,"placement/apples2.out"), err=os.path.join(outdir,"placement/apples2.err")
    shell:
        """
            bash uDance/create_concat_alignment.sh {input.ind} {input.b} {params.o}
            run_apples.py -s {output.d}/backbone.fa -q {output.d}/query.fa \
            -t {output.d}/backbone.tree -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
        """


checkpoint decompose:
    input: j=os.path.join(outdir,"placement.jplace"), ind=alndir
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
    output: expand("%s/{{cluster}}/astral_output.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: o=outdir, method=config["decompose_config"]["infer"]
    shell:
        """
            python run_udance.py refine -c {params.o}/{wildcards.cluster} -m {params.method} -M 1
        """

def aggregate_stitch_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{i}/species.txt"))
    return [f"%s/%s/astral_output.%s.nwk" % (outdir, i, j) for i in wc.i for j in ["incremental", "updates"]]


rule stitch:
    input: aggregate_stitch_input
    output: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: o=outdir
    shell:
        """
           python run_udance.py stitch -o {params.o}
        """

