
import os
import multiprocessing as mp
from uDance.prep_partition_alignments import prep_partition_alignments

# configfile: "config.yaml"

#include: "workflows/decompose.smk"

wdr = config["workdir"]
outdir = os.path.join(wdr, "output")
alndir = os.path.join(wdr, "alignments")
bbone = os.path.join(outdir, "backbone.nwk")

localrules: all, clean

rule all:
    input: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])

rule clean:
    params: outdir
    shell:
        """
            rm -r {params} data2/count.txt
        """

checkpoint prepbackbonegenes:
    input: s=os.path.join(outdir,"backbone/species.txt"), a=alndir
    output: touch(os.path.join(outdir,"backbone/done.txt"))
    resources: cpus=config["prep_config"]["cores"]
    params: ov=config["prep_config"]["overlap"]
    run:
        mp.set_start_method('fork')
        prep_partition_alignments(input.a,
                                  config["chartype"] == "prot",
                                  [input.s],
                                  resources.cpus,
                                  params.ov
                                  )



def aggregate_refine_bb_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.prepbackbonegenes.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{j}/aln.fa"))
    return ["%s/backbone/%s/bestTree.nwk" % (outdir,j) for j in wc.j]


rule refine_bb:
    input: aggregate_refine_bb_input
    output: bbone
    params:
        o=outdir,
        method=config["infer_config"]["method"]
    shell:
        '''
            python run_udance.py refine -c {params.o}/backbone -m {params.method} -M 1
            nw_reroot -d {params.o}/backbone/astral_output.incremental.nwk > {output}
        '''

rule placement:
    input: b = bbone, ind = alndir
    output: j=os.path.join(outdir,"placement.jplace"), d=directory(os.path.join(outdir,"placement"))
    params: o=outdir,
            f=config["apples_config"]["filter"],
            m=config["apples_config"]["method"], b=config["apples_config"]["base"]
    resources: cpus=config["prep_config"]["cores"]
    log: out=os.path.join(outdir,"placement/apples2.out"), err=os.path.join(outdir,"placement/apples2.err")
    shell:
        """
            bash uDance/create_concat_alignment.sh {input.ind} {input.b} {params.o}
            run_apples.py -s {output.d}/backbone.fa -q {output.d}/query.fa -T {resources.cpus} \
            -t {output.d}/backbone.tree -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
        """


checkpoint decompose:
    input: j=os.path.join(outdir,"placement.jplace"), ind=alndir
    output: cst=os.path.join(outdir,"color_spanning_tree.nwk")
    params:
        size=config["prep_config"]["cluster_size"],
        method=config["infer_config"]["method"],
        outd=outdir,
        ov=config["prep_config"]["overlap"]
    resources: cpus=config["prep_config"]["cores"]
    shell:
        """
            python run_udance.py decompose -s {input.ind} -o {params.outd} -f {params.size} -j {input.j} -m {params.method} -T {resources.cpus} -l {params.ov}
        """

# phy inf
rule genetreeinfer:
    input:
        "%s/{cluster}/{gene}/aln.fa" % outdir
    output:
        "%s/{cluster}/{gene}/bestTree.nwk" % outdir
    params:
        c=config["chartype"], s=config["infer_config"]["numstart"]
    shell:
        '''
            bash uDance/process_a_marker.sh {input} {params.c} {params.s} > $(dirname {input})/process.log 
        '''


def aggregate_refine_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "%s/{j}/aln.fa" % wildcards.cluster))
    return ["%s/%s/%s/bestTree.nwk" % (outdir, wildcards.cluster,j) for j in wc.j]


rule refine:
    input: aggregate_refine_input
    output: expand("%s/{{cluster}}/astral_output.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: o=outdir, method=config["infer_config"]["method"]
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
