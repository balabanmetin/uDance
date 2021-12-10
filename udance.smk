
import os
import multiprocessing as mp
from uDance.prep_partition_alignments import prep_partition_alignments

# configfile: "config.yaml"

#include: "workflows/decompose.smk"

wdr = config["workdir"]
outdir = os.path.join(wdr, "output")
alndir = os.path.join(wdr, "alignments")
bbspec = os.path.join(wdr, "species.txt")
bbone = os.path.join(outdir, "backbone.nwk")
trimalndir = os.path.join(outdir, "trimmed")


TRIMMEDGENES = [os.path.join(outdir, "trimdump", f) for f in os.listdir(alndir) if
                  os.path.isfile(os.path.join(alndir, f))]

localrules: all, clean

rule all:
    input: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])

rule clean:
    params: outdir
    shell:
        """
            rm -r {params} data2/count.txt
        """

rule trimtaper:
    input: "%s/{gene}" % alndir
    output: "%s/trimdump/{gene}" % outdir
    params: thr=config["trim_config"]["percent_nongap"]
    shell:
        """
            TFA=`mktemp -t XXXXXX.fa`
            trimal -in {input} -out $TFA -gt {params.thr} 
            julia uDance/correction_multi.jl $TFA > {output}
            rm $TFA
        """

rule trimcollect:
    input: TRIMMEDGENES
    output: directory(os.path.join(outdir, "trimmed"))
    params: o=outdir
    shell:
        """
            mv {params.o}/trimdump {output}
        """

rule mainlines:
    input: trimalndir
    output: bbspec
    params:
            n=config["mainlines_config"]["n"],
            l=config["mainlines_config"]["length"],
            char=config["chartype"]
    shell:
        """
            if [ "{params.char}" == "nuc" ]; then
                python run_udance.py mainlines -s {input} -n {params.n} -l {params.l} > {output}
            else
                python run_udance.py mainlines -s {input} -n {params.n} -l {params.l} -p > {output}
            fi
        """

rule copyspeciesfile:
    input: s=bbspec
    output: os.path.join(outdir, "backbone/0/species.txt")
    shell:
        """
            cp {input} {output}
        """
checkpoint prepbackbonegenes:
    input: s=os.path.join(outdir, "backbone/0/species.txt"), a=trimalndir
    output: touch(os.path.join(outdir,"backbone/0/done.txt"))
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
    return ["%s/backbone/0/%s/bestTree.nwk" % (outdir,j) for j in wc.j]


rule refine_bb:
    input: aggregate_refine_bb_input
    output: bbone
    params:
        o=outdir,
        method=config["infer_config"]["method"],
        c=config["refine_config"]["contract"]
    shell:
        '''
            python run_udance.py refine -p {params.o}/backbone/0 -m {params.method} -M 1 -c {params.c}
            nw_reroot -d {params.o}/backbone/0/astral_output.incremental.nwk > {output}
        '''

rule placement:
    input: b = bbone, ind = trimalndir
    output: j=os.path.join(outdir,"placement.jplace"), d=directory(os.path.join(outdir,"placement"))
    params: o=outdir,
            f=config["apples_config"]["filter"],
            m=config["apples_config"]["method"],
            b=config["apples_config"]["base"],
            char=config["chartype"]
    resources: cpus=config["prep_config"]["cores"]
    log: out=os.path.join(outdir,"placement/apples2.out"), err=os.path.join(outdir,"placement/apples2.err")
    shell:
        """
            bash uDance/create_concat_alignment.sh {input.ind} {input.b} {params.o}
            if [ "{params.char}" == "nuc" ]; then
                run_apples.py -s {output.d}/backbone.fa -q {output.d}/query.fa -T {resources.cpus} \
                -t {output.d}/backbone.tree -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
            else
                run_apples.py -p -s {output.d}/backbone.fa -q {output.d}/query.fa -T {resources.cpus} \
                -t {output.d}/backbone.tree -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
            fi
        """


checkpoint decompose:
    input: j=os.path.join(outdir,"placement.jplace"), ind=trimalndir
    output: cst=os.path.join(outdir,"udance/color_spanning_tree.nwk")
    params:
        size=config["prep_config"]["cluster_size"],
        method=config["infer_config"]["method"],
        outd=outdir,
        ov=config["prep_config"]["overlap"]
    resources: cpus=config["prep_config"]["cores"]
    shell:
        """
            cp {input.j} {params.outd}/udance
            python run_udance.py decompose -s {input.ind} -o {params.outd}/udance -f {params.size} -j {input.j} -m {params.method} -T {resources.cpus} -l {params.ov}
        """

# phy inf
rule genetreeinfer:
    input:
        "%s/{stage}/{cluster}/{gene}/aln.fa" % outdir
    output:
        "%s/{stage}/{cluster}/{gene}/bestTree.nwk" % outdir
    params:
        c=config["chartype"], s=config["infer_config"]["numstart"]
    shell:
        '''
            bash uDance/process_a_marker.sh {input} {params.c} {params.s} > $(dirname {input})/process.log 
        '''


def aggregate_refine_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "%s/{j}/aln.fa" % wildcards.cluster))
    return ["%s/udance/%s/%s/bestTree.nwk" % (outdir, wildcards.cluster,j) for j in wc.j]


rule refine:
    input: aggregate_refine_input
    output: expand("%s/udance/{{cluster}}/astral_output.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: o=outdir,
            method=config["infer_config"]["method"],
            c=config["refine_config"]["contract"]
    shell:
        """
            python run_udance.py refine -p {params.o}/udance/{wildcards.cluster} -m {params.method} -M 1 -c {params.c}
        """

def aggregate_stitch_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{i}/species.txt"))
    return [f"%s/udance/%s/astral_output.%s.nwk" % (outdir, i, j) for i in wc.i for j in ["incremental", "updates"]]


rule stitch:
    input: aggregate_stitch_input
    output: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: o=outdir
    shell:
        """
           python run_udance.py stitch -o {params.o}/udance
           cp {outdir}/udance/udance.*.nwk {outdir}
        """

