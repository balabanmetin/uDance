
import os

#include: "workflows/decompose.smk"
wdr = config["workdir"]
outdir = os.path.join(wdr, "output")
alndir = os.path.join(wdr, "alignments")
bbspec = os.path.join(wdr, "species.txt")
bbone = os.path.join(outdir, "backbone.nwk")
trimalndir = os.path.join(outdir, "trimmed")

TRIMMEDGENES = [os.path.join(outdir, "trimdump", f) for f in os.listdir(alndir) if
                  os.path.isfile(os.path.join(alndir, f))]

udance_logpath = os.path.abspath(os.path.join(wdr, "udance.log"))

localrules: all, clean, copyspeciesfile, trimcollect

rule all:
    input: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])

onstart:
    shell( "if [ -f '{udance_logpath}' ]; then echo '{udance_logpath} already exists."
            "Moving it to {udance_logpath}.bak'; mv {udance_logpath} {udance_logpath}.bak;"
            "else touch {udance_logpath}; fi")

onerror:
    print("Execution failed. Logfile of the execution: ")
    print(udance_logpath)

onsuccess:
    print("Execution is successful. Logfile of the execution: ")
    print(udance_logpath)

rule clean:
    shell:
        """
            rm -r {outdir} 
        """

rule trimtaper:
    input: "%s/{gene}" % alndir
    output: "%s/trimdump/{gene}" % outdir
    params: thr=config["trim_config"]["percent_nongap"]
    shell:
        """
            (
            TFA=`mktemp -t XXXXXX.fa`
            trimal -in {input} -out $TFA -gt {params.thr} 
            julia uDance/correction_multi.jl $TFA > {output}
            rm $TFA
            ) >> {udance_logpath} 2>&1
        """

rule trimcollect:
    input: TRIMMEDGENES
    output: directory(os.path.join(outdir, "trimmed"))
    shell:
        """
            (
            mv {outdir}/trimdump {output}
            ) >> {udance_logpath} 2>&1
        """

rule mainlines:
    input: trimalndir
    output: bbspec
    params:
            n=config["mainlines_config"]["n"],
            l=config["mainlines_config"]["length"],
            char=config["chartype"]
    resources: mem_mb=config["resources"]["large_memory"]
    shell:
        """
            (
            if [ "{params.char}" == "nuc" ]; then
                python run_udance.py mainlines -s {input} -n {params.n} -l {params.l} > {output}
            else
                python run_udance.py mainlines -s {input} -n {params.n} -l {params.l} -p > {output}
            fi
            ) >> {udance_logpath} 2>&1
        """

rule copyspeciesfile:
    input: s=bbspec
    output: os.path.join(outdir, "backbone/0/species.txt")
    shell:
        """
            (
            cp {input} {output}
            ) >> {udance_logpath} 2>&1
        """
checkpoint prepbackbonegenes:
    input: s=os.path.join(outdir, "backbone/0/species.txt"), a=trimalndir
    output: touch(os.path.join(outdir,"backbone/0/done.txt"))
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    params: sub=config["prep_config"]["sublength"],
            frag=config["prep_config"]["fraglength"],
            char=config["chartype"]
    shell:
         # python calling shell calling python. looks terrible but
         # this way we are avoiding forking in snakemake main process
        '''
            (
            python -c  "import multiprocessing as mp; \
                        mp.set_start_method('fork'); \
                        from uDance.prep_partition_alignments import prep_partition_alignments; \
                        prep_partition_alignments('{input.a}', \
                                      '{params.char}' == 'prot', \
                                      ['{input.s}'], \
                                      {resources.cpus}, \
                                      {params.sub}, \
                                      {params.frag})"
            ) >> {udance_logpath} 2>&1
        '''
        # with open(params.logpath, "a") as log_file:
        #     print(params.logpath)



def aggregate_refine_bb_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.prepbackbonegenes.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{j}/aln.fa"))
    return ["%s/backbone/0/%s/bestTree.nwk" % (outdir,j) for j in wc.j]


rule refine_bb:
    input: aggregate_refine_bb_input
    output: bbone = bbone
    params:
        method=config["infer_config"]["method"],
        c=config["refine_config"]["contract"],
        occup = config["refine_config"]["occupancy"]
    resources: mem_mb=config["resources"]["large_memory"]
    shell:
        '''
            (
            python run_udance.py refine -p {outdir}/backbone/0 -m {params.method} -M {resources.mem_mb} -c {params.c} -o {params.occup}
            nw_reroot -d {outdir}/backbone/0/astral_output.incremental.nwk > {output}
            ) >> {udance_logpath} 2>&1
        '''

rule placement_prep:
    input: b = bbone,
           ind = trimalndir
    output: aln = os.path.join(outdir,"placement/backbone.fa"),
            qry = os.path.join(outdir,"placement/query.fa"),
            tre = os.path.join(outdir,"placement/backbone.tree")
    params: char=config["chartype"]
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    shell:
        """
            (
            bash uDance/create_concat_alignment.sh {input.ind} {input.b} {outdir} {params.char} {resources.cpus}
            ) >> {udance_logpath} 2>&1            
        """

rule placement:
    input: aln = os.path.join(outdir,"placement/backbone.fa"),
           qry = os.path.join(outdir,"placement/query.fa"),
           tre = os.path.join(outdir,"placement/backbone.tree")
    output: j=os.path.join(outdir,"placement.jplace")
    params: f=config["apples_config"]["filter"],
            m=config["apples_config"]["method"],
            b=config["apples_config"]["base"],
            char=config["chartype"]
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    log: out=os.path.join(outdir,"placement/apples2.out"), err=os.path.join(outdir,"placement/apples2.err")
    shell:
        """
            (
            if [ "{params.char}" == "nuc" ]; then
                run_apples.py -s {input.aln} -q {input.qry} -T {resources.cpus} \
                -t {input.tre} -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
            else
                run_apples.py -p -s {input.aln} -q {input.qry} -T {resources.cpus} \
                 -t {input.tre} -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
            fi
            ) >> {udance_logpath} 2>&1            
        """


checkpoint decompose:
    input: j=os.path.join(outdir,"placement.jplace"), ind=trimalndir
    output: cst=os.path.join(outdir,"udance/color_spanning_tree.nwk")
    params:
        size=config["prep_config"]["cluster_size"],
        method=config["infer_config"]["method"],
        sub=config["prep_config"]["sublength"],
        frag=config["prep_config"]["fraglength"],
        char=config["chartype"]

    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    shell:
        """
            (
            cp {input.j} {outdir}/udance
            if [ "{params.char}" == "nuc" ]; then
                python run_udance.py decompose -s {input.ind} -o {outdir}/udance -t {params.size} -j {input.j} -m {params.method} -T {resources.cpus} -l {params.sub} -f {params.frag}
            else
                python run_udance.py decompose -p -s {input.ind} -o {outdir}/udance -t {params.size} -j {input.j} -m {params.method} -T {resources.cpus} -l {params.sub} -f {params.frag}
            fi
            ) >> {udance_logpath} 2>&1
        """

# phy inf
rule genetreeinfer:
    input:
        "%s/{stage}/{cluster}/{gene}/aln.fa" % outdir
    output:
        "%s/{stage}/{cluster}/{gene}/bestTree.nwk" % outdir
    params:
          c=config["chartype"],
          s=config["infer_config"]["numstart"],
          t=config["infer_config"]["method"]
    shell:
        '''
            # many instances of this rule may run simultaneously. To reduce the IO overhead, we "sponge" the output
            # before appending to udance_logpath
            source uDance/mysponge.sh
            (
            bash uDance/process_a_marker.sh {input} {params.c} {params.s} {params.t}
            ) 2>&1 | mysponge >> {udance_logpath} 
        '''



def aggregate_refine_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "%s/{j}/aln.fa" % wildcards.cluster))
    return ["%s/udance/%s/%s/bestTree.nwk" % (outdir, wildcards.cluster,j) for j in wc.j]


rule refine:
    input: aggregate_refine_input
    output: expand("%s/udance/{{cluster}}/astral_output.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: method=config["infer_config"]["method"],
            c=config["refine_config"]["contract"],
            occup=config["refine_config"]["occupancy"]
    resources: mem_mb=config["resources"]["large_memory"]
    shell:
        """
            (
            python run_udance.py refine -p {outdir}/udance/{wildcards.cluster} -m {params.method} -M {resources.mem_mb} -c {params.c} -o {params.occup}
            ) >> {udance_logpath} 2>&1
        """

def aggregate_stitch_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{i}/species.txt"))
    return [f"%s/udance/%s/astral_output.%s.nwk" % (outdir, i, j) for i in wc.i for j in ["incremental", "updates"]]


rule stitch:
    input: aggregate_stitch_input
    output: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    shell:
        """
           (
           python run_udance.py stitch -o {outdir}/udance
           cp {outdir}/udance/udance.*.nwk {outdir}
           ) >> {udance_logpath} 2>&1
        """
