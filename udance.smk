
import os

#include: "workflows/decompose.smk"
wdr = config["workdir"]
outdir = os.path.join(wdr, "output")
alndir = os.path.join(wdr, "alignments")
bbspec = os.path.join(wdr, "species.txt")
input_bbone = os.path.join(wdr, "backbone.nwk")
bbone = os.path.join(outdir, "backbone.nwk")
trimalndir = os.path.join(outdir, "trimmed")

TRIMMEDGENES = [os.path.join(outdir, "trimdump", f) for f in os.listdir(alndir) if
                  os.path.isfile(os.path.join(alndir, f))]

udance_logpath = os.path.abspath(os.path.join(wdr, "udance.log"))

localrules: all, clean, trimcollect, copybb

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
    benchmark: "%s/benchmarks/trimtaper.{gene}.txt" % outdir
    shell:
        """
            (
            uDance/trimtaper.sh {input} {params.thr} {output}
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
    output: os.path.join(outdir, "backbone/0/species.txt")
    params:
            n=config["mainlines_config"]["n"],
            l=config["mainlines_config"]["length"],
            char=config["chartype"],
            bck=config["backbone"]
    resources: mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/mainlines.txt" % outdir
    shell:
        """
            (
            if [ "{params.bck}" == "list" ]; then
                cp {bbspec} {output}
            elif [ "{params.bck}" == "tree" ]; then
                nw_labels -I {input_bbone} > {output}
            elif [ "{params.char}" == "nuc" ]; then  # denovo
                python run_udance.py mainlines -s {input} -n {params.n} -l {params.l} > {output}
            else
                python run_udance.py mainlines -s {input} -n {params.n} -l {params.l} -p > {output}
            fi
            ) >> {udance_logpath} 2>&1
        """

checkpoint prepbackbonegenes:
    input: s=os.path.join(outdir, "backbone/0/species.txt"), a=trimalndir
    output: touch(os.path.join(outdir,"backbone/0/done.txt"))
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/prepbackbonegenes.txt" % outdir
    params: sub=config["prep_config"]["sublength"],
            frag=config["prep_config"]["fraglength"],
            char=config["chartype"],
            bck=config["backbone"]
    shell:
         # python calling shell calling python. looks terrible but
         # this way we are avoiding forking in snakemake main process
        '''
            (
            if [ "{params.bck}" == "tree" ]; then
                touch {outdir}/backbone/0/done.txt
            else
                python -c  "import multiprocessing as mp; \
                            mp.set_start_method('fork'); \
                            from uDance.prep_partition_alignments import prep_partition_alignments; \
                            prep_partition_alignments('{input.a}', \
                                          '{params.char}' == 'prot', \
                                          ['{input.s}'], \
                                          {resources.cpus}, \
                                          {params.sub}, \
                                          {params.frag})"
            fi
            ) >> {udance_logpath} 2>&1
        '''
        # with open(params.logpath, "a") as log_file:
        #     print(params.logpath)



def aggregate_refine_bb_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.prepbackbonegenes.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{j}/aln.fa"))
    return ["%s/backbone/0/%s/bestTree.nwk" % (outdir,j) for j in wc.j]

if config["backbone"] != "tree":
    rule refine_bb:
        input: aggregate_refine_bb_input
        output: bbone = bbone
        params:
            method=config["infer_config"]["method"],
            c=config["refine_config"]["contract"],
            occup = config["refine_config"]["occupancy"],
            ol=config["refine_config"]["outlier_sizelimit"],
            od=config["refine_config"]["outlier_difference"]
        resources: cpus=config["resources"]["cores"],
                   mem_mb=config["resources"]["large_memory"]
        benchmark: "%s/benchmarks/refine_copy_bb.txt" % outdir
        shell:
            '''
                (
                python run_udance.py refine -p {outdir}/backbone/0 -m {params.method} -M {resources.mem_mb} -c {params.c} -o {params.occup} -T {resources.cpus} -l {params.ol} -d {params.od}
                nw_reroot -d {outdir}/backbone/0/astral_output.incremental.nwk > {output}
                ) >> {udance_logpath} 2>&1
            '''
else:
    rule copybb:
        input: os.path.join(outdir,"backbone/0/done.txt")
        output: bbone = bbone
        benchmark: "%s/benchmarks/refine_copy_bb.txt" % outdir
        shell:
            '''
                cp {input_bbone} {output.bbone}
            '''

rule placement_prep:
    input: b = bbone,
           ind = trimalndir
    output: aln = os.path.join(outdir,"placement/backbone.fa"),
            qry = os.path.join(outdir,"placement/query.fa"),
            tre = os.path.join(outdir,"placement/backbone.tree")
    params: char=config["chartype"],
            filtering=config["backbone_filtering"],
            f=config["apples_config"]["filter"],
            m=config["apples_config"]["method"],
            b=config["apples_config"]["base"],
            v=config["apples_config"]["overlap"]
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/placement_prep.txt" % outdir
    shell:
        """
            (
            bash uDance/create_concat_alignment.sh {input.ind} {input.b} {outdir} {params.char} {resources.cpus} {params.f} {params.m} {params.b} {params.v} {params.filtering}
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
            v=config["apples_config"]["overlap"],
            char=config["chartype"]
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/placement.txt" % outdir
    log: out=os.path.join(outdir,"placement/apples2.out"), err=os.path.join(outdir,"placement/apples2.err")
    shell:
        """
            (
            export MKL_NUM_THREADS=1
            export NUMEXPR_NUM_THREADS=1
            export OMP_NUM_THREADS=1
            if [ "{params.char}" == "nuc" ]; then
                run_apples.py --exclude -s {input.aln} -q {input.qry} -T {resources.cpus} -V {params.v} \
                -t {input.tre} -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
            else
                run_apples.py --exclude -p -s {input.aln} -q {input.qry} -T {resources.cpus} -V {params.v} \
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
        edg=config["prep_config"]["edge_thr"],
        sub=config["prep_config"]["sublength"],
        frag=config["prep_config"]["fraglength"],
        pra=config["prep_config"]["pruneafter"],
        mps=config["prep_config"]["min_placements"],
        char=config["chartype"]

    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/decompose.txt" % outdir
    shell:
        """
            (
            cp {input.j} {outdir}/udance
            if [ "{params.char}" == "nuc" ]; then
                python run_udance.py decompose -s {input.ind} -o {outdir}/udance -t {params.size} -j {input.j} \
                -m {params.method} -T {resources.cpus} -l {params.sub} -f {params.frag} -e {params.edg} \
                --minplacements {params.mps}
            else
                python run_udance.py decompose -p -s {input.ind} -o {outdir}/udance -t {params.size} -j {input.j} \
                -m {params.method} -T {resources.cpus} -l {params.sub} -f {params.frag} -e {params.edg} \
                --minplacements {params.mps}
            fi
            python prune_similar.py -T {resources.cpus} -o {outdir}/udance -S {params.pra}
            if [  -f {outdir}/udance/dedupe_map.txt ]; then 
                cat {outdir}/udance/dedupe_map.txt > {outdir}/dedupe_map.txt 
            fi 
            if [ -f {outdir}/rm_map.txt ]; then 
                cat {outdir}/rm_map.txt >> {outdir}/dedupe_map.txt 
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
          thrd=config["infer_config"]["numthread"],
          t=config["infer_config"]["method"]
    benchmark: "%s/{stage}/{cluster}/{gene}/benchmark.txt" % outdir
    shell:
        '''
            # many instances of this rule may run simultaneously. To reduce the IO overhead, we "sponge" the output
            # before appending to udance_logpath
            source uDance/mysponge.sh
            (
            bash uDance/process_a_marker.sh {input} {params.c} {params.s} {params.t} {params.thrd}
            ) 2>&1 | mysponge #>> {udance_logpath} 
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
            occup=config["refine_config"]["occupancy"],
            ol=config["refine_config"]["outlier_sizelimit"],
            od=config["refine_config"]["outlier_difference"]
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/refine.{cluster}.txt" % outdir

    shell:
        """
            (
            python run_udance.py refine -p {outdir}/udance/{wildcards.cluster} -m {params.method} -M {resources.mem_mb} -c {params.c} -o {params.occup} -T {resources.cpus} -l {params.ol} -d {params.od}
            ) >> {udance_logpath} 2>&1
        """

rule blinference:
    input: expand("%s/udance/{{cluster}}/astral_output.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    output: expand("%s/udance/{{cluster}}/astral_output.{approach}.nwk.bl" % outdir, approach=["incremental", "updates"])
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/blinference.{cluster}.txt" % outdir
    shell:
        '''
            pwdd=`pwd`
            for approach in incremental updates; do
                if [ -f  {outdir}/udance/{wildcards.cluster}/skip_partition ] ; then
                    cp {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk.bl
                else
                    java -Xmx{resources.mem_mb}M -Djava.library.path=$pwdd/uDance/tools/ASTRAL/lib/ -jar $pwdd/uDance/tools/ASTRAL/astralmp.5.17.2.jar \
                        -q {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk \
                        -i {outdir}/udance/{wildcards.cluster}/astral_input.trees \
                        -o {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk.bl \
                        -C -T {resources.cpus} -u > {outdir}/udance/{wildcards.cluster}/astral.$approach.log.bl 2>&1
                fi
            done
        '''

def aggregate_stitch_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{i}/species.txt"))
    if config["refine_config"]["infer_branchlen"] in [False, "False"]:
        return [f"%s/udance/%s/astral_output.%s.nwk" % (outdir, i, j) for i in wc.i for j in ["incremental", "updates"]]
    else:
        return [f"%s/udance/%s/astral_output.%s.nwk.bl" % (outdir, i, j) for i in wc.i for j in ["incremental", "updates"]]


rule stitch:
    input: aggregate_stitch_input
    output: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: b = config["refine_config"]["infer_branchlen"]
    benchmark: "%s/benchmarks/stitch.txt" % outdir
    shell:
        """
           (
            if [[ "{params.b}" == "False" ]] ; then
                python run_udance.py stitch -o {outdir}/udance
            else
                python run_udance.py stitch -o {outdir}/udance -b
            fi
            cp {outdir}/udance/udance.*.nwk {outdir}
           ) >> {udance_logpath} 2>&1
        """
