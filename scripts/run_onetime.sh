#!/usr/bin/env bash
# $1 backbone tree
# $2 sequences file
# $3 out dir
# $4 run time
# $5 threads
# $6 subset size
# $7 jplace

export SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
bash $SCRIPTS_DIR/dedup.sh ${1} ${2} ${3}/${4}
# nw_reroot -d sometimes fails. So we first pick a random leaf, then reroot at that leaf and finally deroot
nw_reroot ${3}/${4}/backbone.tree `nw_labels -I  ${3}/${4}/backbone.tree | head -n 1 ` | nw_reroot -d - > ${3}/${4}/backbone_deroot.tree
#nw_reroot -d ${3}/${4}/backbone.tree > ${3}/${4}/backbone_deroot.tree
if [ -s ${3}/${4}/backbone_deroot.tree ]; then
        # deroot not fail
	:
else
        # deroot fail
	cp ${3}/${4}/backbone.tree ${3}/${4}/backbone_deroot.tree
fi
if [ -z "$7" ]; then
    # apples commit 3bc50da
	~/apples/run_apples.py -m FM -s ${3}/${4}/backbone.fa -q ${3}/${4}/query.fa -t ${3}/${4}/backbone_deroot.tree -o ${3}/${4}/placement.jplace -V 0.05 -f 0 -b 6 -T ${5} > ${3}/${4}/placement.output
else
	cp ${7} ${3}/${4}/placement.jplace
fi
python $SCRIPTS_DIR/didactic/run_didactic.py decompose -o ${3}/${4} -j ${3}/${4}/placement.jplace -s ${3}/${4}/seqs/ -e 0.02 -f ${6} -T ${5} -m "iqtree"
# we assume 8 cores will be given to each raxml job. PAR denotes number of raxml jobs that should
# will run in parallel
PAR=`python -c "print(max($5//4,1))"`
cat ${3}/${4}/main_raxml_script_0.sh | xargs -n1 -P${PAR} -I% bash -c "%"
python $SCRIPTS_DIR/waiting.py --out-dir ${3}/${4}
python $SCRIPTS_DIR/didactic/run_didactic.py refine -T ${5} -o ${3}/${4} -m "iqtree"
bash ${3}/${4}/main_astral_script.sh 
python $SCRIPTS_DIR/didactic/run_didactic.py stitch -o ${3}/${4}
