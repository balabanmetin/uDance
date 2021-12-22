#!/bin/bash
set -e

# inputs
# $1 alignment path
# $2 protein or nucleotide
# $3 number of starting trees
# $4 method raxml-8, or raxml-ng, TODO or iqtree
# TODO, $5numthreads. currently fixed to 1

export ALN=$(basename $1)
export CHARTYPE=$2
export STARTS=$3
export ITOOL=$4

SHMT=`mktemp -dt processmarkerXXXXXX`
cp $1 $SHMT/$ALN
pushd $SHMT > /dev/null

#pushd $(dirname $1) || (echo "cd to alignment directory failed"; exit 1)



#export PROJ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." &> /dev/null && pwd )"
#
#pushd $PROJ_DIR > /dev/null
#
#mkdir -p $i/genes_processed/$g
#
#/usr/bin/time -o $i/genes_processed/$g/time_UPPFFN.txt -f "%e,%M" run_upp.py -s $i/sequences/ffn/$g.fasta -B 100000 -M -1 -T 0.66 -m dna -d $i/genes_processed/$g -x $thr -o UPPFFN > $i/genes_processed/$g/upp_ffn.out 2> $i/genes_processed/$g/upp_ffn.err
#
#temp=`mktemp -t XXXXXX.fa`
#sitemask=0.95
#seqmask=0.80
#
#echo Masking sites with more than 95% gaps ...
#10kBacGenomes/mask_sites.sh $i/genes_processed/$g/UPPFFN_alignment_masked.fasta $temp $sitemask
#
##echo Masking sequences with more than 66% gaps ...
#echo Masking sequences with more than 80% gaps ...
#10kBacGenomes/mask_sequences.sh $temp $i/genes_processed/$g/UPPFFN_alignment_masked_filtered.fasta $seqmask
#rm $temp


echo "Running Fasttree first time"
export OMP_NUM_THREADS=1

if [[ "$CHARTYPE" == "nuc" ]] ; then
  export FASTMODEL="-nt -gtr"
  export RAXMODEL="GTRCAT"
  export NGIQMODEL="GTR"
else
  export FASTMODEL="-lg"
  export RAXMODEL="PROTCATLG"
  export NGIQMODEL="LG"
fi

fasttree -nopr $FASTMODEL -gamma -seed 12345 -log fasttree.log < $ALN  > fasttree.nwk 2> fasttree.err

echo "Treeshrink'ing"

numsp=`nw_labels -I fasttree.nwk | wc -l`
if [[ "$numsp" -gt 20 ]] ; then
  run_treeshrink.py -t fasttree.nwk -a $ALN -f -o fasttree_shrunk > treeshrink.log 2>&1
  nw_labels -I fasttree_shrunk/output.nwk > remaining_after_shrunk.txt
else
  nw_labels -I fasttree.nwk > remaining_after_shrunk.txt
fi

seqkit grep -f remaining_after_shrunk.txt $ALN -w 0 --quiet -o shrunk.fasta


run_a_start(){
  TREEID=$1
  mkdir -p $TREEID
  pushd $TREEID > /dev/null
    ln -s ../shrunk.fasta
    if [[ "$ITOOL" == "raxml-ng" ]] ; then
      raxml-ng --msa shrunk.fasta --tree pars{1} --model ${NGIQMODEL}+G --threads 1 --seed $TREEID --prefix RUN 2> raxml.err > raxml.log
      nw_topology -bI RUN.raxml.bestTree > RAxML_result.RUN
    else
      fasttree -nopr $FASTMODEL -gamma -seed $TREEID -log fasttree_r2.log < ../shrunk.fasta  > fasttree_r2.nwk 2> fasttree_r2.err
      python -c "import treeswift as ts; t=ts.read_tree_newick(\"fasttree_r2.nwk\"); \
            [c.resolve_polytomies() for c in t.root.children]; print(t)" > fasttree_r2_resolved.nwk
      raxmlHPC -T 1 -m ${RAXMODEL} -F -f D -D -s shrunk.fasta -p $TREEID -n RUN -t fasttree_r2_resolved.nwk 2> raxml.err > raxml.log
    fi
    #if raxmlHPC -T 1 -m ${RAXMODEL}GAMMA -f e -s ../shrunk.fasta -t RAxML_result.RUN -n RUNGAMMA -p 12345 2> raxml_gamma.err > raxml_gamma.log ;
    iqtree -ntmax 1 -abayes -fast -m ${NGIQMODEL}+G -s shrunk.fasta -t RAxML_result.RUN -seed $TREEID > iqtree.out 2> iqtree.err
  popd > /dev/null
}

export -f run_a_start

for i in `seq 1 $STARTS`; do
  run_a_start $i
done

grep -H "BEST SCORE" */shrunk.fasta.log | sort -k5n | tail -n 1 | cut -f1 -d ":"| sed "s/.log/.treefile/g" > bestTreename.txt
cp `cat bestTreename.txt` bestTree_slash.nwk

sed "s/)\//)/g" bestTree_slash.nwk > bestTree.nwk
rm bestTree_slash.nwk

popd > /dev/null

rsync -a "$SHMT"/ $(dirname $1)

rm -rf $SHMT
