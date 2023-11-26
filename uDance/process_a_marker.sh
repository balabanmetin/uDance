#!/bin/bash
set -euo pipefail

# inputs
# $1 alignment path
# $2 protein or nucleotide
# $3 number of starting trees
# $4 method raxml-8, or raxml-ng, TODO or iqtree
# $5 numthreads

export ALN=$(basename $1)
export CHARTYPE=$2
export STARTS=$3
export ITOOL=$4
export NUMTHREADS=$5

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
  export NGMODEL="GTR+G"
  export IQMODEL="GTR+FO+G4"
else
  export FASTMODEL="-lg"
  export RAXMODEL="PROTCATLG"
  export NGMODEL="LG+G"
  export IQMODEL="LG+FO+G4"
fi

fasttree -nopr $FASTMODEL -gamma -seed 12345 -log fasttree.log < $ALN  > fasttree.nwk 2> fasttree.err


numsp=$(nw_labels -I fasttree.nwk | wc -l)
if [[ "$numsp" -gt 20 ]] ; then
  run_treeshrink.py -t fasttree.nwk -a $ALN -f -o fasttree_shrunk > treeshrink.log 2>&1
  nw_labels -I fasttree_shrunk/output.nwk > remaining_after_shrunk.txt
else
  nw_labels -I fasttree.nwk > remaining_after_shrunk.txt
fi

shrinkbefore=$(nw_labels -I fasttree.nwk | wc -l)
shrinkafter=$(cat remaining_after_shrunk.txt | wc -l)

printf "TreeShrink removed %d sequences from $ALN \n" $(python -c "print($shrinkbefore - $shrinkafter)")

seqkit grep -f remaining_after_shrunk.txt $ALN -w 0 --quiet -o shrunk.fasta


run_a_start(){
  TREEID=$1
  mkdir -p $TREEID
  pushd $TREEID > /dev/null
    ln -s ../shrunk.fasta
    if [[ "$ITOOL" == "raxml-ng" ]] ; then
      if raxml-ng --msa shrunk.fasta --tree pars{1} --model ${NGMODEL} --threads 1 --seed $TREEID --prefix RUN --lh-epsilon 0.5 > raxml-ng.log 2>&1; then
        nw_topology -bI RUN.raxml.bestTree | nw_order - > RAxML_result.RUN
        grep "Final" RUN.raxml.log | cut -f3 -d ' ' > likelihood.txt
        iqtree -T 1 -abayes -m ${IQMODEL} -s shrunk.fasta -te RAxML_result.RUN -seed $TREEID --redo > iqtree.out 2> iqtree.err
      else
        grep "ERROR" raxml-ng.log >&2
        echo "WARNING: RAxML-NG failed to infer a tree using $ALN. Continuing using RAxML-8" >&2
        export ITOOL="raxml-8"
        # then run raxml
      fi
    fi
    if [[ "$ITOOL" == "raxml-8" ]] ; then
      # start with iqtree -fast tree
      iqtree -T 1 -fast -m ${IQMODEL} -s shrunk.fasta -seed $TREEID --prefix START> iqtree_start.out 2> iqtree_start.err
      #fasttree -nopr $FASTMODEL -gamma -seed $TREEID -log fasttree_r2.log < ../shrunk.fasta  > fasttree_r2.nwk 2> fasttree_r2.err
      #python -c "import treeswift as ts; t=ts.read_tree_newick(\"fasttree_r2.nwk\"); \
      #      [c.resolve_polytomies() for c in t.root.children]; print(t)" > fasttree_r2_resolved.nwk
      # TODO raxml8 multithreading
      raxmlHPC -T 1 -m ${RAXMODEL} -F -f D -D -s shrunk.fasta -p $TREEID -n RUN -t START.treefile > raxml-8.log 2>&1
      iqtree -T 1 -abayes -m ${IQMODEL} -s shrunk.fasta -te RAxML_result.RUN -seed $TREEID --redo > iqtree.out 2> iqtree.err
      grep "BEST SCORE" shrunk.fasta.log | cut -f5 -d ' ' > likelihood.txt
    elif [[ "$ITOOL" == "iqtree" ]] ; then
      iqtree -T 1 -abayes -m ${IQMODEL} -s shrunk.fasta -seed $TREEID --redo > iqtree.out 2> iqtree.err
      grep "BEST SCORE" shrunk.fasta.log | cut -f5 -d ' ' > likelihood.txt
    else
      echo "ERROR: the provided tool $ITOOL is none of the following options: raxml-8, raxml-ng, or iqtree" >&2
      return 1
    fi
    #if raxmlHPC -T 1 -m ${RAXMODEL}GAMMA -f e -s ../shrunk.fasta -t RAxML_result.RUN -n RUNGAMMA -p 12345 2> raxml_gamma.err > raxml_gamma.log ;
  popd > /dev/null
}

export -f run_a_start

seq 1 $STARTS | xargs -n1 -P${NUMTHREADS} -I% bash -c "run_a_start %"

#for i in `seq 1 $STARTS`; do
#  run_a_start $i
#done

for i in `seq 1 $STARTS`; do
  cat $i/likelihood.txt | sed "s/$/ $i/g"
done | sort -k1n | tail -n 1 | cut -f2 -d ' '| sed "s/$/\/shrunk.fasta.treefile/g" > bestTreename.txt

cp `cat bestTreename.txt` bestTree_slash.nwk

sed "s/)\//)/g" bestTree_slash.nwk > bestTree.nwk
rm bestTree_slash.nwk

popd > /dev/null

rsync -a "$SHMT"/ $(dirname $1)

rm -rf $SHMT
