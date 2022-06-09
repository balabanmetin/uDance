#!/usr/bin/env bash
set -euo pipefail
# $1 concat alignment
# $2 bbone
# $3 char
# $4 number of threads
# $5 all alignments dir

export ALN=$1
export BBONE=$2
export CHARTYPE=$3
export NUMTHR=$4
export ALLALNS=$5
export APF=$6
export APM=$7
export APB=$8
export APV=$9
export MNTMP=${10}

export SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

export OMP_NUM_THREADS=1
# compute backbone tree with branch lengths
if [ "$CHARTYPE" == "nuc" ]; then
  fasttree -nt -nosupport -nopr -nome -noml -intree $BBONE <$ALN >$MNTMP/backbone_me.tree
else
  fasttree -nosupport -nopr -nome -noml -intree $BBONE <$ALN >$MNTMP/backbone_me.tree
fi

$SCRIPTS_DIR/FastLoo $ALN $BBONE $CHARTYPE > $MNTMP/RF.tsv

NSPCS=$(nw_labels -I $BBONE | wc -l)
THRESH=$(python -c "import math; print(math.floor(math.log2($NSPCS)))")
awk -v thr="$THRESH" '$2 > thr' $MNTMP/RF.tsv | cut -f1 >$MNTMP/removedfirststage.tsv
NUMRMFIRST=$(wc -l < $MNTMP/removedfirststage.tsv)

# if there is any large RF placements, remove from the tree and place them one by one
if [ "$NUMRMFIRST" -gt 0 ] ; then
  seqkit grep -vf $MNTMP/removedfirststage.tsv -w 0 --quiet $ALN >$MNTMP/backbone_secondstage.fa

  nw_prune $BBONE $(cat $MNTMP/removedfirststage.tsv) >$MNTMP/backbone_secondstage.tree

  if [ "$CHARTYPE" == "nuc" ]; then
    fasttree -nt -nosupport -nopr -nome -noml -intree $MNTMP/backbone_secondstage.tree <$MNTMP/backbone_secondstage.fa >$MNTMP/backbone_secondstage_reestimated.tree
    build_applesdtb.py -T $NUMTHR -f $APF -s $MNTMP/backbone_secondstage.fa -t $MNTMP/backbone_secondstage_reestimated.tree -D -o $MNTMP/apples_secondstage.dtb
  else
    fasttree -nosupport -nopr -nome -noml -intree $MNTMP/backbone_secondstage.tree <$MNTMP/backbone_secondstage.fa >$MNTMP/backbone_secondstage_reestimated.tree
    build_applesdtb.py -p -T $NUMTHR -f $APF -s $MNTMP/backbone_secondstage.fa -t $MNTMP/backbone_secondstage_reestimated.tree -D -o $MNTMP/apples_secondstage.dtb
  fi

# no need for parallelization using xargs because we expect only a handful of species in removedfirststage.tsv
  while read sp; do
    seqkit grep -f <(echo $sp) -w 0 --quiet $ALN >$MNTMP/query_secondstage.fa
    if [ "$CHARTYPE" == "nuc" ]; then
      run_apples.py -a $MNTMP/apples_secondstage.dtb -q $MNTMP/query_secondstage.fa -f $APF -m $APM -b $APB -V $APV -o $MNTMP/apples.jplace -T 1
    else
      run_apples.py -p -a $MNTMP/apples_secondstage.dtb -q $MNTMP/query_secondstage.fa -f $APF -m $APM -b $APB -V $APV -o $MNTMP/apples.jplace -T 1
    fi
    gappa examine graft --jplace-path=$MNTMP/apples.jplace --out-dir=$MNTMP --allow-file-overwriting >/dev/null 2>/dev/null
    n1=$($SCRIPTS_DIR/tools/compareTrees.missingBranch $BBONE $MNTMP/apples.newick -simplify | awk '{printf $2}')
    printf "$sp\t$n1\n"
  done < $MNTMP/removedfirststage.tsv > $MNTMP/RF2.tsv

  NSPCS=$(nw_labels -I $MNTMP/backbone_secondstage.tree | wc -l)
  THRESH=$(python -c "import math; print(math.floor(math.log2($NSPCS)))")
  awk -v thr="$THRESH" '$2 > thr' $MNTMP/RF2.tsv | cut -f1 >$MNTMP/removedsecondstage.tsv
  NUMRMSECOND=$(wc -l < $MNTMP/removedsecondstage.tsv)

  if [ "$NUMRMSECOND" -gt 0 ] ; then
    nw_prune $MNTMP/backbone_me.tree `cat $MNTMP/removedsecondstage.tsv` > $MNTMP/backbone_thirdstage.tree
  else
    # nothing to be removed
    cp $MNTMP/backbone_me.tree $MNTMP/backbone_thirdstage.tree
  fi

#there is no large RF placements in the first stage. Proceed with occupancy filter.
else
  cp $MNTMP/backbone_me.tree $MNTMP/backbone_thirdstage.tree
  touch $MNTMP/removedsecondstage.tsv
fi

TreeCluster.py -i $MNTMP/backbone_thirdstage.tree -m max -t 0.7 > $MNTMP/clusters.txt

python -c "from uDance.occupancy_outliers import occupancy_outliers; \
           occupancy_outliers(\"$ALLALNS\", \"$MNTMP/clusters.txt\", \"$CHARTYPE\"=='prot')" >$MNTMP/removedthirdstage.tsv
cat $MNTMP/removedsecondstage.tsv $MNTMP/removedthirdstage.tsv

