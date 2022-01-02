#!/usr/bin/env bash
set -euo pipefail
# $1 concat alignment
# $2 bbone
# $3 char
# $4 number of threads
# $5 all alignments dir

export ALN=$1
BBONE=$2
export CHARTYPE=$3
NUMTHR=$4
ALLALNS=$5

export SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

export MNTMP=$(mktemp -dt filterXXXXXX)

#echo $MNTMP

export OMP_NUM_THREADS=1
# compute backbone tree
if [ "$CHARTYPE" == "nuc" ]; then
  fasttree -nt -nosupport -nopr -nome -noml -intree $BBONE <$ALN >$MNTMP/backbone_me.tree
  build_applesdtb.py -T $NUMTHR -f 0.2 -s $ALN -t $MNTMP/backbone_me.tree -D -o $MNTMP/apples.dtb
else
  fasttree -nosupport -nopr -nome -noml -intree $BBONE <$ALN >$MNTMP/backbone_me.tree
  build_applesdtb.py -p -T $NUMTHR -f 0.2 -s $ALN -t $MNTMP/backbone_me.tree -D -o $MNTMP/apples.dtb
fi
# create apples database -D

onequery() {
  # $1 query name
  QUERY=$1

  TMP=$(mktemp -dt onequeryXXXXXX)
  nw_prune $MNTMP/backbone_me.tree $QUERY >$TMP/backbone.nwk
  seqkit grep -f <(echo $QUERY) -w 0 --quiet $ALN >$TMP/query.fa
  if [ "$CHARTYPE" == "nuc" ]; then
    run_apples.py -a $MNTMP/apples.dtb -t $TMP/backbone.nwk -q $TMP/query.fa -f 0.2 -b 25 -D -o $TMP/apples.jplace -T 1
  else
    run_apples.py -p -a $MNTMP/apples.dtb -t $TMP/backbone.nwk -q $TMP/query.fa -f 0.2 -b 25 -D -o $TMP/apples.jplace -T 1
  fi
  gappa examine graft --jplace-path=$TMP/apples.jplace --out-dir=$TMP >/dev/null 2>/dev/null
  n1=$($SCRIPTS_DIR/tools/compareTrees.missingBranch $MNTMP/backbone_me.tree $TMP/apples.newick | awk '{printf $2}')
  printf "$QUERY\t$n1\n"
  rm -r $TMP
}

export -f onequery

nw_labels -I $BBONE | xargs -n1 -P$NUMTHR -I% bash -c "onequery %" | sort -k2n >$MNTMP/RF.tsv

NSPCS=$(nw_labels -I $BBONE | wc -l)
THRESH=$(python -c "import math; print(math.floor(math.log2($NSPCS)))")
awk -v thr="$THRESH" '$2 > thr' $MNTMP/RF.tsv | cut -f1 >$MNTMP/removedfirststage.tsv

seqkit grep -vf $MNTMP/removedfirststage.tsv -w 0 --quiet $ALN >$MNTMP/backbone_secondstage.fa

nw_prune $MNTMP/backbone_me.tree $(cat $MNTMP/removedfirststage.tsv) >$MNTMP/backbone_secondstage.tree

if [ "$CHARTYPE" == "nuc" ]; then
  fasttree -nt -nosupport -nopr -nome -noml -intree $MNTMP/backbone_secondstage.tree <$MNTMP/backbone_secondstage.fa >$MNTMP/backbone_secondstage_reestimated.tree
  build_applesdtb.py -T $NUMTHR -f 0.2 -s $MNTMP/backbone_secondstage.fa -t $MNTMP/backbone_secondstage_reestimated.tree -D -o $MNTMP/apples_secondstage.dtb
else
  fasttree -nosupport -nopr -nome -noml -intree $MNTMP/backbone_secondstage.tree <$MNTMP/backbone_secondstage.fa >$MNTMP/backbone_secondstage_reestimated.tree
  build_applesdtb.py -p -T $NUMTHR -f 0.2 -s $MNTMP/backbone_secondstage.fa -t $MNTMP/backbone_secondstage_reestimated.tree -D -o $MNTMP/apples_secondstage.dtb
fi

while read sp; do
  seqkit grep -f <(echo $sp) -w 0 --quiet $ALN >$MNTMP/query_secondstage.fa
  if [ "$CHARTYPE" == "nuc" ]; then
    run_apples.py -a $MNTMP/apples_secondstage.dtb -q $MNTMP/query_secondstage.fa -f 0.2 -b 25 -o $MNTMP/apples.jplace -T 1
  else
    run_apples.py -p -a $MNTMP/apples_secondstage.dtb -q $MNTMP/query_secondstage.fa -f 0.2 -b 25 -o $MNTMP/apples.jplace -T 1
  fi
  gappa examine graft --jplace-path=$MNTMP/apples.jplace --out-dir=$MNTMP --allow-file-overwriting >/dev/null 2>/dev/null
  n1=$($SCRIPTS_DIR/tools/compareTrees.missingBranch $MNTMP/backbone_me.tree $MNTMP/apples.newick -simplify | awk '{printf $2}')
  printf "$sp\t$n1\n"
done < $MNTMP/removedfirststage.tsv > $MNTMP/RF2.tsv

NSPCS=$(nw_labels -I $MNTMP/backbone_secondstage.tree | wc -l)
THRESH=$(python -c "import math; print(math.floor(math.log2($NSPCS)))")
awk -v thr="$THRESH" '$2 > thr' $MNTMP/RF2.tsv | cut -f1 >$MNTMP/removedsecondstage.tsv

nw_prune $MNTMP/backbone_me.tree `cat $MNTMP/removedsecondstage.tsv` > $MNTMP/backbone_thirdstage.tree
TreeCluster.py -i $MNTMP/backbone_thirdstage.tree -m max -t 0.7 > $MNTMP/clusters.txt

python -c "from uDance.occupancy_outliers import occupancy_outliers; \
           occupancy_outliers(\"$ALLALNS\", \"$MNTMP/clusters.txt\", \"$CHARTYPE\"=='prot')" >$MNTMP/removedthirdstage.tsv
cat $MNTMP/removedsecondstage.tsv $MNTMP/removedthirdstage.tsv

rm -r $MNTMP
