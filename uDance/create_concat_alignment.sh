#!/usr/bin/env bash
set -euo pipefail
# Yueyu Jiang and Metin Balaban

# $1 alignment dir
# $2 bbone
# $3 output dir
# $4 chartype
# $5 num threads
ALNDIR=$1
BBONE=$2
OUTDIR=$3
CHARTYPE=$4
NUMTHREADS=$5
APF=$6
APM=$7
APB=$8
APV=$9


TDR=`mktemp -d $OUTDIR/createconcat`

cat $ALNDIR/* | grep ">" | sed "s/>//g" | sort -u > $TDR/alltaxa.txt
touch $TDR/alnpaths.txt
ls $ALNDIR | while read aln; do
    # awk instead of wc -L because OSX doesn't have wc -L
    sze=$(sed -e "s/>\(.*\)/@>\1@/g" $ALNDIR/$aln |tr -d "\n"|tr "@" "\n"| sed -n '3p' | awk '{print length}')
    grep ">" $ALNDIR/$aln | sed "s/>//g" | sort -u > $TDR/thistaxa.txt
    comm -23 $TDR/alltaxa.txt $TDR/thistaxa.txt > $TDR/missing.txt
    (cat $ALNDIR/$aln; while read tx; do printf ">$tx\n"; dd if=/dev/zero bs=$sze count=1 2>/dev/null | tr '\0' "-" | tr '\0' '-'; printf "\n"; done < $TDR/missing.txt) | gzip -c - > $TDR/$aln.gz
    echo $TDR/$aln.gz >> $TDR/alnpaths.txt
done
# concatenate all alignments
seqkit concat --quiet -w 0 $(cat $TDR/alnpaths.txt) | gzip -c - > $TDR/concat.fa.gz
# order it so that bacbone sequences are the first N records
(seqkit grep -f <(nw_labels -I $BBONE) -w 0 --quiet $TDR/concat.fa.gz; seqkit grep -v -f <(nw_labels -I $BBONE) -w 0 --quiet $TDR/concat.fa.gz) | gzip -c - >  $TDR/concat_ordered.fa.gz
# Duplicate map is written to OUTDIR.
# TODO should we replace N's with dashes?
gzip -cd $TDR/concat_ordered.fa.gz | seqkit rmdup --quiet -s -i -w 0 -D $OUTDIR/rm_map.txt | gzip -c - > $TDR/concat_dedupe.fa.gz
rm $TDR/concat_ordered.fa.gz
# extract backbone and query alignments from deduplicated alignment
seqkit grep -f <(nw_labels -I $BBONE) $TDR/concat_dedupe.fa.gz -w 0 --quiet -o $OUTDIR/placement/backbone.fa
seqkit grep -v -f <(nw_labels -I $BBONE) $TDR/concat_dedupe.fa.gz -w 0 --quiet -o $OUTDIR/placement/query.fa
rm $TDR/concat_dedupe.fa.gz

grep ">" $OUTDIR/placement/backbone.fa | sed "s/>//g" > $TDR/backbone_id_dedup.txt
mapfile -t < $TDR/backbone_id_dedup.txt

# backbone tree with duplicates removed. This tree might potentially have polytomies.
# Polytomies need to be resolved as its required for APPLES-2 (FastTree).
# TODO resolution using raxml
nw_prune -v <(nw_topology -bI $BBONE) "${MAPFILE[@]}" > $TDR/backbone.tree
python -c "import treeswift as ts; t=ts.read_tree_newick(\"$TDR/backbone.tree\"); \
          [c.resolve_polytomies() for c in t.root.children]; print(t)" > $OUTDIR/placement/backbone.tree

# $1 concat alignment
# $2 bbone
# $3 char
# $4 number of threads
# $5 all alignments dir
bash uDance/filter_backbone.sh $OUTDIR/placement/backbone.fa $OUTDIR/placement/backbone.tree \
      $CHARTYPE $NUMTHREADS $ALNDIR $APF $APM $APB $APV $TDR 2> $OUTDIR/placement/filtering.log > $OUTDIR/placement/filtered.txt

# useless cat
NUMFILT=`cat $OUTDIR/placement/filtered.txt | wc -l`

printf "%d low quality sequences are removed from the backbone. The sequences will not be present in the final tree.\n" $NUMFILT >&2

if [[ "$NUMFILT" -gt 0 ]]; then
  seqkit grep -vf $OUTDIR/placement/filtered.txt $OUTDIR/placement/backbone.fa -w 0 --quiet -o $TDR/backbone.fa
  cp $OUTDIR/placement/backbone.fa $OUTDIR/placement/backbone.fa.bak
  # seqkit grep -f $OUTDIR/placement/filtered.txt $OUTDIR/placement/backbone.fa -w 0 --quiet -o $TDR/query.fa

  mv $TDR/backbone.fa $OUTDIR/placement/backbone.fa
  #cat $TDR/query.fa >> $OUTDIR/placement/query.fa
  #cat $TDR/query.fa > $OUTDIR/placement/query.fa

  nw_prune $OUTDIR/placement/backbone.tree `cat $OUTDIR/placement/filtered.txt` > $TDR/backbone.tree
  mv $TDR/backbone.tree $OUTDIR/placement/backbone.tree
fi

#rm -r $TDR
