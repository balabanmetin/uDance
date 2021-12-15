#!/usr/bin/env bash

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


TMPDIR=`mktemp -d placementruleXXXXXX`

cat $ALNDIR/* | grep ">" | sed "s/>//g" | sort -u > $TMPDIR/alltaxa.txt
touch $TMPDIR/alnpaths.txt
ls $ALNDIR | while read aln; do
    # awk instead of wc -L because OSX doesn't have wc -L
    sze=`sed -e "s/>\(.*\)/@>\1@/g" $ALNDIR/$aln|tr -d "\n"|tr "@" "\n"|tail -n+2 | head -n 2 | tail -n 1 | awk '{print length}'`
    grep ">" $ALNDIR/$aln | sed "s/>//g" | sort -u > $TMPDIR/thistaxa.txt
    comm -23 $TMPDIR/alltaxa.txt $TMPDIR/thistaxa.txt > $TMPDIR/missing.txt
    (cat $ALNDIR/$aln; while read tx; do printf ">$tx\n"; dd if=/dev/zero bs=$sze count=1 2>/dev/null | tr '\0' "-" | tr '\0' '-'; printf "\n"; done < $TMPDIR/missing.txt) | gzip -c - > $TMPDIR/$aln.gz
    echo $TMPDIR/$aln.gz >> $TMPDIR/alnpaths.txt
done
# concatenate all alignments
seqkit concat --quiet -w 0 $(cat $TMPDIR/alnpaths.txt) | gzip -c - > $TMPDIR/concat.fa.gz
# order it so that bacbone sequences are the first N records
(seqkit grep -f <(nw_labels -I $BBONE) -w 0 --quiet $TMPDIR/concat.fa.gz; seqkit grep -v -f <(nw_labels -I $BBONE) -w 0 --quiet $TMPDIR/concat.fa.gz) | gzip -c - >  $TMPDIR/concat_ordered.fa.gz
# Duplicate map is written to OUTDIR.
# TODO should we replace N's with dashes?
gzip -cd $TMPDIR/concat_ordered.fa.gz | seqkit rmdup --quiet -s -i -w 0 -D $OUTDIR/rm_map.txt | gzip -c - > $TMPDIR/concat_dedupe.fa.gz
rm $TMPDIR/concat_ordered.fa.gz
# extract backbone and query alignments from deduplicated alignment
seqkit grep -f <(nw_labels -I $BBONE) $TMPDIR/concat_dedupe.fa.gz -w 0 --quiet -o $OUTDIR/placement/backbone.fa
seqkit grep -v -f <(nw_labels -I $BBONE) $TMPDIR/concat_dedupe.fa.gz -w 0 --quiet -o $OUTDIR/placement/query.fa
rm $TMPDIR/concat_dedupe.fa.gz

grep ">" $OUTDIR/placement/backbone.fa | sed "s/>//g" > $TMPDIR/backbone_id_dedup.txt
mapfile -t < $TMPDIR/backbone_id_dedup.txt

# backbone tree with duplicates removed. This tree might potentially have polytomies.
# Polytomies need to be resolved as its required for APPLES-2 (FastTree).
# TODO resolution using raxml
nw_prune -v $BBONE "${MAPFILE[@]}" > $TMPDIR/backbone.tree
python -c "import treeswift as ts; t=ts.read_tree_newick(\"$TMPDIR/backbone.tree\"); \
          [c.resolve_polytomies() for c in t.root.children]; print(t)" > $OUTDIR/placement/backbone.tree

# $1 concat alignment
# $2 bbone
# $3 char
# $4 number of threads
# $5 all alignments dir
bash uDance/filter_backbone.sh $OUTDIR/placement/backbone.fa $OUTDIR/placement/backbone.tree \
      $CHARTYPE $NUMTHREADS $ALNDIR 2> $OUTDIR/placement/filtering.log > $OUTDIR/placement/filtered.txt

NUMFILT=`wc -l $OUTDIR/placement/filtered.txt`

if [[ "$NUMFILT" -gt 0 ]]; then
  seqkit grep -vf $OUTDIR/placement/filtered.txt $OUTDIR/placement/backbone.fa -w 0 --quiet -o $TMPDIR/backbone.fa
  mv $TMPDIR/backbone.fa $OUTDIR/placement/backbone.fa

  seqkit grep -f $OUTDIR/placement/filtered.txt $OUTDIR/placement/backbone.fa -w 0 --quiet -o $TMPDIR/query.fa
  cat $TMPDIR/query.fa >> $OUTDIR/placement/query.fa

  nw_prune $OUTDIR/placement/backbone.tree `cat $OUTDIR/placement/filtered.txt` > $TMPDIR/backbone.tree
  mv $TMPDIR/backbone.tree $OUTDIR/placement/backbone.tree
fi

rm -r $TMPDIR
