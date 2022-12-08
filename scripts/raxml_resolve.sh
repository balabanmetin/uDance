#!/usr/bin/env bash

# $1 backbone_tree
# $2 seqs
# $3 output newick file

TMP=`mktemp -d -t XXXXXX`
echo $TMP
export SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

nw_labels -I $1 > $TMP/bb.txt
seqkit grep -f $TMP/bb.txt -w 0 $2 > $TMP/bb.fna
seqkit rmdup -o $TMP/bb_dedup.fna -i -w 0 -D $TMP/rm_map.txt -s $TMP/bb.fna
grep ">" $TMP/bb_dedup.fna  | sed "s/>//g" > $TMP/bb_dedup.txt
nw_prune -v $1 `cat $TMP/bb_dedup.txt` > $TMP/bb_dedup_constraint.nwk
raxmlHPC-PTHREADS -T 16 -s $TMP/bb_dedup.fna -g $TMP/bb_dedup_constraint.nwk -n file -p 12345 -m GTRCAT -w `realpath $TMP` > $TMP/raxml.out 2> $TMP/raxml.err
sed -i "s/[0-9]\+\t//g" $TMP/rm_map.txt
python $SCRIPTS_DIR/add_back.py --rm-duplicate $TMP/rm_map.txt --tree-file $TMP/RAxML_bestTree.file --out-dir $TMP
cp $TMP/didactic_full.nwk ${3} || ( echo "cannot write output; you can find the output tree at $TMP/didactic_full.nwk. exiting"; exit )

rm -r $TMP

