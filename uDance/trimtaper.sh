#!/usr/bin/env bash
set -euo pipefail

INP=$1
THR=$2
OUT=$3

TFA=`mktemp -dt trimtaperXXXXXX`
trimal -in $INP -out $TFA/shr.fa -gt $THR
cat $TFA/shr.fa | seqkit rmdup -s -o $TFA/seqs.fa -i -w 0 -D $TFA/rm_map.txt
julia uDance/correction_multi.jl $TFA/seqs.fa > $TFA/tapered.fa

if [ -f "$TFA/rm_map.txt" ]; then
    cut -f2 $TFA/rm_map.txt | sed "s/, /,/g" > $TFA/rm_con.txt
    python scripts/expand_dedupe.py $TFA/tapered.fa $TFA/rm_con.txt $OUT
else
    mv $TFA/tapered.fa $OUT
fi

rm -r $TFA
