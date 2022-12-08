#!/usr/bin/env bash
# $1 out_dir
# $2 rm_map

export SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

sed -i "s/[0-9]\+\t//g" ${2}
python $SCRIPTS_DIR/add_back.py --rm-duplicate ${2} --tree-file ${1}/didactic.nwk --out-dir ${1}/../
