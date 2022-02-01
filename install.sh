#!/usr/bin/env bash

export SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# create a conda environment including the tools benchmarked.
if conda info --envs | grep "udance" > /dev/null; then
        echo "conda environment udance exists"
else
        conda create -y -c bioconda -c conda-forge --channel smirarab --name udance python=3.8 \
        pip newick_utils=1.6 setuptools seqkit=2.1.0 scipy dendropy=4.5.2 pandas=1.3.0 snakemake \
        raxml=8.2.12 tqdist iqtree=2.1.2 treeshrink=1.3.9 fasttree=2.1.10 julia=1.7.1 gappa=0.7.1
        source activate udance
        conda activate udance
        pip install apples==2.0.7 kmeans1d==0.3.1

        # julia needs to download packages during the first run. Let's do it now while we are online :)
        julia $SCRIPTS_DIR/uDance/correction_multi.jl datasmall/alignments/p0309.fasta > /dev/null
fi
