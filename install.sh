#!/usr/bin/env bash
# create a conda environment including the tools benchmarked.
if conda info --envs | grep "sarscov2monitor" > /dev/null; then
        echo "conda environment sarscov2monitor exists"
else
        conda create -y -c bioconda --channel smirarab -c conda-forge --name sarscov2monitor python=3.8 pip newick_utils=1.6 setuptools seqkit==0.16.1 scipy dendropy==4.5.2 pandas==1.3.0 iqtree
        source activate sarscov2monitor
        conda activate sarscov2monitor
        pip install git+https://github.com/balabanmetin/apples@dev/hgtremove
fi
