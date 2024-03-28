#!/usr/bin/env bash
export SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# create a conda environment including the tools benchmarked.
if conda info --envs | grep "udance" > /dev/null; then
        echo "conda environment udance exists"
else
	basecomm=$(echo conda create -y -c bioconda -c conda-forge --channel smirarab --name udance python=3.9  pip newick_utils=1.6 setuptools seqkit=2.1.0 scipy dendropy=4.5.2 pandas=1.3.0 snakemake raxml=8.2.12 iqtree=2.1.2 treeshrink=1.3.9 fasttree=2.1.10 julia=1.7.1 gappa=0.7.1 trimal=1.4.1 raxml-ng) 
	condaplat=$(conda info | grep "platform" | awk '{print $3}')
	if [ "$condaplat" == "osx-arm64" ] ; then
		echo "Installing conda environment on Mac OS X Apple Chip platform"
		echo "tqdist will not be installed. Installing tqdist package is optional for uDance."
		CONDA_SUBDIR=osx-64 $basecomm
	elif [ "$condaplat" == "osx-64" ] ; then
		echo "Installing conda environment on Mac OS X platform"
                echo "tqdist will not be installed. Installing tqdist package is optional for uDance."
		$basecomm
	elif [ "$condaplat" == "linux-64" ] ; then
		echo "Installing conda environment on Linux 64-bit platform"
                $basecomm tqdist
	else
		echo "Your platform is not supported. The following conda platforms are supported: osx-arm64, osx-64, and linux-64."
		exit 1
	fi
	
        source activate udance
        conda activate udance
        pip install apples==2.0.11 kmeans1d==0.3.1
	# install raxml-ng from the github release (conda conflict)
	echo "Warning. if you want to use raxml-ng, please install it manually and make the executable "raxml-ng" available in the environment."
        # julia needs to download packages during the first run. Let's do it now while we are online :)
        julia $SCRIPTS_DIR/uDance/correction_multi.jl datasmall/alignments/p0309.fasta > /dev/null
        # install astral
        pushd $SCRIPTS_DIR/uDance/tools/ASTRAL
       	./make.sh > .astralInstallation.log 2>&1
        java -D"java.library.path=lib/" -jar native_library_tester.jar
        popd
fi
