#!/bin/bash
#SBATCH --job-name="16S"
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH -t 144:00:00
#SBATCH --ntasks-per-node=64

source activate sarscov2monitor
#time bash run.sh -b treefile/resolved/RAxML_result.2021-06-19.resolved.treefile -s alignments_2021-07-28_17-22-44-all/passQC/2021-07-28_17-22-44-all_passQC_refs_hist.trimmed.aln -l 50 -o test_pipeline2 -t 128
#time bash run.sh -b treefile/resolved/RAxML_result.2021-06-19.resolved.treefile -s seqs/passQC/2021-07-13-all_passQC_refs_hist.trimmed.aln -l 10 -o test_pipeline_new -t 128
#time bash run.sh -b treefile/2021-06-19_04-41-31_passQC_refs_hist.trimmed.aln.rooted.treefile -s alignments_2021-07-28_17-22-44-all/passQC/2021-07-28_17-22-44-all_passQC_refs_hist.trimmed.aln -l 10 -o test_pipeline_random_resolve -t 128
#time bash run.sh -b gg2_protocols/backbone_16S_nochim_second_stage.tree -s gg2_protocols/single_copy_extended.fa -l 10000000 -j gg2_protocols/05282022-placement-reformatted.jplace -o 16s_result/run8 -t 32 -c 15000
time bash run.sh -b gg2_protocols/selected_seqs_l_dedup.nwk -s gg2_protocols/selected_seqs_l_extended.fa -l 10000000  -o 16s_result/runb -t 32 -c 15000

