# uDance: Updating phylogenetic and phylogenomic trees using divide and conquer

<img src="https://github.com/balabanmetin/uDance/raw/master/logo.png" alt="uDance logo" width="250"/>

### Latest Version

uDance version 1.6.3

## What is uDance?

uDance is highly-scalable end-to-end workflow for inferring phylogenomic trees or updating existing ones. The input to uDance is a backbone tree, a set of DNA xor amino-acid multiple sequence alignments (MSAs) of backbone sequences, and new (query) sequences. Alternatively, when a backbone tree is not available, uDance can select a set of backbone species with high diversity to reconstruct a backbone tree. At a high level, uDance inserts the query sequences on the backbone tree independently and then refines the tree locally in different parts. The backbone is allowed to change based on the new information provided by the query sequences, but uDance also outputs an *incremental* tree with the backbone relationships fixed for users that require consistency between updates in their analyses. Since uDance aims for automatic analyses of large data, it has many built-in quality control and filtering strategies; it may decide a set of query sequences cannot be confidently put in the output tree (i.e., are unplacable) and some backbones need to be removed. When more sequences become available, the output from the previous iteration can be used as the input in the next iteration to incrementally grow the tree.

## Installation

### Prerequisites:

-   Linux or OS X (Intel or Apple chip)
-   Anaconda

### Installation steps:

1.  Clone the repository
2.  `cd uDance`
3.  `bash install.sh`
4.  `conda activate udance`

If you want to use raxml-ng in your workflow, please install raxml-ng manually and make sure the executable `raxml-ng` is available in your path.

## Quick start

You can run uDance on the **test dataset** with 10 gene multiple sequence alignments, 99 backbone sequences, and 146 query sequences using the following command. The command should complete in several minutes using 4 CPU cores:

`snakemake -c 4 --configfile config.yaml --snakefile udance.smk all`

The specification of the run is defined in the config file `config.yaml`. In the `config.yaml`, the working directory is specified using the field `workdir`. Input multiple sequence alignments should be located under `<workdir>/alignments`. For example, for the test dataset, the working directory is `datasmall` and the input alignments are located at `datasmall/alignments`.

If there is an input backbone tree as well, it should be locate it at `<workdir>/backbone.nwk` and the field in the config file `backbone` should be set to `"tree"`.

## Workflow overview

uDance workflow is implemented using Snakemake.

<img src="https://github.com/balabanmetin/uDance/raw/master/rules.png" alt="workflow rules" width="750"/>

## Output

uDance writes its output under the directory `<workdir>/output`. The workflow outputs three phylogenies: `<workdir>/output/udance.maxqs.nwk`, `<workdir>/output/udance.incremental.nwk`, and `<workdir>/output/udance.updates.nwk`. `incremental` tree guarantees that the backbone topology is fixed in the output tree. `maxqs` tree is the "best" one inferred by uDance and the location of backbone sequences might change after insertion of query sequences. We will shortly come back to the exact definition of these three output trees.

uDance workflow also outputs many intermediate files. These intermediate files can be useful for user for debugging as well as to supplement the downstream analysis. All paths below are given in relative to the output directory `<workdir>/output`.

1.  `trimmed` containes trimmed input MSAs.
2.  `backbone.nwk` is the backbone file. If there is an input backbone tree, it's identical to the input backbone tree. Otherwise, backbone sequences are selected using Mainlines algorithm and selected sequence IDs are written to `backbone/0/species.txt`. Gene trees for the selected sequences are found in the directory `backbone/0/<gene>`. This directory contains subdirectories named `1` to `k`, where `k` is the number of starting trees for the inference of the gene tree for this gene. `backbone/0/<gene>/<i>/shrunk.fasta.treefile` is the inferred maximum likelihood tree for the gene using the starting tree `i`. The highest likelihood tree among all `k` starts is `backbone/0/<gene>/bestTree.nwk` and relative-path to the starting tree that yielded the highest likelihood is `backbone/0/<gene>/bestTreename.txt`. All gene trees in newick format are written to `backbone/0/astral_input.trees`, one line per tree. This file is provided to ASTRAL as the input. No constraint tree is used during the estimation of the backbone tree. ASTRAL's output file (a newick tree) is `backbone/0/astral_output.updates.nwk`.
3.  `placement` contains the backbone and query alignments and the backbone tree used at the placement stage. `placement.jplace` is the jplace file output by APPLES-2.
4.  The partitions created by uDance are located under the directory `udance`. This directory contains subdirectories named `0` to `p-1`, where `p` is the number of partitions designated for the uDance run. `atasmall/output/udance/<partition>/species.txt` contains the list of backbone, query, and outgroup sequences in the partition. The organization of the partition is almost identical to the "partition" `backbone/0` given in the item 2 above but there are a few differences. `datasmall/output/udance/<partition>/astral_constraint.nwk` and `datasmall/output/udance/<partition>/raxml_constraint.nwk` are the two kinds of constraint trees used in ASTRAL stage. The former results in an ASTRAL tree (`datasmall/output/udance/<partition>/astral_output.incremental.nwk`) that retains the backbone tree topology and the latter allows topological changes aming the backbone sequences (`datasmall/output/udance/<partition>/astral_output.updates.nwk`).
5.  (continued) The spanning tree of partitions of the placement tree is available at `datasmall/output/udance/color_spanning_tree.nwk`. This can be regarded as the hierarchy or relative positions of the partitions. `datasmall/output/udance/outgroup_map.json` is a dictionary where, for each partition, we list the outgroup sequences (stored in keys `children` and `up`). These two files are used during stitching ASTRAL output trees of the partitions.

## uDance configuration settings

|            Parameter             |                                                Description                                                 |
|:--------------------:|:------------------------------------------------:|
|             chartype             |                                    Amino-acid or nucleotide characters                                     |
|             backbone             |             Three backbone tree source options. (1) de-novo, (2) user tree, and (3) user list              |
|      resources.large_memory      |                                    Large memory jobs memory limit (MB)                                     |
|         resources.cores          |                                     Large memory jobs CPU cores limit                                      |
|    trim_config.percent_nongap    |                           Sites with less non-gap fraction than below is removed                           |
|        mainlines_config.n        |                            Target number of backbone taxa in de novo inference                             |
|     mainlines_config.length      |                                       concatenation alignment length                                       |
|        backbone_filtering        |            Backbone filtering is recommended if backbone contains misplaced or noisy sequences             |
|       apples_config.method       |                                          APPLES-2 placement mode                                           |
|       apples_config.filter       |                                  APPLES-2 -f parameter (filter diameter)                                   |
|        apples_config.base        |                                APPLES-2 -b parameter (minimum observations)                                |
|      apples_config.overlap       |                                APPLES-2 minimum alignment overlap fraction                                 |
|       prep_config.edge_thr       |                                          Partition diameter limit                                          |
|     prep_config.cluster_size     |           Approximate partition size. Options are (1) auto, (2) fast, (3) user defined integer.            |
|      prep_config.sublength       |             Minimum partition alignment length. If not satisfied, the partition is discarded.              |
|      prep_config.pruneafter      |                                           Maximum partition size                                           |
|    prep_config.min_placements    |       The maximum number of placements occurred for the partition to be skipped to save running time       |
|       infer_config.method        |                        Gene tree inference method (RAxML-ng, IQTree-2, or RAxML-8)                         |
|      infer_config.numstart       |                                        The number of starting trees                                        |
|     infer_config.num_threads     |                             The number of threads used in gene tree inference                              |
|      refine_config.contract      |                            Contract low support (IQTree-2 aBayes test) branches                            |
|     refine_config.occupancy      |                              Gene occupancy threshold for sequence inclusion                               |
| refine_config.outlier_sizelimit  | The next two are 1D k-means-based (k=2) outlier gene detection parameter. Limit for outlier size fraction. |
| refine_config.outlier_difference |    Centroid difference must be larger than this value to designate the first cluster as the outlier set    |
|  refine_config.infer_branchlen   |                           Infer branch lengths in substitution unit using ASTRAL                           |

## AUTHORS:

Metin Balaban

Yueyu Jiang

## Changelog

1.6.3

-   Cluster size can be automatically set

-   Bug fix

1.6.2

-   Backbone filtering can be turned on or off

1.6.1

-   Benchmarks directives are used in udance.smk

1.6.0

-   Skip a partition (return backbone) if the number of placements is it less than a desired number

-   Bug fixes

1.5.1:

-   Output tree may have branch lengths in substitution unit if desired.

1.5.0:

-   uDance uses ASTRAL version 5.17.2, which supports both multithreading and constrained search.

1.4.1:

-   Performance improvements in placement_prep

-   Keep placement_prep intermediate files.

1.4.0:

-   Stitching algorithm outputs a new tree named maxqs which picks the best one of the two ASTRAL trees for each partition.

1.3.3:

-   Default decompose edge threshold changed

-   Bug fixes

-   Large cluster pruning strategy is changed to serial search (instead of binary)

1.3.2:

-   Expose backbone selection strategy to the user

-   Expose gene tree filtering parameters at refine stage to the user

1.3.1:

-   Getting rid of pruning thresholds (automated finding)

1.3.0:

-   Pruning Large partitions

-   Changes in TreeCluster logic disallowing formation of very small partitions

1.2.1:

-   APPLES2 excludes sequences that are placed on internal nodes.

-   Filtered backbone sequences are no longer added to the query set.

### Tips/Tricks:

-   Set min_placements: 9999999 to only filter out low quality backbone sequences and return.

-   Set config["mainlines_config"]["n"] to number of species in the dataset and config["backbone"] to "de-novo". Then run snakemake with target {outdir}/backbone.nwk. This will do species tree inference without divide and conquer.
