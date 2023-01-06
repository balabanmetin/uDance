# uDance: Updating phylogenetic and phylogenomic trees using divide and conquer

<img src="https://github.com/balabanmetin/uDance/raw/master/logo.png" alt="uDance logo" width="300"/>

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

## uDance configuration settings

|            Parameter             |                                                Description                                                 |
|:--------------------------------:|:----------------------------------------------------------------------------------------------------------:|
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
