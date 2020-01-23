# wgcna-toolbox

Perslab toolbox for Weighted Gene Co-Expression Network Analysis

## Robust WGCNA pipeline (rwgnca_main_seurat3.0.R)

*Dependencies:*
* optparse
* parallel
* Seurat 3.0
* WGCNA
* dplyr
* reshape
* reshape2
* Matrix
* data.table

### Overview

#### Subset the data and select variable genes for WGCNA
* Split cells into subsets by celltype
* Select genes for the analysis with the arguments `featuresUse` and `nFeatures`. Use either `var.features` for variable features using Seurats `FindVariableFeatures` function (vst), `PCLoading` for top loading PC genes, `JackStrawPCLoading` for top loading PC genes PCs found to be significant using the JackStraw test, or `JackStrawPCSignif` to select genes with best p-values on significant PCs identified by the JackStraw test.

#### Compute adjacency matrix, consensus Topological Overlap Matrix (TOM) and find gene modules

* Compute soft power for adjacency matrix
* Resample expression matrices 
* Compute consensus TOM across resampled expression matrices
* Run hclust to find clusters in the (inverse) TOM 
* Run cutreeHybrid cut out modules in the distance matrices (iterate over different sets of parameters if given)
* Merge close modules iteratively using kIMs or kMEs
* Compute fuzzy module membership for every gene-module pair, either Module Eigengenes (MEs), kMEs ('fuzzy' module membership) and primary kMEs (kME w.r.t. 'own' module); or intramodular k (kIMs), which are average distance between a gene and each gene in the module
* Perform an addition k-means-like reclustering step, reassigning genes to modules if the kME/kIM is more than 1.2 times higher
* Filter out genes without a significant gene-module expression profile correlation using a t-test
* Make module names unique across celltypes

#### Outputs
* table with columns cell type, module, gene and gene weight (primary kM) 
* module gene weights for all genes and modules i.e. kMEs / kIMs as csv
* matrix of original expression data embedded into module space
* table of parameters
* table with run summary statistics
* table with per celltype subset summary statistics

* final session image with all of the above (in \RObjects subdir)

### Usage

time Rscript ./rwgcna_main_seurat3.0.R --pathDatExpr ./expressionData.csv.gz --pathMetadata ./metadata.csv --dirProject /projects/jonatan/mousebrain_8/ --prefixData mb_neurons --prefixRun 81 --dataType sc --colIdents ClusterName--minGeneCells 20 --minCellClusterSize 50 --featuresUse PCLoading --nFeatures 5000 --nRepJackStraw 0 --corFnc cor --networkType "c('signed hybrid')" --nRepTOM 100  --hclustMethod average --minClusterSize 15 --deepSplit "c(2)" --moduleMergeCutHeight 'c(0.15)' --pamStage 'c(T)' --kMReassign F --kMSignifFilter T --fuzzyModMembership kME  

### Help on parameters
Rscript /projects/jonatan/tools/wgcna-src/wgcna-toolbox/rwgcna_main_seurat3.0.R --help

### Issues / bugs

Submit an issue with
* Descriptive title summarising the issue
* What you did
* Expected behaviour
* Actual behaviour
* R version and, if relevant, package versions
* Relevant output
