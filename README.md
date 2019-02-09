# wgcna-toolbox

Perslab toolbox for Weighted Gene Co-Expression Network Analysis

## Robust WGCNA pipeline

*Dependencies:*
* optparse
* parallel
* Seurat 3.0
* WGCNA
* dplyr
* reshape
* reshape2
* Matrix

### Overview

#### Pre-process the expression data: perform QC and select significant genes ([workflow schematic](https://drive.google.com/file/d/1fntPIANPdC5ix1zKf1-mmcSRvIFQ24aB/view?usp=sharing)) 
* Split cells into subsets by celltype and perform basic QC
* Filter genes for the analysis with the arguments featuresUse and nFeatures. Use either 'var.features' for variable features using Seurats FindVariableFeatures function (vst), PCLoading for top loading PC genes (using 50 PCs), JackStrawPCLoading for top loading PC genes on significant PCs after the JackStraw test, or JackStrawPCSignif to select genes with best p-values after JackStraw test.

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
* A single matrix of cell module embeddings using logNormalized data from @data slot in Seurat 
* Run parameters
* Run summary statistics
* Celltype summary statistics
* Dataframe with columns cell type, module, gene and weight (pkM) 
* kMEs / kIM tables as csv
* final session image with all of the above (in \RObjects subdir)


### Usage
export PATH="/usr/local/R-3.5.1/bin/:$PATH" 

time Rscript /projects/jonatan/tools/wgcna-src/wgcna-toolbox/rwgcna_main_seurat3.0.R --pathDatExpr /projects/jonatan/tmp-mousebrain/RObjects/L5_Neurons_sub_20190208.RDS.gz --dirProject /projects/jonatan/mousebrain_8/ --prefixData mb --prefixRun Neurons_sub_ClusterName_8 --dataType sc --colIdents ClusterName --vec_colAnnot 'c(NULL)' --assayUse RNA --slotUse scale.data --varsToRegress 'c(NULL)' --minGeneCells 20 --minCellClusterSize 50 --featuresUse PCLoading --nFeatures 5000 --nRepJackStraw 0 --corFnc cor --networkType "c('signed hybrid')" --nRepTOM 100  --hclustMethod average --minClusterSize 15 --deepSplit "c(2)" --moduleMergeCutHeight 'c(0.15)' --pamStage 'c(T)' --kMReassign F --kMSignifFilter T --fuzzyModMembership kME  

### Help
Rscript /projects/jonatan/tools/wgcna-src/wgcna-toolbox/rwgcna_main_seurat3.0.R --help
