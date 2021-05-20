# wgcna-toolbox

Perslab toolbox for Weighted Gene Co-Expression Network Analysis

## Robust WGCNA pipeline (rwgnca_main_seurat3.0.R)

*Dependencies:*
* optparse
* parallel
* Seurat >= 3.0
* WGCNA >- 1.68
* magrittr
* Matrix
* data.table
* Biobase
* reshape
* reshape2

### Quickstart

Make sure you have log-normalized expression data and metadata with cell names in the first column. Then call

```
time Rscript ./rwgcna_main_seurat3.0.R --pathDatExpr ./expressionData.csv.gz --pathMetadata ./metadata.csv --dirProject /my/project/dir/ --prefixData mousebrain_neurons --prefixRun 1 --dataType sc --colIdents ClusterName 
```

### Parameter hints

```
Rscript ./rwgcna_main_seurat3.0.R --help
```

### Workflow details

#### Subset the data and select variable genes for WGCNA
* Split cells into subsets by annotation (i.e. celltype for single-cell, tissue-type for bulk RNA-seq)
* Select genes for the analysis with the arguments `featuresUse` and `nFeatures`. Use either `var.features` for variable features using Seurats `FindVariableFeatures` function (vst), `PCLoading` for top loading PC genes, `JackStrawPCLoading` for top loading PC genes PCs found to be significant using the JackStraw test, or `JackStrawPCSignif` to select genes with best p-values on significant PCs identified by the JackStraw test.

#### Compute adjacency matrix, consensus Topological Overlap Matrix (TOM) and find gene modules

* Compute soft power for adjacency matrix
* Resample expression matrices (.66, without replacement)
* Compute consensus TOM across resampled expression matrices
* Run hclust to find clusters in the (inverse) TOM 
* Run cutreeHybrid to cut out modules in the distance matrices (iterate over different sets of parameters if given)
* Merge close modules iteratively using kIMs or kMEs
* Compute fuzzy module membership for every gene-module pair, either kMEs ('fuzzy' module membership) or intramodular k (kIMs), which are average distance between a gene and each gene in the module
* Perform an additional k-means-like reclustering step, reassigning genes to modules if the kME/kIM is more than 1.2 times higher
* Filter out genes without a significant gene-module expression profile correlation using a t-test
* Make module names unique across all celltypes/tissues

#### Outputs

##### Results:
* table with columns cell type, module, gene and gene weight/centrality (kME or kIM)
* module gene weights for all genes with respect to all modules (kME / kIM)
* matrix of original expression data embedded into module space

##### Metadata:
* table of parameters
* table with run summary statistics
* table with per celltype subset summary statistics
* final session image with all of the above (in \RObjects subdir)

### Issues / bugs

Submit an issue with
* Descriptive title summarising the issue
* What you did
* Expected behaviour
* Actual behaviour
* R version and, if relevant, package versions
* Relevant output
