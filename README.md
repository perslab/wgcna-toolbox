# wgcna-toolbox

Perslab toolbox for Weighted Gene Co-Expression Network Analysis

## Robust WGCNA pipeline

*Dependencies:*
* optparse
* parallel
* irlba
* Seurat
* WGCNA
* dplyr
* reshape
* reshape2
* Matrix
* STRINGdb
* boot
* liger


### Overview

#### Pre-process the expression data: perform QC and select significant genes ([workflow schematic](https://drive.google.com/file/d/1fntPIANPdC5ix1zKf1-mmcSRvIFQ24aB/view?usp=sharing)) 

* Map genes from symbol to ensembl ID
* Convert factor metadata to one column per level with dummy variables (similar to model.matrix function)
* Save user parameters and built-in parameters to tables; keep track of n cells and n genes in each celltype subset
* Filter out genes expressed in very few cells
* Scale the data and regress out selected potential confounders, e.g. "nUMI", "percent.mito", "percent.ribo,  if present in obj@meta.data
* Find variable genes
* Perform Principal Component Analysis. Use Jackstraw, i.e. repeated PCA on resampled datasets, to determine significant PCs and genes; or to save time use top 5000 PCA loading genes or just variable genes 

#### Compute consensus Topological Overlap Matrix (TOM) and find gene modules

* Compute soft power for adjacency matrix
* Resample expression matrices 
* Compute consensus TOM across resampled expression matrices
* Run hclust to find clusters in the (inverse) TOM 
* Run cutreeHybrid cut out modules in the distance matrices (iterate over different sets of parameters if given)
* Merge close modules iteratively using kIMs or kMEs
* Compute fuzzy module membership for every gene-module pair, either Module Eigengenes (MEs), kMEs ('fuzzy' module membership) and primary kMEs (kME w.r.t. 'own' module); or intramodular k (kIMs), which are average distance between a gene and each gene in the module
* Perform an addition k-means-like reclustering step, reassigning genes to modules if the kME/kIM is more than 1.05 times higher
* Filter out genes without a significant gene-module expression profile correlation using a t-test
* Match modules between different parameter settings

#### Find Protein-Protein Interactions (PPI)
* Run gene modules against STRINGdb PPI enrichment database
* If desired, filter out modules which are not significant and any cell types which have no significant modules
* For each celltype, select the set of parameters with highest number of genes assigned to modules after checking the modules for PPI enrichment.

#### Find genetic enrichment
* Make colors unique across celltypes
* Compute rare variant and mendelian gene enrichment using Gene Set Enrichment Analysis
* Compute MAGMA gwas enrichment using a paired samples t-test on module gene versus all MAGMA gene p-values
* Re-compute MEs and kMEs after MAGMA GWAS filter; remove any celltypes without significant modules if desired

#### Compute eigengene - metadata correlations
* Compute correlation between metadata and eigengene expression, per cell type
* Compute p value for the correlations, adjusted for multiple testing across all cell types
* If desired, filter modules and cell types to retain only those significantly correlated with metadata

#### Compute and output tables
* Run parameters
* Run summary statistics
* Per cell subset summary statistics
* Module gene lists ordered by kME/kIM
* Module-module correlations and clustering 
* Dataframe with columns cell type, module, gene ensembl id and, if available, gene symbol 
* kMEs / kIM tables as csv
* Module eigenvectors ('u' of Singular Value Decomposition), of length n genes of each module, for computing cell embeddings later

#### Plotting - separate notebook
* Module assignment under different parameters
* Module assignment after each round of filtering (PPI, gwas, metadata correlations)
* Module GWAS enrichment with MAGNA
* Module geneset enrichment analysis 
* Cell x eigengene cell embedding heatmap across whole dataset
* Module-metadata correlations 
* Cell t-SNE plot
* Eigengene expression featureplots in t-SNE space
* module-module correlation matrix across all celltypes

### Usage

e.g.

`time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_clust_all.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell-1/  --data_prefix campbell  --run_index 1 --regress_out "c('percent.mito', 'percent.ribo', 'nUMI')" --min.cells 5  --genes_use PCA --pca_genes pca.genes --corFnc cor --networkType signed --hclustMethod average --minClusterSize "c(15, 25)" --deepSplit "c(2,4)" --moduleMergeCutHeight "c(0.2)" --fuzzyModMembership kIM --jackstrawnReplicate 0 --TOMnReplicate 50 --kM_reassign T --kM_signif_filter T  --PPI_filter T --data_organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits "c('t1d', 't2d', 'BMI')" --metadata_corr_col "c('SEX', 'AGE')" --metadata_corr_filter_vals "c('male', 'p22-25')" --n_cores 10`

### Args

* `seurat_path`: Path to Seurat object, preferably with some QC based on nUMI and percent.mito
* `project_dir`: Optional. Project directory. Directory and subdirectories RObjects, plots, tables and log will be created if they do not exist. Defaults to the directory one level up from input data dir.
* `data_prefix`: Dataset prefix for output files.
* `run_prefix`: Run prefix for output files. 
* `autosave`: Save at four checkpoints during the script to enable resuming the session in case of problems? Defaults to `T`.
* `resume`: Resume from a previous session image? Must have same path and data_prefix. Options are `checkpoint_1` - `checkpoint_4`. Defaults to `NULL`.
* `quit_session`: Specify a checkpoint (`checkpoint_1` to `checkpoint_4`) after which to quit the script. Defaults to `NULL`.
* `metadata_subset_col`: Specify a `seurat@meta.data$` column to use for subsetting the Seurat object. If NULL (default) uses the `@ident` slot. 
* `metadata_corr_col`: Specify `seurat@meta.data$` column(s) for which to compute correlations with gene modules. Takes a character with a vector of `meta.data` column names e.g. `nUMI` or `c('nUMI', 'Age')`. For factor or character metadata, each levels is analysed as a dummy variable, so exercise caution. Defaults to `NULL`.
* `metadata_corr_filter_vals`: Specify one or more values within the `seurat@meta.data$` column(s). Retain only modules which are significantly (anti-) correlated (at present, there is no threshold value for the correlation). Takes a character with a vector of `meta.data` column names e.g. `Female` or `c('fasted', 'HFD')`. Case-insensitive. Defaults to `NULL`.
* `use.imputed`: Use data in the `obj@imputed` slot for the computations to replace the `@data` slot? If the `@imputed` slot is empty, will revert to the default (`FALSE`).
* `regress_out`: e.g. `"c('percent.mito','percent.ribo', 'nUMI')"`. Defaults to `NULL`
* `min.cells`: What is the minimum number of cells in each subset in the data in which a gene should be detected to not be filtered out? Integer, defaults to 5. 
* `genes_remove_dir`: "Path to a directory of gene lists saved in .csv, .tab, .txt, .RData, or .RDS, to remove before the analysis".
* `genes_use`: One of `"all"`, `"var.genes"` for seurat var.genes, or `"PCA"` for genes that load significantly on at least one significant PC. Defaults to `"PCA"`
* `pca_genes`: If `genes_use` == `"PCA"`, use `var.genes` or `all` genes to perform PCA to select genes based on PC loadings? If `num.replicate` is zero, select 5000 genes loading highly on the top PCs for downstream analysis. If non-zero, use JackStraw to identify significant PCs. If using `all` genes for the PCA, select significant genes based on the JackStraw; otherwise just use loadings.
* `corFnc`: Correlation function: either `"cor"` (Pearson) or `"bicor"` - biweighted midcorrelation. Defaults to `"cor"`
* `networkType`: `"signed"`, `"signed hybrid"` or `"unsigned"`. '"signed"' scales correlations to [0:1]; '"unsigned"' takes the absolute value (but the TOM can still be '"signed"'); '"signed hybrid"' sets negative correlations to zero. Defaults to `"signed"`. `"signed hybrid"` may be used together with `anti_cor_action == "kME_reassign"` to reassign genes to a new module if they are more anticorrelated with the module eigengene than they are positively correlated with their current module eigengene.
* `hclustMethod`: Hierarchical clustering agglomeration method. See `hclust()` documentation. For bulk RNAseq data, consider `average`. Defaults to `complete``
* `minClusterSize`: Minimum genes needed to form a module. Takes a vector with one or more values, given as a string, e.g. `"c(15,20)"`. Recommended range 5-25. Defaults to `"c(15)"`
* `deepSplit`: Controls the sensitivity of the `cutreeDynamic` algorithm. Takes a vector with one or more values, given as a string, e.g. `"c(2,3)"`. Takes integer values 0-4, defaults to `"c(2)"`. 
* `moduleMergeCutHeight`: Cut-off level for the variable (1-correlation) for merging eigengenes. Takes a vector with one or more values, given as a string, e.g. `"c(0.20, 0.25)"`. Recommended value range 0.05-0.25
* `pamStage`: For `cutreeHybrid`. Perform additional Partition Around Medroids step? Takes a vector with one or two values, given as a string, e.g. `"c(TRUE,FALSE)"`, default `"c(TRUE)"`
* `kM_reassign`: reassign genes according to their correlation with the module Eigengenes (if fuzzyModMembership=='kME') or mean distance from each module's genes? Default 'TRUE')
* `kM_signif_filter`: remove genes with a non-significant fdr-corrected correlation p-value with their module expression?
* `jackstrawnReplicate`: Number of times to re-run PCA after permuting a small proportion of genes to perform empirical significance tests, i.e. the `JackStraw` procedure (see `pca_genes` above). Integer, defaults to 500. 
* `TOMnReplicate`: Number of times to permute the dataset, defaults to 100
* `fuzzyModMembership: `kME`: gene-eigengene correlation or `kIM`: average distance from gene to each module gene. Also used when computing cell embeddings on gene modules
* `scale_MEs_by_kIMs`: Disactivated 
* `PPI_filter`: Use Protein-Protein Interactions to validate modules? Defaults to TRUE
* `data_organism`: `hsapiens` or `mmusculus`.
* `magma_gwas_dir`: MAGMA input GWAS data directory as a character. Outputs results per subdirectory. Defaults to `/projects/jonatan/tmp-bmi-brain/data/magma/`.
* `gwas_filter_traits`: Filter out modules not significantly correlated with matching gwas studies within the magma_gwas_dir. Takes a character with a vector of character names to match within the filename of the GWAS , e.g. `body_BMI_Locke2015` or `c('BMI', 'T1D', 'T2D')`. Case-insensitive. Defaults to `NULL`
* `n_cores`: Number of cores to use for parallelization. Defaults to 5
