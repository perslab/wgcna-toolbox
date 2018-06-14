# wgcna-toolbox

Perslab toolbox for Weighted Gene Co-Expression Network Analysis

## Robust WGCNA pipeline

*Dependencies:*
* optparse
* parallel
* Seurat
* WGCNA
* dplyr
* reshape
* reshape2
* Matrix
* Stringdb
* dplyr
* ggplot2
* liger

### Overview

#### Pre-process the expression data: perform QC and select significant genes ([workflow schematic](https://drive.google.com/file/d/1fntPIANPdC5ix1zKf1-mmcSRvIFQ24aB/view?usp=sharing)) 

* Map genes from symbol to ensembl ID
* Convert factor metadata to one column per level
* Save user parameters and built-in parameters to tables
* Filter out genes expressed in very few cells
* Scale the data and regress out "nUMI", "percent.mito", "percent.ribo" if present in meta.data
* Find variable genes
* Perform Principal Component Analysis
* Use Jackstraw, i.e. repeated PCA on resampled datasets, to determine significant PCs and genes

#### Compute consensus Topological Overlap Matrix (TOM) and find gene modules

* Compute soft power for adjacency matrix
* Resample expression matrices 
* Compute consensus TOM
* Run hclust to find clusters in the (inverse) TOM 
* Run cutreeHybrid cut out modules in the distance matrices (iterate over different sets of parameters if given)
* merge close modules
* filter out runs with only grey modules
* compute Module Eigengenes (MEs), kMEs ('fuzzy' module membership) and primary kMEs (kME w.r.t. 'own' module)
* Match modules between different parameter settings

#### Filter modules and parameterisations using STRINGdb's Protein-Protein Interactions (PPI)
* Check STRINGdb PPI
* Filter out modules which are not significant and any cell types which have no significant modules
* For each celltype, select the set of parameters with highest number of genes assigned to modules after checking the modules for PPI enrichment.

#### Filter modules on genetic enrichment
* Compute MAGMA gwas correlations
* Filter modules based on gwas enrichment
* _TODO: add rare variants test_
* Assign new unique module names across all cell clusters to avoid identical or similarly named modules
* Re-compute MEs and kMEs after MAGMA GWAS filter; remove any celltypes without significant modules

#### Compute eigengene - metadata correlations
* Compute correlation between metadata and eigengene expression, per cell type
* Compute p value for the correlations, adjusted for multiple testing across all cell types
* Filter modules and cell types to retain only those significantly correlated with metadata

#### Perform Gene Set Enrichment Analysis (GSEA) and compute eigengene matrix
* Use LIGER to perform GSEA
* Filter GSEA results for fdr p-value significance
* _TODO: Filter modules for significant gene set enrichment_
* _TODO: use semantic correlation to cluster GSE terms_

#### Compute and output tables
* Prepare module gene lists ordered by kME
* _TODO: Compute module-module correlations; merge highly correlated modules_
* Compute cell x eigengene cell embedding matrix across all celltypes and modules
* Prepare dataframe with columns cell type, module, gene ensembl id and, if available, gene symbol 
* Output kMEs tables as csv
* Output GSEA tables as csv

### Plotting - separate notebook
* _TODO: Module assignment after each round of filtering_
* Module GWAS enrichment 
* _TODO: Cell x eigengene cell embedding matrix with genetic and metadata annotations_
* Module-metadata correlations _TODO: incorporate into eigenmatrix plot above_
* Cell t-SNE plot
* Eigengene expression featureplots
* _TODO: module-module correlation matrix across all celltypes_
* _TODO: module preservation across different tissues_

### Usage

e.g.

`time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-22/ --meta.data_corr_col 'c("X2.group", "X3.batches", "X4.sex", "X5.Diet", "X6.FvF")' --meta.data_corr_filter_vals "c('X2.group', 'X3.batches', 'X4.sex', 'X5.Diet', 'X6.FvF')" --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits 'c("BMI", "T1D", "T2D")' --data_prefix campbell-AgRP-22 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod complete --minClusterSize "c(15)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 400 --TOM.num.replicate 100 --replace T --organism mmusculus --plot_permuted F --n_cores 5`

### Args

* `seurat_path`: Path to Seurat object, preferably with some QC based on nUMI and percent.mito
* `project_dir`: Optional. Project directory. Directory and subdirectories RObjects, plots, tables and log will be created if they do not exist. Defaults to the directory one level up from input data dir.
* `data_prefix`: Dataset prefix for output files. Defaults to today's date.
* `meta.data_ID`: Specify the name of a `seurat@meta.data$...` column to use for subsetting the Seurat object. If NULL (default) uses the `@ident` slot.
* `resume`: Resume from a previous session image? Must have same path and data_prefix. Options are `checkpoint_1` - `checkpoint_4`. Defaults to `NULL`.
* `meta.data_subset_col`: Specify a `seurat@meta.data$` column to use for subsetting the Seurat object. If NULL (default) uses the `@ident` slot. 
* `meta.data_corr_col`: Specify `seurat@meta.data$` column(s) for which to compute correlations with gene modules. Takes a character with a vector of `meta.data` column names e.g. `nUMI` or `c('nUMI', 'Age')`. For factor or character metadata, each levels is analysed as a dummy variable, so exercise caution. Defaults to `NULL`.
* `meta.data_corr_filter_vals`: Specify one or more values within the `seurat@meta.data$` column(s). Retain only modules which are significantly (anti-) correlated (at present, there is no threshold value for the correlation). Takes a character with a vector of `meta.data` column names e.g. `Female` or `c('fasted', 'HFD')`. Case-insensitive. Defaults to `NULL`.
* `use.imputed`: Use data in the `obj@imputed` slot for the computations to replace the `@data` slot? If the `@imputed` slot is empty, will revert to the default (`FALSE`).
* `min.cells`: What is the minimum number of cells in each subset in the data in which a gene should be detected to not be filtered out? Integer, defaults to 5. 
* `do.center`: Use centered data? In either case data is scaled and nUMI and mitochrondrial genes are regressed out. Default to `TRUE`
* `genes_use`: One of `"all"`, `"var.genes"` for seurat var.genes, or `"PCA"` for genes that load significantly on at least one significant PC. Defaults to `"PCA"`
* `pca_genes`: If `genes_use` == `"PCA"`, use `var.genes` or `all` genes to perform PCA to select genes based on PC loadings? If `num.replicate` is zero, select 5000 genes loading highly on the top PCs for downstream analysis. If non-zero, use JackStraw to identify significant PCs. If using `all` genes for the PCA, select significant genes based on the JackStraw; otherwise just use loadings.
* `corFnc`: Correlation function: either `"cor"` (Pearson) or `"bicor"` - biweighted midcorrelation. Defaults to `"cor"`
* `networkType`: `"signed"`, `"signed hybrid"` or `"unsigned"`. '"signed"' scales correlations to [0:1]; '"unsigned"' takes the absolute value (but the TOM can still be '"signed"'); '"signed hybrid"' sets negative correlations to zero. Defaults to `"signed"`. `"signed hybrid"` may be used together with `anti_cor_action == "kME_reassign"` to reassign genes to a new module if they are more anticorrelated with the module eigengene than they are positively correlated with their current module eigengene.
* ~~`anti_cor_action`: Optional. '"kME_reassign"' reassigns genes with a negative kME more than 1.25 the kME w.r.t. their own (primary) module. Should be used only with 'networkType == "signed hybrid"'~~
* `hclustMethod`: Hierarchical clustering agglomeration method. See `hclust()` documentation. For bulk RNAseq data, consider `average`. Defaults to `complete``
* `minClusterSize`: Minimum genes needed to form a module. Takes a vector with one or more values, given as a string, e.g. `"c(15,20)"`. Recommended range 5-25. Defaults to `"c(15)"`
* `deepSplit`: Controls the sensitivity of the `cutreeDynamic` algorithm. Takes a vector with one or more values, given as a string, e.g. `"c(2,3)"`. Takes integer values 0-4, defaults to `"c(2)"`. 
* `moduleMergeCutHeight`: Cut-off level for the variable (1-correlation) for merging eigengenes. Takes a vector with one or more values, given as a string, e.g. `"c(0.20, 0.25)"`. Recommended value range 0.05-0.25
* `pamStage`: For `cutreeHybrid`. Perform additional Partition Around Medroids step? Takes a vector with one or two values, given as a string, e.g. `"c(TRUE,FALSE)"`, default `"c(TRUE)"`
* `jackstraw.num.replicate`: Number of times to resample to make null distributions for empirical significance tests in `JackStraw` (see `pca_genes` above). Integer, defaults to 500. 
* `TOM.num.replicate`: Number of times to permute the dataset, defaults to 100
* `replace`: Sample with replacement? If `TRUE`, uses all samples, if `FALSE`, uses 66% each time. Defaults to `TRUE`.
* `organism`: `hsapiens` or `mmusculus`. 
* `magma_gwas_dir`: MAGMA input GWAS data directory as a character. Outputs results per subdirectory. Defaults to `/projects/jonatan/tmp-bmi-brain/data/magma/`.
* `gwas_filter_traits`: Filter out modules not significantly correlated with matching gwas studies within the magma_gwas_dir. Takes a character with a vector of character names to match within the filename of the GWAS , e.g. `body_BMI_Locke2015` or `c('BMI', 'T1D', 'T2D')`. Case-insensitive. Defaults to `NULL`
* `plot_permuted`: Compute and plot the modules on each resampled dataset? Good for visually inspecting how robust modules are to resampling, but computationally expensive. Defaults to `FALSE`
* `n_cores`: Number of cores to use for parallelization. Defaults to 5
      
### Returns _OUT OF DATE_

`/tables`

INFO_user_parameters: Parameters provided by the user when running the script from terminal

INFO_subsets_n_cells: numbers of cells in each subset within the seurat objected provided

INFO_built_in_parameters: a subset of built-in parameters (sourced from a separate script)

kMEs / kMEs_PPI: One csv file per subsetName, saved to dir_tables, each with the following columns:
* colors:       module assignment
* genes:        gene name
* kMEs:         The signed kMEs for the gene to every module, one column per module

pkMEs / pkMEs_PPI: One csv file per cell type, saved to dir_tables, each  with the following columns:
* genes:        gene name
* pkMEs:        the signed kME (gene expression correlation with the module 'eigengene', i.e. first PC)
                for the gene to its own module

ensembl_out / ensembl_PPI_out: one file per cell type, saved to `dir_tables`, each with
* module:       module colors
* ensemble_IDs: ensemble gene IDs

ensembleID_proportion_not_mapped / PPI: a file whose title gives the prop not mapped

`/RObjects`

session image for each cell type, saved to `dir_RObjects`, containing amongst others
  
1. datExprFilter: expression matrix filtered at the `consensusTOM` step
2. metaData:   metadata for the subset
3. geneTree:   output of `hclust`
4. cutree:     output of `cutreeHybrid` or `cutreeDynamic`
5. merged:     output of `mergeCloseModules`
6. colors:     colors assigning each gene to a module
7. ensembl_out: as for csv above
8. ensembl_PPI_out: as for csv above
9. MEs:        module eigengenes
10. kMEs:       all module kMEs
11. pkMEs:      primary module kMEs
12. module_PPI:  database with null hypothesis PPI interactions and p-value for each module
              used to filter modules
13. colors_PPI: filtered colors
14. MEs_PPI:    filtered Module Eigengenes
15. kMEs_PPI:   filtered kMEs
16. pkMEs_PPI:  filtered principal kMEs
17. fitIndices: output of `pickSoftThreshold`
18. softPower:  power selected from options given to `pickSoftThreshold` based
              on scale-free topology fit
19. clust_qual_params: if `--compare_params TRUE`, the parameters tested
20. cutree_params_final:
              final parameters used for `cutreeHybrid` and `mergeCloseModules`

`/plots`

* pickSoftThresholdSFTfit: pickSoftThreshold scale free topology plot
* pickSoftThresholdMeanCon: pickSoftThreshold mean connectivity plot
* diffpermutedColors: if plot_permuted == T, the plot of the modules found in each permuted dataset
* diffParams_colors_PPI_order: if the user gave several sets of parameters, a plot of the modules, ordered by the number of assigned genes after PPI enrichment test.
* final_colors: the final module assignment before and after PPI enrichment test
`/projects/jonatan/tmp-bmi-brain/`
* MAGMA outputs 

