# wgcna-toolbox

Perslab toolbox for Weighted Gene Co-Expression Network Analysis

## Robust WGCNA pipeline

### Overview

Finds 'robust' gene modules in each `obj@ident` class in a Seurat format dataset:

1. Do seurat pre-processing [workflow](https://drive.google.com/file/d/1fntPIANPdC5ix1zKf1-mmcSRvIFQ24aB/view?usp=sharing) 
2. Permute the dataset and compute a consensus Topological Overlap Matrix (TOM)
3. Cluster on the consensus TOM to find modules. Identify eigengenes. Merge modules with highly correlated eigengenes.
      If `--compare_params TRUE`, plot the modules found with different parameters. Select a single set of modules, corresponding to a single set of parameters, based on a module quality statistic (experimental).
4. Use the eigengenes to assign a kME value (pairwise correlation between gene and eigengene expression, a fuzzy module membership score) for each gene-module pair.
5. If `anti_cor_action == kME_reassign`, reassign genes with a negative kME to another module more than 1.25 times the magnitude of the kME to their own module.
6. Filter modules for PPI enrichment.
7. Perform MAGMA GWAS analysis - outputs to `/projects/jonatan/tmp-bmi-brain/`

### Usage

e.g.

`time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --data_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_neurons.RData --dir_project /projects/jonatan/wgcna/campbell-neurons-1/ --data_prefix campbell-neurons-1 --compare_params FALSE --scale_data FALSE --genes_use PCA_5000 --corFnc cor --networkType signed --anti_cor_action NULL --minClusterSize 20 --deepSplit 2 --moduleMergeCutHeight 0.2 --nPermutations 100 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --save_plots TRUE --plot_permuted F --n_cores 5`

### Args

* `data_path`: Path to Seurat object with cell type identities in the `obj@ident` slot stores as .RData
* `dir_project`: Optional. Project directory. Directory and subdirectories RObjects, plots, tables will be created if they do not exist. Defaults to the directory one level up from input data dir.
* `data_prefix`: Dataset prefix for output files. Defaults to today's date.
* `compare_params` = Compare many different cutreeHybrid parameters, plots different colorings and select best using module statistics. Experimental feature. Defaults to `FALSE`
* `do.center`: Use centered data? In either case data is scaled and nUMI and mitochrondrial genes are regressed out. Default to `TRUE`
* `genes_use`: One of `"all"`, `"var.genes`", "`hvg_ _number of highly variable genes_"`, `"PCA_ _number of high loading genes_"` Defaults to `"PCA_5000"`
* `corFnc`: Correlation function: either `"cor"` (Pearson) or `"bicor"` - biweighted midcorrelation. Defaults to `"cor"`
* `networkType`: `"signed"`, `"signed hybrid"` or `"unsigned"`. '"signed"' scales correlations to [0:1]; '"unsigned"' takes the absolute value (but the TOM can still be '"signed"'); '"signed hybrid"' sets negative correlations to zero. Defaults to `"signed"`. `"signed hybrid"` may be used together with `anti_cor_action == "kME_reassign"` to reassign genes to a new module if they are more anticorrelated with the module eigengene than they are positively correlated with their current module eigengene.
* `anti_cor_action`: Optional. '"kME_reassign"' reassigns genes with a negative kME more than 1.25 the kME w.r.t. their own (primary) module. Should be used only with 'networkType == "signed hybrid"'
* `minClusterSize`: Minimum genes needed to form a module, recommended range 5-25. Defaults to 15
* `deepSplit`: Controls the sensitivity of the `cutreeDynamic` algorithm. Takes integer values 0-4, defaults to 2. Only applicable if `test_params == FALSE`
* `moduleMergeCutHeight`: Cut-off level for the variable (1-correlation) for merging eigengenes. Recommended value range 0.1-0.2
* `replace`: Sample with replacement? If `TRUE`, uses all samples, if `FALSE`, uses 66% each time. Defaults to `TRUE`.
* `nPermutations`: Number of times to permute the dataset, defaults to 100
* `STRINGdb_species`: Species for which to retrieve protein data from `STRINGdb` to validate clusters. Defaults to 10090, which is mus musculus. To skip this step set to `NULL`
* `ensembl_dataset`: Dataset for ensemblIDs for outputting colors for LD score regression. Defaults to `"mmusculus_gene_ensembl"`. To skip this step set to `NULL`
* `save_plots`: Save plots? Defaults to `TRUE`
* `plot_permuted`: Compute and plot the modules on each resampled dataset? Good for visually inspecting how robust modules are to resampling, but computationally intensive. Defaults to `FALSE`
* `n_cores`: Number of cores to use for parallelization. Defaults to 5
      
### Returns

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
* colors: the final colors before and after filtering for PPI enrichment
* compareParams: if compare_params == T, colors found using each set of parameters

`/projects/jonatan/tmp-bmi-brain/`
* MAGMA outputs 

