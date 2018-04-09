# wgcna-toolbox
Perslab toolbox for Weighted Gene Co-Expression Network Analysis

robust_wgcna_pipeline

######################################################################
############################# OVERVIEW ###############################
######################################################################

Finds 'robust' gene modules in each obj@ident class in a Seurat format dataset:

 1. Permute the dataset and compute a consensus Topological Overlap Matrix (TOM)
 2. Cluster on the consensus TOM to find modules. Identify eigengenes. Merge modules with highly correlated eigengenes.
      If compare_params == TRUE,
        Plot the modules found with different parameters.
        Select a single set of modules, corresponding to a single set of parameters, based on a module quality statistic.
 3. Use the eigengenes to assign a kME value (pairwise correlation between gene and eigengene expression,
      a fuzzy module membership score) for each gene-module pair.
 4. If anti_cor_action = "kME_reassign", reassign genes with a negative kME to another module more than
      1.5 times the magnitude of the kME to their own module
 5. Filter module for PPI enrichment to retain significant modules.
 6. Export as Robjects and .csv

######################################################################
############################## USAGE #################################
######################################################################

e.g.

time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/robust_wgcna_pipeline_0.5.R --data_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --dir_project /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-n12n13-8/ --data_prefix campbell-n12n13-8 --compare_params FALSE --scale_data TRUE --genes_use PCA_5000 --corFnc cor --networkType signed --anti_cor_action NULL --minClusterSize 20 --deepSplit 2 --moduleMergeCutHeight 0.2 --nPermutations 20 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --save_plots TRUE --plot_permuted F --n_cores 5

######################################################################
############################### IN ###################################
######################################################################

Seurat object with cell type identities in the obj@ident slot

######################################################################
############################## OUT ###################################
######################################################################

############################## /TABLES ###############################

INFO_user_parameters

INFO_subsets_n_cells

INFO_built_in_parameters

kMEs / kMEs_PPI: One csv file per subsetName, saved to dir_tables, each with the following columns:
  colors:       module assignment
  genes:        gene name
  kMEs:         *One column per module*. The signed kMEs for the gene to every module.

pkMEs / pkMEs_PPI: One csv file per cell type, saved to dir_tables, each  with the following columns:
  genes:        gene name
  pkMEs:        the signed kME (gene expression correlation with the module 'eigengene', i.e. first PC)
                for the gene to its own module

ensembl_out / ensembl_PPI_out: one file per cell type, saved to dir_tables, each with
  module:       module colors
  ensemble_IDs: ensemble gene IDs

ensembleID_proportion_not_mapped / PPI: a file whose title gives the prop not mapped

############################### #ROBJECTS ##############################

results: a list per cell type, saved to dir_RObjects, containing
  metaData:   metadata for the subset
  datExprFilter: expression matrix filtered at the consensusTOM step
  geneTree:   output of hclust
  cutree:     output of cutreeHybrid or cutreeDynamic
  merged:     output of mergeCloseModules
  colors:     colors assigning each gene to a module
  ensembl_out: as for csv above
  ensembl_PPI_out:
              as for csv above
  MEs:        module eigengenes
  kMEs:       all module kMEs
  pkMEs:      primary module kMEs
  module_PPI:  database with null hypothesis PPI interactions and p-value for each module
              used to filter modules
  colors_PPI: filtered colors
  MEs_PPI:    filtered Module Eigengenes
  kMEs_PPI:   filtered kMEs
  pkMEs_PPI:  filtered principal kMEs
  fitIndices: output of pickSoftThreshold
  softPower:  power selected from options given to pickSoftThreshold based
              on scale-free topology fit
  consTomDS:  consensus Topological Overlap Matrix after permuting data
              and merging individual TOMs
  clust_qual_params: if "compare_params"=TRUE, the parameters tested
  cutree_params_final:
              final parameters used for cutreeHybrid and mergeCloseModules

############################### /PLOTS #################################

pickSoftThresholdSFTfit: pickSoftThreshold scale free topology plot
pickSoftThresholdMeanCon: pickSoftThreshold mean connectivity plot
diffpermutedColors: if plot_permuted == T, the plot of the modules found in each permuted dataset
colors: the final colors before and after filtering for PPI enrichment
compareParams: if compare_params == T, colors found using each set of parameters
