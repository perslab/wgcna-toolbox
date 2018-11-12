# multi-analysis module preservation pipeline
# overview: TODO

# loomR docs:
## Introduction to loomR: https://satijalab.org/loomR/loomR_tutorial.html
## loomR MCA tutorial: https://satijalab.org/seurat/mca_loom.html
## loomR full specs: https://satijalab.org/loomR/loomR.pdf
## loom file specification: http://linnarssonlab.org/loompy/format/index.html


# usage (loomR requires R/3.4.3)
## module use /tools/modules
## module unload R
## module load R/3.4.3  
## export R_MAX_NUM_DLLS=999
## time Rscript /projects/jonatan/tools/wgcna-src/wgcna-toolbox/WGCNA_multi_compare.R --dir_project_WGCNA /projects/jonatan/tmp-mousebrain/ --dir_project_data /projects/jonatan/tmp-mousebrain/ --path_data /projects/jonatan/tmp-mousebrain/RObjects/L5.RDS --prefixes_WGCNA_run 'c("Neurons_ClusterName_2")' --prefix_out mousebrain_neuronMods_1  --scale_center_data F --n_cores 20

######################################################################
########################## DEFINE FUNCTIONS ##########################
######################################################################

LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}

current.dir = paste0(LocationOfThisScript(), "/")
#current.dir = "/projects/jonatan/tools/wgcna-src/wgcna-toolbox/"

source(file=paste0(current.dir,"rwgcna_functions.R"))

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("optparse", "Seurat", "WGCNA", "dplyr", "Biobase", "Matrix", "parallel", "readr"))#, "NetRep"))

######################################################################
########################### OPTPARSE #################################
######################################################################

option_list <- list(
  # make_option("--dirs_project_nwa", type="character", default=NULL,
  #             help = "Folders of gene network analyses, each must contain /tables subdirectory."),
  # make_option("--prefixes_nwa_run", type="character",
  #             help = "Quoted character vector of gene network analysis run prefixes, e.g. ''c('mousebrain_Neurons_1','mousebrain_Neurons_2')''"),
  make_option("--path_nwa_df", type="character",
              help = "path to dataframe from a gene network analysis run, in long format (i.e. one row per gene per celltype), containing 'cell_cluster', 'module', one or two gene name columns, and a column of numeric scores. I.e. one row per gene, e.g. rwgcna cell_cluster_module_genes.csv files"),  
  # make_option("--patterns_gene_module_loading_df", type="character",
  #             help = "Pattern to match alongside prefixes_nwa_run in each of dirs_project_nwa/tables to find generalised gene x module membership score dataframes, either one per run or one per celltype."),
  # # make_option("--gene_module_score_colnames", type="character", default=NULL,
  #             help = "list of named character vectors, one vector per module_df giving column names for 'gene', 'module' and 'score', optionally plus additional metadata. E.g. ''list(c(gene='score',module='blue',weight='pkIMs',tissue='tissue'), c())''. The order should match the other related arguments. If not provided, the script greps for 'module', 'gene', and looks for a numerical column"),
  make_option("--paths_datExpr", type="character",
              help = "Paths to expression data in loomR, seurat or other RObject format"),  
  make_option("--IDs_datExpr", type="character",
              help = "quoted vector of IDs for each of paths_datExpr"),
  # make_option("--paths_cellModEmbed", type="character", default = NULL,
  #             help = "Paths to precomputed cell x module embedding matrices [default %default]"),  
  make_option("--paths_datExpr_metadata", type="character", default = NULL,
              help = "Path to metadata in one of the standard (compressed) character separated formats. Should correspond to files in paths_datExpr; use NA_character_ to skip a datExpr. The script will always (also) use metadata stored within the path_datExpr objects if they are of the Seurat class. [default %default]"),  
  make_option("--metadata_ident_cols", type="character",
              help = "vector of characters to identify columns in metadata by which to split the expression data, named by the IDs_datExpr, e.g. ''c(datExpr_ID1='celltype', datExpr_ID2='clust.res.1')''. If only one is provided and dataset_ID names are omitted and there is one one IDs_datExpr, the script assigns that ID as name to all."),
  make_option("--nwa_cell_clusters", type="character", default=NULL,
              help = "a quoted vector of cell subsets for which to evaluate colorings. If left as NULL, uses all levels of the nwa df cell_cluster column in alphabetical order"),
  make_option("--ref_sets", type="character", default=NULL,
              help="a quoted list with a single character vector of identities found in metadata_ident_cols. The list component name gives the expression data ID. If dataset_ID names are omitted and there is one one IDs_datExpr, the script assigns that ID as name to all. The vector should be ordered by nwa_cell_clusters E.g. ''list(datExpr_ID1=c('n01','n02','n03'))''. If left as NULL, will grep the metadata_ident_cols for nwa_cell_clusters."),
  make_option("--test_sets", type="character",
              help="a quoted list of character vector of identities found in metadata_ident_cols, one vector per nwa_cell_clusters. List components should have corresponding reference sets as names; if not, the lists are assumed to be in the same order. Names of vector components give the dataset. E.g. ''list('n01'=c(datExpr_ID1='n02',datExpr_ID1='n03'), 'n04'=c(datExpr_ID1='n01',datExpr_ID1='n03'),'n03'=c(datExpr_ID1='n01',datExpr_ID1='n02'))''. If dataset ID names are omitted and there is one one IDs_datExpr, the script assigns that ID as name to all. If the argument is left as NULL, will grep the metadata_ident_cols for nwa_cell_clusters."),       #help = "quoted list of named character vectors with named components. List names are labelling levels (e.g. 'tissue', 'celltype'), vector names are column in datExpr metadata, vector values are regex to match levels. Use NA in a vector to skip a dataset for that labelling. E.g. ''list('tissue' = c(tissue='.*', NA), 'sub_celltype'=c(tissue_cell_type = '.*', ClusterName = '.*'))''"),
  make_option("--minCellClusterSize", type="integer", default=20L,
              help="What is the minimum number of cells in a subset to continue? Integer, [default %default]."),
              #help = "Similar to ident_groups_datExpr but for gene networks. List names should match those of ident_groups_datExpr (e.g. 'tissue', 'celltype'). Each vector should have an entry for each df in paths_nwa_df. E.g. ''list('tissue' = c(cell_type='.*', NA), 'celltype'=c(cell_type = '.*', cell_type='.*'))''  To bypass a dataset, place a 'NA_character_' in the vector. "),
  # make_option("--ident_groups_combs", type="character", default="c('all_to_all')",
  #             help = "Quoted vector with either a single component or one component for each labelling, i.e., vector in ident_groups, to indicate whether to evaluate all combinations ''c('all_to_all')'' or just pairwise matches ''c('pairwise')''; if the latter, the number of levels must be equal across datasets [default %default]"),
  make_option("--gene_names", type="character", default="hgnc|symbol|gene_name_optimal",
              help ="string or regex for grepping nwa_df and gene mapping dataframe and for output, e.g. 'hgnc|symbol|gene_name' or 'ensembl', [default %default]"),
  make_option("--dir_out", type="character",
              help = "Outputs go to /tables and /RObjects subdirectories"),  
  make_option("--prefix_out", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"),
  make_option("--dir_scratch", type="character", default = "/scratch/tmp-wgcna/",
              help = "Outputs go to /tables and /RObjects subdirectories"),
  make_option("--networkType", type="character", default = "signed hybrid",
              help = "for WGCNA modulePreservation"),
  make_option("--corFnc", type="character", default = "cor",
              help = "for WGCNA modulePreservation"),

  make_option("--data_organism", type="character", default="mmusculus",
              help = "'hsapiens' or 'mmusculus', [default %default]"),
  make_option("--scale_center_regress_data", type="logical", default = "TRUE",
              help = "[default %default]. If TRUE, and path_data points to an RObject, the script will try to use the @scale.data slot and, if missing, will normalize and scale the data and regress out nUMI, percent.mito and percent.ribo. If FALSE, for RObjects will use the @data slot, or, if empty, the @raw.data slot"),  
  # make_option("--scale_center_mods", type="logical", default = "TRUE",
  #             help = "scale expression of the modules? Recommended if the gene modules originate from different analyses with different scales. Set to F if desiring to scale at a later stage, e.g. after merging several embeddings matrices [default %default]."),  
  make_option("--do_plot", type="logical", default = "FALSE",
              help = "Output plots? [default %default]."),
  make_option("--RAM_Gb_max", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
  # make_option("--n_cores", type="integer", default = 15L,
  #             help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]")
  # make_option("--use_loom", type="logical", default = F,
  #             help = "Use loom HF5 format? Less RAM-heavy but error-prone [default %default]")
)

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

#if (!is.null(opt$dirs_project_nwa)) dirs_project_nwa <- eval(parse(text=opt$dirs_project_nwa))
#prefixes_nwa_run <- eval(parse(text=opt$prefixes_nwa_run))
path_nwa_df <- opt$path_nwa_df 
#patterns_gene_module_df <- eval(parse(text=opt$patterns_gene_module_df))
#patterns_gene_module_k_df <- eval(parse(text=opt$patterns_gene_module_k_df))
paths_datExpr <- eval(parse(text=opt$paths_datExpr))
IDs_datExpr <-  eval(parse(text=opt$IDs_datExpr))
if (!is.null(opt$paths_datExpr_metadata)) paths_datExpr_metadata <- eval(parse(text=opt$paths_datExpr_metadata))
metadata_ident_cols <-eval(parse(test=opt$metadata_ident_cols))
if (!is.null(opt$nwa_cell_clusters)) nwa_cell_clusters <- eval(parse(text=opt$nwa_cell_clusters)) 
ref_sets <- eval(parse(text=opt$ref_sets))
test_sets <- eval(parse(text=opt$test_sets))
minCellClusterSize <- opt$minCellClusterSize
#ident_groups_datExpr <- eval(parse(text=opt$ident_groups_datExpr))
#ident_groups_nwa <- eval(parse(text=opt$ident_groups_nwa))
gene_names <- opt$gene_names
dir_out <- opt$dir_out
prefix_out <- opt$prefix_out
dir_scratch <- opt$dir_scratch
networkType <- opt$networkType
corFnc <- opt$corFnc
data_organism <- opt$data_organism
scale_center_regress_data <- opt$scale_center_regress_data
#scale_center_mods <- opt$scale_center_mods
do_plot <- opt$do_plot
#n_cores <- opt$n_cores
RAM_Gb_max <- opt$RAM_Gb_max

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F, use="pairwise.complete.obs")

######################################################################
############################ CONSTANTS ###############################
######################################################################

# if specified output directory doesn't exist, create it 
if (!file.exists(dir_out)) {
  dir.create(dir_out) 
  message("dir_out not found, new one created")
}

dir_plots = paste0(dir_out,"plots/")
if (!file.exists(dir_plots)) dir.create(dir_plots) 

dir_tables = paste0(dir_out,"tables/")
if (!file.exists(dir_tables)) dir.create(dir_tables)

dir_RObjects = paste0(dir_out,"RObjects/")
if (!file.exists(dir_RObjects)) dir.create(dir_RObjects)

dir_log = paste0(dir_out,"log/")
if (!file.exists(dir_log)) dir.create(dir_log)

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

######################################################################
########################### VERIFY INPUT #############################
######################################################################

if (is.null(names(ref_sets)) & length(IDs_datExpr)==1) { names(ref_sets) <- IDs_datExpr }

if (all(sapply(test_sets, function(test_set) is.null(names(test_set)))) & length(IDs_datExpr)==1) { 
 test_sets <- lapply(test_sets, function(test_set) {
    names(test_set) <- rep(IDs_datExpr, times=length(test_set)) 
    return(test_set)
  })
}

if (is.null(names(test_sets))) names(test_sets) <- ref_sets[[1]]

if (is.null(names(metadata_ident_cols)) & length(IDs_datExpr)==1 & length(metadata_ident_cols)==1) {names(metadata_ident_cols) <- IDs_datExpr} 
# TODO
#if(!is.null(paths_metadata)) stopifnot(length(paths_metadata)==length(paths_datExpr))
# stopifnot(all.equal(sapply(ident_groups, length)))
# stopifnot(length(ident_groups_combs) %in% c(1, length(ident_groups)))

######################################################################
######################## LOAD MODULE DATA ############################
######################################################################

message("Loading gene module data")

gene_module_df <- load_obj(path_nwa_df)

# If there are duplicate modules between WGCNA runs, prefix module names with the WGCNA run 
# if (sapply(X=list_gene_module_df , function(df) unique(df[["module"]]), simplify = T) %>% unlist(use.names=F) %>% duplicated %>% any) {
#   for (i in 1:length(list_gene_module_df)) {
#     list_gene_module_df[[i]][["module"]] <- paste0(prefixes_nwa_run[i], "_",list_gene_module_df[[i]][["module"]]) # ensure modules from different runs remain distinct
#   }
# }

# If there are duplicate levels for celltype or tissue between WGCNA runs, prefix level names with the WGCNA run 

# for (ident_group in ident_groups_nwa) {
#   for (colname in names(ident_group)) {
#     if (sapply(X=list_gene_module_df , function(df) unique(df[[colname]]), simplify = T) %>% unlist(use.names = F) %>% duplicated %>% any) {
#       for (i in 1:length(list_gene_module_df)) {
#         list_gene_module_df[[i]][[colname]] <- paste0(prefixes_nwa_run[i], "_",list_gene_module_df[[i]][[colname]]) # ensure modules from different runs remain distinct
#       }
#     }
#   }
# }

# Merge into a single df
# gene_module_df <- Reduce(x = list_gene_module_df,f = rbind)
# rm(list_gene_module_df)

# Filter out tiny modules
mods_too_small <- names(table(gene_module_df[["module"]]))[table(gene_module_df[["module"]])<5]
gene_module_df <- gene_module_df[!gene_module_df[["module"]] %in% mods_too_small,]

# Identify gene and score columns
gene_col <- grep(gene_names, colnames(gene_module_df), value=T)
loading_col <- colnames(gene_module_df)[which(sapply(X=colnames(gene_module_df), FUN =function(colname) class(gene_module_df[[colname]]))=="numeric")]

######################################################################
####################### LOAD DATEXPR METADATA ########################
######################################################################
list_datExpr_metadata <- NULL
if (!is.null(paths_datExpr_metadata)){ 
  message("Loading metadata") 
  list_datExpr_metadata <- lapply(paths_datExpr_metadata, load_obj)
} 

######################################################################
######################## LOAD EXPRESSION DATA ########################
######################################################################

message("Loading expression matrix")

# load/open connection to expression matrix
# NB: may be huge
list_seurat_obj <- lapply(paths_datExpr, function(path_data) {
  if (grepl(pattern = "\\.loom", path_data)) {
    data_obj <-  connect(filename = path_data, mode = "r+")
  } else {
    data_obj <- load_obj(path_data)
  } 
  if (!any(c("seurat","loom") %in% class(data_obj))) { # not a loom object, nor a seurat object
    data_obj <- CreateSeuratObject(raw.data=data_obj, project= prefix_out, min.cells = -Inf, min.genes = -Inf)
  } else if ("loom" %in% class(data_obj)) {
    data_obj <- Seurat::Convert(from=data_obj, to="seurat")
  } # if already a seurat object, do nothing
  return(data_obj)
})

names(list_seurat_obj) <- IDs_datExpr
  
# Add any further metadata coming from another file
if (!is.null(list_datExpr_metadata)) {
  list_seurat_obj<- mapply(function(seurat_obj, datExpr_metadata) {
    if (all_equal(rownames(datExpr_metadata), colnames(seurat_obj@raw.data))) {
      seurat_obj <- AddMetaData(object = seurat_obj, metadata = data.frame(datExpr_metadata))
    }
  }, seurat_obj=list_seurat_obj, datExpr_metadata=list_datExpr_metadata, SIMPLIFY=F)
}

######################################################################
####################### MAP GENES TO COMMON NAMES ####################
######################################################################

# TODO: make orthologue mapping function
  
mapping = if (data_organism=="hsapiens") {
  load_obj("/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz")
} else if (data_organism=="mmusculus") {
  load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
}  

args <- list("X"=list_seurat_obj)
fun <- function(seurat_obj) {
  if (any(grepl("ENSG|ENSMUSG",rownames(seurat_obj@raw.data)))) {
    current_ensembl <- T
    tmp <- gene_map(df = seurat_obj@raw.data,
                    idx_gene_column = NULL,
                    mapping=mapping,
                    from="ensembl",
                    to="gene_name_optimal",
                    replace = T,
                    na.rm = T)
  } else {
    current_ensembl <- F
    tmp <- seurat_obj@raw.data
  }
  message("Counting percent ribo and percent mito")
  idx_mito.genes <- grepl(pattern = "^mt-", x = rownames(tmp), ignore.case=T)
  idx_ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(tmp), ignore.case=T)
  percent_mito <- Matrix::colSums(tmp[idx_mito.genes, ])/Matrix::colSums(tmp)
  percent_ribo <- Matrix::colSums(tmp[idx_ribo.genes, ])/Matrix::colSums(tmp)
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent_mito, col.name = "percent_mito")
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent_ribo, col.name = "percent_ribo")
  if (grepl("ensembl", gene_names) & !current_ensembl) {
    seurat_obj@raw.data <- gene_map(df = seurat_obj@raw.data,
                                    idx_gene_column = NULL,
                                    mapping=mapping,
                                    from="gene_name_optimal",
                                    to="ensembl",
                                    replace = T,
                                    na.rm = T)
  } else if (grepl("symbol|hgcn", gene_names) & current_ensembl) {
    seurat_obj@raw.data <- tmp
  }
  return(seurat_obj)
}
# Call the function
list_seurat_obj <- safeParallel(fun=fun, args=args, mapping=mapping)
rm(mapping)

######################################################################
######################## NORMALIZE EXPRESSION DATA ###################
######################################################################
  
message("Normalizing the expression data")

args <- list("X"=list_seurat_obj)
fun <- function(seurat_obj) NormalizeData(seurat_obj)
list_seurat_obj <- safeParallel(fun=fun, args=args)

######################################################################
################# SCALE AND REGRESS EXPRESSION DATA ##################
######################################################################

if (scale_center_regress_data) {
  message("scaling and regressing out confounders from the data")
  
  list_seurat_obj<-lapply(X = list_seurat_obj, FUN = function(seurat_obj) ScaleData(seurat_obj,
                        vars.to.regress = c("percent_mito", "nUMI", "percent_ribo")[c("percent_mito", "nUMI", "percent_ribo") %in% colnames(seurat_obj@meta.data)], 
                        do.scale=scale_center_regress_data,
                        do.center=scale_center_regress_data,
                        display.progress= T, 
                        do.par=scale_center_regress_data, 
                        num.cores=10))
}

######################################################################
############## EXTRACT SCALED DATA AND DELETE SEURAT OBJ #############
######################################################################

list_datExpr <- lapply(list_seurat_obj, function(seurat_obj) {
  out <- if (scale_center_regress_data) seurat_obj@scale.data else seurat_obj@data
  return(t(out))
  })
list_metadata <- lapply(list_seurat_obj, function(seurat_obj) seurat_obj@meta.data)
names(list_metadata) <- names(list_datExpr) <- IDs_datExpr
rm(list_seurat_obj)

######################################################################
################### COMPUTE MODULE PRESERVATION ######################
######################################################################

list_modulePreservation_out <- list()
for (i in 1:length(unique(gene_module_df[["cell_cluster"]]))) {
  # find the appropriate expression data and subset it 

  cell_cluster <- sort(unique(gene_module_df[["cell_cluster"]]))[i]
  message(paste0("starting cell cluster ", i, ": ", cell_cluster))
  # reference sets
  # format: list(datExpr_ID1=c('n01','n02','n03'))
  ID_datExpr_ref <- names(ref_sets)
  colname_subset_ref <- metadata_ident_cols[ID_datExpr_ref]
  pattern_subset_ref <- ref_sets[[1]][i]

  metadata_ref <- list_metadata[[ID_datExpr_ref]]
  idx_subset_ref <- which(metadata_ref[[colname_subset_ref]] == unique(grep(paste0("^", pattern_subset_ref, "$"), metadata_ref[[colname_subset_ref]], value=T)))
  datExpr_ref <- list_datExpr[[ID_datExpr_ref]][idx_subset_ref,]

  names(datExpr_ref) <- cell_cluster
  # Check that there are sufficient cells
  if (nrow(datExpr_ref) < minCellClusterSize) {
    list_modulePreservation_out[[cell_cluster]] <- NA_character_
    message(paste0("cellcluster ", i, ": ", cell_cluster, " skipped due to having fewer than ", minCellClusterSize, " cells"))
    next
  }

  # Prepare test sets
  # format: list('n01'=c(datExpr_ID1='n02',datExpr_ID1='n03'), 'n04'=c(datExpr_ID1='n01',datExpr_ID1='n03'),'n03'=c(datExpr_ID1='n01',datExpr_ID1='n02'))
  IDs_datExpr_test <- names(test_sets[[i]])
  colnames_subset_test <- metadata_ident_cols[IDs_datExpr_test]
  patterns_subset_test <- test_sets[[i]]
  
  list_metadata_test <- list_metadata[IDs_datExpr_test]
  list_idx_subset_test <- mapply(function(metadata_test, colname_subset_test, pattern_subset_test) {
    which(metadata_test[[colname_subset_test]] == unique(grep(paste0("^", pattern_subset_test, "$"), metadata_test[[colname_subset_test]], value=T)))
  }, metadata_test = list_metadata_test, colname_subset_test=colnames_subset_test, pattern_subset_test=patterns_subset_test, SIMPLIFY=F)
  
  list_datExpr_test <- list_datExpr[IDs_datExpr_test]
  
  list_datExpr_test<- mapply(function(datExpr_test, idx_subset_test) {
    datExpr_test[idx_subset_test,]
  }, datExpr_test= list_datExpr_test, idx_subset_test = list_idx_subset_test, SIMPLIFY=F)
  
  names(list_datExpr_test) <- patterns_subset_test
  # Check that there are sufficient cells
  list_datExpr_test <- lapply(list_datExpr_test,function(datExpr_test) {
    return( if (nrow(datExpr_test) >= minCellClusterSize) datExpr_test else NULL)
  })
  
  # check test sets
  idx_test_NULL <- sapply(X=list_datExpr_test, FUN=is.null,simplify=T)
  if (sum(idx_test_NULL)==length(list_datExpr_test)) {
    list_modulePreservation_out[[cell_cluster]] <- NA_character_
    next
  } else if (sum(idx_test_NULL)>1 & sum(idx_test_NULL)<length(list_datExpr_test)) {
    list_datExpr_test <- list_datExpr_test[!idx_test_NULL]
    message(paste0(patterns_subset_test[idx_test_NULL], " dropped as test sets due to having fewer than ", minCellClusterSize, " cells"))
  }
  
  ####################################################
  #### Get gene network coloring for cell_cluster ####
  ####################################################

  coloring <- gene_module_df[["module"]][gene_module_df[["cell_cluster"]] == cell_cluster]

  names(coloring) <- make.unique(gene_module_df[[gene_col]][gene_module_df[["cell_cluster"]] == cell_cluster])
  
  # Make the coloring vector the same length as the number of genes in the dataset
  colors_tmp <- rep("grey", length=ncol(datExpr_ref))
  names(colors_tmp) <- colnames(datExpr_ref)
  colors_tmp[as.integer(na.omit(match(names(coloring), names(colors_tmp))))] <- coloring[as.integer(na.omit(match(names(colors_tmp), names(coloring))))]
  coloring <- colors_tmp
  
  ####################################################
  ############# Run modulePreservation ###############
  ####################################################

  # Prepare multidata
  multiData <- vector(mode="list",length = 1+length(list_datExpr_test)) # 1 for reference network
  multiData[[1]] <- list("data"=datExpr_ref)
  for (i in 1:(length(multiData)-1)){
    multiData[[i+1]] <- list("data"=list_datExpr_test[[i]])
  }
  names(multiData)[1] <- cell_cluster
  names(multiData)[2:length(multiData)] <- names(list_datExpr_test)
    
  # Check the multidata samples 
  goodSamplesGenesMS_out <- goodSamplesGenesMS(
                              multiExpr=multiData,
                              multiWeights = NULL,
                              minFraction = 1/500,
                              minNSamples = 5,
                              minNGenes = 100,
                              tol = NULL,
                              minRelativeWeight = 0.01,
                              verbose = 3,
                              indent = 0)

  multiData <- mapply(function(dataset, idx_goodSamples) {
    dataset$data <- dataset$data[idx_goodSamples, goodSamplesGenesMS_out$goodGenes]
    }, dataset=multiData, idx_goodSamples=goodSamplesGenesMS_out$goodSamples)
  
  
  # Prepare multiColor
  multiColor <- list(coloring)
  names(multiColor) <- cell_cluster
    
  # Prepare for parallel computation (WGCNA multi threads)
  
  require("doParallel")
  additional_Gb = max(as.numeric(sapply(multiData, FUN = function(x) object.size(x), simplify = T)))/1024^3
  obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
  n_cores <- max(1, min(detectCores() %/% 3, RAM_Gb_max %/% (obj_size_Gb + additional_Gb))-1)
  n_cores <- min(n_cores, 40)

  enableWGCNAThreads(nThreads = n_cores)
  #WGCNA::disableWGCNAThreads()
  setwd(dir = dir_scratch)

  list_modulePreservation_out[[cell_cluster]] <- modulePreservation(multiData=multiData,
                                                                multiColor=multiColor,
                                                                dataIsExpr = TRUE,
                                                                networkType = networkType, 
                                                                corFnc = corFnc,
                                                                corOptions =NULL, #if (corFnc == "cor") list(use = 'p') else NULL,
                                                                referenceNetworks = 1, 
                                                                testNetworks = NULL,
                                                                nPermutations = 100, 
                                                                includekMEallInSummary = FALSE,
                                                                restrictSummaryForGeneralNetworks = TRUE,
                                                                calculateQvalue = FALSE,
                                                                randomSeed = 12345, 
                                                                maxGoldModuleSize = 1000, 
                                                                maxModuleSize = 1000, 
                                                                quickCor = 0, # set to one to speed up/risk errors and inaccuracies 
                                                                ccTupletSize = 2, 
                                                                calculateCor.kIMall = FALSE,
                                                                calculateClusterCoeff = T,
                                                                useInterpolation = FALSE, 
                                                                checkData = TRUE, 
                                                                greyName = NULL, 
                                                                savePermutedStatistics = TRUE, 
                                                                loadPermutedStatistics = FALSE, 
                                                                permutedStatisticsFile = "permutedStats-actualModules.RData",#if (useInterpolation) "permutedStats-intrModules.RData" else "permutedStats-actualModules.RData", 
                                                                plotInterpolation = TRUE, 
                                                                interpolationPlotFile = "modulePreservationInterpolationPlots.pdf", 
                                                                discardInvalidOutput = TRUE,
                                                                parallelCalculation = TRUE,
                                                                verbose = 3, 
                                                                indent = 0)
  
}

######################################################################
################### COMPUTE MODULE CORRELATIONS ######################
######################################################################

 # Get kMs


  
  
  
  
  
  
  
  