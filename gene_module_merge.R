# Cross-dataset gene network correlation and preservation analysis pipeline
# Usage: e.g.
#
#

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  # inputs
  make_option("--dirs_project_WGCNA", type="character",
              help = "Directories of WGCNA projects, each must contain /tables and /RObjects subdirectories."), 
  make_option("--prefixes_WGCNA_run", type="character", default = NULL,
              help = "Search for matching files within the dirs_project_WGCNA, if left as NULL will take all runs, [default %default]"),
  make_option("--dir_out", type="character", 
              help = "Directory for outputs, must contain /tables and /RObjects subdirectories"),
  make_option("--prefix_out", type="character", 
              help = "Prefix for output files"),
  make_option("--corr_space", type="character", default="cell",
              help = "Measure correlation in cell embedding ('cell') or gene loading ('gene') space, where the correlation takes the union of genes? [default %default]"),
  make_option("--path_data", type="character", default = NULL,
              help = "If corr_space == cell, path to dataset with which to find module cell embeddings to measure module correlation in cell space."), 
  make_option("--gene_weights_recompute", type="character", default="kIM",
              help = "How should the script recompute gene weights in merged modules? One of 'union_normalise': take the union of gene weights while keeping the highest weights within the intersect, then re-normalize, 'kIM': compute kIMs in the dataset given as path_data, 'kME': compute module kMEs in in the dataset [default %default]"),
  make_option("--moduleMergeCutHeight", type="numeric", default=0.1,
              help = "Cut-off level for merging modules based on pairwise distance. [default %default]"),
  make_option("--n_cores", type="integer", default = 10L,
              help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]")
)

######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

if (FALSE) { 
  dirs_project_WGCNA = c("/projects/jonatan/tmp-maca/", "/projects/jonatan/tmp-mousebrain/")
  prefixes_WGCNA_run = c("")
  dir_out = "/projects/jonatan/tmp-rwgcna-tests/module_merge_test/"
  prefix_out = "mod_merge_test"
  corr_space = "cell"
  path_data = "/projects/jonatan/tmp-mousebrain/mousebrain_hypoth.RDS"
  gene_weights_recompute = "kIM"
  moduleMergeCutHeight = 0.1
  n_cores = 10
}

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

dirs_project_WGCNA <- opt$dirs_project_WGCNA 
if (!is.null(dirs_project_WGCNA)) dirs_project_WGCNA <- eval(parse(text=dirs_project_WGCNA))
prefixes_WGCNA_run <- opt$prefixes_WGCNA_run
if (!is.null(prefixes_WGCNA_run)) prefixes_WGCNA_run <- eval(parse(text=prefixes_WGCNA_run))
path_data <- opt$path_data

######################################################################
############################## PACKAGES #############################
######################################################################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))

######################################################################
########################## SOURCE FUNCTIONS ##########################
######################################################################

source(file = "/projects/jonatan/wgcna-src/wgcna-toolbox/rwgcna_functions.R")

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F)

######################################################################
###################### SET AND VERIFY PATHS ##########################
######################################################################

if (!file.exists(dirs_project_WGCNA)) stop("dir_project_WGCNA not found")
if (!file.exists(path_data)) stop("path_data not found")

# set up and verify WGCNA project directory paths

dirs_RObjects_WGCNA = paste0(dirs_project_WGCNA,"RObjects/")
if (any(sapply(dirs_RObjects_WGCNA, function(dir_RObjects_WGCNA) !file.exists(dir_RObjects_WGCNA), simplify = T))) stop("Each dirs_RObjects_WGCNA must contain a /RObjects subdir")

# set up output directories 

dir_out_plots = paste0(dir_out,"plots/")
if (!file.exists(dir_out_plots)) dir.create(dir_out_plots) 

dir_out_tables = paste0(dir_out,"tables/")
if (!file.exists(dir_out_tables)) dir.create(dir_out_tables)

dir_out_RObjects = paste0(dir_out,"RObjects/")
if (!file.exists(dir_out_RObjects)) dir.create(dir_out_RObjects)

dir_out_log = paste0(dir_out, "log/")
if (!file.exists(dir_out_log)) dir.create(dir_out_log)

dir_scratch = "/scratch/tmp-wgcna/"

WGCNA_projects <- gsub(".*/", "", sapply(dirs_project_WGCNA, function(dir_proj) substr(dir_proj, start=1, stop = nchar(dir_proj)-1)))

######################################################################
######################### LOAD GENE LOADINGS VECTORS #################
######################################################################

message("Loading module gene loading vectors")

WGCNAproj_run_cellType_module_u <- list()

for (i in 1:length(dirs_WGCNA_RObjects)) {
  run_cellType_module_u_paths <- dir(path = dirs_RObjects_WGCNA[i], pattern = paste0(prefixes_WGCNA_run, "_list_list_module_u.RDS", collapse = "|"), full.names = T)
  for (cellType_module_u_path in run_cellType_module_u_paths) {
    cellType_module_u_name <- gsub(".*/|_list_list_module_u.RDS", "", cellType_module_u_path)
    WGCNAproj_cellType_module_u[[WGCNA_projects[i]]][[cellType_module_u_name]] <- load_obj(f = cellType_module_u_path)
  }
  names(WGCNAproj_run_cellType_module_u[[i]]) <- gsub(".*/|_list_list_module_u.RDS", "", run_cellType_module_u_paths)
}

names(WGCNAproj_run_cellType_module_u) = WGCNA_projects

######################################################################
######################### LOAD GENE MAPPING ##########################
######################################################################

message("Loading gene mapping files")

WGCNAproj_run_cellType_mapping <- list()

for (i in 1:length(dirs_WGCNA_RObjects)) {
  run_cellType_mapping_paths <- dir(path = dirs_tables_WGCNA[i], pattern = paste0(prefixes_WGCNA_run, ".*_hgnc_to_ensembl_mapping_df.csv", collapse = "|"), full.names = T)
  for (cellType_module_u_path in run_cellType_module_u_paths) {
    cellType_module_u_name <- gsub(".*/|_hgnc_to_ensembl_mapping_df.csv", "", cellType_module_u_path)
    WGCNAproj_cellType_module_u[[WGCNA_projects[i]]][[cellType_module_u_name]] <- load_obj(f = cellType_module_u_path)
  }
  names(WGCNAproj_run_cellType_mapping[[i]]) <- gsub(".*/|_hgnc_to_ensembl_mapping_df.csv", "", run_cellType_mapping_paths)
}

names(WGCNAproj_run_cellType_mapping) <- WGCNA_projects

######################################################################
####################### NORMALIZE GENE LOADING VECTORS #################
######################################################################

if (gene_weights_recompute %in% c("kIM", "kME")) {
  
}
######################################################################
####################### RESCALE GENE LOADING VECTORS #################
######################################################################


# 1. correlation and merging networks
#     In: set of vectors of gene weights. May be overlapping.
#     Usage:  
#       i. scale each vector to sum to 1
#       ii. correlation
#       iii. hiearchical clustering and cutting
#       iv. Iterate: merge close modules
#           Take union of genes.
#            Either: 
#             Recompute gene weights in dataset using kME or kIM
#            Or:
#             Re-normalise the union so they still sum to 1.
#     Out: reduced set of vectors of gene weights