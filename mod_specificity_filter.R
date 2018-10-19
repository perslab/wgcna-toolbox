# module specificity filter
# 

# Usage: 
# export R_MAX_NUM_DLLS=999
# time RScript..

# TODO

######################################################################
######################### UTILITY FUNCTIONS ##########################
######################################################################

source(file="/projects/jonatan/tools/functions-src/utility_functions.R")

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("optparse", "Matrix", "ggplot2", "dplyr", "parallel", "pSI"))#, "loomR", "doubletFinder")

######################################################################
########################### OptParse #################################
######################################################################

option_list <- list(
  make_option("--path_typeModEmbedMat", type="character", 
              help = "Path to celltype module embedding matrix. Embedding may be expression or preservation or any other measure"),  
  make_option("--dir_out", type="character",
              help = "Project directory to which to write files. Should include subdirectories /tables, /RObjects, /plots, /log"),
  make_option("--n_cores", type="integer", default = 5L,
              help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]")
)

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))
path_typeModEmbedMat <- opt$path_typeModEmbedMat
dir_out <- opt$dir_out
n_cores <- opt$n_cores

######################################################################
########################### SET OPTIONS (DEV) ########################
######################################################################

if (FALSE) {
  path_typeModEmbedMat = ""
  dir_out = ""
  n_cores = 5
}

######################################################################
########################### SET OPTIONS ##############################
######################################################################

options(stringsAsFactors = F)

######################################################################
############################ CONSTANTS ###############################
######################################################################

# if specified output directory doesn't exist, create it 
if (!file.exists(dir_out)) {
  dir.create(dir_out) 
  message("Project directory not found, new one created")
}

dir_plots = paste0(dir_out,"plots/")
if (!file.exists(dir_plots)) dir.create(dir_plots) 

dir_tables = paste0(dir_out,"tables/")
if (!file.exists(dir_tables)) dir.create(dir_tables)

dir_RObjects = paste0(dir_out,"RObjects/")
if (!file.exists(dir_RObjects)) dir.create(dir_RObjects)

dir_log = paste0(dir_out,"log/")
if (!file.exists(dir_log)) dir.create(dir_log)

dir_scratch = "/scratch/tmp-wgcna/"

dir_current = paste0(LocationOfThisScript(), "/")

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

######################################################################
############################ LOAD DATA ###############################
######################################################################

typeModEmbedMat <- load_obj(path_typeModEmbedMat)

######################################################################
############################ LOAD DATA ###############################
######################################################################

# TODO
spec_index <- specificity.index()

