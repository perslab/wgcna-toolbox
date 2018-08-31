# Make cellType x module embedding matrix
# Use gene loadings from one data (sub)set and possibly expression matrix from another

######################################################################
######################### LOAD LIBRARIES #############################
######################################################################

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(parallel))

######################################################################
########################## SET PARAMS ################################
######################################################################

source(file = "/projects/jonatan/wgcna-src/wgcna-toolbox/rwgcna_params.R")
source(file = "/projects/jonatan/wgcna-src/wgcna-toolbox/rwgcna_functions.R")
options(stringsAsFactors = F)

######################################################################
######################### SET CONSTANTS ##############################
######################################################################
seurat_path <- "/projects/jonatan/tmp-mousebrain/RObjects/L5.RDS"
project_dir <- "/projects/jonatan/tmp-mousebrain/"
data_prefix <- "mousebrain"
vec_wgcna_run_prefix <- c("Neurons_ClusterName_2")
n_cores = 20
metadata_subset_col <- "ClusterName"

# if specified output directory doesn't exist, create it 
if (!file.exists(project_dir)) {
  dir.create(project_dir) 
  message("Project directory not found, new one created")
}

plots_dir = paste0(project_dir,"plots/")
if (!file.exists(plots_dir)) dir.create(plots_dir) 

tables_dir = paste0(project_dir,"tables/")
if (!file.exists(tables_dir)) dir.create(tables_dir)

RObjects_dir = paste0(project_dir,"RObjects/")
if (!file.exists(RObjects_dir)) dir.create(RObjects_dir)

log_dir = paste0(project_dir,"log/")
if (!file.exists(log_dir)) dir.create(log_dir)

scratch_dir = "/scratch/tmp-wgcna/"

######################################################################
######################### LOAD GENE LOADINGS VECTORS #################
######################################################################

run_cellType_module_u <- list()
for (run in vec_wgcna_run_prefix) {
  run_cellType_module_u[[run]] <- load_obj(f = sprintf("%s%s_%s_list_list_module_u.RDS", RObjects_dir, data_prefix, run))
}

cellType_module_u <- unlist(x = run_cellType_module_u, recursive = F, use.names = T) # TODO: Does this give names as desired?
module_u <- unlist(x = cellType_module_u, recursive = F, use.names = T)
names(module_u) <- paste0(data_prefix, ".", names(module_u))

######################################################################
######################### LOAD EXPRESSION MATRIX #####################
######################################################################

# load/open connection to expression matrix
if (grepl(pattern = "\\.loom", seurat_path)) {
  stop("Stop: at present script only takes as input a Seurat object in the form of an R Object with ScaleData done")
  # datExpr <- connect(filename = datExpr_mat_path, mode = "r+")
  # fileType <- "loom"
} else if (grepl(pattern = "\\.RData|\\.Rdata|\\.RDS|\\.rds", seurat_path)) {
  data_obj <- load_obj(seurat_path)
  fileType = "RObject"
}

######################################################################
######################### PREPARE DATA ###############################
######################################################################

if (fileType=="RObject") {
  
  datExpr <- data_obj@scale.data 
  datExpr_colnames <- colnames(datExpr)
  
  if (!is.null(metadata_subset_col)) {
    ident = as.character(data_obj@meta.data[[metadata_subset_col]])
    names(ident) = rownames(data_obj@meta.data)
    ident <- ident[match(names(ident), colnames(datExpr))]
  } else {
    ident <- as.character(data_obj@ident)
    ident <- ident[match(names(ident), colnames(datExpr))]
  }
  rm(data_obj)
} else {
  # TODO for loom
}

######################################################################
################ COMPUTE CELL MODULE EMBEDDINGS ######################
######################################################################

cell_chunk_idx <- c(seq(from=0, to = ncol(datExpr), by=1000), ncol(datExpr))

cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(log_dir, data_prefix, "_", run_prefix, "_", "cell_x_module_matrix.txt"))

for (i in 1:(length(cell_chunk_idx)-1)) {
  submat <- datExpr[, (cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]
  #Todo: does submat inherit names?
  subEmbed <- parSapply(cl, X=module_u, FUN=function(u) {
    submat[,match(names(u), colnames(submat))] %*% u
  }, 
  simplify=T) # should give a 1000 * n_module embedding matrix
  
  if (i==1) cellModEmbed <- subEmbed else cellModEmbed <- rbind(cellModEmbed, subEmbed)
}

rownames(cellModEmbed) = datExpr_colnames; colnames(cellModEmbed) <- names(module_u)

saveRDS(cellModEmbed, file=paste0(RObjects_dir, data_prefix, "_cellModEmbed_mat_full.RDS"))

######################################################################
################## AVERAGE EMBEDDINGS FOR EACH CELLTYPE ##############
######################################################################

for (i in 1:length(sort(unique(ident)))) {
  subEmbed <- cellModEmbed[ident==sort(unique(ident))[i],]
  mean <- parSapply(cl, X=subEmbed, FUN=function(column) {
    mean(column)
  }, simplify=T) # it might be ok to do this w/o parallelising
  if (i==1) cellTypeModEmbed <- subEmbed else cellTypeModEmbed <- rbind(cellTypeModEmbed, subEmbed) 
}

stopCluster(cl)

rownames(cellTypeModEmbed) = sort(unique(ident)); colnames(cellTypeModEmbed) <- names(module_u)

saveRDS(cellTypeModEmbed, file=paste0(RObjects_dir, data_prefix, "_cellTypeModEmbed_mat_full.RDS"))

######################################################################
######################## PLOT BOTH HEATMAPS ##########################
######################################################################
# TODO
