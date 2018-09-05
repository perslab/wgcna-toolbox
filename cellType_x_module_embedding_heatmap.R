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
suppressPackageStartupMessages(library(loomR))

######################################################################
########################## SOURCE FUNCTIONS ##########################
######################################################################

source(file = "/projects/jonatan/wgcna-src/wgcna-toolbox/rwgcna_functions.R")

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F)

######################################################################
######################### SET CONSTANTS ##############################
######################################################################
data_path <- "/projects/jonatan/tmp-mousebrain/RObjects/L5.RDS"
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

######################################################################
######################### LOAD GENE MAPPING ##########################
######################################################################

run_cellType_mapping <- list()
for (run in vec_wgcna_run_prefix) {
  path_mapping <- dir(pattern = paste0(data_prefix, "_", run, ".*_hgnc_to_ensembl_mapping"), path = tables_dir, full.names = T)
  run_cellType_mapping[[run]] <- read.csv(file = path_mapping, header = T, quote = "")
}

######################################################################
################### MAP U VECTOR GENE NAMES TO SYMBOL ################
######################################################################

cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(log_dir, data_prefix,"_", "map_u_from_ensembl_to_hgnc.txt"))

run_cellType_module_u <- mapply(function(cellType_module_u, mapping) {
  parLapplyLB(cl, cellType_module_u, function(module_u) {
    lapply(module_u, function(u) {
      names(u) <- mapping$symbol[match(names(u), mapping$ensembl, nomatch = 0)]
      u  
    })
  })
},  cellType_module_u = run_cellType_module_u,
mapping = run_cellType_mapping,
SIMPLIFY=F)

stopCluster(cl)

cellType_module_u <- unlist(x = run_cellType_module_u, recursive = F, use.names = T) # TODO: Does this give names as desired?
module_u <- unlist(x = cellType_module_u, recursive = F, use.names = T)
names(module_u) <- paste0(data_prefix, ".", names(module_u))

######################################################################
######################### LOAD EXPRESSION MATRIX #####################
######################################################################
# Leave this as late as possible to avoid hogging memory..

# load/open connection to expression matrix
if (grepl(pattern = "\\.loom", data_path)) {
  data_obj <-  connect(filename = data_path, mode = "r+")
  fileType <- "loom"
} else if (grepl(pattern = "\\.RData|\\.Rdata|\\.RDS|\\.rds", data_path)) {
  data_obj <- load_obj(data_path)
  fileType = "RObject"
}

######################################################################
####################### COUNT PERCENT RIBO, MITO #####################
######################################################################
# todo: make sure row/column is correct: matrix is transposed
if (fileType =="loom") {
  mito.genes_idx <- grepl(pattern = "^mt-", x =data_obj$row.attrs$gene_names[], ignore.case=T)
  ribo.genes_idx <- grepl(pattern = "^Rp[sl][[:digit:]]", data_obj$row.attrs$gene_names[], ignore.case=T)
  
  # TODO
  
  nUMI_map <- data_obj$apply(FUN = function(x) {
    rowSums(x[,mito.genes_idx])
  }, MARGIN = 2, 
  chunk.size = 500, 
  dataset.use = "matrix", 
  display.progress = FALSE)
  
  
  percent.mito <- Matrix::rowSums(data_obj[["matrix"]][,mito.genes_idx])/Matrix::rowSums(data_obj[["matrix"]][,])
  percent.ribo <- Matrix::rowSums(data_obj[["matrix"]][,ribo.genes_idx])/Matrix::rowSums(data_obj[["matrix"]][,])
  data_obj <- lfile$add.col.attribute(list(percent.mito=percent.mito), overwrite = TRUE)
  data_obj <- lfile$add.col.attribute(list(percent.ribo=percent.ribo), overwrite = TRUE)
  
} else if (fileType == "RObject") {
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = data_obj@data), value = TRUE, ignore.case=T)
  ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = data_obj@data), value = TRUE)
  percent.mito <- Matrix::colSums(data_obj@raw.data[mito.genes, ])/Matrix::colSums(data_obj@raw.data)
  percent.ribo <- Matrix::colSums(data_obj@raw.data[ribo.genes, ])/Matrix::colSums(data_obj@raw.data)
  data_obj <- AddMetaData(object = data_obj, metadata = percent.mito, col.name = "percent.mito")
  data_obj <- AddMetaData(object = data_obj, metadata = percent.ribo, col.name = "percent.ribo")
}

######################################################################
######################### PREPARE DATA ###############################
######################################################################

if (fileType=="RObject") {

  ###################### NORMALIZE AND SCALE DATA ######################

  data_obj <- NormalizeData(object=data_obj, display.progress = T)
  data_obj <- ScaleData(object = data_obj, vars.to.regress = c("percent.mito", "percent.ribo", "nUMI"), do.scale = T, do.center=T, display.progress = T, do.par = T, num.cores = 40)
  datExpr <- data_obj@scale.data 
  datExpr_colnames <- colnames(datExpr)
  saveRDS(datExpr, file=paste0(RObjects_dir, data_prefix, "full_scale_data.RDS"))
  
  if (!is.null(metadata_subset_col)) {
    ident = as.character(data_obj@meta.data[[metadata_subset_col]])
    names(ident) = rownames(data_obj@meta.data)
    ident <- ident[match(names(ident), colnames(datExpr))]
  } else {
    ident <- as.character(data_obj@ident)
    ident <- ident[match(names(ident), colnames(datExpr))]
  }
  rm(data_obj)
} else if (fileType=="loom") {
  
  NormalizeData(object = data_obj, overwrite = TRUE, display.progress = FALSE)
  data_obj$apply(name = "col_attrs/umi_apply", FUN = rowSums, MARGIN = 2, chunk.size = 500, 
              dataset.use = "matrix", display.progress = FALSE, overwrite = TRUE)
  
  # TODO for loom
}


######################################################################
################ COMPUTE CELL MODULE EMBEDDINGS ######################
######################################################################

cell_chunk_idx <- c(seq(from=0, to = ncol(datExpr), by=10000), ncol(datExpr))

cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(log_dir, data_prefix,"_", "cell_x_module_matrix.txt"))

for (i in 1:(length(cell_chunk_idx)-1)) { # loop over chunks of cell
  submat <- datExpr[, (cell_chunk_idx[i]+1):cell_chunk_idx[i+1]] #inherits big matrix names
  subEmbed <- parSapply(cl, X=module_u, FUN=function(u) {
    t(submat[match(names(u), rownames(submat)),]) %*% matrix(u, nrow = length(u))
  }, 
  simplify=T) # should give a 1000 * n_module embedding matrix
  
  if (i==1) cellModEmbed <- subEmbed else cellModEmbed <- rbind(cellModEmbed, subEmbed)
}

rownames(cellModEmbed) = datExpr_colnames; colnames(cellModEmbed) <- names(module_u)

saveRDS(cellModEmbed, file=paste0(scratch_dir, data_prefix, "_cellModEmbed_mat_full.RDS"))

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

######################################################################
############################# WRAP UP ################################
######################################################################

if (fileType == "loom") data_obj$close_all()