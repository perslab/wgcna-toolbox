# Make cellType x module embedding matrix
# Use gene loadings from one data (sub)set and possibly expression matrix from another

# loomR:
# Introduction to loomR: https://satijalab.org/loomR/loomR_tutorial.html
# loomR MCA tutorial: https://satijalab.org/seurat/mca_loom.html
# loomR full specs: https://satijalab.org/loomR/loomR.pdf
# loom file specification: http://linnarssonlab.org/loompy/format/index.html

# usage: 
# export R_MAX_NUM_DLLS=999
# time Rscript /projects/jonatan/tools/wgcna-src/wgcna-toolbox/cell_module_heatmap.R --dir_project_WGCNA /projects/jonatan/tmp-mousebrain/ --dir_project_data /projects/jonatan/tmp-mousebrain/ --path_data /projects/jonatan/tmp-mousebrain/RObjects/L5.RDS --prefixes_WGCNA_run 'c("Neurons_ClusterName_2")' --prefix_out mousebrain_neuronMods_1  --scale_data F --n_cores 20

######################################################################
########################## DEFINE FUNCTIONS ##########################
######################################################################

source(file="/projects/jonatan/tools/functions-src/utility_functions.R")

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("optparse", "Seurat", "dplyr", "Biobase", "Matrix", "parallel", "loomR", "ComplexHeatmap"))

######################################################################
########################### OPTPARSE #################################
######################################################################

option_list <- list(
  make_option("--dir_project_WGCNA", type="character",
              help = "Folder of WGCNA project, must contain /tables and /RObjects subdirectories."),  
  make_option("--dir_project_data", type="character",
              help = "Folder of expression data project. The script sends output here. must contain /tables and /RObjects subdirectories, may be same as dir_project_WGCNA"),  
  make_option("--path_data", type="character",
              help = "Path to expression data in Seurat loomR or RObject format to score on WGCNA modules."),  
  make_option("--path_metadata", type="character", default = NULL,
              help = "Path to metadata in one of the standard (compressed) character separated formats. If not given, takes any metadata stored within the path_data object. [default %default]"),  
  make_option("--prefixes_WGCNA_run", type="character",
              help = "WGCNA run prefixes given as vector in double (single) quotes of characters in single (double) quotes. Used to search for WGCNA outputs from dir_project_WGCNA subdirs"), 
  make_option("--prefix_out", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"), 
  make_option("--datExpr_ident_cols", type="character", default=NULL,
              help = "Specify one or more columns in metadata with cell identity labels. Used to chunk the analysis and annotate the heatmap. Format as character vector, with additional outer quotation marks and no whitespace ''c('tissue','celltype')''. If NULL use the seurat_object@ident slot. [default %default]"),
  make_option("--WGCNA_ident_cats", type="character", default=NULL,
              help = "Specify one or more level of cell identity labels for the WGCNA analyses to identify in the filenames. This is in addition to the WGCNA run and cell subset that the analysis was done on. Used to annotate the heatmap. Takes a named list, in additional outer quotation marks, of identity levels, with no whitespace, e.g. ''list(tissue=c('cortex','hindbrain','midbrain'),cellType=c('neuron','OPC','microglia','astrocyte'))''. (case insensitive) Unlike for the expression data which is embedded in the WGCNA modules, all levels must be specified [default %default]"),
  make_option("--scale_data", type="logical", default = "TRUE",
              help = "[default %default]. If TRUE, and path_data points to an RObject, the script will try to use the @scale.data slot and, if missing, will normalize and scale the data and regress out nUMI, percent.mito and percent.ribo. If FALSE, for RObjects will use the @data slot, or, if empty, the @raw.data slot"),  
  make_option("--do_plot", type="logical", default = "TRUE",
              help = "[default %default]. If FALSE, script outputs embedding matrices without plotting them as heatmaps"),  
  make_option("--plot_PPI_enrichment", type="logical", default = "F",
              help = "Add module PPI enrichment to plots, using the rWGCNA pipeline outputs [default %default]"), 
  make_option("--plot_magma_enrichment", type="character", default = NULL,
              help = "Add module gwas enrichments to plots? Takes a vector, in quotes and without whitespace, of (parts of) rWGCNA gwas column names (case insensitive). [default %default]"), 
  make_option("--plot_LDscore_enrichment", type="logical", default = "F",
              help = "Add LDscore enrichment to plots, [default %default]"), 
  make_option("--plot_metadata_corr", type="character", default = NULL,
              help = "Add metadata correlation coefficients to plots? Takes a vector, in quotes and without whitespace, of (parts of) rWGCNA metadata correlation traits (case insensitive). [default %default]"), 
  make_option("--n_cores", type="integer", default = 15L,
              help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]")
  # make_option("--use_loom", type="logical", default = F,
  #             help = "Use loom HF5 format? Less RAM-heavy but error-prone [default %default]")
)

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

dir_project_WGCNA <- opt$dir_project_WGCNA 
dir_project_data <- opt$dir_project_data
path_data <- opt$path_data
path_metadata <- opt$path_metadata
prefixes_WGCNA_run <- eval(parse(text=opt$prefixes_WGCNA_run))
prefix_out <- opt$prefix_out 
if (is.null(opt$datExpr_ident_cols)) datExpr_ident_cols <- eval(parse(text=opt$datExpr_ident_cols))
if (is.null(opt$WGCNA_ident_cats)) WGCNA_ident_cats <- eval(parse(text=opt$WGCNA_ident_cats))
scale_data <- opt$scale_data
do_plot <- opt$do_plot
plot_PPI_enrichment <- opt$plot_PPI_enrichment
if (is.null(opt$plot_LSscore_enrichment)) plot_LDscore_enrichment <- eval(parse(text=opt$plot_LSscore_enrichment))
if (is.null(opt$plot_magma_enrichment)) plot_magma_enrichment <- eval(parse(text=opt$plot_magma_enrichment))
if (is.null(opt$plot_metadata_corr)) plot_metadata_corr <- eval(parse(text=opt$plot_metadata_corr))
n_cores <- opt$n_cores
use_loom <- F#opt$use_loom # For the initial Seurat processing, before computing embeddings. The script will still use loomR to compute the embedding matrix

######################################################################
######################### SET OPTIONS (DEV) ##########################
######################################################################

if (FALSE) { 
  dir_project_WGCNA <- "/projects/jonatan/tmp-epilepsy/" 
  path_data <- "/projects/jonatan/tmp-epilepsy/RObjects/EX_so2_QC.RDS" #"/data/pub-others/zeisel-biorxiv-2018/data/l5_all.loom"
  path_metadata <- NULL
  dir_project_data <- "/projects/jonatan/tmp-epilepsy/"
  prefixes_WGCNA_run <- c("EP_UPF_5", "EP_lake_5")
  prefix_out <- "EP_heatmap"
  datExpr_ident_cols <- c("dataset","res.0.6") #c("tissue","cell_ontology_class")
  WGCNA_ident_cats <- NULL 
  scale_data <- T
  do_plot = T
  plot_PPI_enrichment = F
  plot_variant_enrichment = F
  plot_magma_enrichment = NULL
  plot_LDscore_enrichment = F
  plot_metadata_corr = F 
  n_cores = 5
  #use_loom = T
}

if (FALSE) { 
  dir_project_WGCNA <- "/projects/jonatan/tmp-mousebrain/" 
  path_data <- "/projects/jonatan/tmp-mousebrain/RObjects/L5.RDS" #"/data/pub-others/zeisel-biorxiv-2018/data/l5_all.loom"
  path_metadata <- NULL
  dir_project_data <- "/projects/jonatan/tmp-mousebrain/"
  prefixes_WGCNA_run <- c("Astrocytes_ClusterName_2","Ependymal_ClusterName_2","Immune_ClusterName_2","Neurons_ClusterName_2","Oligos_ClusterName_2","PeripheralGlia_ClusterName_2","Vascular_ClusterName_2")
  prefix_out <- "mousebrain_mod_cell_kME_1"
  datExpr_ident_cols <- #c("tissue","cell_ontology_class")
  WGCNA_ident_cats <- NULL 
  scale_data <- T
  do_plot = T
  plot_PPI_enrichment = F
  plot_variant_enrichment = F
  plot_magma_enrichment = NULL
  plot_LDscore_enrichment = F
  plot_metadata_corr = F 
  n_cores = 30
  #use_loom = T
}

if (FALSE) { 
  dir_project_WGCNA <-"/projects/jonatan/tools/tmp-rwgcna-tests/maca-pancreas-test/"
  path_data <- "/projects/jonatan/tmp-maca/RObjects/maca_seurat_pancreas.Rdata" # "/data/pub-others/zeisel-biorxiv-2018/data/l5_all.loom"
  path_metadata <- NULL
  dir_project_data <- "/projects/jonatan/tmp-maca/"
  prefixes_WGCNA_run <- c("test_corr_f_v")
  prefix_out <- "maca_pancreas_heatmap_test_1"
  datExpr_ident_cols <- #c("tissue","cell_ontology_class")
  WGCNA_ident_cats <- NULL 
  datExpr_ident_cols <- c("tissue","cell_ontology_class")
  scale_data <- T
  do_plot = T
  plot_PPI_enrichment = T
  plot_magma_enrichment = NULL
  plot_LDscore_enrichment = NULL 
  plot_metadata_corr = NULL 
  n_cores = 7
  #use_loom = T
}

if (FALSE) { 
  dir_project_WGCNA <-"/projects/jonatan/tmp-maca/"
  path_data <- "/data/pub-others/tabula_muris/figshare/180126-facs/maca.seurat_obj.facs.figshare_180126.RData" # "/data/pub-others/zeisel-biorxiv-2018/data/l5_all.loom"
  path_metadata <- NULL
  dir_project_data <- "/projects/jonatan/tmp-maca/"
  prefixes_WGCNA_run <- c("tissue_cell_type_kIM")
  prefix_out <- "maca_kIM_embed_1"
  datExpr_ident_cols <- c("tissue","tissue_cell_type")
  WGCNA_ident_cats <- c("tissue", "tissue_cell_type")
  scale_data <- T
  do_plot = T
  plot_PPI_enrichment = T
  plot_magma_enrichment = T
  plot_LDscore_enrichment = F 
  plot_metadata_corr = F 
  n_cores = 15
  #use_loom = T
}

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F)

######################################################################
###################### SET AND VERIFY PATHS ##########################
######################################################################

if (!file.exists(dir_project_WGCNA)) stop("dir_project_WGCNA not found")
if (!file.exists(dir_project_data)) stop("dir_project_data not found")
if (!file.exists(path_data)) stop("path_data not found")
if (!is.null(path_metadata)) if (!file.exists(path_metadata)) stop("path_data not found")

# set up and verify WGCNA directory paths

dir_tables_WGCNA = paste0(dir_project_WGCNA,"tables/")
if (!file.exists(dir_tables_WGCNA)) stop("dir_project WGCNA must contain a /tables subdir")

dir_RObjects_WGCNA = paste0(dir_project_WGCNA,"RObjects/")
if (!file.exists(dir_RObjects_WGCNA)) stop("dir_project_WGCNA must contain a /RObjects subdir")

# set up expression data / output directories 

dir_plots_data = paste0(dir_project_data,"plots/")
if (!file.exists(dir_plots_data)) dir.create(dir_plots_data) 

dir_tables_data = paste0(dir_project_data,"tables/")
if (!file.exists(dir_tables_data)) dir.create(dir_tables_data)

dir_RObjects_data = paste0(dir_project_data,"RObjects/")
if (!file.exists(dir_RObjects_data)) dir.create(dir_RObjects_data)

dir_log_data = paste0(dir_project_data,"log/")
if (!file.exists(dir_log_data)) dir.create(dir_log_data)

dir_scratch = "/scratch/tmp-wgcna/"

######################################################################
######################### LOAD GENE LOADINGS VECTORS #################
######################################################################

run_cellType_module_u <- list()
WGCNA_cellTypes = c()
for (prefix_run in prefixes_WGCNA_run) {
  cellType_module_u_path <- dir(path = dir_RObjects_WGCNA, pattern = paste0(prefix_run, "_list_list_module_u.RDS"), full.names = T)
  run_cellType_module_u[[prefix_run]] <- load_obj(f = cellType_module_u_path)
  WGCNA_cellTypes <- c(WGCNA_cellTypes, names(run_cellType_module_u[[prefix_run]]))
}

names(run_cellType_module_u) <- prefixes_WGCNA_run

######################################################################
######################### LOAD METADATA ############################
######################################################################

if (!is.null(path_metadata)){ 
  message("Loading metadata")
  if (grepl(pattern = "\\.gz", x = path_metadata)) {
    if (grepl(pattern = "\\.csv", x = path_metadata)) {
      metadata <- read.csv(file=gzfile(path_metadata), stringsAsFactors = F, quote="", header=T) 
    } else if (grepl(pattern = "\\.tab", x = path_metadata)) {
      metadata <- read.table(file=gzfile(path_metadata), sep="\t", stringsAsFactors = F, header = T, quote="")
    } else if (grepl(pattern = "\\.txt", x = path_metadata)) {
      metadata <- read.delim(file=gzfile(path_metadata), stringsAsFactors = F, quote = "", header=T)
    }
  } else if (grepl(pattern = "\\.csv", x = path_metadata)) {
    metadata <- read.csv(file=path_metadata, stringsAsFactors = F, quote="", header=T)
  } else if (grepl(pattern = "\\.tab", x = path_metadata)) {
    metadata <- read.table(file=path_metadata,sep="\t", stringsAsFactors = F, quote="", header=T)
  } else if (grepl(pattern = "\\.txt", x = path_metadata)) {
    metadata <- read.delim(file=path_metadata, stringsAsFactors = F, quote="", header=T)
  }
}

######################################################################
######################### LOAD GENE MAPPING ##########################
######################################################################

message("Loading gene mapping files")

run_cellType_mapping <- list()

for (run in prefixes_WGCNA_run) {
  path_mapping <- dir(pattern = paste0(run, ".*_hgnc_to_ensembl_mapping_df.csv$"), path = dir_tables_WGCNA, full.names = T)
  if (length(path_mapping)>1) warning(paste0(run, " WGCNA run: more than one hgcn to ensembl mapping file found, selecting the first one"))
  path_mapping <- path_mapping[1]
  run_cellType_mapping[[run]] <- read.csv(file = path_mapping, header = T, quote = "")
}

names(run_cellType_mapping) <- prefixes_WGCNA_run


######################################################################
######################### LOAD EXPRESSION MATRIX #####################
######################################################################

message("Loading expression matrix")

# load/open connection to expression matrix
# NB: may be huge

if (grepl(pattern = "\\.loom", path_data)) {
  data_obj <-  connect(filename = path_data, mode = "r+")
} else if (grepl(pattern = "\\.RData|\\.RDS|\\.csv|\\.tab|\\.txt|\\.tsv", path_data, ignore.case = T)) {
  data_obj <- load_obj(path_data)
} 

if (!"loom" %in% class(data_obj) & class(data_obj)!="seurat") { # not a loom object, nor a seurat object
  data_obj <- CreateSeuratObject(raw.data=data_obj, project= prefix_out, min.cells = -Inf, min.genes = -Inf)
} else if ("loom" %in% class(data_obj)) {
  data_obj <- Seurat::Convert(from=data_obj, to="seurat")
} # if already a seurat object, do nothing

# Add any further metadata coming from another file
if (!is.null(path_metadata)) {
  if (all_equal(rownames(metadata_tmp), colnames(data_obj@raw.data))) {
    data_obj <- AddMetaData(object = data_obj, metadata = data.frame(metadata))
  } 
}

######################################################################
#################### MAP ENSEMBL TO SYMBOL ###########################
######################################################################

if (!any(grepl("ENS|ENSMUSG", rownames(data_obj@raw.data)))) {
  
  message("Mapping u vector gene names to hgnc symbol")

  cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "map_u_from_ensembl_to_hgnc.txt"))
  
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
}

######################################################################
################ UNLIST MODULES TO ONE LAYER #########################
######################################################################

cellType_module_u <- unlist(x = run_cellType_module_u, recursive = F, use.names = T)
module_u <- unlist(x = cellType_module_u, recursive = F, use.names = T)
names(module_u) <- paste0(prefix_out, ".", names(module_u))

######################################################################
####################### COUNT PERCENT RIBO, MITO #####################
######################################################################

message("Counting percent ribo and percent mito")

idx_mito.genes <- grepl(pattern = "^mt-", x = rownames(data_obj@raw.data), ignore.case=T)
idx_ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(data_obj@raw.data), ignore.case=T)

percent_mito <- Matrix::colSums(data_obj@raw.data[idx_mito.genes, ])/Matrix::colSums(data_obj@raw.data)
percent_ribo <- Matrix::colSums(data_obj@raw.data[idx_ribo.genes, ])/Matrix::colSums(data_obj@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
data_obj <- AddMetaData(object = data_obj, metadata = percent_mito, col.name = "percent_mito")
data_obj <- AddMetaData(object = data_obj, metadata = percent_ribo, col.name = "percent_ribo")

# mito.genes <- grepl(pattern = "^mt-", x = loom_obj[["row_attrs/gene_names"]][], ignore.case=T)
# ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = loom_obj[["row_attrs/gene_names"]][], ignore.case=T)
# 
# loom_obj$apply(name = "col_attrs/percent_mito", FUN = function(mat) {
#   return(rowSums(x = mat[, mito.genes])/rowSums(x = mat))
# }, MARGIN = 2, dataset.use = "matrix")
# 
# loom_obj$apply(name = "col_attrs/percent_ribo", FUN = function(mat) {
#   return(rowSums(x = mat[, ribo.genes])/rowSums(x = mat))
# }, MARGIN = 2, dataset.use = "matrix")

######################################################################
######################### PREPARE DATA ###############################
######################################################################

message("Normalising the data")

data_obj <- NormalizeData(data_obj)
# loom_obj$add.layer(layers=list(data=Seurat::LogNormalize(t(loom_obj[["matrix"]][,]))))
# dimnames(loom_obj$add.layers)
# if (class(x = data.use) == "dgCMatrix" 
#                       "dgTMatrix") { Seurat::FastSparseRowScale(mat = loom_obj[["layers/data"]][,], drop = F], scale = do.scale, center = do.center, 
#                                    scale_max = scale.max, display_progress = FALSE)
# }
# else {
#   data.scale <- FastRowScale(mat = as.matrix(x = data.use[genes.use[my.inds], 
#                                                           , drop = F]), scale = do.scale, center = do.center, 
#                              scale_max = scale.max, display_progress = FALSE)
# }
message("Scaling the data and regressing out confounders percent.mito, percent.ribo and nUMI")

data_obj<-ScaleData(data_obj, 
          vars.to.regress = c("percent_mito", "nUMI", "percent_ribo")[c("percent_mito", "nUMI", "percent_ribo") %in% colnames(data_obj@meta.data)], 
          display.progress= T, 
          do.par=T, 
          num.cores=n_cores)

ident <- if (!is.null(datExpr_ident_cols)) data_obj@meta.data[[datExpr_ident_cols[1]]] else data_obj@ident

datExpr <- data_obj@scale.data
metadata <- data_obj@meta.data
rm(data_obj)

invisible(gc()); invisible(R.utils::gcDLLs())

######################################################################
########################## SET UP LOOM OBJECT ########################
######################################################################

# # Convert (back) to loom
# if (file.exists(paste0(dir_scratch, prefix_out, "_heatmap_datExpr_obj.loom"))) file.remove(paste0(dir_scratch, prefix_out, "_heatmap_datExpr_obj.loom"))
# loom_obj <- Seurat::Convert(from = data_obj,
#                             to = "loom",
#                             filename=paste0(dir_scratch, prefix_out, "_heatmap_datExpr_obj.loom"),
#                            # overwrite=T,
#                             #chunk.size=min(ncol(data_obj@raw.data), 10000),
#                             display.progress=T)
# 
# loom_obj <- loomR::create(filename=paste0(dir_scratch, prefix_out, "_heatmap_datExpr_obj.loom"),
#                           data= data_obj@scale.data,
#                           gene.attrs = NULL,
#                           cell.attrs = data_obj@meta.data,
#                           layers = NULL,
#                           chunk.dims = "auto",
#                           chunk.size = min(ncol(data_obj@scale.data), 10000),
#                           do.transpose = TRUE,
#                           calc.numi = T,
#                           overwrite = T,
#                           display.progress = TRUE)
# 
# rm(data_obj)
# # NormalizeData
# if (is.null(data_obj@data)) {
#   message("Normalising the data")
#   NormalizeData(object=data_obj, display.progress = T)
# } else {
#   message("Normalised data detected")
# }

# ScaleData
# if (scale_data) {
#   if (is.null(data_obj@scale.data)) {
#     message("Scaling the data and regressing out confounders percent.mito, percent.ribo and nUMI")
#     data_obj <- ScaleData(object = data_obj, 
#                           vars.to.regress = c("percent.mito", "percent.ribo", "nUMI"), 
#                           do.scale = T, do.center=T, 
#                           display.progress = T, 
#                           do.par = T, 
#                           num.cores = n_cores)
#     
#     invisible(gc())
#     invisible(R.utils::gcDLLs())
#     
#     datExpr <- data_obj@scale.data 
#   } else {
#     message("scale data detected")
#     datExpr <- data_obj@scale.data 
#   }
# } else {
#   datExpr <- data_obj@data
# }

# Save as loomR
#if (file.exists(paste0(dir_scratch, prefix_out, "_data_obj.loom"))) file.remove(paste0(dir_scratch, prefix_out, "_data_obj.loom"))
# data_obj <- create(filename = paste0(dir_scratch, prefix_out, "_data_obj.loom"), 
#                    data = datExpr, 
#                    #gene.attrs = list(gene_names = rownames(datExpr)),
#                    #cell.attrs = list(cell_names = colnames(datExpr)),
#                    display.progress = T, 
#                    calc.numi = if (scale_data) T else F,
#                    overwrite=T)
#data_obj$add.col.attribute(attribute = list(cell_names = colnames(datExpr)))
#data_obj$add.col.attribute(attribute = list(ident = ident), overwrite = TRUE)
#data_obj$add.row.attribute(attribute = list(gene_names = rownames(datExpr)))
#rm(datExpr)

######################################################################
################ COMPUTE CELL MODULE EMBEDDINGS ######################
######################################################################
message("Computing cell module embeddings")

if (!use_loom) {

  ncell <- ncol(datExpr) 
  
  cell_chunk_idx <- c(seq(from=0, to = ncell, by=10000), ncell)
  
  cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "cell_x_module_matrix.txt"))
  
  for (i in 1:(length(cell_chunk_idx)-1)) { # loop over chunks of cell
    
    submat <- datExpr[,(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]] 
    
    rownames(submat) <- rownames(datExpr)
    
    subEmbed <- sapply(X=module_u, FUN=function(u) {
      t(submat[match(names(u), rownames(submat)),]) %*% matrix(u, nrow = length(u))
    }, 
    simplify=T)
    
    rownames(subEmbed) <- rownames(datExpr)[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]
    colnames(subEmbed) <- names(module_u)
    
    if (i==1) {
      
      if (file.exists(paste0(dir_scratch, 
                             prefix_out, "_cellModEmbed.loom"))) invisible(file.remove(paste0(dir_scratch, prefix_out, "_cellModEmbed.loom")))
      
      cellModEmbed_loom <- loomR::create(filename = paste0(dir_scratch, 
                                                           prefix_out, "_cellModEmbed.loom"), 
                                         data = t(subEmbed), 
                                         #gene.attrs = list(gene_names = rownames(subEmbed), ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]),
                                         #cell.attrs = list(cell_names = colnames(subEmbed)),
                                         display.progress = T, 
                                         calc.numi = F,
                                         overwrite=T)
      
      # NB: We flit dimensions again since the embeddings are cell * module
      #cellModEmbed_loom_chunk$add.row.attribute(list(cell_names = rownames(subEmbed)), overwrite = TRUE)
      
      cellModEmbed_loom$add.col.attribute(list(ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]), overwrite = TRUE)
      #cellModEmbed_loom_chunk$add.row.attribute(list(gene_names = colnames(subEmbed)), overwrite = TRUE)
      cellModEmbed_loom$add.row.attribute(list(module = colnames(subEmbed)), overwrite = TRUE)
      
    } else {
      
      cellModEmbed_loom_chunk <- loomR::create(filename = paste0(dir_scratch, 
                                                                 prefix_out, "_cellModEmbed_chunk.loom"), 
                                               data = t(subEmbed), 
                                               #gene.attrs = list(gene_names = rownames(subEmbed), ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]),
                                               #cell.attrs = list(cell_names = colnames(subEmbed)),
                                               display.progress = T, 
                                               calc.numi = F,
                                               overwrite=T)
      
      cellModEmbed_loom_chunk$add.col.attribute(list(ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]), overwrite = TRUE)
      cellModEmbed_loom_chunk$add.row.attribute(list(module = colnames(subEmbed)), overwrite = TRUE)
      
      try({
        cellModEmbed_loom$close_all()
        cellModEmbed_loom_chunk$close_all()
      })
      
      cellModEmbed_comb <- loomR::combine(looms=list(paste0(dir_scratch, 
                                                            prefix_out, "_cellModEmbed.loom"), 
                                                     paste0(dir_scratch, 
                                                            prefix_out, "_cellModEmbed_chunk.loom")), 
                                          filename=paste0(dir_scratch, prefix_out, "_cellModEmbed_comb.loom"), 
                                          chunk.size=1000, 
                                          overwrite=T, 
                                          display.progress=T)
      
      
      file.remove(paste0(dir_scratch, 
                         prefix_out, "_cellModEmbed.loom"), 
                  paste0(dir_scratch, 
                         prefix_out, "_cellModEmbed_chunk.loom"))
      
      cellModEmbed_loom <- loomR::create(filename=paste0(dir_scratch, 
                                                         prefix_out, "_cellModEmbed.loom"), 
                                         data = t(cellModEmbed_comb[["matrix"]][,])
                                         #gene.attrs = list(modules= cellModEmbed_comb[["row_attrs/gene_names"]]),
                                         #cell.attrs = list(ident= cellModEmbed_comb[["col_attrs/ident"]])
      )
      cellModEmbed_loom$add.row.attribute(list(gene_names=cellModEmbed_comb[["row_attrs/gene_names"]][]), overwrite = TRUE)
      cellModEmbed_loom$add.row.attribute(list(module=cellModEmbed_comb[["row_attrs/module"]][]), overwrite = TRUE)
      cellModEmbed_loom$add.col.attribute(list(cell_names=cellModEmbed_comb[["col_attrs/cell_names"]][]), overwrite = TRUE)
      cellModEmbed_loom$add.col.attribute(list(ident=cellModEmbed_comb[["col_attrs/ident"]][]), overwrite = TRUE)
      
      try({
        cellModEmbed_comb$close_all()
        invisible(file.remove(paste0(dir_scratch, prefix_out, "_cellModEmbed_comb.loom")))
      })
      
    }
    
    rm(submat, subEmbed)
  }
  
  stopCluster(cl)

} else if (use_loom) {

  ncell <- dim(data_obj[["matrix"]][,])[1] # NB: the matrix is stored transposed
  
  cell_chunk_idx <- c(seq(from=0, to = ncell, by=3000), ncell)
  
  cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "cell_x_module_matrix.txt"))
  
  for (i in 1:(length(cell_chunk_idx)-1)) { # loop over chunks of cell
  
    submat <- t(data_obj[["matrix"]][(cell_chunk_idx[i]+1):cell_chunk_idx[i+1],]) # transpose to have cells in rows
  
    rownames(submat) <- data_obj[["row_attrs/gene_names"]][]
  
    subEmbed <- sapply(X=module_u, FUN=function(u) {
      t(submat[match(names(u), rownames(submat)),]) %*% matrix(u, nrow = length(u))
    },
    simplify=T)
  
    rownames(subEmbed) <- data_obj[["col_attrs/cell_names"]][(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]
    colnames(subEmbed) <- names(module_u)
  
  
    if (i==1) {
  
      if (file.exists(paste0(dir_scratch,
                             prefix_out, "_cellModEmbed.loom"))) invisible(file.remove(paste0(dir_scratch, prefix_out, "_cellModEmbed.loom")))
  
      cellModEmbed_loom <- loomR::create(filename = paste0(dir_scratch,
                                                                 prefix_out, "_cellModEmbed.loom"),
                                               data = t(subEmbed),
                                               #gene.attrs = list(gene_names = rownames(subEmbed), ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]),
                                               #cell.attrs = list(cell_names = colnames(subEmbed)),
                                               display.progress = T,
                                               calc.numi = F,
                                               overwrite=T)
  
      # NB: We flit dimensions again since the embeddings are cell * module
      #cellModEmbed_loom_chunk$add.row.attribute(list(cell_names = rownames(subEmbed)), overwrite = TRUE)
  
      cellModEmbed_loom$add.col.attribute(list(ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]), overwrite = TRUE)
      #cellModEmbed_loom_chunk$add.row.attribute(list(gene_names = colnames(subEmbed)), overwrite = TRUE)
      cellModEmbed_loom$add.row.attribute(list(module = colnames(subEmbed)), overwrite = TRUE)
  
    } else {
  
      cellModEmbed_loom_chunk <- loomR::create(filename = paste0(dir_scratch,
                                                                 prefix_out, "_cellModEmbed_chunk.loom"),
                                               data = t(subEmbed),
                                               #gene.attrs = list(gene_names = rownames(subEmbed), ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]),
                                               #cell.attrs = list(cell_names = colnames(subEmbed)),
                                               display.progress = T,
                                               calc.numi = F,
                                               overwrite=T)
  
      cellModEmbed_loom_chunk$add.col.attribute(list(ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]), overwrite = TRUE)
      cellModEmbed_loom_chunk$add.row.attribute(list(module = colnames(subEmbed)), overwrite = TRUE)
  
      try({
        cellModEmbed_loom$close_all()
        cellModEmbed_loom_chunk$close_all()
      })
  
      cellModEmbed_comb <- loomR::combine(looms=list(paste0(dir_scratch,
                                                            prefix_out, "_cellModEmbed.loom"),
                                                     paste0(dir_scratch,
                                                            prefix_out, "_cellModEmbed_chunk.loom")),
                                          filename=paste0(dir_scratch, prefix_out, "_cellModEmbed_comb.loom"),
                                          chunk.size=1000,
                                          overwrite=T,
                                          display.progress=T)
  
  
      file.remove(paste0(dir_scratch,
                         prefix_out, "_cellModEmbed.loom"),
                  paste0(dir_scratch,
                         prefix_out, "_cellModEmbed_chunk.loom"))
  
        cellModEmbed_loom <- loomR::create(filename=paste0(dir_scratch,
                                                           prefix_out, "_cellModEmbed.loom"),
                                           data = t(cellModEmbed_comb[["matrix"]][,])
                                           #gene.attrs = list(modules= cellModEmbed_comb[["row_attrs/gene_names"]]),
                                           #cell.attrs = list(ident= cellModEmbed_comb[["col_attrs/ident"]])
        )
        cellModEmbed_loom$add.row.attribute(list(gene_names=cellModEmbed_comb[["row_attrs/gene_names"]][]), overwrite = TRUE)
        cellModEmbed_loom$add.row.attribute(list(module=cellModEmbed_comb[["row_attrs/module"]][]), overwrite = TRUE)
        cellModEmbed_loom$add.col.attribute(list(cell_names=cellModEmbed_comb[["col_attrs/cell_names"]][]), overwrite = TRUE)
        cellModEmbed_loom$add.col.attribute(list(ident=cellModEmbed_comb[["col_attrs/ident"]][]), overwrite = TRUE)
  
        try({
          cellModEmbed_comb$close_all()
          invisible(file.remove(paste0(dir_scratch, prefix_out, "_cellModEmbed_comb.loom")))
        })
  
    }
  
  
  
  
    rm(submat, subEmbed)
  }
  
  stopCluster(cl)

}

# AVERAGE EMBEDDINGS FOR datExpr_ident_cols

message(paste0("Computing expression averaged over ",paste0(datExpr_ident_cols,collapse=", ")))

# remove NAs from ident columns

for (id in datExpr_ident_cols) {
  metadata[[ident_col]][is.na(metadata[[ident_col]])] <- "Not_attributed"
}

list_subsetModEmbed <- lapply(datExpr_ident_cols, function(ident_col) {

  cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", ident_col, "_mean_mod_expr.txt"))
  
  for (i in 1:length(sort(unique(ident)))) {
    # cellModEmbed_loom$apply(name = paste0("row_attrs/",unique(ident)[i]), FUN = colMeans, MARGIN = 1, index.use = which(ident==unique(ident)[i]),
    #             dataset.use = "matrix", display.progress = FALSE, overwrite = TRUE)
    ident <- metadata[[ident_col]]
    subEmbed <- cellModEmbed_loom[["matrix"]][ident==sort(unique(ident))[i],]
    meanExpr <- matrix(parApply(cl, X=subEmbed, FUN=mean, MARGIN=2), nrow=1)
    rownames(meanExpr) <- sort(unique(ident))[i]
    colnames(meanExpr) <- cellModEmbed_loom[["row_attrs/module"]][]
    # it might be ok to do this w/o parallelising
    if (i==1) subsetModEmbed <- meanExpr else subsetModEmbed <- rbind(subsetModEmbed, meanExpr)
  }
  stopCluster(cl)
  
  subsetModEmbed

})

names(list_subsetModEmbed) <- datExpr_ident_cols

# cellTypeModEmbed <- t(parSapply(cl, X=paste0("row_attrs/",unique(ident)), FUN = function(index) cellModEmbed_loom[[index]][], simplify = T))
# rownames(cellTypeModEmbed) <- gsub("col_attrs/", "", rownames(cellTypeModEmbed))
# colnames(cellTypeModEmbed) <- 

#rownames(cellTypeModEmbed) = sort(unique(ident)); colnames(cellTypeModEmbed) <- names(module_u)

######################################################################
######################## CLUSTER CELLS, MODULES  #####################
######################################################################

modDist <- dist(x=t(cellModEmbed_loom[["matrix"]][,]), method="euclidean", diag=F, upper=F) 
modDendro <- hclust(d=modDist, method="average")
modClust_order <- modDendro$order

# cellDist <- dist(x=cellModEmbed_loom[["matrix"]][,], method="euclidean", diag=F, upper=F) 
# cellDendro <- hclust(d=cellDist, method="average")
# cellClust_order <- cellDendro$order

cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", ident_col, "_mean_mod_expr_clust.txt"))

list_subsetModEmbed_modClustOrder <- parLapplyLB(cl, list_subsetModEmbed, function(subsetModEmbed) {
  cellType_modDist <- dist(x=t(subsetModEmbed), method="euclidean", diag=F, upper=F) 
  cellType_modDendro <- hclust(d=cellType_modDist, method="average")
  cellType_modDendro$order
})

list_subsetModEmbed_cellClustOrder <- parLapplyLB(cl, list_subsetModEmbed, function(subsetModEmbed) {
  cellType_cellDist <- dist(x=subsetModEmbed, method="euclidean", diag=F, upper=F) 
  cellType_cellDendro <- hclust(d=cellType_modDist, method="average")
  cellType_cellDendro$order
})

stopCluster(cl)

######################################################################
#################### PLOT CELL*MODULE MATRIX #########################
######################################################################

if (do_plot) {
  # generate colors for annotation
  colors_uniq <- gsub("\\d", "", colors()) %>% unique
  
  ######################## ROW ANNOTATION ##############################
  
  list_htra <- lapply(datExpr_ident_cols, function(ident_col) {
    ident <- metadata[[ident_col]]
    htra_colors <-   sample(if(length(unique(ident))<length(colors_uniq)) colors_uniq else colors(), size=length(unique(ident)), replace=F) 
    names(htra_colors) <- sort(unique(ident), decreasing=F)  
    
    htra <- rowAnnotation(ident_col = sort(ident, decreasing=F),
                          #annotation_legend_param = list(cellCluster = list(nrow = length(unique(ident)), title = "Cell cluster", title_position = "topcenter")),
                          col=list(ident_col = htra_colors),
                          width = unit(5, "mm"),
                          show_legend=F,
                          #annotation_legend_param = list(show_legend = F), 
                          show_annotation_name = F)#,
    htra
  })
    
  # make row annotation labels as a separate row annotation
    list_htra_link <- lapply(datExpr_ident_cols, function(ident_col)  {
      
      ident <- metadata[[ident_col]]
      htra_labels <- sort(unique(ident), decreasing=F)
      
      fontsize <- min(20, max(8, as.integer(500/length(unique(ident)))))
  
      ident <- metadata[[ident_col]]
      htra_label_pos <- cumsum(table(ident)) - as.integer(round(table(ident)/2,0))
      
      rowAnnotation(link = anno_link(which= "row",
                    at = htra_label_pos, 
                    labels = htra_labels,
                    link_width = unit(3, "mm"),
                    labels_gp=gpar(fontsize=fontsize), padding=0.5), 
                    width = unit(5, "mm") + max_text_width(labels, gp=gpar(fontsize=fontsize)))
      })
  
  ######################## COLUMN ANNOTATION ###########################
  
  #  WGCNA run column annotation
  WGCNA_run <- unlist(mapply(function(names_run, cellType_module_u) rep(names_run, times= length(unlist(cellType_module_u, recursive = F))), names_run=prefixes_WGCNA_run, cellType_module_u=run_cellType_module_u, SIMPLIFY=T))
  htca_WGCNA_run_colors <- sample(if(length(prefixes_WGCNA_run)<length(colors_uniq)) colors_uniq else colors(), size=length(prefixes_WGCNA_run), replace=F) 
  names(htca_WGCNA_run_colors) <- sort(unique(prefixes_WGCNA_run, decreasing=F))
  
  # Prepare WGCNA celltype column annotation
  #WGCNA_cellCluster <- rep(gsub(pattern = paste0(paste0(names(run_cellType_module_u), collapse = "|"),"\\."), "", names(cellType_module_u)), sapply(cellType_module_u, length))
  WGCNA_cellCluster <- rep(WGCNA_cellTypes,times= sapply(cellType_module_u, length))
  
 # WGCNA_cellCluster <- gsub("\\.","_",WGCNA_cellCluster)
  htca_WGCNA_cellCluster_colors <- sample(if (length(WGCNA_cellTypes)<length(colors_uniq)) colors_uniq else colors(), size=length(WGCNA_cellTypes), replace=F) 
  names(htca_WGCNA_cellCluster_colors) <- sort(WGCNA_cellTypes, decreasing=F) # no need to sort?
  
  # Prepare WGCNA celltype column annotation labels as separate column annotation

  #htca_labels <- sort(unique(WGCNA_cellCluster), decreasing=F)
  htca_labels <- sort(WGCNA_cellTypes)
  htca_label_pos <- cumsum(table(WGCNA_cellCluster)) - as.integer(round(table(WGCNA_cellCluster)/2,0))

  # Add PPI enrichment column annotation
  if (plot_PPI_enrichment) {
    run_cellType_STRINGdb_output_all <- list()
    for (run in prefixes_WGCNA_run) {
      cellType_STRINGdb_output_all_paths <- dir(path = dir_tables_WGCNA, pattern = paste0(run, ".*STRINGsb_output_all\\.csv$"), full.names = T)
      run_cellType_STRINGdb_output_all[[run]] <- lapply(cellType_STRINGdb_output_all_paths, load_obj) # list of list of data.frame
    }
    names(run_cellType_STRINGdb_output_all) <- prefixes_WGCNA_run
  }
  
  cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "make_PPI_column_annotation.txt"))
  
  run_cellType_modPPIenrichment <- clusterMap(cl, function(cellType_mod_u, cellType_STRINGdb_output_all) {
    mapply(function(modname, STRINGdb_output_all) {
        STRINGdb_output_all[["q.value"]][grep(modname, STRINGdb_output_all[["colors"]])] %>% -log10
      },
      modname = names(cellType_mod_u), STRINGdb_output_all=cellType_STRINGdb_output_all, SIMPLIFY=F)
  }, cellType_mod_u = run_cellType_module_u, cellType_STRINGdb_output_all= run_cellType_STRINGdb_output_all, SIMPLIFY=F, .scheduling=c("dynamic"))
  
  stopCluster(cl)
  
  cellType_modPPIenrichment <- unlist(run_cellType_modPPIenrichment, recursive = F, use.names = T)
  modPPIenrichment <- unlist(cellType_modPPIenrichment, recursive = F, use.names = T)
  htca_modPPIenrichment_colors <- colorRamp2(c(0, 20), c("white", "red"))
  
  # Add MAGMA enrichment column annotation
  if (!is.null(plot_magma_enrichment)) {
    run_magma <- list()
    for (run in prefixes_WGCNA_run) {
      path_magma <- dir(pattern = paste0(run, "_magma\\.fdr\\.log\\.csv$"), path = dir_tables_WGCNA, full.names = T)
      run_magma[[run]] <- read.csv(file = path_magma, header = T, quote = "")
    }
    names(run_magma) <- prefixes_WGCNA_run
    
    cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "make_magma_column_annotation.txt"))
    
    run_cellType_modMagmaEnrichment <- clusterMap(cl, function(cellType_mod_u, magma) {
      lapply(names(cellType_mod_u), function(modname){
        gwas_columns <- unique(sapply(plot_magma_enrichment, function(gwas_pattern) {
          grep(pattern=gwas_pattern, x = colnames(magma))
          }
          , SIMPLIFY=T))
        out <- magma[[grep(modname, magma[["module"]]), c(grep("module", colnames(magma)),gwas_columns), drop=T]] 
        out
      })
    }, cellType_mod_u = run_cellType_module_u, magma = run_magma, SIMPLIFY=F, .scheduling=c("dynamic"))
    
    stopCluster(cl)
  }
  
  cellType_modMagmaEnrichment <- unlist(run_cellType_modMagmaEnrichment, recursive = F, use.names = T)
  
  # deal with multiple GWAS magma columns 
  # also check module order 
  list_modMagmaEnrichment <- list()
  for (cname in colnames(cellType_modMagmaEnrichment[[1]])[-grep("module", colnames(cellType_modMagmaEnrichment[[1]]))]) { #assumes the GWAS are always the same
    list_modMagmaEnrichment[[cname]] <- sapply(cellType_modMagmaEnrichment, function(modMagmaEnrichment) modMagmaEnrichment[[cname]])
    names(list_modMagmaEnrichment[[cname]])  <- modMagmaEnrichment[["module"]] # list of named annotation vectors
  }
  
  
  # TODO: check that the magma vectors (of length module) are ordered the same way as the rest of the annotaiton
  # TODO: because we omit the neext unlist step, does the name scheme work out?
  #modMagmaEnrichment <- unlist(list_modMagmaEnrichment, recursive = F, use.names = T)
  # TODO: generate a vector of 'up' colors, one per gwas phenotype. Then generate a list of colours scheme vectors
  #htca_modMagmaEnrichment_colors <- colorRamp2(c(0, 20), c("white", "blue"))
  
  # Add metadata correlation column annotation
  #TODO 
  
  # Add WGCNA ident group column annotation
  htca_WGCNA_ident_group_df  <- NULL 
  if (!is.null(WGCNA_ident_cats)) {
    
    sapply(X=WGCNA_ident_cats, FUN = function(ident_cat) {
      level_vec <- rep(NA_character_, length=length(module_u))
      for (lvl in unique(ident_cat)){
        idx <- grepl(lvl, names(module_u))
        level_vec[idx] <- lvl
      }
    }, simplify = T) %>% data.frame -> htca_WGCNA_ident_df_tmp
    
    # Prepare WGCNA ident group colors
    list_htcs_WGCNA_ident_group_colors <- apply(X = htca_WGCNA_ident_df_tmp, MARGIN = 2, FUN = function(level_vec){
      out <- sample(if (length(unique(level_vec))<length(colors_uniq)) colors_uniq else colors(), size=length(unique(level_vec)), replace=F) 
      names(out) <- sort(unique(level_vec), decreasing=F)
      out
    })
    
  }

  # Finally gathercolumn annotations in a dataframe
  # TODO htca_df <- data.frame()
  
  # If there are common celltypes between the data and the WGCNA module celltypes, they get the same heatmap annotation colors
  if (!TODO) { 
    match_cellType_idx <- match(names(htca_WGCNA_cellCluster_colors), names(htra_cellCluster_colors), nomatch = NA)
    new_colors_tmp <- htra_cellCluster_colors[na.omit(match_cellType_idx)]
    if (length(new_colors_tmp) > 0) htca_WGCNA_cellCluster_colors[names(htca_WGCNA_cellCluster_colors) %in% names(htra_cellCluster_colors)] <- new_colors_tmp
  }
  
  # Make column annotation objects
  htca <- columnAnnotation(#WGCNA_run=WGCNA_run,
                           WGCNA_cellCluster=WGCNA_cellCluster,
 # htca <- columnAnnotation(WGCNA_run=sort(WGCNA_run,  1:length(WGCNA_run), decreasing=F),
 #                          WGCNA_cellCluster=sort(WGCNA_cellCluster, decreasing=F),
                           show_legend=F,
                           link = anno_link(which = c("column"),
                                            at = htca_label_pos, 
                                            labels = htca_labels,
                                            link_width = unit(0, "mm"),
                                            labels_gp=gpar(fontsize=14), 
                                            padding=0.3),  
                           #col = list(WGCNA_run = htca_WGCNA_run_colors, WGCNA_cellCluster = htca_WGCNA_cellCluster_colors),
                           col = list(WGCNA_cellCluster = htca_WGCNA_cellCluster_colors),
                           annotation_height = unit.c(unit(1, "cm"), unit(1,"cm"), unit(1, "cm")))
  
  # When clustering module columns omit the rowannotation labels
  
  htca_noLabel <- columnAnnotation(WGCNA_run=sort(WGCNA_run,  1:length(WGCNA_run), decreasing=F),
                           WGCNA_cellCluster=sort(WGCNA_cellCluster, decreasing=F),
                           show_legend=F,
                           col = list(WGCNA_run = htca_WGCNA_run_colors, WGCNA_cellCluster = htca_WGCNA_cellCluster_colors),
                           annotation_height = unit.c(unit(1, "cm"), unit(1,"cm")))
  
  # Make plots
  ht1 <- Heatmap(cellModEmbed_loom[["matrix"]][order(ident, 1:nrow(cellModEmbed_loom[["matrix"]][,]), decreasing=F),order(WGCNA_run,  1:length(WGCNA_run), decreasing=F)], 
                 cluster_rows = F,
                 cluster_columns = F, 
                 show_heatmap_legend = T,
                 show_row_names = F, 
                 show_column_names = F,
                 bottom_annotation = htca,
                 use_raster=T,
                 raster_device = c("png"),
                 raster_quality = 1,
                 heatmap_legend_param = list(title = "Expression"))
  pdf(sprintf("%s%s_cellModEmbed.pdf", dir_plots_data, prefix_out), h=max(25, min(40, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2000)), w=max(30, min(40, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 50)))
  draw(ht1+list_htra[[1]]+list_htra_link[[1]])
  dev.off()
  
  ht1_cellClust <- Heatmap(cellModEmbed_loom[["matrix"]][order(ident, 1:nrow(cellModEmbed_loom[["matrix"]][,]), decreasing=F),order(WGCNA_run,  1:length(WGCNA_run), decreasing=F)], 
                 cluster_rows = T,
                 cluster_columns = F, 
                 show_row_dend = F, 
                 show_heatmap_legend = T,
                 show_row_names = F, 
                 show_column_names = F,
                 bottom_annotation = htca,
                 use_raster=T,
                 raster_device = c("png"),
                 raster_quality = 1,
                 heatmap_legend_param = list(title = "Expression"))
  pdf(sprintf("%s%s_cellModEmbed_cell_clust.pdf", dir_plots_data, prefix_out), h=max(25, min(40, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2000)), w=max(25, min(40, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 50)))
  draw(ht1_cellClust+htra)#+htra_link)
  dev.off()
  
  ht1_modClust <- Heatmap(cellModEmbed_loom[["matrix"]][order(ident, 1:nrow(cellModEmbed_loom[["matrix"]][,]), decreasing=F),order(WGCNA_run,  1:length(WGCNA_run), decreasing=F)], 
                           cluster_rows = F,
                           cluster_columns = T, 
                           show_column_dend = F, 
                           show_heatmap_legend = T,
                           show_row_names = F, 
                           show_column_names = F,
                           bottom_annotation = htca_noLabel,
                           use_raster=T,
                           raster_device = c("png"),
                           raster_quality = 1,
                           heatmap_legend_param = list(title = "Expression"))#,
  pdf(sprintf("%s%s_cellModEmbed_mod_clust.pdf", dir_plots_data, prefix_out), h=max(25, min(40, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2000)), w=max(25, min(40, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 50)))
  draw(ht1_modClust+htra+htra_link)
  dev.off()
  
  ht1_cellModClust <- Heatmap(cellModEmbed_loom[["matrix"]][order(ident, 1:nrow(cellModEmbed_loom[["matrix"]][,]), decreasing=F),order(WGCNA_run,  1:length(WGCNA_run), decreasing=F)], 
                          cluster_rows = T,
                          cluster_columns = T, 
                          show_row_dend = F, 
                          show_column_dend = F, 
                          show_heatmap_legend = T,
                          show_row_names = F, 
                          show_column_names = F,
                          bottom_annotation = htca_noLabel,
                          use_raster=T,
                          raster_device = c("png"),
                          raster_quality = 1,
                          heatmap_legend_param = list(title = "Expression"))
  pdf(sprintf("%s%s_cellModEmbed_cell_mod_clust.pdf", dir_plots_data, prefix_out), h=max(25, min(40, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2000)), w=max(25, min(40, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 50)))
  draw(ht1_cellModClust+htra)#+htra_link)
  dev.off()
  
  ht2 <- Heatmap(cellTypeModEmbed[,], 
                 cluster_rows = F,
                 cluster_columns = F, 
                 show_heatmap_legend = T, 
                 show_row_names = T, 
                 show_column_names = F,
                 bottom_annotation = htca,
                 use_raster=T,
                 raster_device = c("png"),
                 raster_quality = 2,
                 heatmap_legend_param = list(title = "Expression"))
  pdf(sprintf("%s%s_cellTypeModEmbed.pdf", dir_plots_data, prefix_out), h=max(10, min(40, nrow(cellTypeModEmbed) %/% 4)), w=max(10, min(40, ncol(cellTypeModEmbed) %/% 50)))
  draw(ht2)
  dev.off()
  
  ht2_cellTypeClust <- Heatmap(cellTypeModEmbed[,], 
                 cluster_rows = T,
                 cluster_columns = F, 
                 show_row_dend = F, 
                 show_heatmap_legend = T, 
                 show_row_names = T, 
                 show_column_names = F,
                 bottom_annotation = htca,
                 use_raster=T,
                 raster_device = c("png"),
                 raster_quality = 2,
                 heatmap_legend_param = list(title = "Expression"))#
  pdf(sprintf("%s%s_cellTypeModEmbed_cellType_clust.pdf", dir_plots_data, prefix_out), h=max(10, min(40, nrow(cellTypeModEmbed) %/% 4)), w=max(10, min(40, ncol(cellTypeModEmbed) %/% 50)))
  draw(ht2_cellTypeClust)
  dev.off()
  
  ht2_modClust <- Heatmap(cellTypeModEmbed[,], 
                           cluster_rows = F,
                           cluster_columns = T, 
                           show_column_dend = F, 
                           show_heatmap_legend = T, 
                           show_row_names = T, 
                           show_column_names = F,
                           bottom_annotation = htca_noLabel,
                           use_raster=T,
                           raster_device = c("png"),
                           raster_quality = 2,
                           heatmap_legend_param = list(title = "Expression"))
  
  pdf(sprintf("%s%s_cellTypeModEmbed_mod_clust.pdf", dir_plots_data, prefix_out), h=max(10, min(40, nrow(cellTypeModEmbed) %/% 4)), w=max(10, min(40, ncol(cellTypeModEmbed) %/% 50)))
  draw(ht2_modClust)
  dev.off()
  
  ht2_cellModClust <- Heatmap(cellTypeModEmbed[,], 
                          cluster_rows = T,
                          cluster_columns = T, 
                          show_row_dend = F, 
                          show_column_dend = F, 
                          show_heatmap_legend = T, 
                          show_row_names = T, 
                          show_column_names = F,
                          bottom_annotation = htca_noLabel,
                          use_raster=T,
                          raster_device = c("png"),
                          raster_quality = 2,
                          heatmap_legend_param = list(title = "Expression"))#
  pdf(sprintf("%s%s_cellTypeModEmbed_cellType_mod_clust.pdf", dir_plots_data, prefix_out), h=max(10, min(40, nrow(cellTypeModEmbed) %/% 4)), w=max(10, min(40, ncol(cellTypeModEmbed) %/% 50)))
  draw(ht2_cellModClust)
  dev.off()
}
######################################################################
############################# WRAP UP ################################
######################################################################

try(data_obj$close_all())
if (file.exists(paste0(dir_scratch, prefix_out, "_data_obj.loom"))) try(invisible(file.remove(paste0(dir_scratch, prefix_out, "_data_obj.loom"))))

try(cellModEmbed_loom_chunk$close_all())
if (file.exists(paste0(dir_scratch, 
                       prefix_out, "_cellModEmbed_chunk.loom"))) try(invisible(file.remove(paste0(dir_scratch, 
                                                                                                  prefix_out, "_cellModEmbed_chunk.loom"))))                  


cellModEmbed <- cellModEmbed_loom[["matrix"]][,]
rownames(cellModEmbed) <- cellModEmbed_loom[["col_attrs/cell_names"]][]
colnames(cellModEmbed) <- cellModEmbed_loom[["row_attrs/module"]][]  

try(cellModEmbed_loom$close_all())
try(invisible(file.remove(paste0(dir_scratch, prefix_out, "_cellModEmbed.loom"))))

saveRDS(cellModEmbed, file = paste0(dir_RObjects_data, prefix_out, "cellModEmbed.RDS"), compress = "gzip")
# Save the embeddings averaged over tissue, celltype etc
#date = substr(gsub("-","",as.character(Sys.Date())),3,1000)
saveRDS(list_subsetModEmbed, file = paste0(dir_RObjects_data, prefix_out, "list_subsetModEmbed.RDS"), compress = "gzip")
#save.image(paste0(dir_scratch, prefix_out, "_", date, "_session_image.RData"))               

message("DONE!")
