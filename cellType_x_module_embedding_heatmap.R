# Make cellType x module embedding matrix
# Use gene loadings from one data (sub)set and possibly expression matrix from another

# loomR:
# Introduction to loomR: https://satijalab.org/loomR/loomR_tutorial.html
# MCA tutorial: https://satijalab.org/seurat/mca_loom.html
# loom file specification: http://linnarssonlab.org/loompy/format/index.html

# usage: 
# export R_MAX_NUM_DLLS=999
# time Rscript /projects/jonatan/wgcna-src/wgcna-toolbox/cellType_x_module_embedding_heatmap.R --dir_project_WGCNA /projects/jonatan/tmp-mousebrain/ --dir_project_data /projects/jonatan/tmp-mousebrain/ --path_data /projects/jonatan/tmp-mousebrain/RObjects/L5.RDS --prefixes_WGCNA_run 'c("Neurons_ClusterName_2")' --prefix_out mousebrain_neuronMods_1 --metadata_subset_col ClusterName --scale_data F --n_cores 20
# time Rscript /projects/jonatan/wgcna-src/wgcna-toolbox/cellType_x_module_embedding_heatmap.R --dir_project_WGCNA /projects/jonatan/tmp-maca/ --dir_project_data /projects/jonatan/tmp-maca/ --path_data /data/pub-others/tabula_muris/figshare/180126-facs/maca.seurat_obj.facs.figshare_180126.RData --prefixes_WGCNA_run 'c("tissue_cell_type_kME")' --prefix_out maca_mods_1 --metadata_subset_col tissue_cell_type --scale_data F --n_cores 30
# time Rscript /projects/jonatan/wgcna-src/wgcna-toolbox/cellType_x_module_embedding_heatmap.R --dir_project_WGCNA /projects/jonatan/tmp-mousebrain/ --dir_project_data /projects/jonatan/tmp-mousebrain/ --path_data /projects/jonatan/tmp-mousebrain/RObjects/L5.RDS --prefixes_WGCNA_run 'c("Neurons_ClusterName_2")' --prefix_out mousebrain_neuronMods_scale_1 --metadata_subset_col ClusterName --scale_data T --n_cores 20
# time Rscript /projects/jonatan/wgcna-src/wgcna-toolbox/cellType_x_module_embedding_heatmap.R --dir_project_WGCNA /projects/jonatan/tmp-maca/ --dir_project_data /projects/jonatan/tmp-maca/ --path_data /data/pub-others/tabula_muris/figshare/180126-facs/maca.seurat_obj.facs.figshare_180126.RData --prefixes_WGCNA_run 'c("tissue_cell_type_kME")' --prefix_out maca_mods_scaled_1 --metadata_subset_col tissue_cell_type --scale_data T --n_cores 20
# time Rscript /projects/jonatan/wgcna-src/wgcna-toolbox/cellType_x_module_embedding_heatmap.R --dir_project_WGCNA /projects/jonatan/tmp-maca/ --dir_project_data /projects/jonatan/tmp-maca/ --path_data /projects/jonatan/tmp-maca/RObjects/maca_seurat_pancreas.Rdata --prefixes_WGCNA_run 'c("tissue_cell_type_kME")' --prefix_out test_heatmap_1 --scale_data T --n_cores 7

# TODO: Pull metadata from loom object?
# TODO: Normalize and scale data in loom object?
# TODO: Use sqroot transform to bring out colors?
# TODO: Add column annotation (get it to work..)

suppressPackageStartupMessages(library(optparse))

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
              help = "WGCNA run prefixes given as vector in double (single) quotes of characters in single (double) quotes."), 
  make_option("--prefix_out", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"), 
  make_option("--metadata_subset_col", type="character", default=NULL,
              help = "Specify a seurat@meta.data$... column to use for subsetting the Seurat object. If NULL uses the @ident slot. [default %default]"),
  make_option("--scale_data", type="logical", default = "TRUE",
              help = "[default %default]. If TRUE, and path_data points to an RObject, the script will try to use the @scale.data slot and, if missing, will normalize and scale the data and regress out nUMI, percent.mito and percent.ribo. If FALSE, for RObjects will use the @data slot, or, if empty, the @raw.data slot"),  
  make_option("--add_PPI_enrichment", type="logical", default = "F",
              help = "Add module PPI enrichment to plots, [default %default]"), 
  make_option("--add_variant_enrichment", type="logical", default = "F",
              help = "Add BMI mendelian and rare variant enrichment to plots, [default %default]"), 
  make_option("--add_magma_enrichment", type="logical", default = "F",
              help = "Add module gwas enrichments to plots, [default %default]"), 
  make_option("--add_metadata_corr", type="logical", default = "F",
              help = "Add any significant metadata correlation coefficients to plots, [default %default]"), 
  make_option("--n_cores", type="integer", default = 15L,
              help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]")
)

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
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

dir_project_WGCNA <- opt$dir_project_WGCNA 
path_data <- opt$path_data
path_metadata <- opt$path_metadata
dir_project_data <- opt$dir_project_data
tmp <- opt$prefixes_WGCNA_run
prefixes_WGCNA_run <- eval(parse(text=tmp))
prefix_out <- opt$prefix_out 
metadata_subset_col <- opt$metadata_subset_col
scale_data <- opt$scale_data
add_PPI_enrichment <- opt$add_PPI_enrichment
add_variant_enrichment <- opt$add_variant_enrichment
add_magma_enrichment <- opt$add_magma_enrichment
add_metadata_corr <- opt$add_metadata_corr
n_cores <- opt$n_cores
 

######################################################################
######################### SET OPTIONS (DEV) ##########################
######################################################################

if (FALSE) { 
  dir_project_WGCNA <- "/projects/jonatan/tmp-mousebrain/" 
  path_data <- "/projects/jonatan/tmp-mousebrain/RObjects/L5.RDS"
  path_metadata <- NULL
  dir_project_data <- "/projects/jonatan/tmp-mousebrain/"
  prefixes_WGCNA_run <- c("Neurons_ClusterName_1b")

  #prefixes_WGCNA_run <- c("Vascular_ClusterName_1", "PeripheralGlia_ClusterName_1", "Oligos_ClusterName_1", "Neurons_ClusterName_1b", "Immune_ClusterName_1", "Ependymal_ClusterName_1", "Astrocytes_ClusterName_1") 
  prefix_out <- "mousebrain_mods"
  metadata_subset_col <- "ClusterName" 
  scale_data <- T
  add_PPI_enrichment = F
  add_variant_enrichment = F
  add_magma_enrichment = F
  add_metadata_corr = F
  n_cores = 30
}

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
for (run in prefixes_WGCNA_run) {
  cellType_module_u_path <- dir(path = dir_RObjects_WGCNA, pattern = paste0(run, "_list_list_module_u.RDS"), full.names = T)
  run_cellType_module_u[[run]] <- load_obj(f = cellType_module_u_path)
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
  path_mapping <- dir(pattern = paste0(run, ".*_hgnc_to_ensembl_mapping"), path = dir_tables_WGCNA, full.names = T)
  run_cellType_mapping[[run]] <- read.csv(file = path_mapping, header = T, quote = "")
}

names(run_cellType_mapping) <- prefixes_WGCNA_run

######################################################################
################### MAP U VECTOR GENE NAMES TO SYMBOL ################
######################################################################

message("Mapping u vector gene names to symbol")

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

cellType_module_u <- unlist(x = run_cellType_module_u, recursive = F, use.names = T)
module_u <- unlist(x = cellType_module_u, recursive = F, use.names = T)
names(module_u) <- paste0(prefix_out, ".", names(module_u))

######################################################################
######################### LOAD EXPRESSION MATRIX #####################
######################################################################

message("Loading expression matrix")

# load/open connection to expression matrix
# NB: may be huge

if (grepl(pattern = "\\.loom", path_data)) {
  data_obj <-  connect(filename = path_data, mode = "r+")
  fileType <- "loom"
} else if (grepl(pattern = "\\.RData|\\.Rdata|\\.RDS|\\.rds", path_data)) {
  data_obj <- load_obj(path_data)
  fileType = "RObject"
}

if (fileType == "loom" & scale_data) warning("Currently unable to scale and regress data from a loom file")
# Check that metadata and data cell names match

metadata <- NULL 

if (!is.null(path_metadata)) {
  if (fileType=="RObject") {
    if (!is.null(data_obj@data)) if (any(rownames(metadata) != rownames(data_obj@data))) stop("path_data and path_metadata object rownames do not match")
  } else if (fileType=="loom") {
    if (any(rownames(metadata) != lfile[["col_attrs/cell_names"]][])) stop("path_data and path_metadata object rownames do not match")
  }
}

######################################################################
####################### COUNT PERCENT RIBO, MITO #####################
######################################################################
# todo: make sure row/column is correct: matrix is transposed

if (scale_data) {
  message("Counting percent ribo and percent mito")
  if (fileType =="loom") {
    if (FALSE) { # TODO
      mito.genes <- grepl(pattern = "^mt-", x = data_obj[["row_attrs/gene_names"]][], ignore.case=T)
      ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = data_obj[["row_attrs/gene_names"]][], ignore.case=T)
      
      
      data_obj$apply(name = "col_attrs/percent_mito", FUN = function(mat) {
        return(rowSums(x = mat[, mito.genes])/rowSums(x = mat))
      }, MARGIN = 2, dataset.use = "matrix")
      
  
      data_obj$apply(name = "col_attrs/percent_ribo", FUN = function(mat) {
        return(rowSums(x = mat[, ribo.genes])/rowSums(x = mat))
      }, MARGIN = 2, dataset.use = "matrix")
    }
    # sum.mito <- data_obj$map(FUN = function(x) rowSums(x[,mito.genes_idx]), MARGIN = 2, chunk.size = 500, dataset.use = "matrix", 
    #                       display.progress = FALSE)
    # sum.ribo <- data_obj$map(FUN = function(x) rowSums(x[,ribo.genes_idx]), MARGIN = 2, chunk.size = 500, dataset.use = "matrix", 
    #                          display.progress = FALSE)
    # sum.all <- data_obj$map(FUN = function(x) rowSums, MARGIN = 2, chunk.size = 500, dataset.use = "matrix", 
    #                         display.progress = FALSE)
    # percent.mito <- sum.mito/sum.all
    # percent.ribo <- sum.ribo/sum.all
    # data_obj <- lfile$add.col.attribute(list(percent.mito=percent.mito), overwrite = TRUE)
    # data_obj <- lfile$add.col.attribute(list(percent.ribo=percent.ribo), overwrite = TRUE)
  } else if (fileType == "RObject") {
    mito.genes <- grep(pattern = "^mt-", x = rownames(x = data_obj@data), value = TRUE, ignore.case=T)
    ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = data_obj@data), value = TRUE)
    percent.mito <- Matrix::colSums(data_obj@raw.data[mito.genes, ])/Matrix::colSums(data_obj@raw.data)
    percent.ribo <- Matrix::colSums(data_obj@raw.data[ribo.genes, ])/Matrix::colSums(data_obj@raw.data)
    data_obj <- AddMetaData(object = data_obj, metadata = percent.mito, col.name = "percent.mito")
    data_obj <- AddMetaData(object = data_obj, metadata = percent.ribo, col.name = "percent.ribo")
  }
}
######################################################################
######################### PREPARE DATA ###############################
######################################################################
# If scale_data == T, use scale.data or run ScaleData

if (file.exists(paste0(dir_scratch, prefix_out, "_datExpr.loom"))) file.remove(paste0(dir_scratch, prefix_out, "_datExpr.loom"))
  # message(paste0("Existing loom file found under ", paste0(dir_scratch, prefix_out, "_datExpr.loom"), ". Using that data"))
  # data_obj <- loomR::connect(filename = paste0(dir_scratch, prefix_out, "_datExpr.loom"), mode = "r+")
  # 
#} else {
  
if (fileType=="RObject") {
  # get ident
  ident <- if (!is.null(metadata_subset_col)) as.character(data_obj@meta.data[[metadata_subset_col]]) else as.character(data_obj@ident)
  
  # get metadata
  if (!is.null(data_obj@meta.data)) {
    metadata <- if (!is.null(path_metadata)) data.frame(metadata, data_obj@meta.data) else data_obj@meta.data
  }
  
  # NormalizeData
  #if (is.null(data_obj@data)) data_obj <- NormalizeData(object=data_obj, display.progress = T)
  # If needed, normalize
  message("Normalising the data")
  
  data_obj <- NormalizeData(object=data_obj, display.progress = T)
  
  # ScaleData
  if (scale_data) {
    message("Scaling the data and regressing out confounders")
    if (is.null(data_obj@scale.data)) {
      
      data_obj <- ScaleData(object = data_obj, 
                            vars.to.regress = c("percent.mito", "percent.ribo", "nUMI"), 
                            do.scale = T, do.center=T, 
                            display.progress = T, 
                            do.par = T, 
                            num.cores = n_cores)
      datExpr <- data_obj@scale.data 
    
    } else {
      datExpr <- if(is.null(data_obj@data)) data_obj@data else data_obj@raw.data
    }
  } else {
    datExpr <- if(is.null(data_obj@data)) data_obj@data else data_obj@raw.data
  }
  # Save as loomR
  data_obj <- create(filename = paste0(dir_scratch, prefix_out, "_datExpr.loom"), 
                     data = datExpr, 
                     #gene.attrs = list(gene_names = rownames(datExpr)),
                     #cell.attrs = list(cell_names = colnames(datExpr)),
                     display.progress = T, 
                     calc.numi = if (scale_data) T else F)

  #data_obj$add.col.attribute(attribute = list(cell_names = colnames(datExpr)))
  data_obj$add.col.attribute(attribute = list(ident = ident), overwrite = TRUE)
  #data_obj$add.row.attribute(attribute = list(gene_names = rownames(datExpr)))

  rm(datExpr)

} else if (fileType=="loom") {
  
  ident <- as.character(data_obj[[paste0("col_attrs/", metadata_subset_col)]][])
  data_obj$add.row.attribute(list(ident = ident), overwrite = TRUE)  
  # NB: cannot pull all metadata from loom object without specific names?
   
  #TODO: use add.layer?
  message("Normalising the data")
  tryCatch({NormalizeData(data_obj)
    }, error = function(err) warning("NormalizeData call on loom file failed")) # may fail if data is already normalized. #TODO does the loom object hold different versions of the data? E.g. layers?
  
  if (scale_data) {

    tryCatch({
    FindVariableGenes(object = data_obj)
    }, error = function(err) warning("FindVariableGenes call on loom file failed"))
    tryCatch({
    ScaleData(object = data_obj,
              genes.use = hv.genes,
              chunk.size = 500,
              display.progress = FALSE,
              vars.to.regress = "percent_mito")
    }, error = function(err) warning("ScaleData call on loom file failed"))
  }

}
#}


######################################################################
################ COMPUTE CELL MODULE EMBEDDINGS ######################
######################################################################

message("Computing cell module embeddings")

ncell <- dim(data_obj[["matrix"]][,])[1] # NB: the matrix is stored transposed

cell_chunk_idx <- c(seq(from=0, to = ncell, by=10000), ncell)

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
  
  if (file.exists(paste0(dir_RObjects_data, 
                         prefix_out, "_cellModEmbed_chunk.loom"))) file.remove(paste0(dir_RObjects_data, 
                                                                                      prefix_out, "_cellModEmbed_chunk.loom"))
  
  cellModEmbed_loom_chunk <- loomR::create(filename = paste0(dir_RObjects_data, 
                                                             prefix_out, "_cellModEmbed_chunk.loom"), 
                                           data = t(subEmbed), 
                                           #gene.attrs = list(gene_names = rownames(subEmbed), ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]),
                                           #cell.attrs = list(cell_names = colnames(subEmbed)),
                                           display.progress = T, 
                                           calc.numi = F)
  
  # NB: We flit dimensions again since the embeddings are cell * module
  #cellModEmbed_loom_chunk$add.row.attribute(list(cell_names = rownames(subEmbed)), overwrite = TRUE)
  #cellModEmbed_loom_chunk$add.row.attribute(list(ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]), overwrite = TRUE)
  #cellModEmbed_loom_chunk$add.col.attribute(list(module = colnames(subEmbed)), overwrite = TRUE)
  
  if (i==1) {
    if (file.exists(paste0(dir_RObjects_data, 
                           prefix_out, "_cellModEmbed.loom"))) file.remove(paste0(dir_RObjects_data, 
                                                                                        prefix_out, "_cellModEmbed.loom"))
    
    cellModEmbed_loom <- create(filename=paste0(dir_RObjects_data, 
                                                prefix_out, "_cellModEmbed.loom"), 
                                data=cellModEmbed_loom_chunk[["matrix"]][,], 
                                display.progress=T, 
                                calc.numi=F)
  } else {
    cellModEmbed_loom <- loomR::combine(list(cellModEmbed_loom, 
                                             cellModEmbed_loom_chunk), 
                                        filename=paste0(dir_RObjects_data, prefix_out, "_cellModEmbed.loom"), 
                                        chunk.size=1000, 
                                        overwrite=T, 
                                        display.progress=T)
  }
  
  cellModEmbed_loom$add.col.attribute(list(module = colnames(subEmbed)), overwrite = TRUE)
  cellModEmbed_loom$add.row.attribute(list(cell_name = rownames(subEmbed)), overwrite = TRUE)
  
  try({
    cellModEmbed_loom_chunk$close_all()
    invisible(file.remove(paste0(dir_RObjects_data, prefix_out, "_cellModEmbed_chunk.loom")))
    })
  

  rm(submat, subEmbed)
}

# AVERAGE EMBEDDINGS FOR EACH CELLTYPE 

for (i in 1:length(sort(unique(ident)))) {

  # cellModEmbed_loom$apply(name = paste0("row_attrs/",unique(ident)[i]), FUN = colMeans, MARGIN = 1, index.use = which(ident==unique(ident)[i]),
  #             dataset.use = "matrix", display.progress = FALSE, overwrite = TRUE)
  # 
  subEmbed <- t(cellModEmbed_loom[["matrix"]][,ident==sort(unique(ident))[i]])
  meanExpr <- matrix(parApply(cl, X=subEmbed, FUN=mean, MARGIN=2), nrow=1)
  rownames(meanExpr) <- sort(unique(ident))[i]
  colnames(meanExpr) <- cellModEmbed_loom[["col_attrs/module"]][]
  # it might be ok to do this w/o parallelising
  if (i==1) cellTypeModEmbed <- meanExpr else cellTypeModEmbed <- rbind(cellTypeModEmbed, meanExpr)
}

# cellTypeModEmbed <- t(parSapply(cl, X=paste0("row_attrs/",unique(ident)), FUN = function(index) cellModEmbed_loom[[index]][], simplify = T))
# rownames(cellTypeModEmbed) <- gsub("col_attrs/", "", rownames(cellTypeModEmbed))
# colnames(cellTypeModEmbed) <- 
stopCluster(cl)
#rownames(cellTypeModEmbed) = sort(unique(ident)); colnames(cellTypeModEmbed) <- names(module_u)

saveRDS(cellTypeModEmbed, file=paste0(dir_RObjects_data, prefix_out, "_cellTypeModEmbed.RDS"))

######################################################################
######################## PLOT CELL*MODULE MATRIX #####################
######################################################################
message("Plotting")
# make row annotation: random colors for celltypes
colors_uniq <- gsub("\\d", "", colors()) %>% unique

htra_colors <-   sample(if (length(ident)<length(colors_uniq)) colors_uniq else colors(), size=length(unique(ident)), replace=F) 

names(htra_colors) <- sort(unique(ident), decreasing=T)

htra <- rowAnnotation(cellcluster = sort(ident, decreasing=T),
                    annotation_legend_param = list(cellcluster = list(nrow = 7, title = "Cell cluster", title_position = "topcenter")),
                    col=list(cellcluster = htra_colors),
                    width = unit(5, "mm"))

# htca_colors <- sample( if (length(cellType_module_u)<length(colors_uniq)) colors_uniq else colors(), size=length(cellType_module_u), replace=F) 
# 
# names(htca_colors) <- names(cellType_module_u)
#   
# htca <- columnAnnotation(df = data.frame(WGCNA_cluster = names(cellType_module_u), row.names = htca_colors, stringsAsFactors = F),
#                          annotation_legend_param = list(WGCNA_cluster = list(nrow=7, title="WGCNA_cluster", title_position = "topcenter")),
#                          col = list(WGCNA_cluster=htca_colors),
#                          height = unit(5,"mm"))
  
ht1 <- Heatmap(t(cellModEmbed_loom[["matrix"]][,order(ident, 1:ncol(cellModEmbed_loom[["matrix"]][,]), decreasing=T)]), 
               cluster_rows = F,
               #row_order = ,
               cluster_columns = F, 
               #show_row_dend = F, 
               show_column_dend = F, 
               show_heatmap_legend = T, 
               show_row_names = F, 
               show_column_names = F,
               heatmap_legend_param = list(title = "Expression"), 
               #bottom_annotation = htca,
               #bottom_annotation_height = unit(5, "mm"),
               use_raster = T)
#top_annotation = eval(parse(text = paste0("list_htca_fdr", "[[", 1:length(list_htca_fdr),"]]"))))
pdf(sprintf("%s%s_cellModEmbed_cellOrder.pdf", dir_plots_data, prefix_out), h=max(10, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 5000), w=min(20, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2))
draw(ht1+htra)#+htca)
dev.off()

ht2 <- Heatmap(cellTypeModEmbed[,], 
               cluster_rows = F,
               #row_order = ,
               cluster_columns = F, 
               #show_row_dend = F, 
               show_column_dend = F, 
               show_heatmap_legend = T, 
               show_row_names = T, 
               show_column_names = F,
               heatmap_legend_param = list(title = "Expression"), 
               use_raster = T)#
               #top_annotation = htca)
#top_annotation = eval(parse(text = paste0("list_htca_fdr", "[[", 1:length(list_htca_fdr),"]]"))))
pdf(sprintf("%s%s_cellTypeModEmbed.pdf", dir_plots_data, prefix_out), h=min(6, nrow(cellTypeModEmbed) %/% 2), w=min(20, ncol(cellTypeModEmbed) %/% 2))
draw(ht2)
dev.off()

######################################################################
############################# WRAP UP ################################
######################################################################
rm(cellModEmbed)
data_obj$close_all()
cellModEmbed_loom$close_all()

message("DONE!")
