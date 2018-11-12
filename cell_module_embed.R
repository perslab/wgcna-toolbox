# Make cellType x module embedding matrix
# Use gene loadings from one data (sub)set and possibly expression matrix from another

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

ipak(c("optparse", "Seurat", "dplyr", "Biobase", "Matrix", "parallel", "loomR", "readr"))

######################################################################
########################### OPTPARSE #################################
######################################################################

option_list <- list(
  make_option("--dir_project_WGCNA", type="character",
              help = "Folder of WGCNA project, must contain /tables and /RObjects subdirectories."),  
  make_option("--paths_gene_module_score", type="character", default=NULL,
              help = "vector of paths to dataframes containing WGCNA module genes and scores in long format, i.e. one row per gene e.g. rwgcna cell_cluster_module_genes.csv files, e.g.''c('mousebrain='/projects/jonatan/tmp-mousebrain/tables/example1.csv','/projects/jonatan/tmp-maca/example2.csv)''. The order should match prefixes_WGCNA_run. If not provided, the script searches the WGCNA tables subdir for '_cell_cluster_module_genes.csv', [default %default]"),  
  make_option("--dir_project_data", type="character",
              help = "Folder of expression data project. The script sends output here. must contain /tables and /RObjects subdirectories, may be same as dir_project_WGCNA"),  
  make_option("--dir_scratch", type="character", default="/scratch/tmp-wgcna/",
              help = "Directory for temporary files, [default %default]"),
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
  make_option("--data_organism", type="character", default="mmusculus",
              help = "'hsapiens' or 'mmusculus', [default %default]"),
  make_option("--gene_names", type="character", default="hgnc|symbol|gene_name_optimal",
              help ="string or regex for grepping WGCNA outputs and gene mapping dataframe,s e.g. 'hgnc|symbol|gene_name' or 'ensembl', [default %default]"),
  make_option("--scale_center_regress_data", type="logical", default = "TRUE",
              help = "[default %default]. If TRUE, and path_data points to an RObject, the script will try to use the @scale.data slot and, if missing, will normalize and scale the data and regress out nUMI, percent.mito and percent.ribo. If FALSE, for RObjects will use the @data slot, or, if empty, the @raw.data slot"),  
  make_option("--scale_center_mods", type="logical", default = "TRUE",
              help = "scale expression of the modules? Recommended if the gene modules originate from different analyses with different scales. Set to F if desiring to scale at a later stage, e.g. after merging several embeddings matrices [default %default]."),  
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
              help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]"),
  make_option("--RAM_Gb_max", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
  
)

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

dir_project_WGCNA <- opt$dir_project_WGCNA 
dir_project_data <- opt$dir_project_data
paths_gene_module_score <- opt$paths_gene_module_score
if (!is.null(paths_gene_module_score)) paths_gene_module_score <- eval(parse(text=paths_gene_module_score))
dir_scratch <- opt$dir_scratch
path_data <- opt$path_data
path_metadata <- opt$path_metadata
prefixes_WGCNA_run <- eval(parse(text=opt$prefixes_WGCNA_run))
prefix_out <- opt$prefix_out 
if (!is.null(opt$datExpr_ident_cols)) datExpr_ident_cols <- eval(parse(text=opt$datExpr_ident_cols))
if (!is.null(opt$WGCNA_ident_cats)) WGCNA_ident_cats <- eval(parse(text=opt$WGCNA_ident_cats))
data_organism <- opt$data_organism
gene_names <- opt$gene_names
scale_center_regress_data <- opt$scale_center_regress_data
scale_center_mods <- opt$scale_center_mods
do_plot <- opt$do_plot
plot_PPI_enrichment <- opt$plot_PPI_enrichment
if (!is.null(opt$plot_LSscore_enrichment)) plot_LDscore_enrichment <- eval(parse(text=opt$plot_LSscore_enrichment))
if (!is.null(opt$plot_magma_enrichment)) plot_magma_enrichment <- eval(parse(text=opt$plot_magma_enrichment))
if (!is.null(opt$plot_metadata_corr)) plot_metadata_corr <- eval(parse(text=opt$plot_metadata_corr))
n_cores <- opt$n_cores
RAM_Gb_max <- opt$RAM_Gb_max
#use_loom <- F#opt$use_loom # For the initial Seurat processing, before computing embeddings. The script will still use loomR to compute the embedding matrix

if (do_plot) ipak("ComplexHeatmap")

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F)
use_loom <- F

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

######################################################################
######################### LOAD GENE LOADINGS VECTORS #################
######################################################################

run_cellClusterModuleGenes <- list()
for (i in 1:length(prefixes_WGCNA_run)) {
  if (is.null(paths_gene_module_score)) {
    cellClusterModuleGenes_path <- dir(path = dir_tables_WGCNA, pattern = paste0(prefixes_WGCNA_run[i], ".*_cell_cluster_module_genes\\.csv"), full.names = T)
  } else {
    cellClusterModuleGenes_path <- paths_gene_module_score[i]
  }
  run_cellClusterModuleGenes[[prefixes_WGCNA_run[i]]] <- load_obj(f = cellClusterModuleGenes_path)
}

# If there are duplicate modules between WGCNA runs, prefix module names with the WGCNA run 
if (sapply(X=run_cellClusterModuleGenes, function(df) unique(df[["module"]]), simplify = T) %>% duplicated %>% any) {
  for (i in 1:length(run_cellClusterModuleGenes)) {
    run_cellClusterModuleGenes[[i]][["module"]] <- paste0(prefixes_WGCNA_run[i], "_",run_cellClusterModuleGenes[[i]][["module"]]) # ensure modules from different runs remain distinct
  }
}

# merge into one long data frame
runCellClusterModuleGenes <- Reduce(x = run_cellClusterModuleGenes,f = rbind)
rm(run_cellClusterModuleGenes)

# filter out very small modules. This should have been done in WGCNA main script already
mods_too_small <- names(table(runCellClusterModuleGenes[["module"]]))[table(runCellClusterModuleGenes[["module"]])<5]
runCellClusterModuleGenes <- runCellClusterModuleGenes[!runCellClusterModuleGenes[["module"]] %in% mods_too_small,]
###########

gene_col <- grep(gene_names, colnames(runCellClusterModuleGenes), value=T)
pkM_col <- grep("^p*kI*ME*$", colnames(runCellClusterModuleGenes), ignore.case=T, value=T)

# run_cellType_module_u <- list()
# WGCNA_cellTypes = c()
# for (prefix_run in prefixes_WGCNA_run) {
#   cellType_module_u_path <- dir(path = dir_RObjects_WGCNA, pattern = paste0(prefix_run, "_list_list_module_u.RDS"), full.names = T)
#   run_cellType_module_u[[prefix_run]] <- load_obj(f = cellType_module_u_path)
#   WGCNA_cellTypes <- c(WGCNA_cellTypes, names(run_cellType_module_u[[prefix_run]]))
# }
# 
# names(run_cellType_module_u) <- prefixes_WGCNA_run

######################################################################
######################### LOAD METADATA ############################
######################################################################

if (!is.null(path_metadata)){ 
  message("Loading metadata") 
  metadata <- load_obj(f=path_metadata)
}

######################################################################
######################### LOAD EXPRESSION MATRIX #####################
######################################################################

message("Loading expression matrix")

# load/open connection to expression matrix
# NB: may be huge

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

# Add any further metadata coming from another file
if (!is.null(path_metadata)) {
  if (all_equal(rownames(metadata_tmp), colnames(data_obj@raw.data))) {
    data_obj <- AddMetaData(object = data_obj, metadata = data.frame(metadata))
  } 
}

######################################################################
#################### MAP ENSEMBL TO SYMBOL ###########################
######################################################################

mapping = if (data_organism=="hsapiens") {
  load_obj("/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz")
} else if (data_organism=="mmusculus") {
  load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
}  

if (any(grepl("ENSG|ENSMUSG",rownames(data_obj@raw.data)))) {
  
  current_ensembl <- T
  
  tmp <- gene_map(df = data_obj@raw.data,
                  idx_gene_column = NULL,
                  mapping=mapping,
                  from="ensembl",
                  to="gene_name_optimal",
                  replace = T,
                  na.rm = T)
} else {
  current_ensembl <- F
  tmp <- data_obj@raw.data
}

message("Counting percent ribo and percent mito")

idx_mito.genes <- grepl(pattern = "^mt-", x = rownames(tmp), ignore.case=T)
idx_ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(tmp), ignore.case=T)

percent_mito <- Matrix::colSums(tmp[idx_mito.genes, ])/Matrix::colSums(tmp)
percent_ribo <- Matrix::colSums(tmp[idx_ribo.genes, ])/Matrix::colSums(tmp)

data_obj <- AddMetaData(object = data_obj, metadata = percent_mito, col.name = "percent_mito")
data_obj <- AddMetaData(object = data_obj, metadata = percent_ribo, col.name = "percent_ribo")

if (grepl("ensembl", gene_names) & !current_ensembl) {
  data_obj@raw.data <- gene_map(df = data_obj@raw.data,
                                idx_gene_column = NULL,
                                mapping=mapping,
                                from="gene_name_optimal",
                                to="ensembl",
                                replace = T,
                                na.rm = T)
} else if (grepl("symbol|hgcn", gene_names) & current_ensembl) {
  data_obj@raw.data <- tmp
}

rm(tmp, mapping, idx_mito.genes, idx_ribo.genes, percent_ribo, percent_mito)

######################################################################
######################### PREPARE DATA ###############################
######################################################################

# TODO: move the datExpr processing into loom
# will need to normalize and regress out ourselves but shouldn't be an issue.

message("Normalising the data")

data_obj <- NormalizeData(data_obj)

if (scale_center_regress_data){
  
  message("Scaling the data and regressing out confounders percent.mito, percent.ribo and nUMI")
  
  data_obj<-ScaleData(data_obj,
                      vars.to.regress = c("percent_mito", "nUMI", "percent_ribo")[c("percent_mito", "nUMI", "percent_ribo") %in% colnames(data_obj@meta.data)], 
                      do.scale=scale_center_regress_data,
                      do.center=scale_center_regress_data,
                      display.progress= T, 
                      do.par=scale_center_regress_data, 
                      num.cores=n_cores)
}
ident <- if (!is.null(datExpr_ident_cols)) data_obj@meta.data[[datExpr_ident_cols[1]]] else data_obj@ident

datExpr <- if (scale_center_regress_data) data_obj@scale.data else  data_obj@data

metadata <- data_obj@meta.data

rm(data_obj)

invisible(gc()); invisible(R.utils::gcDLLs())

######################################################################
################ COMPUTE CELL MODULE EMBEDDINGS ######################
######################################################################

message("Computing cell module embeddings")

if (!use_loom) {
  
  ncell <- ncol(datExpr) 
  
  cell_chunk_idx <- c(seq.int(from=0, to = ncell, by=min(5000, ncell %/% 1.9)))
  
  #cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "cell_x_module_matrix.txt"), timeout=30)
  
  for (i in 1:(length(cell_chunk_idx)-1)) { # loop over chunks of cell
    
    submat <- datExpr[,(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]] %>% as.matrix
    
    rownames(submat) <- rownames(datExpr) # all genes
    
    invisible(gc())
    # Compute embeddings
    subEmbed <- sapply(X=unique(runCellClusterModuleGenes[["module"]]), FUN=function(module) {
      idx_mod <- which(runCellClusterModuleGenes[["module"]]==module)
      idx_match <- match(runCellClusterModuleGenes[[gene_col]][idx_mod], rownames(submat)) 
      t(submat[idx_match[!is.na(idx_match)],]) %*% matrix(runCellClusterModuleGenes[[pkM_col]][idx_mod[!is.na(idx_match)]], nrow = length(idx_mod[!is.na(idx_match)]))
    },
    simplify=T)
    
    rownames(subEmbed) <- colnames(datExpr)[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]] #fixed from <- rownames... 181029
    colnames(subEmbed) <- unique(runCellClusterModuleGenes[["module"]])
    
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
  
  if (scale_center_mods) {
    message("Scaling the module embeddings")
    # scale the embeddings by subtracting the mean and dividing by standard deviation and 
    cellModEmbed_loom$apply(name = "layers/scale.data", FUN = function(module_vec) {
      module_vec %>% '-'(mean(.)) %>% '/'(stats::sd(.)) -> out
      out
    }, MARGIN = 2, chunk.size = 5000, 
    dataset.use = "matrix", display.progress = T, overwrite = T)
  }    
  #try(stopCluster(cl))
  
} #else if (use_loom) {

# AVERAGE EMBEDDINGS FOR datExpr_ident_cols

if (FALSE) {
  if (!is.null(datExpr_ident_cols)) {
    
    message(paste0("Computing expression averaged over ",paste0(datExpr_ident_cols, collapse=", ")))
    
    # remove NAs from ident columns
    
    for (ident_col in datExpr_ident_cols) {
      metadata[[ident_col]][is.na(metadata[[ident_col]])] <- "Not_attributed"
    }
    
    list_subsetModEmbed <- lapply(datExpr_ident_cols, function(ident_col) {
      
      cl <- makeCluster(spec=min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = .5)-1), type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", ident_col, "_mean_mod_expr.txt"))
      
      for (i in 1:length(sort(unique(ident)))) {
        # cellModEmbed_loom$apply(name = paste0("row_attrs/",unique(ident)[i]), FUN = colMeans, MARGIN = 1, index.use = which(ident==unique(ident)[i]),
        #             dataset.use = "matrix", display.progress = FALSE, overwrite = TRUE)
        ident <- metadata[[ident_col]]
        subEmbed <- if (scale_center_mods) cellModEmbed_loom[["layers/scale.data"]][ident==sort(unique(ident))[i],] else cellModEmbed_loom[["matrix"]][ident==sort(unique(ident))[i],] 
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
  }
}
# cellTypeModEmbed <- t(parSapply(cl, X=paste0("row_attrs/",unique(ident)), FUN = function(index) cellModEmbed_loom[[index]][], simplify = T))
# rownames(cellTypeModEmbed) <- gsub("col_attrs/", "", rownames(cellTypeModEmbed))
# colnames(cellTypeModEmbed) <- 

#rownames(cellTypeModEmbed) = sort(unique(ident)); colnames(cellTypeModEmbed) <- names(module_u)

######################################################################
######################## CLUSTER CELLS, MODULES  #####################
######################################################################
if (do_plot) {
  # TODO: plot these with scaled data from cellModEmbed_loom[["layers/scale.data"]]
  
  list_subsetModEmbed <- NULL
  if (FALSE){
    modDist <- dist(x=t(cellModEmbed_loom[["matrix"]][,]), method="euclidean", diag=F, upper=F) 
    modDendro <- hclust(d=modDist, method="average")
    modClust_order <- modDendro$order
    
    # cellDist <- dist(x=cellModEmbed_loom[["matrix"]][,], method="euclidean", diag=F, upper=F) 
    # cellDendro <- hclust(d=cellDist, method="average")
    # cellClust_order <- cellDendro$order
    
    cl <- makeCluster(spec=min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = as.numeric(max(sapply(list_subsetModEmbed, object.size())))/1024^3)-1), type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", ident_col, "_mean_mod_expr_clust.txt"))
    
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
  }
}
######################################################################
#################### PLOT CELL*MODULE MATRIX #########################
######################################################################

if (do_plot) {
  # generate colors for annotation
  colors_uniq <- gsub("\\d", "", colors()) %>% unique
  
  ######################## ROW ANNOTATION ##############################
  list_htra <- list_htra_link <- NULL 
  
  if (!is.null(datExpr_ident_cols)) {
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
  }
  ######################## COLUMN ANNOTATION ###########################
  
  #  WGCNA run column annotation
  WGCNA_run <- unlist(mapply(function(names_run, cellType_module_u) rep(names_run, times= length(unlist(cellType_module_u, recursive = F))), names_run=prefixes_WGCNA_run, cellType_module_u=run_cellType_module_u, SIMPLIFY=T))
  htca_WGCNA_run_colors <- sample(if(length(prefixes_WGCNA_run)<length(colors_uniq)) colors_uniq else colors(), size=length(prefixes_WGCNA_run), replace=F) 
  names(htca_WGCNA_run_colors) <- sort(unique(prefixes_WGCNA_run, decreasing=F))
  
  # WGCNA celltype column annotation
  WGCNA_cellCluster <- rep(WGCNA_cellTypes,times= sapply(cellType_module_u, length))
  
  htca_WGCNA_cellCluster_colors <- sample(if (length(WGCNA_cellTypes)<length(colors_uniq)) colors_uniq else colors(), size=length(WGCNA_cellTypes), replace=F) 
  names(htca_WGCNA_cellCluster_colors) <- sort(WGCNA_cellTypes, decreasing=F) # no need to sort?
  
  # Prepare WGCNA celltype column annotation labels as separate column annotation
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
  
  cl <- makeCluster(spec=min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = 0.2)-1), type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "make_PPI_column_annotation.txt"))
  
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
    
    cl <- makeCluster(spec=min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = 0.2)-1), type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "make_magma_column_annotation.txt"))
    
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
message("Writing files to disk")
# Close and delete intermediate loom files
#try(data_obj$close_all())
#if (file.exists(paste0(dir_scratch, prefix_out, "_data_obj.loom"))) try(invisible(file.remove(paste0(dir_scratch, prefix_out, "_data_obj.loom"))))

# try(cellModEmbed_loom_chunk$close_all())
# if (file.exists(paste0(dir_scratch, 
#                        prefix_out, "_cellModEmbed_chunk.loom"))) try(invisible(file.remove(paste0(dir_scratch, 
#                                                                                                   prefix_out, "_cellModEmbed_chunk.loom"))))                  
cellModEmbed <- data.frame(data=rep(prefix_out, times = length(cellModEmbed_loom[["col_attrs/ident"]][])),
                           ident=cellModEmbed_loom[["col_attrs/ident"]][], 
                           cell_names=cellModEmbed_loom[["col_attrs/cell_names"]][], 
                           if (scale_center_mods) cellModEmbed_loom[["layers/scale.data"]][,] else cellModEmbed_loom[["matrix"]][,])
#cellModEmbed[,2:length(cellModEmbed_loom[["row_attrs/gene_names"]][])+2] <- if (scale_center_mods) cellModEmbed_loom[["layers/scale.data"]][,] else cellModEmbed_loom[["matrix"]][,] 
#cellModEmbed <- cbind(cellModEmbed_loom[["col_attrs/ident"]][], cellModEmbed_loom[["col_attrs/cell_names"]][], cellModEmbed)
#rownames(cellModEmbed) <- cellModEmbed_loom[["col_attrs/cell_names"]][]
colnames(cellModEmbed)[4:ncol(cellModEmbed)] <- cellModEmbed_loom[["row_attrs/module"]][]  

write_csv(x=as.data.frame(cellModEmbed), path= paste0(dir_tables_data, prefix_out, "_cellModEmbed.csv"))
rm(cellModEmbed)

# Write out metadatato separate file to cbind in unix (using paste -d',' file1 file2)
# TODO - better way might be to add to loom matrix first?

# write_csv(x=data.frame(ident=cellModEmbed_loom[["col_attrs/ident"]][], cell_names=cellModEmbed_loom[["col_attrs/cell_names"]][]),
#           path= paste0(dir_tables_data, prefix_out, "_full_cellModEmbed_cell_annot.csv"))

# write_csv(x=data.frame(ident=cellModEmbed_loom[["col_attrs/ident"]][], cell_names=cellModEmbed_loom[["col_attrs/cell_names"]][]), 
#           path= paste0(dir_tables_data, prefix_out, "_full_cellModEmbed_cell_annot.csv"))

# Save the embeddings averaged over tissue, celltype etc
if (FALSE) {
  if (!is.null(datExpr_ident_cols)) {
    for (datExpr_ident in datExpr_ident_cols) {
      write_csv(list_subsetModEmbed[[datExpr_ident]], file =  paste0(dir_tables_data, prefix_out, "_", datExpr_ident, "_cellModEmbed.csv") )
    }
  }
}

try(cellModEmbed_loom$close_all())
try(invisible(file.remove(paste0(dir_scratch, prefix_out, "_cellModEmbed.loom"))))

#date = substr(gsub("-","",as.character(Sys.Date())),3,1000)
#if (!is.null(list_subsetModEmbed)) saveRDS(list_subsetModEmbed, file = paste0(dir_RObjects_data, prefix_out, "list_subsetModEmbed.RDS"), compress = "gzip")
#save.image(paste0(dir_scratch, prefix_out, "_", date, "_session_image.RData"))               

message("DONE!")
