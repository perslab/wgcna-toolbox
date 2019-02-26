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
#current.dir = "/projects/jonatan/tools/wgcna-src/wgcna-toolbox/postWGCNAscripts/"
source(file=paste0(gsub("postWGCNAscripts/","WGCNApipeline/",current.dir),"rwgcna_functions_seurat3.0.R"))

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("Seurat", "optparse", "dplyr", "Biobase", "Matrix", "parallel", "loomR", "readr"))

#stopifnot(as.character(packageVersion("Seurat"))=='3.0.0.9000')

######################################################################
########################### OPTPARSE #################################
######################################################################

option_list <- list(
  make_option("--path_datExpr", type="character",
              help = "Path to expression data in Seurat loomR or RObject format to score on WGCNA modules."),  
  make_option("--vec_paths_df_NWA", type="character", default=NULL,
              help = "named vector of paths to dataframes containing WGCNA module genes and scores in long format, i.e. one row per gene e.g. rwgcna cell_cluster_module_genes.csv files, e.g.''c('mousebrain='/projects/jonatan/tmp-mousebrain/tables/example1.csv', 'tab_muris' = '/projects/jonatan/tmp-maca/example2.csv)''"),  
  make_option("--colGeneNames", type="character",
              help = "nwa_df column with gene names"), 
  make_option("--colGeneWeights", type="character",
              help = "nwa_df column with gene weights"), 
  make_option("--colModule", type="character", default="module_merged",
              help = "nwa_df column with module assignment [default %default]"), 
  make_option("--dirOut", type="character",
              help = "Folder of expression data project. The script sends output here. must contain /tables and /RObjects subdirectories, may be same as dir_project_WGCNA"),  
  make_option("--prefixOut", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"), 
  make_option("--minGeneClusterSize", type="integer", default=10L,
              help="What is the minimum number of genes in a module to continue? Integer, [default %default]."),
  make_option("--normalizeGeneWeights", type="logical", default=F,
              help="normalize gene weights to sum to 1?, [default %default]."),
  make_option("--scaleToZScore", type="logical", default = TRUE,
              help = "scale and center module expression to make it into a Z score? Recommended if using kIMs as gene weights as their scale depends on the WGCNA softPower parameter [default %default]."),  
  make_option("--assaySlot", type="character", default = NULL,
              help = "If datExpr is stored in a Seurat object, indicate whether to use raw.data/counts, data or scale.data [default %default]."),  
  make_option("--dirScratch", type="character", default="/scratch/tmp-wgcna/",
              help = "Directory for temporary files, [default %default]"),
  make_option("--RAMGbMax", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]"),
  make_option("--path_runLog", type="character", default=NULL,
              help = "Path to file to log the run and the git commit. If left as NULL, write to a file called runLog.text in the dirLog [default %default]")
  
)

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

path_datExpr <- opt$path_datExpr
vec_paths_df_NWA <- opt$vec_paths_df_NWA
if (!is.null(vec_paths_df_NWA)) vec_paths_df_NWA <- eval(parse(text=vec_paths_df_NWA))
colGeneNames <- opt$colGeneNames
colGeneWeights <- opt$colGeneWeights
colModule <- opt$colModule
dirOut <- opt$dirOut
prefixOut <- opt$prefixOut
minGeneClusterSize <- opt$minGeneClusterSize
normalizeGeneWeights <- opt$normalizeGeneWeights
scaleToZScore <- opt$scaleToZScore
assaySlot <- opt$assaySlot
dirScratch <- opt$dirScratch
RAMGbMax <- opt$RAMGbMax
path_runLog <- opt$path_runLog

use_loom = F # NB!!
######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F, na.action = "na.omit")

######################################################################
###################### SET AND VERIFY PATHS ##########################
######################################################################

#if (!file.exists(dir_project_WGCNA)) stop("dir_project_WGCNA not found")
if (!file.exists(dirOut)) stop("dirOut not found")
if (!file.exists(path_datExpr)) stop("path_datExpr not found")
#if (!is.null(path_metadata)) if (!file.exists(path_metadata)) stop("path_datExpr not found")

# set up and verify WGCNA directory paths

#dir_tables_WGCNA = paste0(dir_project_WGCNA,"tables/")
#if (!file.exists(dir_tables_WGCNA)) stop("dir_project WGCNA must contain a /tables subdir")

#dir_RObjects_WGCNA = paste0(dir_project_WGCNA,"RObjects/")
#if (!file.exists(dir_RObjects_WGCNA)) stop("dir_project_WGCNA must contain a /RObjects subdir")

# set up expression data / output directories 

dirPlots = paste0(dirOut,"plots/")
if (!file.exists(dirPlots)) dir.create(dirPlots) 

dirTables = paste0(dirOut,"tables/")
if (!file.exists(dirTables)) dir.create(dirTables)

dirRObjects = paste0(dirOut,"RObjects/")
if (!file.exists(dirRObjects)) dir.create(dirRObjects)

dirLog = paste0(dirOut,"log/")
if (!file.exists(dirLog)) dir.create(dirLog)

######################################################################
######################### LOAD GENE LOADINGS VECTORS #################
######################################################################

run_cellClusterModuleGenes <- lapply(vec_paths_df_NWA, read.csv)

# If there are duplicate modules between WGCNA runs, prefix module names with the WGCNA run 
if (length(run_cellClusterModuleGenes)>1) {
  if (sapply(X=run_cellClusterModuleGenes, function(df) unique(df[[colModule]]), simplify = T) %>% duplicated %>% any) {
    for (i in 1:length(run_cellClusterModuleGenes)) {
      run_cellClusterModuleGenes[[i]][[colModule]] <- paste0(names(vec_paths_df_NWA)[i], "_",run_cellClusterModuleGenes[[i]][[colModule]]) # ensure modules from different runs remain distinct
    }
  }
}

# merge into one long data frame
df_NWA <- if (length(run_cellClusterModuleGenes)>1) Reduce(x = run_cellClusterModuleGenes,f = rbind) else run_cellClusterModuleGenes[[1]]
rm(run_cellClusterModuleGenes)

# filter out very small modules. This should have been done in WGCNA main script already
mods_too_small <- names(table(df_NWA[[colModule]]))[table(df_NWA[[colModule]])<minGeneClusterSize]
df_NWA <- df_NWA[!df_NWA[[colModule]] %in% mods_too_small,]

df_NWA <- df_NWA[!is.na(df_NWA[[colModule]]),]

if (normalizeGeneWeights) {
  for (module in unique(df_NWA[[colModule]])) {
    normFactor <- sum(df_NWA[[colGeneWeights]][df_NWA[[colModule]]==module])
    df_NWA[[colGeneWeights]][df_NWA[[colModule]]==module] <- df_NWA[[colGeneWeights]][df_NWA[[colModule]]==module]/normFactor
  }
}
###########


# colGeneNames <- grep(gene_names, colnames(df_NWA), value=T)
# colGeneWeights <- grep("^p*kI*ME*$", colnames(df_NWA), ignore.case=T, value=T)

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

# if (!is.null(path_metadata)){ 
#   message("Loading metadata") 
#   metadata <- load_obj(f=path_metadata)
# }

######################################################################
######################### LOAD EXPRESSION MATRIX #####################
######################################################################

message("Loading expression matrix")

# load/open connection to expression matrix
# NB: may be huge

if (grepl(pattern = "\\.loom", path_datExpr)) {
  dataObj <-  connect(filename = path_datExpr, mode = "r+")
} else {
  dataObj <- load_obj(path_datExpr)
} 

if (!any(c("seurat","loom") %in% class(dataObj))) { # not a loom object, nor a seurat object
  #dataObj <- CreateSeuratObject(raw.data=dataObj, project= prefixOut, min.cells = -Inf, min.genes = -Inf)
  datExpr <- dataObj
} else if ("seurat" %in% class(dataObj)) {
  if (dataObj@version=='3.0.0.9000') {
    datExpr <- GetAssayData(dataObj, assay.type="RNA", slot=assaySlot) 
    } else {
      if (assaySlot=="counts") assaySlot <- "raw.data"
      datExpr <- eval(parse(text=paste0("dataObj@", assaySlot))) %>% as.matrix
    }
} else if ("loom" %in% class(dataObj)) {
  datExpr <- dataObj[["matrix"]][,]
  #dataObj <- Seurat::Convert(from=dataObj, to="seurat")
} else {
  datExpr <- dataObj[,-1] 
  rownames(datExpr) <- dataObj[,1]
}
rm (dataObj)

# Add any further metadata coming from another file
# if (!is.null(path_metadata)) {
#   if (all_equal(rownames(metadata_tmp), colnames(dataObj@raw.data))) {
#     dataObj <- AddMetaData(object = dataObj, metadata = data.frame(metadata))
#   } 
# }

######################################################################
#################### MAP ENSEMBL TO SYMBOL ###########################
######################################################################

# todo: can this be optional? (if we have a Seurat object with percent.mito and ribo already counted?)
# mapping = if (data_organism=="hsapiens") {
#   load_obj("/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz")
# } else if (data_organism=="mmusculus") {
#   load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
# }  
# 
# if (any(grepl("ENSG|ENSMUSG",rownames(dataObj@raw.data)))) {
#   
#   current_ensembl <- T
#   
#   tmp <- gene_map(df = dataObj@raw.data,
#                   idx_colGeneNamesumn = NULL,
#                   mapping=mapping,
#                   from="ensembl",
#                   to="gene_name_optimal",
#                   replace = T,
#                   na.rm = T)
# } else {
#   current_ensembl <- F
#   tmp <- dataObj@raw.data
# }
# 
# message("Counting percent ribo and percent mito")
# 
# idx_mito.genes <- grepl(pattern = "^mt-", x = rownames(tmp), ignore.case=T)
# idx_ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(tmp), ignore.case=T)
# 
# percent_mito <- Matrix::colSums(tmp[idx_mito.genes, ])/Matrix::colSums(tmp)
# percent_ribo <- Matrix::colSums(tmp[idx_ribo.genes, ])/Matrix::colSums(tmp)
# 
# dataObj <- AddMetaData(object = dataObj, metadata = percent_mito, col.name = "percent_mito")
# dataObj <- AddMetaData(object = dataObj, metadata = percent_ribo, col.name = "percent_ribo")
# 
# if (grepl("ensembl", gene_names) & !current_ensembl) {
#   dataObj@raw.data <- gene_map(df = dataObj@raw.data,
#                                 idx_colGeneNamesumn = NULL,
#                                 mapping=mapping,
#                                 from="gene_name_optimal",
#                                 to="ensembl",
#                                 replace = T,
#                                 na.rm = T)
# } else if (grepl("symbol|hgcn", gene_names) & current_ensembl) {
#   dataObj@raw.data <- tmp
# }
# 
# rm(tmp, mapping, idx_mito.genes, idx_ribo.genes, percent_ribo, percent_mito)
######################################################################
######################### PREPARE DATA ###############################
######################################################################

# TODO: move the datExpr processing into loom
# will need to normalize and regress out ourselves but shouldn't be an issue.

# message("Normalising the data")
# 
# dataObj <- NormalizeData(dataObj)
# 
# if (scale_center_regress_data){
#   
#   message("Scaling the data and regressing out confounders percent.mito, percent.ribo and nUMI")
#   
#   dataObj<-ScaleData(dataObj,
#                       vars.to.regress = c("percent_mito", "nUMI", "percent_ribo")[c("percent_mito", "nUMI", "percent_ribo") %in% colnames(dataObj@meta.data)], 
#                       do.scale=scale_center_regress_data,
#                       do.center=scale_center_regress_data,
#                       display.progress= T, 
#                       do.par=scale_center_regress_data, 
#                       num.cores=n_cores)
# }
# ident <- if (!is.null(datExpr_ident_cols)) dataObj@meta.data[[datExpr_ident_cols[1]]] else dataObj@ident
# 
# datExpr <- if (scale_center_regress_data) dataObj@scale.data else  dataObj@data
# 
# metadata <- dataObj@meta.data
# 
# rm(dataObj)
# 
# invisible(gc()); invisible(R.utils::gcDLLs())


######################################################################
################ COMPUTE CELL MODULE EMBEDDINGS ######################
######################################################################

message("Computing cell module embeddings")

if (!use_loom) {
  
  ncell <- ncol(datExpr) 
  
  cell_chunk_idx <- c(seq.int(from=0, to = ncell, by=min(5000, ncell %/% 1.9)))
  
  #cl <- makeCluster(spec=n_cores, type="FORK", outfile = paste0(dirLog, prefixOut,"_", "cell_x_module_matrix.txt"), timeout=30)
  
  for (i in 1:(length(cell_chunk_idx)-1)) { # loop over chunks of cell
    
    submat <- datExpr[,(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]] %>% as.matrix
    
    rownames(submat) <- rownames(datExpr) # all genes
    
    invisible(gc())
    # Compute embeddings
    subEmbed <- sapply(X=unique(df_NWA[[colModule]]), FUN=function(module) {
      idx_mod <- which(df_NWA[[colModule]]==module)
      idx_match <- match(df_NWA[[colGeneNames]][idx_mod], rownames(submat)) 
      t(submat[idx_match[!is.na(idx_match)],]) %*% matrix(df_NWA[[colGeneWeights]][idx_mod[!is.na(idx_match)]], nrow = length(idx_mod[!is.na(idx_match)]))
    },
    simplify=T)
    
    rownames(subEmbed) <- colnames(datExpr)[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]] #fixed from <- rownames... 181029
    colnames(subEmbed) <- unique(df_NWA[[colModule]])
    
    if (i==1) {
      
      if (file.exists(paste0(dirScratch, 
                             prefixOut, "_cellModEmbed.loom"))) invisible(file.remove(paste0(dirScratch, prefixOut, "_cellModEmbed.loom")))
      
      cellModEmbed_loom <- loomR::create(filename = paste0(dirScratch, 
                                                           prefixOut, "_cellModEmbed.loom"), 
                                         data = t(subEmbed), 
                                         #gene.attrs = list(gene_names = rownames(subEmbed), ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]),
                                         #cell.attrs = list(cell_names = colnames(subEmbed)),
                                         display.progress = T, 
                                         calc.numi = F,
                                         overwrite=T)
      
      # NB: We flit dimensions again since the embeddings are cell * module
      #cellModEmbed_loom_chunk$add.row.attribute(list(cell_names = rownames(subEmbed)), overwrite = TRUE)
      
      #cellModEmbed_loom$add.col.attribute(list(ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]), overwrite = TRUE)
      cellModEmbed_loom$add.row.attribute(list(gene_names = colnames(subEmbed)), overwrite = TRUE)
      cellModEmbed_loom$add.row.attribute(list(module = colnames(subEmbed)), overwrite = TRUE)
      
    } else {
      
      cellModEmbed_loom_chunk <- loomR::create(filename = paste0(dirScratch, 
                                                                 prefixOut, "_cellModEmbed_chunk.loom"), 
                                               data = t(subEmbed), 
                                               #gene.attrs = list(gene_names = rownames(subEmbed), ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]),
                                               #cell.attrs = list(cell_names = colnames(subEmbed)),
                                               display.progress = T, 
                                               calc.numi = F,
                                               overwrite=T)
      
      #cellModEmbed_loom_chunk$add.col.attribute(list(ident = ident[(cell_chunk_idx[i]+1):cell_chunk_idx[i+1]]), overwrite = TRUE)
      cellModEmbed_loom_chunk$add.row.attribute(list(module = colnames(subEmbed)), overwrite = TRUE)
      
      try({
        cellModEmbed_loom$close_all()
        cellModEmbed_loom_chunk$close_all()
      })
      
      cellModEmbed_comb <- loomR::combine(looms=list(paste0(dirScratch, 
                                                            prefixOut, "_cellModEmbed.loom"), 
                                                     paste0(dirScratch, 
                                                            prefixOut, "_cellModEmbed_chunk.loom")), 
                                          filename=paste0(dirScratch, prefixOut, "_cellModEmbed_comb.loom"), 
                                          chunk.size=1000, 
                                          overwrite=T, 
                                          display.progress=T)
      
      
      file.remove(paste0(dirScratch, 
                         prefixOut, "_cellModEmbed.loom"), 
                  paste0(dirScratch, 
                         prefixOut, "_cellModEmbed_chunk.loom"))
      
      cellModEmbed_loom <- loomR::create(filename=paste0(dirScratch, 
                                                         prefixOut, "_cellModEmbed.loom"), 
                                         data = t(cellModEmbed_comb[["matrix"]][,])
                                         #gene.attrs = list(modules= cellModEmbed_comb[["row_attrs/gene_names"]]),
                                         #cell.attrs = list(ident= cellModEmbed_comb[["col_attrs/ident"]])
      )
      cellModEmbed_loom$add.row.attribute(list(gene_names=cellModEmbed_comb[["row_attrs/gene_names"]][]), overwrite = TRUE)
      cellModEmbed_loom$add.row.attribute(list(module=cellModEmbed_comb[["row_attrs/module"]][]), overwrite = TRUE)
      cellModEmbed_loom$add.col.attribute(list(cell_names=cellModEmbed_comb[["col_attrs/cell_names"]][]), overwrite = TRUE)
      #cellModEmbed_loom$add.col.attribute(list(ident=cellModEmbed_comb[["col_attrs/ident"]][]), overwrite = TRUE)
      
      try({
        cellModEmbed_comb$close_all()
        invisible(file.remove(paste0(dirScratch, prefixOut, "_cellModEmbed_comb.loom")))
      })
      
    }
    
    rm(submat, subEmbed)
  }
  
  if (scaleToZScore) {
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

######################################################################
######################### WRITE FILES TO DISK ########################
######################################################################

message("Writing files to disk")
# Close and delete intermediate loom files
#try(dataObj$close_all())
#if (file.exists(paste0(dirScratch, prefixOut, "_dataObj.loom"))) try(invisible(file.remove(paste0(dirScratch, prefixOut, "_dataObj.loom"))))

# try(cellModEmbed_loom_chunk$close_all())
# if (file.exists(paste0(dirScratch, 
#                        prefixOut, "_cellModEmbed_chunk.loom"))) try(invisible(file.remove(paste0(dirScratch, 
#                                                                                                   prefixOut, "_cellModEmbed_chunk.loom"))))                  
cellModEmbed <- if (scaleToZScore) cellModEmbed_loom[["layers/scale.data"]][,] else cellModEmbed_loom[["matrix"]][,]
#data.frame(#data=rep(prefixOut, times = length(cellModEmbed_loom[["col_attrs/ident"]][])),
                           #ident=cellModEmbed_loom[["col_attrs/ident"]][], 
                           #cell_names=cellModEmbed_loom[["col_attrs/cell_names"]][], 
                           #if (scaleToZScore) cellModEmbed_loom[["layers/scale.data"]][,] else cellModEmbed_loom[["matrix"]][,])
rownames(cellModEmbed) <- cellModEmbed_loom[["col_attrs/cell_names"]][]
#cellModEmbed[,2:length(cellModEmbed_loom[["row_attrs/gene_names"]][])+2] <- if (scaleToZScore) cellModEmbed_loom[["layers/scale.data"]][,] else cellModEmbed_loom[["matrix"]][,] 
#cellModEmbed <- cbind(cellModEmbed_loom[["col_attrs/ident"]][], cellModEmbed_loom[["col_attrs/cell_names"]][], cellModEmbed)
#rownames(cellModEmbed) <- cellModEmbed_loom[["col_attrs/cell_names"]][]
#colnames(cellModEmbed)[4:ncol(cellModEmbed)] <- cellModEmbed_loom[["row_attrs/module"]][]  
colnames(cellModEmbed) <- cellModEmbed_loom[["row_attrs/module"]][]  

write_csv(x=as.data.frame(cellModEmbed), path= paste0(dirTables, prefixOut, "_cellModEmbed.csv"))

rm(cellModEmbed)

# Write out metadatato separate file to cbind in unix (using paste -d',' file1 file2)
# TODO - better way might be to add to loom matrix first?

# write_csv(x=data.frame(ident=cellModEmbed_loom[["col_attrs/ident"]][], cell_names=cellModEmbed_loom[["col_attrs/cell_names"]][]),
#           path= paste0(dirTables, prefixOut, "_full_cellModEmbed_cell_annot.csv"))

# write_csv(x=data.frame(ident=cellModEmbed_loom[["col_attrs/ident"]][], cell_names=cellModEmbed_loom[["col_attrs/cell_names"]][]), 
#           path= paste0(dirTables, prefixOut, "_full_cellModEmbed_cell_annot.csv"))

# Save the embeddings averaged over tissue, celltype etc

try(cellModEmbed_loom$close_all())
try(invisible(file.remove(paste0(dirScratch, prefixOut, "_cellModEmbed.loom"))))

#date = substr(gsub("-","",as.character(Sys.Date())),3,1000)
#if (!is.null(list_subsetModEmbed)) saveRDS(list_subsetModEmbed, file = paste0(dirRObjects, prefixOut, "list_subsetModEmbed.RDS"), compress = "gzip")
#save.image(paste0(dirScratch, prefixOut, "_", date, "_session_image.RData"))               


######################################################################
############### LOG PARAMETERS AND FILE VERSION ######################
######################################################################

as.character(Sys.time()) %>% gsub("\\ ", "_",.) %>% gsub("\\:", ".", .) ->tStop

if (is.null(path_runLog)) path_runLog <- paste0(dirLog, "_preservation_runLog.txt")

dirCurrent = paste0(LocationOfThisScript(), "/") # need to have this function defined
setwd(dirCurrent) # this should be a git directory

# get the latest git commit
gitCommitEntry <- try(system2(command="git", args=c("log", "-n 1 --oneline"), stdout=TRUE))

# Write to text file
cat(text = "\n" , file =  path_runLog, append=T, sep = "\n")
cat(text = "##########################" , file =  path_runLog, append=T, sep = "\n")

cat(text = prefixOut , file =  path_runLog, append=T, sep = "\n")
cat(sessionInfo()[[1]]$version.string, file=path_runLog, append=T, sep="\n")
if (!"try-error" %in% class(gitCommitEntry)) cat(text = paste0("git commit: ", gitCommitEntry) , file =  path_runLog, append=T, sep = "\n")
cat(text = paste0("tStart: ", tStart) , file =  path_runLog, append=T, sep = "\n")
cat(text = paste0("tStop: ", tStop) , file =  path_runLog, append=T, sep = "\n")

# output parameters (assumping use of optparse package)
cat(text = "\nPARAMETERS: " , file =  path_runLog, append=T, sep = "\n")
for (i in 1:length(opt)) {
  cat(text = paste0(names(opt)[i], "= ", opt[[i]]) , file =  path_runLog, append=T, sep = "\n")
}

cat(text = "##########################" , file =  path_runLog, append=T, sep = "\n")


############################ WRAP UP ################################

message("DONE!")
