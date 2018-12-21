# multi-analysis module preservation pipeline
# overview
# *inputs: 
#   * WGCNA outputs
#   * one or more expression datasets
# *outputs:
#   * Preservation Z-scores
#   * Preservation p-values

# usage: e.g. TODO: UPDATE
## time Rscript /projects/jonatan/tools/wgcna-src/wgcna-toolbox/wgcna_multianalysis.R --df_pathNWA /projects/jonatan/tmp-epilepsy/tables/ep_ex_4_10000_genes_cell_cluster_module_genes.csv --dirOut /projects/jonatan/tmp-epilepsy/ --prefixOut ep_in_4_multi_2 --vec_pathsDatExpr "c('/projects/jonatan/tmp-epilepsy/RObjects/EP_ex_seuratObj_filtered.RDS.gz')" --IDs_datExpr 'c("EP_ex")' --networkType "c('signed hybrid')" --corFnc cor --minCellClusterSize 50 --vec_metadataIdentCols "c'subtypesBoth')" --dirScratch /scratch/tmp-wgcna/ --dataOrganism hsapiens --scaleCenterRegress T --do_plot F --RAMGbMax 200 --geneNames hgnc|symbol|gene_name_optimal 

# TODO:

# * filter coloring so module have genes present in all test sets for a given ref? Or in any?
# * modulePreservation only uses genes in common between the reference and test set. Color vector must
# match ref set. No need to filter test set genes (although fine to do so)
# * to evaluate network connectivity preservation, we 
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

dir_current = paste0(LocationOfThisScript(), "/")
#dir_current = "/projects/jonatan/tools/wgcna-src/wgcna-toolbox/"

source(file=paste0(dir_current,"rwgcna_functions.R"))

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("optparse", "Seurat", "WGCNA", "dplyr", "Biobase", "Matrix", "parallel", "readr"))#, "NetRep"))

######################################################################
########################### OPTPARSE #################################
######################################################################

option_list <- list(
  make_option("--df_pathNWA", type="character",
              help = "path to dataframe from a gene network analysis run, in long format (i.e. one row per gene per celltype), containing 'cell_cluster', 'module', one or two gene name columns, and a column of numeric scores. I.e. one row per gene, e.g. rwgcna cell_cluster_module_genes.csv files"),  
  make_option("--geneWeightsCol", type="character",
              help = "nwa_df column with gene weights"),  
  make_option("--vec_pathsDatExpr", type="character",
              help = "Quoted vector of named strings with paths to expression datasets in seurat or matrix dataframe datatype saved as (compressed) RObject, e.g. ''c('datExpr1'='/projects/mydata/datExpr1.RDS.gz', 'datExpr2'='/projects/mydata2/datExpr2.csv'). The first dataset is assumed to contain the reference subsets, but subsets of it may also be included among test datasets''"),  
  # make_option("--IDs_datExpr", type="character",
  #             help = "Quoted vector of IDs for each of vec_pathsDatExpr"),
  make_option("--vec_pathsDatExprMetadata", type="character", default = NULL,
              help = "Optional: Quoted vector of full paths to metadata in one of the standard (compressed) character separated formats. Should correspond to files in vec_pathsDatExpr; use NA_character_ to skip a datExpr. The script will always (also) use metadata stored within the path_datExpr objects if they are of the Seurat class. [default %default]"),  
  make_option("--vec_metadataIdentCols", type="character",
              help = "vector of characters to identify columns in metadata by which to split the expression data, named by the datExpr names, e.g. ''c(datExpr1='celltype', datExpr2='clust.res.1')''."),
  # make_option("--list_metadata_corr_col", type="character", default=NULL, # currently does nothing
  #             help = "quoted list of vectors of characters to identify columns in metadata with group identities for which to compute association with preserved modules. List components are named by datasets, else assumed to correspond to vec_pathsDatExprMetadata. Use NA to skip a dataset. Takes a character with a vector in single (double) quotes of seuratObj@meta.data column names in double (single) quotes, without whitspace, e.g. 'nUMI' or 'c('Sex','Age')'. For factor or character metadata, each levels is analysed as a dummy variable, so exercise caution.  [default %default]"),
  # make_option("--nwa_cell_clusters", type="character", default=NULL,
  #             help = "Optional: Quoted vector of levels in the nwa_df$cell_cluster column for which to evaluate module preservation, e.g. ''c('neurons1','neurons2','neurons3')'' If left as NULL, uses all levels of the nwa df cell_cluster column in alphabetical order, [default %default]"),
  # make_option("--lvls_ref", type="character", default=NULL,
  #             help="a quoted list of vectors, named by expression dataset, giving reference cell identities to find in metadata[[vec_metadataIdentCols]]. The script assumes that the levels also exist in the nwa_df[['cell_cluster']] column and if not stops with an error. The list names must correspond to reference expression sets, e.g. ''list('datExpr1'=c('n01','n02','n03'))''. If left as NULL, will use all levels in the metadata, [default %default]"),
  make_option("--list_list_vec_identLvlsMatch", type="character", default=NULL,
              help="a quoted list of named lists, named by 'reference' celltype, of vectors of components designating 'test' celltype(s). All , e.g. ''list('datExpr2'=list('n01'=c('microgl1',migrogl2','migrogl3'), 's02'=c('oligo1', 'oligo2','oligo3')), 'datExpr3'=list('microglia'=c('s01','s02'), 'oligo'=c('s03','s04')))''. If the argument is left as NULL, will match all levels to all others (TODO).  [default %default]"),       #help = "quoted list of named character vectors with named components. List names are labelling levels (e.g. 'tissue', 'celltype'), vector names are column in datExpr metadata, vector values are regex to match levels. Use NA in a vector to skip a dataset for that labelling. E.g. ''list('tissue' = c(tissue='.*', NA), 'sub_celltype'=c(tissue_cell_type = '.*', ClusterName = '.*'))''"),
  make_option("--minCellClusterSize", type="integer", default=100L,
              help="What is the minimum number of cells in a cell_cluster to continue? Integer, [default %default]."),
  make_option("--minGeneClusterSize", type="integer", default=10L,
              help="What is the minimum number of genes in a module to continue? Integer, [default %default]."),
  make_option("--geneNames", type="character", default="hgnc|symbol|gene_name_optimal",
              help ="string or regex for grepping nwa_df and gene mapping dataframe and for output, e.g. 'hgnc|symbol|gene_name' or 'ensembl', [default %default]"),
  make_option("--dirOut", type="character",
              help = "Outputs go to /tables and /RObjects subdirectories"),  
  make_option("--prefixOut", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"),
  make_option("--dirScratch", type="character", default = "/scratch/tmp-wgcna/",
              help = "Outputs go to /tables and /RObjects subdirectories"),
  make_option("--networkType", type="character", default = "c('signed hybrid')",
              help = "for WGCNA modulePreservation function: one of 'signed', 'unsigned', ''c('signed hybrid')'',  [default %default]"),
  make_option("--corFnc", type="character", default = "cor",
              help = "for WGCNA modulePreservation function, 'cor' or 'bicor'"),
  make_option("--dataOrganism", type="character", default="mmusculus",
              help = "'hsapiens' or 'mmusculus', [default %default]"),
  make_option("--scaleCenterRegress", type="logical", default = "TRUE",
              help = "[default %default]. If TRUE, and vec_pathsDatExpr points to an RObject, the script will try to use the @scale.data slot and, if missing, will normalize and scale the data and regress out nUMI, percent.mito and percent.ribo. If FALSE, for RObjects will use the @data slot, or, if empty, the @raw.data slot"),  
  make_option("--RAMGbMax", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

df_pathNWA <- opt$df_pathNWA 
geneWeightsCol <- opt$geneWeightsCol
vec_pathsDatExpr <- eval(parse(text=opt$vec_pathsDatExpr))

vec_pathsDatExprMetadata <- opt$vec_pathsDatExprMetadata
if (!is.null(vec_pathsDatExprMetadata)) vec_pathsDatExprMetadata <- eval(parse(text=opt$vec_pathsDatExprMetadata)) # can be null

vec_metadataIdentCols <- eval(parse(text=opt$vec_metadataIdentCols)) # cannot be null

list_list_vec_identLvlsMatch <- opt$list_list_vec_identLvlsMatch
if (!is.null(list_list_vec_identLvlsMatch)) list_list_vec_identLvlsMatch <- eval(parse(text=list_list_vec_identLvlsMatch))

minCellClusterSize <- opt$minCellClusterSize
minGeneClusterSize <- opt$minGeneClusterSize
geneNames <- opt$geneNames
dirOut <- opt$dirOut
prefixOut <- opt$prefixOut
dirScratch <- opt$dirScratch
networkType <- opt$networkType
corFnc <- opt$corFnc
dataOrganism <- opt$dataOrganism
scaleCenterRegress <- opt$scaleCenterRegress
RAMGbMax <- opt$RAMGbMax

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F, use="pairwise.complete.obs")

######################################################################
############################ CONSTANTS ###############################
######################################################################
  
# if specified output directory doesn't exist, create it 
if (!file.exists(dirOut)) {
  dir.create(dirOut) 
  message("dirOut not found, new one created")
}

dirPlots = paste0(dirOut,"plots/")
if (!file.exists(dirPlots)) dir.create(dirPlots) 

dirTables = paste0(dirOut,"tables/")
if (!file.exists(dirTables)) dir.create(dirTables)

dirRObjects = paste0(dirOut,"RObjects/")
if (!file.exists(dirRObjects)) dir.create(dirRObjects)

dirLog = paste0(dirOut,"log/")
if (!file.exists(dirLog)) dir.create(dirLog)

flagDate = substr(gsub("-","",as.character(Sys.Date())),3,1000)

randomSeed = 12345
######################################################################
##################### VERIFY AND NAME INPUT ##########################
######################################################################

# if (!is.null(lvls_ref)) {
#   if (is.null(names(lvls_ref)) & length(IDs_datExpr)==1) { names(lvls_ref) <- IDs_datExpr }
# }

# if (all(sapply(lvls_test, function(test_lvl) is.null(names(test_lvl)))) & length(IDs_datExpr)==1) { 
#  lvls_test <- lapply(lvls_test, function(test_lvl) {
#     names(test_lvl) <- rep(IDs_datExpr, times=length(test_lvl)) 
#     return(test_lvl)
#   })
# }

#if (is.null(names(lvls_test))) names(lvls_test) <- lvls_ref[[1]]

#if (is.null(names(vec_metadataIdentCols)) & length(IDs_datExpr)==1 & length(vec_metadataIdentCols)==1) {names(vec_metadataIdentCols) <- IDs_datExpr} 
# TODO
#if(!is.null(paths_metadata)) stopifnot(length(paths_metadata)==length(vec_pathsDatExpr))
# stopifnot(all.equal(sapply(ident_groups, length)))
# stopifnot(length(ident_groups_combs) %in% c(1, length(ident_groups)))
#names(vec_metadataIdentCols) <- IDs_datExpr
######################################################################
######################## LOAD MODULE DATA ############################
######################################################################

message("Loading gene module data")

df_geneModule <- load_obj(df_pathNWA)

# If there are duplicate modules between WGCNA runs, prefix module names with the WGCNA run 
# if (sapply(X=list_df_geneModule , function(df) unique(df[["module"]]), simplify = T) %>% unlist(use.names=F) %>% duplicated %>% any) {
#   for (i in 1:length(list_df_geneModule)) {
#     list_df_geneModule[[i]][["module"]] <- paste0(prefixes_nwa_run[i], "_",list_df_geneModule[[i]][["module"]]) # ensure modules from different runs remain distinct
#   }
# }

# If there are duplicate levels for celltype or tissue between WGCNA runs, prefix level names with the WGCNA run 

# for (ident_group in ident_groups_nwa) {
#   for (colname in names(ident_group)) {
#     if (sapply(X=list_df_geneModule , function(df) unique(df[[colname]]), simplify = T) %>% unlist(use.names = F) %>% duplicated %>% any) {
#       for (i in 1:length(list_df_geneModule)) {
#         list_df_geneModule[[i]][[colname]] <- paste0(prefixes_nwa_run[i], "_",list_df_geneModule[[i]][[colname]]) # ensure modules from different runs remain distinct
#       }
#     }
#   }
# }

# Identify gene and score columns
colnameGene <- grep(geneNames, colnames(df_geneModule), value=T)
colnameLoading <- grep(geneWeightsCol, colnames(df_geneModule), value=T)#colnames(df_geneModule)[which(sapply(X=colnames(df_geneModule), FUN =function(colname) class(df_geneModule[[colname]]))=="numeric")]

######################################################################
####################### LOAD DATEXPR METADATA ########################
######################################################################

# check if the specified reference levels exist in the gene module df

if (!all(names(list_list_vec_identLvlsMatch) %in% df_geneModule[["module"]]))  {
  missing = names(list_list_vec_identLvlsMatch)[!names(list_list_vec_identLvlsMatch) %in% df_geneModule[["cell_cluster"]]] 
  if (length(missing)>0) stop(paste0(missing, " not found in gene module df ", collapse=" "))
}

######################################################################
####################### LOAD DATEXPR METADATA ########################
######################################################################

list_datExprMetadata <- NULL
if (!is.null(vec_pathsDatExprMetadata)){ 
  message("Loading metadata") 
  list_datExprMetadata <- lapply(vec_pathsDatExprMetadata, load_obj)
} 

######################################################################
######################## LOAD EXPRESSION DATA ########################
######################################################################

message("Loading expression matrix")

# load/open connection to expression matrix
# NB: may be huge
list_seuratObj <- lapply(vec_pathsDatExpr, function(vec_pathsDatExpr) {
  if (grepl(pattern = "\\.loom", vec_pathsDatExpr)) {
    dataObj <-  connect(filename = vec_pathsDatExpr, mode = "r+")
  } else {
    dataObj <- load_obj(vec_pathsDatExpr)
  } 
  if (!any(c("seurat","loom") %in% class(dataObj))) { # not a loom object, nor a seurat object
    dataObj <- CreateSeuratObject(raw.data=dataObj, project= prefixOut, min.cells = -Inf, min.genes = -Inf)
  } else if ("loom" %in% class(dataObj)) {
    dataObj <- Seurat::Convert(from=dataObj, to="seurat")
  } # if already a seurat object, do nothing
  return(dataObj)
})

names(list_seuratObj) <- names(vec_pathsDatExpr)
  
# Add any further metadata coming from another file
if (!is.null(list_datExprMetadata)) {
  list_seuratObj<- mapply(function(seuratObj, datExpr_metadata) {
    if (all_equal(rownames(datExpr_metadata), colnames(seuratObj@raw.data))) {
      seuratObj <- AddMetaData(object = seuratObj, metadata = data.frame(datExpr_metadata))
    }
  }, seuratObj=list_seuratObj, datExpr_metadata=list_datExprMetadata, SIMPLIFY=F)
}

######################################################################
###################### EXTRACT METADATA ##############################
######################################################################

list_metadata <- lapply(list_seuratObj, function(seuratObj) seuratObj@meta.data)
names(list_metadata) <- names(vec_pathsDatExpr)

######################################################################
######### VERIFY THAT lvls_ref AND lvls_test EXIST IN METADATA #######
######################################################################

# check if the specified reference levels exist in the metadata annotation

if (!all(names(list_list_vec_identLvlsMatch) %in% list_metadata[[1]][[vec_metadataIdentCols[1]]]))  {
  missing = names(list_list_vec_identLvlsMatch)[!names(list_list_vec_identLvlsMatch) %in% list_metadata[[1]][[vec_metadataIdentCols[1]]]] 
  stop(paste0(missing, " not found in ", names(list_metadata)[1], " metadata", collapse=" "))
}

# check if the specified test levels exist in the metadata annotation

for (lvlTest in names(list_list_vec_identLvlsMatch)) {
  for (datExprTestName in names(list_list_vec_identLvlsMatch[[lvlTest]])) {
    if(!all(list_list_vec_identLvlsMatch[[lvlTest]][[datExprTestName]] %in% list_metadata[[datExprTestName]][[vec_metadataIdentCols[datExprTestName]]])) {
      missing = list_list_vec_identLvlsMatch[[lvlTest]][[datExprTestName]][!list_list_vec_identLvlsMatch[[lvlTest]][[datExprTestName]] %in% list_metadata[[datExprTestName]][vec_metadataIdentCols[datExprTestName]]]
      stop(paste0(missing, " not found in ", datExprTestName, " metadata", collapse=" "))
    }
  }
}


######################################################################
######################## MAP GENES TO SAME NAMES #####################
######################################################################

# TODO: make orthologue mapping function

if (FALSE) {  
  mapping = if (dataOrganism=="hsapiens") {
    load_obj("/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz")
  } else if (dataOrganism=="mmusculus") {
    load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
  }  
  
  args <- list("X"=list_seuratObj)
  fun <- function(seuratObj) {
    if (any(grepl("ENSG|ENSMUSG",rownames(seuratObj@raw.data)))) {
      current_ensembl <- T
      tmp <- gene_map(df = seuratObj@raw.data,
                      idx_colnameGeneumn = NULL,
                      mapping=mapping,
                      from="ensembl",
                      to="gene_name_optimal",
                      replace = T,
                      na.rm = T)
    } else {
      current_ensembl <- F
      tmp <- seuratObj@raw.data
    }
    message("Counting percent ribo and percent mito")
    idx_mito.genes <- grepl(pattern = "^mt-", x = rownames(tmp), ignore.case=T)
    idx_ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(tmp), ignore.case=T)
    percent_mito <- Matrix::colSums(tmp[idx_mito.genes, ])/Matrix::colSums(tmp)
    percent_ribo <- Matrix::colSums(tmp[idx_ribo.genes, ])/Matrix::colSums(tmp)
    seuratObj <- AddMetaData(object = seuratObj, metadata = percent_mito, col.name = "percent_mito")
    seuratObj <- AddMetaData(object = seuratObj, metadata = percent_ribo, col.name = "percent_ribo")
    if (grepl("ensembl", geneNames) & !current_ensembl) {
      seuratObj@raw.data <- gene_map(df = seuratObj@raw.data,
                                      idx_colnameGeneumn = NULL,
                                      mapping=mapping,
                                      from="gene_name_optimal",
                                      to="ensembl",
                                      replace = T,
                                      na.rm = T)
    } else if (grepl("symbol|hgcn", geneNames) & current_ensembl) {
      seuratObj@raw.data <- tmp
    }
    return(seuratObj)
  }
  # Call the function
  list_seuratObj <- safeParallel(fun=fun, args=args, mapping=mapping)
  rm(mapping)
}
######################################################################
######################## NORMALIZE EXPRESSION DATA ###################
######################################################################
  
message("Normalizing the expression data")

args <- list("X"=list_seuratObj)
fun <- function(seuratObj) NormalizeData(seuratObj)
list_seuratObj <- safeParallel(fun=fun, args=args)

######################################################################
################# SCALE AND REGRESS EXPRESSION DATA ##################
######################################################################

if (scaleCenterRegress) {
  message("scaling and regressing out confounders from the data")
  
  list_seuratObj<-lapply(X = list_seuratObj, FUN = function(seuratObj) ScaleData(seuratObj,
                        vars.to.regress = c("percent_mito", "nUMI", "percent_ribo")[c("percent_mito", "nUMI", "percent_ribo") %in% colnames(seuratObj@meta.data)], 
                        do.scale=scaleCenterRegress,
                        do.center=scaleCenterRegress,
                        display.progress= T, 
                        do.par=scaleCenterRegress, 
                        num.cores=10))
}

######################################################################
############## EXTRACT SCALED DATA AND DELETE SEURAT OBJ #############
######################################################################
  
list_datExpr <- lapply(list_seuratObj, function(seuratObj) {
  out <- if (scaleCenterRegress) seuratObj@scale.data else seuratObj@data
  return(t(out))
  })
names(list_datExpr) <- names(vec_pathsDatExpr)
rm(list_seuratObj)

######################################################################
######## EXTRACT METADATA AND CONVERT FACTORS TO MODEL MATRIX ########
######################################################################
#TODO: Reintroduce correlation analysis?
if (FALSE) {
  list_metadat_names <- lapply(list_metadata, colnames)
  # Convert any character or factor meta.data to numeric dummy variables each level with its own numeric column
  if (!is.null(list_metadata_corr_col)) {
    list_metadata <- mapply(function(metadata_corr_col, metadat_names) {
      if (any(colnames(meta.data) %in% metadata_corr_col)) {
        metadata <- matrix(NA, nrow=nrow(seuratObj@meta.data), ncol=1)
        include <- seuratObj@meta.data[,colnames(seuratObj@meta.data) %in% metadata_corr_col, drop=F]
        for (i in 1:ncol(include)) {
          if (class(include[,i]) %in% c("factor", "character")) {
            metadata <- cbind(metadata, factorToIndicator(include[,i, drop=T]))        
          } else {
            metadata <- cbind(metadata, include[,i, drop=F])
          }
        }
        metadata <- metadata[,-1, drop=F]
        if (!is.null(metadata_corr_filter_vals)) metadata <- metadata[, toupper(colnames(metadata)) %in% toupper(metadata_corr_filter_vals), drop = F]
        metadata <- as.data.frame(metadata)
        # Filter out any metadata columns where all the values are identical
        metadata <- metadata[apply(metadata, MARGIN=2, FUN = function(x) length(unique(x))>1)]
        if (ncol(metadata) == 0) metadata <- NULL    
        rownames(metadata) = rownames(seuratObj@meta.data)
        return(metadata)
      } else {
        NULL
      }
    },metadata_corr_col = list_metadata_corr_col, 
    metadat_names = list_metadat_names, 
    SIMPLIFY=F)
  } else {
    list_metadata <- NULL
  }
}

######################################################################
################## RUN PRESERVATION ANALYSIS #########################
######################################################################

outfile = paste0(dirLog, prefixOut, "_modulePreservation_log.txt")

lvlRef_datExprTest_presModsOut <- list()
lvlRef_datExprTest_presNwOut <- list()

#WGCNA::disableWGCNAThreads()
setwd(dir = dirScratch)

for (lvlRef in names(list_list_vec_identLvlsMatch)) {
  
  lvlRef_datExprTest_presModsOut[[lvlRef]] <- list()
  lvlRef_datExprTest_presNwOut[[lvlRef]] <- list()
  
  ####################################################
  #### Get gene network coloring for cell_cluster ####
  ####################################################
  
  idxDuplicateGenes <- duplicated(df_geneModule[[colnameGene]][df_geneModule[["cell_cluster"]] == lvlRef])
  coloring <- df_geneModule[["module"]][df_geneModule[["cell_cluster"]] == lvlRef][!idxDuplicateGenes]
  names(coloring) <- df_geneModule[[colnameGene]][df_geneModule[["cell_cluster"]] == lvlRef][!idxDuplicateGenes]
  
  ####################################################
  ############ Prepare reference dataset #############
  ####################################################
  
  # Get metadata df
  metadataRef <- list_metadata[[1]]
    metadataColnameRef <- vec_metadataIdentCols[1]
  
  idxLvlRefCells <- as.integer(na.omit(match(rownames(metadataRef)[(metadataRef[[metadataColnameRef]] == lvlRef)], 
                                              rownames(list_datExpr[[datExprTestName]]))))
  
  if (length(idxLvlRefCells) < minCellClusterSize) {
    lvlRef_datExprTest_presModsOut[[lvlRef]] <- NA_character_
    lvlRef_datExprTest_presNwOut[[lvlRef]] <- NA_character_
    message(paste0(lvlRef, " skipped because it has < ", minCellClusterSize, " cells"))
    next
  } 
  
  idxGenes <- as.integer(na.omit(match(names(coloring), colnames(list_datExpr[[1]]))))
  datExprRefLvl <- list_datExpr[[1]][idxLvlRefCells, idxGenes]
  
  ######################################################################
  #### Filter out coloring genes which are missing in datExpr_ref ######
  ######################################################################
  
  coloring <- coloring[names(coloring) %in% colnames(datExprRefLvl)]
  
  ####################################################
  ############ Filter out small modules ##############
  ####################################################
  
  modsTooSmall <- names(table(coloring))[table(coloring)<minGeneClusterSize]
  coloring <- coloring[!coloring %in% modsTooSmall]
  
  ####################################################
  ############## Prepare test datasets ###############
  ####################################################
  
  for (datExprTestName in names(list_list_vec_identLvlsMatch[[lvlRef]])) {

    #lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]] <- list()
    #lvlRef_datExprTest_presNwOut[[lvlRef]][[datExprTestName]] <- list()

    message(paste0("Reference: ", lvlRef, ", testing preservation in subsets of ", datExprTestName))
    
    metadataTest <- list_metadata[[datExprTestName]]
    metadataColnameTest <- vec_metadataIdentCols[datExprTestName]
    
    idxGenes <- as.integer(na.omit(match(names(coloring), colnames(list_datExpr[[datExprTestName]]))))
      
    # filter out modules missing completely in datExprTest
    mods_absent <- names(table(coloring))[!names(table(coloring)) %in% names(table(coloring[names(coloring) %in% colnames(list_datExpr[[datExprTestName]])[idxGenes]]))]
    coloring_f <- coloring
    coloring_f[coloring %in% mods_absent] <- "grey"
    
    # Filter out modules where fewer than half the genes present in datExprTest
    vec_modPropGenesInTest <- table(coloring_f)/table(coloring[coloring %in% names(table(coloring_f))])
    vec_modPropGenesInTestTooSmall<- names(vec_modPropGenesInTest[vec_modPropGenesInTest<0.5])
    coloring_f <- coloring_f[!coloring_f %in% vec_modPropGenesInTestTooSmall]

    list_datExprTest <- list()
    
    for (lvlTest in list_list_vec_identLvlsMatch[[lvlRef]][[datExprTestName]]) {
      
      idxCells <- as.integer(na.omit(match(rownames(metadataTest)[(metadataTest[[metadataColnameTest]] == lvlTest)], 
                         rownames(list_datExpr[[datExprTestName]]))))
      datExprLvlTest <- list_datExpr[[datExprTestName]][idxCells, idxGenes]
      
      ######################################################################
      ################ Check that test subset has sufficient cells #########
      ######################################################################
      
      if (nrow(datExprLvlTest) < minCellClusterSize) {
        message(paste0(lvlRef, " preservation in ", lvlTest, " skipped because ", lvlTest, " has < ", minCellClusterSize, " cells"))
        next
      }
      list_datExprTest[[lvlTest]] <- datExprTest
    }
    
    if (length(list_datExprTest)==0) {
      lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]] <- NA_character_
      lvlRef_datExprTest_presNwOut[[lvlRef]][[datExprTestName]] <- NA_character_
      message(paste0(datExprTestName, " skipped because no test sets have < ", minCellClusterSize, " cells"))
      next
    } 
    # list_idx_subset_test <- mapply(function(metadataTest, colname_subset_test, pattern_subset_test) {
    #   which(metadataTest[[colname_subset_test]] == unique(grep(paste0("^", pattern_subset_test, "$"), metadataTest[[colname_subset_test]], value=T)))
    # }, metadataTest = list_metadataTest, colname_subset_test=colnames_subset_test, pattern_subset_test=patterns_subset_test, SIMPLIFY=F)
    # 
    # list_datExprTest <- list_datExpr[IDs_datExprTest]
    # 
    # list_datExprTest<- mapply(function(datExprTest, idx_subset_test) {
    #   idxGenes <- colnames(datExprTest) %in% names(coloring)
    #   datExprTest <- datExprTest[idx_subset_test,idxGenes]
    # }, datExprTest = list_datExprTest, idx_subset_test = list_idx_subset_test, SIMPLIFY=F)
    # 
    # names(list_datExprTest) <- patterns_subset_test
    # # Check that there are sufficient cells
    # list_datExprTest <- lapply(list_datExprTest,function(datExprTest) {
    #   return( if (nrow(datExprTest) >= minCellClusterSize) datExprTest else NULL)
    # })
    # 
    # check test sets
    # idx_test_NULL <- sapply(X=list_datExprTest, FUN=is.null,simplify=T)
    # if (sum(idx_test_NULL)==length(list_datExprTest)) {
    #   list_modulePreservation_out[[cell_cluster]] <- NA_character_
    #   next
    # } else if (sum(idx_test_NULL)>1 & sum(idx_test_NULL)<length(list_datExprTest)) {
    #   list_datExprTest <- list_datExprTest[!idx_test_NULL]
    #   message(paste0(patterns_subset_test[idx_test_NULL], " dropped as test sets due to having fewer than ", minCellClusterSize, " cells"))
    # }
  
    ################################################################################
    ################## CHECK MODULE GENES IN TEST SET ##############################
    ################################################################################
    
    # Make the coloring vector the same length as the number of genes in the ref set
    # coloring <- coloring[names(coloring) %in% colnames(datExpr_ref)]
    # colors_tmp <- rep("grey", length=ncol(datExpr_ref))
    # names(colors_tmp) <- colnames(datExpr_ref)
    # colors_tmp[as.integer(na.omit(match(names(coloring), names(colors_tmp))))] <- coloring[as.integer(na.omit(match(names(colors_tmp), names(coloring))))]
    # coloring <- colors_tmp
    
    ####################################################
    ################ Prepare multidata #################
    ####################################################
  
    # Prepare multidata
    multiData <- vector(mode="list",length = 1+length(list_datExprTest)) # 1 for reference network
    multiData[[1]] <- list("data"=datExprRefLvl)
    for (k in 1:(length(multiData)-1)){
      multiData[[k+1]] <- list("data"=list_datExprTest[[k]])
    }
    names(multiData)[1] <- lvlRef
    names(multiData)[2:length(multiData)] <- names(list_datExprTest)
    
    rm(list_datExprTest)
    # if (FALSE) {   
    #   # Check the multidata samples 
    #   goodSamplesGenesMS_out <- try(goodSamplesGenesMS(
    #                               multiExpr=multiData,
    #                               multiWeights = NULL,
    #                               minFraction = 1/2000,
    #                               minNSamples = 1,
    #                               minNGenes = 50,
    #                               tol = NULL,
    #                               minRelativeWeight = -Inf,
    #                               verbose = 3,
    #                               indent = 0))
    #   if ("try-error" %in% class(goodSamplesGenesMS_out)) {
    #     message(paste0(cell_cluster, ": goodSamplesGenesMS failed, skipping to next ref set"))
    #     list_modulePreservation_out[[cell_cluster]] <- NA_character_
    #     next
    #   }
    # 
    #   multiData <- mapply(function(dataset, idx_goodSamples) {
    #     dataset$data <- dataset$data[idx_goodSamples, goodSamplesGenesMS_out$goodGenes]
    #     }, dataset=multiData, idx_goodSamples=goodSamplesGenesMS_out$goodSamples)
    # }
    
    ####################################################
    ####### Run preservationNetworkConnectivity ########
    ####################################################
    
    # This should test preservation in each test_lvl datExpr
    lvlRef_datExprTest_presNwOut[[lvlRef]][[datExprTestName]] <- 
      WGCNA::preservationNetworkConnectivity(multiExpr = multiData, 
                                             corFnc = corFnc, 
                                             corOptions = NULL, 
                                             networkType = networkType, 
                                             power=8, 
                                             sampleLinks=T, 
                                             nLinks=50000, 
                                             blockSize=5000, 
                                             setSeed=randomSeed, 
                                             weightPower=2, 
                                             verbose=3, 
                                             indent=0)
    
    ####################################################
    ############### Prepare multicolor #################
    ####################################################
    
    multiColor <- list(coloring)
    names(multiColor) <- lvlRef
      
    # Prepare for parallel computation (WGCNA multi threads)
    
    require("doParallel")
    additionalGb = max(as.numeric(sapply(multiData, FUN = function(x) object.size(x), simplify = T)))/1024^3
    objSizeGb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
    nCores <- max(1, min(detectCores() %/% 3, RAMGbMax %/% (objSizeGb + additionalGb))-1)
    nCores <- min(nCores, 40)
    enableWGCNAThreads(nThreads = nCores)
    
      lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]] <- tryCatch({
      modulePreservation(multiData=multiData,
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
                        randomSeed = randomSeed, 
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
                        indent = 0)}, error = function(err) {
                          message(paste0(lvlRef, ": modulePreservation failed with the following error: ", err))
                          write(x = paste0(lvlRef, ": modulePreservation failed with the following error: ", err), file = outfile,
                                ncolumns = 1, append = T)
                          return(NA_character_)
                          })
    
  }
}


######################################################################
################# REFORMAT MODULE PRESERVATION SCORES ################
######################################################################

lvlRef_df_datExprTestStats <- lapply(X=names(list_list_vec_identLvlsMatch), FUN=function(lvlRef) {
  
  if (all(is.na(lvlRef_datExprTest_presModsOut[[lvlRef]]))) return(NA_character_)
  
  list_datExprTestStats <- lapply(names(list_list_vec_identLvlsMatch[[lvlRef]]), function(datExprTestName) {
    
    if (all(is.na(lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]])) |
        all(is.null(lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]]))) return(NA_character_)

    modulePreservationOut = lvlRef_datExprTest_presModsOut[[lvlRef]][[datExprTestName]]
    # we should get a matrix with a column for each test set
    list_df_Zsummary.pres <- lapply("X"=names(modulePreservationOut[["preservation"]][["Z"]][[1]][-1]), 
                      FUN=function(lvlTest) {
                        df_presStats <- modulePreservationOut[["preservation"]][["Z"]][[1]][-1][[lvlTest]]                 
                        if(!all(is.na(df_presStats))) {
                          lvlTestCrop <- gsub("inColumnsAlsoPresentIn\\.", "", lvlTest)
                          df_Zsummary.pres <- data.frame(
                            "ref_dataset"=names(vec_pathsDatExpr)[1],
                            "ref_lvl"=lvlRef,
                            "test_dataset"=datExprTestName,
                            "test_lvl"=rep(lvlTestCrop, times= nrow(df_presStats)), 
                            "module" = rownames(df_presStats), 
                            "Z_stat"=df_presStats[["Zsummary.pres"]])
                 
                          return(df_Zsummary.pres)
                          } else {
                            return(NA_character_)
                          }
    })

    list_df_Zsummary.pres <- list_df_Zsummary.pres[!sapply(list_df_Zsummary.pres, function(df) all(is.na(df)))]
    df_Zsummary.pres <- if (length(list_df_Zsummary.pres)==1) list_df_Zsummary.pres[[1]] else Reduce(f=rbind, x=list_df_Zsummary.pres)
    
    #Z_stats_long <- reshape2::melt(data=Z_stats)
    #colnames(Z_stats_long) <- c("module","test_lvl","Zsummary.pres")
    # we should get a matrix with a column for each test set
    list_df_pBonf.pres <-lapply(names(modulePreservationOut[["preservation"]][["log.pBonf"]][[1]][-1]), FUN=function(lvlTest) {
      df_presStats <- modulePreservationOut[["preservation"]][["log.pBonf"]][[1]][-1][[lvlTest]] 
      if(!all(is.na(df_presStats))) {
        lvlTestCrop <- gsub("inColumnsAlsoPresentIn\\.", "", lvlTest)
        df_pBonf.pres <- data.frame(
          "module" = rownames(df_presStats), 
          "log.p.Bonfsummary.pres"=df_presStats[["log.p.Bonfsummary.pres"]])
        return(df_pBonf.pres)
      } else {
        return(NA_character_)
      }
    })
  
    list_df_pBonf.pres <- list_df_pBonf.pres[!sapply(list_df_pBonf.pres, function(df) all(is.na(df)))]
    df_pBonf.pres <- if (length(list_df_pBonf.pres)==1) list_df_pBonf.pres[[1]] else Reduce(f=rbind, x=list_df_pBonf.pres)
    
    # colnames(logpBonf_stats) <- gsub("inColumnsAlsoPresentIn\\.","",colnames(logpBonf_stats))
    # logpBonf_stats_long <- reshape2::melt(logpBonf_stats)[-c(1,2)]
    # colnames(logpBonf_stats_long) <- c("log.p.Bonfsummary.pres")

    df_presStats <- dplyr::full_join(df_Zsummary.pres, df_pBonf.pres)
      
      
      # data.frame("ref_dataset" = names(vec_pathsDatExpr)[1],
      #                        "ref_lvl"= rep(lvlRef, times=nrow(Z_stats_long)),
      #                        "test_dataset"= rep(datExprTestName, times=nrow(Z_stats_long)),
      #                        Z_stats_long,
      #                        logpBonf_stats_long)
    return(df_presStats)
  })
  df_datExprTestStats <- if (length(list_datExprTestStats)==1) list_datExprTestStats[[1]] else Reduce(f=rbind, x=list_datExprTestStats)
  return(df_datExprTestStats)
})

names(lvlRef_df_datExprTestStats) <-names(list_list_vec_identLvlsMatch)
lvlRef_df_datExprTestStats <- lvlRef_df_datExprTestStats[!sapply(lvlRef_df_datExprTestStats, function(df) all(is.na(df)))]
df_modPresStats <- if (length(lvlRef_df_datExprTestStats)==1) lvlRef_df_datExprTestStats[[1]] else Reduce(f=rbind, x=lvlRef_df_datExprTestStats)
df_modPresStats <- df_modPresStats[!df_modPresStats[["module"]]=="gold",]

######################################################################
########### GET NETWORK CONNECTIVITY PRESERVATION STATS ##############
######################################################################
df_nwPresStats <- data.frame(ref_dataset=NA_character_, 
                         ref_lvl=NA_character_, 
                         test_dataset=NA_character_, 
                         test_lvl =NA_character_)


k=0
for (lvlRef in names(lvlRef_datExprTest_presNwOut)) { 
  datExprTest_presNwOut <- lvlRef_datExprTest_presNwOut[[lvlRef]]
  if (all(is.na(datExprTest_presNwOut))) next
    for (datExprTestName in names(datExprTest_presNwOut)) {
      presNwOut <- datExprTest_presNwOut[[datExprTestName]]
      if (all(is.na(presNwOut))) next
      for (i in 1:ncol(presNwOut[["pairwise"]])) {
        k=k+1
        df_nwPresStats[k,"ref_dataset"] <- names(vec_pathsDatExpr)[1]
        df_nwPresStats[k,"ref_lvl"] <- lvlRef
        df_nwPresStats[k,"test_dataset"] <- datExprTestName
        df_nwPresStats[k,"test_lvl"] <- list_list_vec_identLvlsMatch[[lvlRef]][[datExprTestName]][i]
        df_nwPresStats[k,"network_pairwise_min"] <- min(presNwOut[["pairwise"]][,i])
        df_nwPresStats[k,"network_pairwise_median"] <- median(presNwOut[["pairwise"]][,i])
        df_nwPresStats[k,"network_pairwise_max"] <- max(presNwOut[["pairwise"]][,i])
        df_nwPresStats[k,"network_pairwiseWeighted_min"] <- min(presNwOut[["pairwiseWeighted"]][,i])
        df_nwPresStats[k,"network_pairwiseWeighted_median"] <- median(presNwOut[["pairwiseWeighted"]][,i])
        df_nwPresStats[k,"network_pairwiseWeighted_max"] <- max(presNwOut[["pairwiseWeighted"]][,i])
    }
    # will come out as character vector
  }
}

######################################################################
########################### JOIN RESULTS #############################
######################################################################

df_presStats <- dplyr::left_join(df_modPresStats, df_nwPresStats, by=c("ref_dataset","ref_lvl", "test_dataset", "test_lvl"))

######################################################################
############################ OUTPUTS #################################
######################################################################

write.csv(x=df_presStats, file = paste0(dirTables, prefixOut, "_df_preservationStats.csv"), quote = F, row.names = F)
# write.csv(x=df_nwPresStats, file = paste0(dirTables, prefixOut, "_df_nwPresStats.csv"), quote = F, row.names = F)
#   
########################### SAVE LOG OF PARAMS #######################

paramsRun <- unlist(opt, recursive=F)# all the user options/defaults
matrix(unlist(paramsRun), nrow=1) %>% as.data.frame(row.names=NULL, stringsAsFactors=F) -> df_paramsRun
colnames(df_paramsRun) <- names(paramsRun)
paramsRunPath = sprintf("%s%s_paramsRun.tab", dirLog, prefixOut)
append = F#file.exists(params_run_path) # NB: if left as NULL, they won't get included!
write.table(df_paramsRun, 
            file=paramsRunPath, 
            quote = F, 
            sep = "\t", 
            row.names=F, 
            append = append, 
            col.names = !append)

############################ FILTERED ################################

message("Script done!")
