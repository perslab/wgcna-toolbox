# Title: Script to find robust WGCNA modules

######################################################################
############################## USAGE #################################
#######################o###############################################

# e.g.
# use R OpenBLAS for faster linear algebra. The following command achieves this on yggdrasil (Pers lab).
# see https://stackoverflow.com/questions/26897335/how-can-i-load-a-specific-version-of-r-in-linux and https://gist.github.com/cheuerde/8fb9fd0dc8c0eca17c16#file-r_openblas-L64
# export PATH="/usr/local/R-3.5.1/bin/:$PATH"  
# export R_MAX_NUM_DLLS=999
# 
######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--pathDatExpr", type="character",
              help = "Provide full path to raw counts expression data in delimited (with gene names in the first row), seurat or seurat loom object format, saved as .RData, .RDS or .loom file, optionally gzip compressed"),
  make_option("--pathMetadata", type="character", default = NULL,
              help = "If datExpr is not provided as a Seurat object, provide path to metadata containing celltype identities under the colIdents. File should be in one of the standard (compressed) character separated formats. First column should have cellnames matching column names in datExpr [default %default]"),  
  make_option("--dirProject", type="character", default=NULL,
              help = "Optional. Provide project directory. Must have subdirs RObjects, plots, tables. If not provided, assumed to be dir one level up from input data dir. [default %default]"),
  make_option("--dirTmp", type="character", default=NULL,
              help = "Directory for temporary files, defaults to creating a 'tmp' folder within the project directory, [default %default]"),
  make_option("--prefixData", type="character", 
              help = "Dataset prefix for the project, [default %default]"),
  make_option("--prefixRun", type="character", default="run1",
              help = "Run prefix to distinguish runs e.g. with different parameters, [default %default]"),
  make_option("--dataType", type="character", default="sc",
              help = "Expression data type, one of 'sc' or 'bulk', [default %default]"),
  make_option("--autosave", type="logical", default=T,
              help = "Autosave session images at regular intervals to resume later? [default %default]."),
  make_option("--resume", type="character", default=NULL,
              help = "Resume from a checkpoint? Must have same path and prefixData as provided. Options are 'checkpoint_1' - '5' [default %default]"),
  make_option("--colIdents", type="character", default=NULL,
              help = "Specify a metadata column to use for subsetting the data by cell cluster. If NULL uses the @ident slot. [default %default]"),
  make_option("--vec_colAnnot", type="character", default=NULL,
              help="Provide column names for higher levels of sample hierarchy above colIdents, e.g. tissue, to add to output. Cannot cross-cut colIdents. Does not affect the analysis. [default %default]"),
  make_option("--assayUse", type="character", default="RNA",
              help="If datExpr is a Seurat object, indicate whether to use 'RNA', 'SCT' or 'integrated' assay, [default %default]"),
  make_option("--slotUse", type="character", default="scale.data",
              help="Use 'counts', logNormalized ('data') or Z-scored ('scale.data') for the analysis? The script carries out the scaling. scale.data is unavoidable if we wish to regress out confounders [default %default]"),
  make_option("--varsToRegress", type="character", default='c("nCount_RNA", "percent.mito", "percent.ribo")',
              help="Provide arguments to Seurat's ScaleData function in the form of a vector in quotes, defaults to c('nUMI', 'percent.mito', 'percent.ribo') [default %default]"),
  make_option("--minGeneCells", type="integer", default=20L,
              help="What is the minimum number of cells in each subset in the data in which a gene should be detected to not be filtered out? Integer, [default %default]."),
  make_option("--minCellClusterSize", type="integer", default=50L,
              help="What is the minimum number of cells in a subset to continue? Integer, [default %default]."),
  make_option("--featuresUse", type="character", default="PCLoading",
              help="One of 'var.features', 'PCLoading', 'JackStrawPCLoading' or 'JackStrawPCSignif'. PCLoading uses the top nFeatures PC loading genes; JackStrawPCLoading does the same but only on PCs retained after JackStraw significance testing; JackStrawPCSignif uses gene p-values instead of loadings [default %default]"), 
  make_option("--nFeatures", type="numeric", default=5000,
              help="How many genes to select under the criteria set with the featuresUse argument, set to Inf to take all [default %default]"), 
  make_option("--nPC", type="integer", default=100L,
              help = "Number of principal components to compute for selecting genes based on their loadings, [default %default]."),
  make_option("--nRepJackStraw", type="integer", default=250L,
              help = "Number of times to re-run PCA after permuting a small proportion of genes to perform empirical significance tests, i.e. the `JackStraw` procedure (see `featuresUse` above), [default %default]."),
  make_option("--corFnc", type="character", default="cor",
              help="Use 'cor' for Pearson or 'bicor' for midweighted bicorrelation function (https://en.wikipedia.org/wiki/Biweight_midcorrelation). [default %default]"), 
  make_option("--networkType", type="character", default = "signed hybrid",
              help="'signed' scales correlations to [0:1]; 'unsigned' takes the absolute value (but the TOM can still be 'signed'); ''c('signed hybrid')'' (quoted vector) sets negative correlations to zero. [default %default]"),
  make_option("--nRepTOM", type="integer", default=100L,
              help = "Number of times to resample the dataset when finding the consensus TOM [default %default]"),
  make_option("--hclustMethod", type="character", default="average",
              help = "Hierarchical clustering agglomeration method. One of 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA), 'median' (= WPGMC) or 'centroid' (= UPGMC). See hclust() documentation for further information. [default %default]"),
  make_option("--minClusterSize", type="character", default="15L",
              help = "Minimum features needed to form a module, or an initial cluster before the Partitioning Around Medoids-like step. WGCNA authors recommend decreasing the minimum cluster size when using higher settings of deepSplit. Takes a character with a vector, without whitespace, of integer values to try 'c(15,20,25)' [default %default]"),
  make_option("--deepSplit", type="character", default="2L",
              help = "Controls the sensitivity of the cutreeDynamic/cutreeHybrid algorithm. Takes a character with a vector, without whitespace, of integer values between 0-4 to try, e.g. 'c(1,2,3)', [default %default]"),
  make_option("--moduleMergeCutHeight", type="character", default="c(0.15)",
              help = "Cut-off level for 1-cor(eigen-gene, eigen-gene) for merging modules. Takes a character with a vector, without whitespace, of double values to try, e.g. 'c(0.1, 0.2)', [default %default]"),
  make_option("--pamStage", type="character", default="c(FALSE)",
              help = "cutreeHybrid. Perform additional Partition Around Medroids step? Users favoring specificity over sensitivity (that is, tight clusters with few mis-classifications) should set pamStage = FALSE, or at least specify the maximum allowable object-cluster dissimilarity to be zero (in which case objects are assigned to clusters only if the object-cluster dissimilarity is smaller than the radius of the cluster). Takes a character with a vector, without whitespace, of logicals to try, e.g. 'c(TRUE,FALSE)', [default %default]"),
  make_option("--kMReassign", type="logical", default="TRUE",
              help = "Following hierarchical clustering, do additional k-means clustering to reassign features whose proximity to another module is >= 1.2 times the proximity to their own? [default %default]"),
  make_option("--kMSignifFilter", type="logical", default="TRUE",
              help = "Do a t test to filter out insignificant features from modules? If fuzzyModMembership=='kME', carries out a correlation t-test; if kIM, does a two-sample t-test between IntraModular and ExtraModular connectivities"),
  make_option("--fuzzyModMembership", type="character", default="kIM",
              help="Which 'fuzzy' measure of gene membership in a module should be used? Options are 'kME' (correlation between gene expression and module PC1 expression) and 'kIM' (sum of edges between a gene and genes in the module, normalized by number of genes in the module; i.e. average distance in the TOM between a gene and genes in the module [default %default]."),  
  make_option("--RAMGbMax", type="integer", default=400,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)

######################################################################
########################## GET CURRENT DIR ###########################
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

######################################################################
############################## FUNCTIONS #############################
######################################################################

source(file = paste0(current.dir, "rwgcna_functions_seurat3.0.R"))

######################################################################
############################## PACKAGES ##############################
######################################################################

# install and require packages
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("reshape"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("WGCNA"))

#ipak(c("dplyr", "Biobase", "Matrix", "Seurat", "parallel", "reshape", "reshape2", "WGCNA"))

message("Packages loaded")

# verify Seurat version
#stopifnot(as.character(packageVersion("Seurat"))=='3.0.0.9000')

######################################################################
################### GET COMMAND LINE OPTIONS #########################
######################################################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 

opt <- parse_args(OptionParser(option_list=option_list))

resume <- opt$resume
prefixData <- opt$prefixData 
dirTmp <- opt$dirTmp
prefixRun <- opt$prefixRun

# load saved image?
if (!is.null(resume)) {
  message(paste0("loading from ", resume))
  tryCatch({load(file=sprintf("%s%s_%s_%s_image.RData.gz", dirTmp, prefixData, prefixRun, resume))},
           error = function(x) {stop(paste0(resume, " session image file not found in ", dirTmp))})
  
}

# get opt again because loading from checkpoint will have overwritten it..
opt <- parse_args(OptionParser(option_list=option_list))

pathDatExpr <- opt$pathDatExpr 
pathMetadata <- opt$pathMetadata
dirProject <- opt$dirProject
dirTmp <- opt$dirTmp
dataType = opt$dataType
autosave <- opt$autosave
assayUse <- opt$assayUse
slotUse <- opt$slotUse
colIdents <- opt$colIdents
vec_colAnnot <- opt$vec_colAnnot
if (!is.null(vec_colAnnot)) vec_colAnnot <- eval(parse(text=vec_colAnnot))
varsToRegress <- opt$varsToRegress
if (!is.null(varsToRegress)) varsToRegress <- eval(parse(text=varsToRegress))
minGeneCells <- opt$minGeneCells
minCellClusterSize <- opt$minCellClusterSize
featuresUse <- opt$featuresUse
#genesPCA <- opt$genesPCA
nFeatures <- opt$nFeatures
corFnc <- opt$corFnc 
networkType <- opt$networkType
if (grepl("hybrid", networkType)) networkType <- eval(parse(text=opt$networkType))
nRepTOM <- opt$nRepTOM
hclustMethod <- opt$hclustMethod
minClusterSize <- eval(parse(text=opt$minClusterSize))
deepSplit <- eval(parse(text=opt$deepSplit))
moduleMergeCutHeight <- eval(parse(text=opt$moduleMergeCutHeight))
pamStage <- eval(parse(text=opt$pamStage))
kMReassign <- opt$kMReassign
kMSignifFilter <- opt$kMSignifFilter
nPC <- opt$nPC
nRepJackStraw <- opt$nRepJackStraw
fuzzyModMembership <- opt$fuzzyModMembership
RAMGbMax <- opt$RAMGbMax

######################################################################
########################## FUNCTIONS (AGAIN) #########################
######################################################################

if (!is.null(resume)) source(file = paste0(current.dir, "rwgcna_functions_seurat3.0.R"))

######################################################################
############################## CONSTANTS #############################
######################################################################

if (!file.exists(pathDatExpr) & is.null(resume)) stop("Input data path not found and no previous session to resume")

# if no output directory provided, use the dir one level above that of the input file.

if (is.null(dirProject)) {
  pos <- tail(gregexpr("/", pathDatExpr)[[1]],2)[1]
  dirProject = substr(pathDatExpr, 1, pos)
}

# if specified output directory doesn't exist, create it 
if (!file.exists(dirProject)) {
  dir.create(dirProject) 
  message("Project directory not found, new one created")
}

dirPlots = paste0(dirProject,"plots/")
if (!file.exists(dirPlots)) dir.create(dirPlots) 

dirTables = paste0(dirProject,"tables/")
if (!file.exists(dirTables)) dir.create(dirTables)

dirRObjects = paste0(dirProject,"RObjects/")
if (!file.exists(dirRObjects)) dir.create(dirRObjects)

dirLog = paste0(dirProject,"log/")
if (!file.exists(dirLog)) dir.create(dirLog)

if (is.null(dirTmp)) dirTmp <- paste0(dirProject,"tmp/")
if (!file.exists(dirTmp)) dir.create(dirTmp)

flagDate = substr(gsub("-","",as.character(Sys.Date())),3,1000)
tStart <- as.character(Sys.time())

# source parameter values
source(file = paste0(current.dir, "rwgcna_params_seurat3.0.R"))

set.seed(randomSeed)

######################################################################
########################## PREPROCESS DATA ###########################
######################################################################

if (is.null(resume)) {
  
  ######################################################################
  ############################ CHECK INPUT #############################
  ######################################################################
  
  if (minGeneCells < 0) stop("minGeneCells must be a non-negative integer") 
  
  if (!featuresUse %in% c("var.features", "PCLoading", "JackStrawPCLoading", "JackStrawPCSignif")) stop('featuresUse must be one of "var.features", "PCLoading", "JackStrawPCLoading", "JackStrawPCSignif"')
  
  if (nRepJackStraw < 1 & featuresUse %in% c("JackStrawPCLoading", "JackStrawPCSignif")) stop("nRepJackStraw must be 0 or higher") 
    
  if (!corFnc %in% c("cor", "bicor")) stop("corFnc must be one of 'cor' for Pearson's or 'bicor' for biweighted midcorrelation")
  
  if (!networkType %in% c('signed', 'unsigned', 'signed hybrid')) stop("networkType must be one of 'signed', 'unsigned' or 'signed hybrid' (not 'signed_hybrid')")
  
  if (! (nRepTOM >= 0)) stop("nRepTOM must be non-negative")
  
  if (length(c(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight))>4) warning("Comparing different parameters increases processing time")
  
  if (min(minClusterSize) < 5) stop("minClusterSize must be a vector of integers over 5")
  
  if (max(deepSplit) > 4 | min(deepSplit) < 0) stop("deepSplit must be a vector of integers between 0 and 4")
  
  if (min(moduleMergeCutHeight) < 0 | max(moduleMergeCutHeight) > 1) stop("moduleMergeCutHeight must be a vector of doubles between 0 and 1, recommended range is between 0.1 and 0.2")
  
  if (!fuzzyModMembership %in% c("kME", "kIM")) {
    warning("Invalid fuzzeModuleMembership value -  reset to kIM (default)")
    fuzzyModMembership <- "kIM"
  }
  
  if (!dataType %in% c("sc", "bulk")) stop("dataType must be 'sc' or 'bulk'")
  
  if (!RAMGbMax >= 1 | !RAMGbMax <= 1000) stop("verify RAMGbMax value")
  
  ######################################################################
  ######################## CREATE ERROR LOG ############################
  ######################################################################
  
  pathLogCellClustersDropped <- paste0(dirLog, prefixData, "_", prefixRun, "_cellClustersDropped.txt")

  ######################################################################
  ################# LOAD AND SUBSET EXPRESSSION DATA ###################
  ######################################################################
  
  message("Loading expression data..")
  
  if (grepl(pattern = "\\.loom", pathDatExpr)) {
    tmp <-  connect(filename = pathDatExpr, mode = "r+")
  } else {
    tmp <- load_obj(pathDatExpr)
  } 
  
  ######################################################################
  ########################### LOAD METADATA ############################
  ######################################################################
    
  if (!is.null(pathMetadata)) {
    metadata <- load_obj(pathMetadata) 
    if (any(duplicated(metadata[,1]))) metadata[,1] <- make.unique(metadata[,1]) 
    rownames(metadata) <- metadata[,1]
    metadata <- metadata[,-1]
  } else {
    metadata <- NULL
  }
  
  # Convert to suitable new Seurat format
  if (any(c("seurat", "Seurat") %in% class(tmp))) {
    seurat_obj <- tmp#if (tmp@version!='3.0.0.9000') {
    #   tryCatch({
    #     Seurat::UpdateSeuratObject(object = tmp)}, error=function(err) {
    #       counts <- tmp@raw.data
    #       if (any(duplicated(rownames(counts)))) rownames(counts) <- make.unique(rownames(counts)) # should not be possible
    #       CreateSeuratObject(counts=counts, 
    #                          project = prefixData,
    #                          assay="RNA", 
    #                          min.cells=0,
    #                          min.features=0, 
    #                          meta.data = tmp@meta.data)
    #   })
    # }  else  {
    #   tmp
    # }
  } else if ("loom" %in% class(tmp)) {
    seurat_obj <- Seurat::Convert(from=tmp, to="seurat")
    if (!is.null(metadata)) seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
  } else { # not a loom object, nor a seurat object
    # Get genes, make them unique, make them into rownames
    if (any(duplicated(tmp[,1]))) tmp[,1] <- make.unique(tmp[,1]) 
    rownames(tmp) <- tmp[,1]
    tmp <- tmp[,-1]
    
    seurat_obj <- CreateSeuratObject(counts= tmp, 
                                     project= paste0(prefixData, "_", prefixRun),
                                     assay = "RNA",
                                     meta.data = metadata) #todo is this
  }
  
  rm(tmp)

  ######################################################################
  ############################## ANNOTATION ############################
  ######################################################################
  
  # replace inappropriate characters
  for (colMetadata in c(colIdents, vec_colAnnot)) {
    if (!colMetadata %in% colnames(seurat_obj@meta.data)) stop(paste0(colMetadata), " not found in metadata")
    seurat_obj@meta.data[[colMetadata]] <- gsub("[^a-z|0-9|_]", "_", seurat_obj@meta.data[[colMetadata]], ignore.case = T)
    seurat_obj@meta.data[[colMetadata]] <- gsub("__", "_", seurat_obj@meta.data[[colMetadata]], ignore.case=T)
  }

  ######################################################################
  ########################### SET DEFAULT ASSAY ########################
  ######################################################################
  
  DefaultAssay(seurat_obj) <- assayUse
  
  ######################################################################
  #################### SET SEURAT IDENT AND SUBSET NAMES ###############
  ######################################################################
  
  # Set seurat object ident
  if (!is.null(colIdents)) {
    Idents(seurat_obj) <- seurat_obj@meta.data[[colIdents]]
  }
  
  vec_idents <- Idents(seurat_obj) %>% as.character
  sNames_0 <- vec_idents %>% unique %>% sort
  
  mat_addAnnot <- if (!is.null(vec_colAnnot)) seurat_obj@meta.data[,vec_colAnnot] else NULL

  ######################################################################
  ################ determine organism and gene names ###################
  ######################################################################
  
  geneNames <- dataOrganism <- NULL
  
  if (sum(grepl("ENSMUSG",rownames(seurat_obj), ignore.case=F))/nrow(seurat_obj) > 0.9) {
    geneNames <- "ensembl"
    dataOrganism <- "mmusculus"
  } else if (sum(grepl("ENSG",rownames(seurat_obj), ignore.case=F))/nrow(seurat_obj)>0.9) {
    geneNames <- "ensembl"
    dataOrganism <- "hsapiens"
  } else if (sum(grepl("[a-z]", rownames(seurat_obj), ignore.case=F))/nrow(seurat_obj)>0.9) {
    geneNames <- "symbol"
    dataOrganism <- "mmusculus"
  } else if (sum(grepl("[A-Z]", rownames(seurat_obj), ignore.case=F))/nrow(seurat_obj) > 0.9) {
    geneNames <- "symbol"
    dataOrganism <- "hsapiens"
  }
  
  ######################################################################
  ################# IF SYMBOL, REMAP TO ENSEMBL ID #####################
  ######################################################################
  
  # mapTo <- if (geneNames == "ensembl") "symbol" else "ensembl"
  # df_genes <- if (geneNames=="ensembl") data.frame("ensembl"=rownames(seurat_obj)) else data.frame("symbol"=rownames(seurat_obj))
  # 
  # if (!is.null(path_df_mapping)) {
  # 
  #   df_ensemblMapping <- load_obj(path_df_ensemblMapping[[dataOrganism]])
  # 
  #   df_genes <- gene_map(dataIn=df_genes,
  #                            colGene = geneNames,
  #                            mapping=df_ensemblMapping,
  #                            from=geneNames,
  #                            to=mapTo,
  #                            replace = F,
  #                            na.rm = T)
  # 
  # }
  
  ######################################################################
  ################### COMPUTE PERCENT MITO, PERCENT RIBO ###############
  ######################################################################
  
  riboMitoMeta <- NULL
  
  if (any(c("percent.ribo", "percent.mito") %in% varsToRegress) & geneNames=="symbol" & !all(c("percent.ribo", "percent.mito") %in% names(seurat_obj@meta.data))) {
    if (dataType=="sc") {
      
      vec_logicalMito.genes <- grepl(pattern = "^mt-", x = df_genes[["symbol"]], ignore.case=T)
      vec_logicalRibo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = df_genes[["symbol"]], ignore.case=T)
      colSums_tmp <- GetAssayData(object=seurat_obj, slot = "counts") %>% colSums
      
      riboMitoMeta <- data.frame(percent.mito=colSums(x = GetAssayData(object=seurat_obj, slot = "counts")[vec_logicalMito.genes,])/colSums_tmp,
                                 percent.ribo = colSums(x = GetAssayData(object=seurat_obj, slot = "counts")[vec_logicalRibo.genes,])/colSums_tmp,
                                 row.names=colnames(seurat_obj))
      seurat_obj <- AddMetaData(seurat_obj, metadata =riboMitoMeta)
    }
  } else { 
    if (any(c("percent.ribo", "percent.mito") %in% varsToRegress)) varsToRegress <- varsToRegress[!varsToRegress %in% c("percent.ribo", "percent.mito")]
  }
  ######################################################################
  ######## DO SEURAT PROCESSING ON FULL EXPRESSION MATRIX ##############
  ######################################################################
  
  # Normalise
  if (assayUse == "RNA")  {
    if (all(dim(GetAssayData(seurat_obj, slot="data")))==0) seurat_obj <- NormalizeData(object = seurat_obj)#, display.progress = T)
  }
  # Scale and regress out confounders
  #if (!is.null(varsToRegress)) varsToRegress <- varsToRegress[varsToRegress %in% names(seurat_obj@meta.data)]
  
  # if (is.null(seurat_obj@assays[[assayUse]]@scale.data)) {
  #   seurat_obj <- ScaleData(object = seurat_obj,
  #                         vars.to.regress = vars.to.regress,
  #                         model.use="linear",
  #                         do.par=T,
  #                         num.cores = min(10, detectCores_plus(Gb_max = RAMGbMax, 
  #                                                              additional_Gb = as.numeric(object.size(seurat_obj@data))/1024^3)-1),
  #                         do.scale=T,
  #                         do.center=T,
  #                         display.progress = T)
   
  # # Make sure we close socket workers
  # invisible(gc()); invisible(R.utils::gcDLLs())
  # 
  # scale_data <- seurat_obj@scale.data
  #ident <- seurat_obj@ident
  #}
  
  datExprNormFull <- if (assayUse=="RNA") {
    GetAssayData(object=seurat_obj, slot="data") %>% as.matrix 
  } else if (assayUse = "integrated") {
    GetAssayData(object=seurat_obj) %>% as.matrix
  }  else if (assayUse =="SCT") {
    GetAssayData(object=seurat_obj, slot="counts") %>% as.matrix
  }
  # Save scale and regressed whole expression matrix with ensembl rownames for later use
  message("Saving full expression matrix")
  
  saveRDS(datExprNormFull, file = sprintf("%s%s_%s_datExprNormFull.RDS.gz", 
                                  dirTmp, 
                                  prefixData, 
                                  prefixRun), 
          compress = "gzip")
  #save(ident, file = sprintf("%s%s_%s_ident.RData", dirTmp, prefixData, prefixRun))
  rm(datExprNormFull)

  ######################################################################
  ######## DO SEURAT PROCESSING ON SUBSETTED EXPRESSION MATRICES #######
  ######################################################################
  
  message("Subsetting the dataset")
  subsets <- lapply(sNames_0, function(name) {
    # using createSeuratObject instead of subset allows us to filter on min.cells and min.features
    CreateSeuratObject(counts = GetAssayData(object=seurat_obj, 
                                             assay=assayUse, 
                                             slot="counts")[,Idents(seurat_obj)==name], 
                       project = prefixData,
                       assay = assayUse,
                       min.cells = minGeneCells,
                       #min.features = 1000,
                       meta.data = seurat_obj@meta.data[Idents(seurat_obj)==name,]
                       )})
  
  # Save stats for outputting at the end
  #n_cells_subsets = table(Idents(seurat_obj)) #warning: this is in the wrong order because Idents is a factor..
  n_cells_subsets = sapply(subsets, ncol)
  vec_logicalCellClustOK <- if (dataType=="sc") n_cells_subsets > minCellClusterSize else !logical(length = length(subsets))
  #sapply(subsets, function(seuratObj) ncol(seuratObj), simplify = T)
  n_genes_subsets = sapply(subsets, nrow, simplify=T)
  
  # Free up space 
  rm(seurat_obj)
  
  ### Filter out cell types that have too few cells
  # We do this do avoid downstream problems with Seurat or WGCNA. 
  # E.g. Seurat will ScaleData will fail if regressing out variables when there are only 2 cells in the data.
  
  # if (dataType=="sc") {
  #   message(paste0("Filtering out cell clusters with fewer than ", minCellClusterSize, " cells"))
  #   vec_logicalCellClustOK <- n_cells_subsets>=minCellClusterSize
  # } else if (dataType=="bulk") { # no need to filter out any subsets
  #   vec_logicalCellClustOK = !logical(length = length(subsets))
  # }
  
  if (!all(vec_logicalCellClustOK)) {
    log_entry <- paste0(sNames_0[!vec_logicalCellClustOK], ": had <", minCellClusterSize," cells and was therefore dropped")
    warning(log_entry)
    cat(text = log_entry, file =  pathLogCellClustersDropped, append=T, sep = "\n")
    # Filter
    #subsets <- subsets[vec_logicalCellClustOK]
  }
  
  sNames_1 <- sNames_0[vec_logicalCellClustOK]
  names(subsets) <- sNames_1 
  
  # # Filter genes expressed in fewer than min.cells in a subset as these will also lead to spurious associations and computational difficulties
  # if (dataType=="sc") {
  #   
  #   message(paste0("Filtering out genes expressed in fewer than ", minGeneCells, " cells"))
  #   # see https://github.com/satijalab/seurat/issues/147
  #   # The authors prefer that you do it before creating a seurat object
  #   
  #   outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "FilterGenes.txt")
  #   list_iterable = list("X"=subsets)
  #   fun = function(seurat_obj) {
  #     num.cells <- rowSums(GetAssayData(seurat_obj, assay= assayUse, slot="counts") > 0)
  #     vec_logicalfeaturesUse <- which(x = num.cells >= minGeneCells)
  #     if (sum(genes.use)>0) {
  #       SetAssayData(object= seurat_obj, assay=assayUse, slot="counts", new.data = GetAssayData(seurat_obj, assay=assayUse, slot="counts")[vec_logicalfeaturesUse, ]) # will this work?
  #       SetAssayData(object= seurat_obj, assay=assayUse, slot="data", new.data = GetAssayData(seurat_obj, assay=assayUse, slot="data")[vec_logicalfeaturesUse, ]) # will this work?
  #       seurat_obj@assays$RNA@meta.features <- seurat_obj@assays$RNA@meta.features[vec_logicalfeaturesUse,]
  #       return(seurat_obj)
  #     } else {
  #       return(NULL)
  #     }
  #   }
  #   subsets <- safeParallel(fun=fun, 
  #                           list_iterable=list_iterable, 
  #                           outfile=outfile, 
  #                           minGeneCells=minGeneCells)
  # }
  # 
  # # Filter
  # vec_logicalCellClustOK <- sapply(subsets, function(subset) {!is.null(rownames(subset))}, simplify=T)
  # 
  # if (!all(vec_logicalCellClustOK)) {
  #   log_entry <-paste0(sNames_1[!vec_logicalCellClustOK], ": had no genes left after filtering out genes expressed in fewer than ", minGeneCells, " cells, therefore dropped")
  #   warning(log_entry)
  #   cat(log_entry, file = pathLogCellClustersDropped, append=T, sep = "\n")
  #   # filter
  #   subsets <- subsets[vec_logicalCellClustOK]
  # }
  # 
  # sNames_1 <- sNames_1[vec_logicalCellClustOK]
  # 
  # subsets <- lapply(subsets, function(seuSubset) { 
  #   idx <- rownames(GetAssayData(object=seuSubset, assay=assayUse, slot="scale.data")) %in% rownames(seuSubset)
  #   SetAssayData(object=seuSubset, assay=assayUse, slot="scale.data", 
  #                new.data=GetAssayData(object=seuSubset, assay=assayUse, slot="scale.data")[idx,])
  #   seuSubset@assays$RNA@meta.features <- seuSubset@assays$RNA@meta.features[idx,]
  # })
  # Need to renormalise after removing genes
  # if (dataType=="sc") {
  #   outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "NormalizeData.txt")
  #   list_iterable = list("X"=subsets)
  #   subsets <- safeParallel(fun=NormalizeData,
  #                           list_iterable=list_iterable,
  #                           outfile=outfile)
  # }

  message("Finding variable features")
  fun <- function(seurat_obj, seuName) {
    tryCatch({
      
      FindVariableFeatures(object = seurat_obj,
                        selection.method = "vst",
                        do.plot=F,
                        nfeatures = if (featuresUse=="var.features") {
                          min(nFeatures, nrow(seurat_obj))
                          } else if (featuresUse %in% c('JackStrawPCLoading', 'JackStrawPCSignif')) {
                            nrow(seurat_obj)
                          } else 1000)

      }, error = function(err) {
        warning(paste0(seuName, ": FindVariableFeatures failed with selection.method = 'vst'. Trying selection.method = 'mean.var.plot' instead."))
        FindVariableFeatures(object = seurat_obj,
                           selection.method = "mean.var.plot",
                           do.plot=F)
      })
  }
  
  list_iterable = list("seurat_obj"=subsets, "seuName" = names(subsets))
  outfile <- paste0(dirLog, prefixData, "_", prefixRun, "_FindVariableFeatures_log.txt") 
  subsets <- safeParallel(fun=fun,
                          list_iterable=list_iterable,
                          outfile=outfile)
  
  # Scale and regress
  if (featuresUse %in% c('PCLoading', 'JackStrawPCLoading', 'JackStrawPCSignif') | 
      slotUse =="scale.data") {
   
    message("Scaling data subsets")
    
    if (!is.null(varsToRegress)) varsToRegress = varsToRegress[varsToRegress %in% colnames(subsets[[1]]@meta.data)] 
    
    # Scale data and regress out confounders
    fun<- function(seurat_obj, name) {
      tryCatch({
        ScaleData(object = seurat_obj,
                  # choice of features to scale constrains the genes for which we can get PCA loadings downstream
                  features = if (featuresUse=="var.features") VariableFeatures(object=seurat_obj, assay = assayUse) else rownames(seurat_obj), # 
                  vars.to.regress = if (dataType=="sc" & length(varsToRegress)>0) varsToRegress else NULL,
                  #model.use="linear",
                  do.scale=T,
                  #model.use = "negbinom",
                  do.center=do.center,
                  block.size=15000,
                  assay=assayUse,
                  min.cells.to.block=5000)},
        error = function(err) {
          stop(paste0(name,": ScaleData failed with the error: ", err))
        })
    }
    list_iterable = list(seurat_obj = subsets, name = names(subsets))
    outfile <- paste0(dirLog, prefixData, "_", prefixRun, "_ScaleData_log.txt") 
    subsets <- mapply(FUN = fun,seurat_obj=subsets, name=names(subsets))
    #subsets <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
  }

  # PCA 
  if (featuresUse %in% c('PCLoading', 'JackStrawPCLoading', 'JackStrawPCSignif')) {
    
    message("Performing Principal Component Analysis")
    
    start_time <- Sys.time()
    
    fun = function(seurat_obj, name) {
      seurat_tmp <- tryCatch({ 
        out <- RunPCA(object = seurat_obj,
               assay = assayUse,
               features = if (featuresUse %in% c('JackStrawPCLoading', 'JackStrawPCSignif')) rownames(seurat_obj) else VariableFeatures(seurat_obj),
               #npcs = nPC,
               npcs = min(nPC, min(length(VariableFeatures(seurat_obj))%/% 2, ncol(seurat_obj) %/% 2)),
               #rev.pca = F,
               weight.by.var = T, # weighs cell embeddings 
               do.print = F,
               seed.use = randomSeed,
               maxit = maxit, # set to 500 as default
               fastpath = fastpath, 
               verbose=T) 
        }, 
        error = function(err) {
          message(paste0(name, ": RunPCA's IRLBA algorithm failed with the error: ", err))
          message("Trying RunPCA with var.features, fewer components, double max iterations and fastpath == F")
          tryCatch({
            out <- RunPCA(object = seurat_obj,
                   assay = assayUse,
                   features = if (featuresUse %in% c('JackStrawPCLoading', 'JackStrawPCSignif')) rownames(seurat_obj) else VariableFeatures(seurat_obj),
                   #npcs = nPC %/% 1.5,#min(nPC, length(seurat_obj@assays[[assayUse]] %>% '@'('var.features')) %/% 3),
                   npcs = min(nPC, min(length(VariableFeatures(seurat_obj))%/% 3, ncol(seurat_obj) %/% 3)),
                   weight.by.var = T,
                   seed.use = randomSeed,
                   maxit = maxit*2, # set to 500 as default
                   fastpath = F,
                   verbose=T)
            
            
          }, error = function(err1) {
            message(paste0(name, ": RunPCA's IRLBA algorithm failed again with error: ", err1))
            message("Returning the original Seurat object with empty dr$pca slot")
            return(seurat_obj)
            }) 
        })
      # if (genesPCA == "all") {
      #   seurat_tmp@reductions$pca@feature.loadings.full <- seurat_tmp@dr$pca@gene.loadings
      # }
      return(seurat_tmp)
    }
    
    list_iterable <- list("seurat_obj" = subsets, "name" = names(subsets))
    
    outfile <- paste0(dirLog, prefixData, "_", prefixRun, "_PCA_log.txt")
    
    subsets <- safeParallel(fun=fun, 
                            list_iterable=list_iterable, 
                            outfile=outfile)#, 
                            #nPC=nPC,
                            #randomSeed=randomSeed,
                            #maxit=maxit,
                            #fastpath=fastpath)
    
    end_time <- Sys.time()
    
    message(sprintf("PCA done, time elapsed: %s seconds", round(end_time - start_time,2)))
    
  }
  
  if (featuresUse=='PCLoading') {
    
    list_datExpr <- lapply(subsets, function(seurat_obj) {
      seurat_obj <- ProjectDim(object=seurat_obj, 
                               reduction = "pca", 
                               verbose=F,
                               overwrite=F)
      loadingsAbs <- seurat_obj@reductions$pca@feature.loadings.projected %>% abs 
      apply(loadingsAbs, MARGIN=1, FUN=max) %>% sort(., decreasing=T) %>% 
        '['(1:(min(nFeatures, nrow(loadingsAbs)))) %>% names -> namesFeaturesUse
      GetAssayData(object=seurat_obj,slot = slotUse, assay = assayUse)[namesFeaturesUse,] %>% t
      })
    
    } else if (featuresUse %in% c('JackStrawPCLoading', 'JackStrawPCSignif')) {
    
    # Select significant PCs using empirical p-value based on JackStraw resampling
    # Source: https://rdrr.io/cran/Seurat/src/R/plotting.R
    # score the PCs using JackStraw resampling to get an empirical null distribution to get p-values for the PCs 
    # based on the p-values of gene loadings. Use these to select genes
    
    message(sprintf("Performing JackStraw with %s replications to select genes that load on significant PCs", nRepJackStraw))
    fun = function(seurat_obj, name) {
      # seurat_obj <- ProjectDim(object=seurat_obj, 
      #                          reduction = "pca", 
      #                          verbose=F,
      #                          overwrite=T)
      wrapJackStraw(seurat_obj_sub = seurat_obj, 
                    nRepJackStraw = nRepJackStraw, 
                    pvalThreshold = pvalThreshold,
                    featuresUse=featuresUse,
                    nFeatures = nFeatures,
                    nPC = nPC)

    }
    
    list_iterable = list(seurat_obj = subsets, name = names(subsets))
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_JackStraw.txt")
    list_datExpr<- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
    
  } else if (featuresUse == "var.features") {
    list_datExpr <- lapply(subsets, function(seurat_obj) GetAssayData(seurat_obj, assay=assayUse, slot=slotUse)[rownames(seurat_obj) %in% VariableFeatures(seurat_obj),] %>% t)
  } 
  
  n_genes_subsets_used <- sapply(list_datExpr, ncol) 
  
  # Clear up
  rm(subsets) 
  invisible(gc())
  
  ######################################################################
  ########################### CHECKPOINT ###############################
  ######################################################################
  
  # Save or load session image 
  resume = "checkpoint_1"
  if (autosave) save.image(file=sprintf("%s%s_%s_checkpoint_1_image.RData.gz", dirTmp, prefixData, prefixRun), compress = "gzip")
  
} 

if (resume == "checkpoint_1") {
  
  ######################################################################
  ###################### (UN) LOAD PACKAGES ############################
  ######################################################################
  
  # Unload packages
  try(detach("package:Seurat", unload=TRUE))
  
  # Free up DLLs
  invisible(R.utils::gcDLLs())
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("WGCNA"))
  
  #pkgs <- c("dplyr", "Biobase", "Matrix", "parallel", "reshape", "reshape2", "WGCNA")
  
  #ipak(pkgs)
  
  disableWGCNAThreads()
  
  ######################################################################
  ####### PICK SOFT THRESHOLD POWER FOR ADJACENCY MATRIX ###############
  ######################################################################
  message("Computing soft threshold powers to maximise the fit of a scale free topology to the adjacency matrix")
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  softPower <- 8 # Set a default value as fall back
  
  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_compute_softPowers.txt")
  list_iterable = list("datExpr" = list_datExpr, "subsetName"=sNames_1)
  fun = function(datExpr, subsetName) {
    powers = c(1:30)
    sft = pickSoftThreshold(data=datExpr,
                            powerVector = powers,
                            blockSize = min(maxBlockSize , ncol(datExpr)), #try to prevent crashing
                            corFnc = corFnc,
                            corOptions =  corOptions,
                            networkType = networkType,
                            verbose = 1)
    pdf(sprintf("%s%s_%s_%s_pickSoftThresholdSFTFit.pdf", dirPlots, prefixData, prefixRun, subsetName),width=10,height=5)
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1],
         -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",
         type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of 0.9
    abline(h=0.90,col="red")
    dev.off()
    # Mean connectivity as a function of the soft-thresholding power
    pdf(sprintf("%s%s_%s_%s_pickSoftThresholdMeanCon.pdf", dirPlots, prefixData, prefixRun, subsetName),width=10,height=5)
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",
         ylab="Mean Connectivity",
         type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    
    # Select softPower: of lower .95 percentile connectivity, if several softPowers achieve 0.9 R.sq, 
    # take smallest softPower; else take best fitting one
    
    fitIndices <- as.data.frame(sft$fitIndices)
    fitIndices %>% dplyr::filter(median.k. <= quantile(median.k.,0.95, na.rm=T)) -> fitIndices_filter 
    
    if (sum(fitIndices_filter$SFT.R.sq >= 0.93) > 1) {
      fitIndices_filter_2 <- fitIndices_filter[fitIndices_filter$SFT.R.sq>=0.93,] 
      fitIndices_filter_2[which.min(fitIndices_filter_2$Power), c(1,2,6)] -> df_sft
    } else {
      (fitIndices_filter %>% dplyr::arrange(desc(SFT.R.sq)))[1,c(1,2,6)]  -> df_sft
    }
    return(df_sft)
  }
  
  list_sft = safeParallel(fun=fun, 
                          list_iterable=list_iterable, 
                          outfile=outfile)#, 
                          #maxBlockSize=maxBlockSize, 
                          #corFnc=corFnc, 
                          #corOptions=corOptions, 
                          #networkType=networkType,
                          #dirPlots=dirPlots, 
                          #prefixData=prefixData, 
                          #prefixRun=prefixRun)
  
  names(list_sft) <- sNames_1
  
  # TOM files will be saved here by consensusTOM 
  setwd(dirTmp)
  
  if (nRepTOM > 0) {
    
    ######################################################################
    ############### RESAMPLE THE DATA FOR ROBUSTNESS #####################
    ######################################################################
    
    message(sprintf("Resampling the expression data and computing the consensus Topological Overlap Matrix with fraction %s with replace = %s and %s replicates", fraction, replace, nRepTOM))
    #message(sprintf("Computing the consensus Topological Overlap Matrix with %s permutations", nRepTOM))
    
    fun <- function(datExpr,name) {
      #tryCatch({
          
        multiExpr <- bootstrap(datExpr=datExpr, 
                      nPermutations = nRepTOM,
                      replace = replace,
                      fraction = fraction,
                      randomSeed = randomSeed)
          
        consensusTOM(multiExpr = multiExpr, 
                               checkMissingData = checkMissingData,
                               maxBlockSize = maxBlockSize, 
                               blockSizePenaltyPower = blockSizePenaltyPower, 
                               randomSeed = randomSeed,
                               corType = corType,
                               maxPOutliers = maxPOutliers,
                               quickCor = quickCor,
                               pearsonFallback = pearsonFallback,
                               cosineCorrelation = cosineCorrelation,
                               replaceMissingAdjacencies = replaceMissingAdjacencies,
                               power = list_sft[[name]][["Power"]],
                               networkType = networkType,
                               TOMDenom = TOMDenom,
                               saveIndividualTOMs = saveIndividualTOMs,
                               individualTOMFileNames = paste0(prefixData, "_", prefixRun, "_", name, "_individualTOM-Set%s-Block%b.RData"),
                               networkCalibration = networkCalibration,
                               sampleForCalibration = sampleForCalibration,
                               sampleForCalibrationFactor = sampleForCalibrationFactor,
                               getNetworkCalibrationSamples = getNetworkCalibrationSamples,
                               consensusQuantile = consensusQuantile,
                               useMean = useMean,
                               saveConsensusTOMs = saveConsensusTOMs,
                               consensusTOMFilePattern = paste0(prefixData,"_", prefixRun, "_", name,"_consensusTOM-block.%b.RData"),
                               returnTOMs = F,
                               useDiskCache = T,
                               cacheDir = dirTmp,
                               cacheBase = ".blockConsModsCache",
                               verbose = verbose,
                               indent = indent)

      #}, error = function(err) {
        # message(paste0(name, ": consensusTOM failed, computing normal TOM"))
        # adjacency = adjacency(datExpr=list_datExpr[[name]], 
        #                       type=type, 
        #                       power = list_sft[[name]]$Power, 
        #                       corFnc = corFnc, 
        #                       corOptions = corOptions)
        #                       
        # consTomDS = TOMsimilarity(adjMat=adjacency,
        #                           TOMType=TOMType,
        #                           TOMDenom=TOMDenom,
        #                           verbose=verbose,
        #                           indent = indent)
        # colnames(consTomDS) <- rownames(consTomDS) <- colnames(list_datExpr[[name]])
        # save(consTomDS, file=sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", dirTmp, prefixData, prefixRun, name)) # Save TOM the way consensusTOM would have done
     # })
    }
    
    list_iterable <- list("datExpr"=list_datExpr, "name"=sNames_1)
    outfile = outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_consensusTOM.txt")

    additional_Gb = max(as.numeric(sapply(list_iterable, FUN = function(x) {object.size(x)*nRepTOM}, simplify = T)))/1024^3 
    obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
    n_cores <- min(max(sapply(list_iterable, length)), min(detectCores()%/%3, RAMGbMax %/% (obj_size_Gb + additional_Gb))-1)
    
    list_consensus <- safeParallel(fun=fun, 
                                   list_iterable=list_iterable,
                                   outfile=outfile#,
                                   #n_cores=n_cores,
                                   #nPermutations=nRepTOM,
                                   #replace=replace,
                                   #fraction=fraction,
                                   #randomSeed=randomSeed,
                                   #checkMissingData = checkMissingData,
                                   #maxBlockSize = maxBlockSize, 
                                   #blockSizePenaltyPower = blockSizePenaltyPower, 
                                   #randomSeed = randomSeed,
                                   #corType = corType,
                                   #corFnc = corFnc, 
                                   #corOptions = corOptions,
                                   #maxPOutliers = maxPOutliers,
                                   #quickCor = quickCor,
                                   #pearsonFallback = pearsonFallback,
                                   #cosineCorrelation = cosineCorrelation,
                                   #replaceMissingAdjacencies = replaceMissingAdjacencies,
                                   #list_sft = list_sft,
                                   #networkType = networkType,
                                   #type = networkType,
                                   #TOMType=TOMType,
                                   #TOMDenom = TOMDenom,
                                   #saveIndividualTOMs = saveIndividualTOMs,
                                   #prefixData=prefixData, 
                                   #prefixRun=prefixRun,
                                   #networkCalibration = networkCalibration,
                                   #sampleForCalibration = sampleForCalibration,
                                   #sampleForCalibrationFactor = sampleForCalibrationFactor,
                                   #getNetworkCalibrationSamples = getNetworkCalibrationSamples,
                                   #consensusQuantile = consensusQuantile,
                                   #useMean = useMean,
                                   #saveConsensusTOMs = saveConsensusTOMs,
                                   #returnTOMs = F,
                                   #useDiskCache = T,
                                   #cacheDir = dirTmp,
                                   #cacheBase = ".blockConsModsCache",
                                   #verbose = verbose,
                                   #indent = indent,
                                   #list_datExpr=list_datExpr,
                                   #dirTmp = dirTmp
                                   )

  } else if (nRepTOM==0) {
    
    message("Computing the Topological Overlap Matrix")
    
    ######################################################################
    ##### COMPUTE THE ADJACENCY MATRIX AND TOM WITHOUT RESAMPLING ########
    ######################################################################
    
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_TOM_for_par.txt")
    
    fun <- function(datExpr,name) { 
      adjacency = adjacency(datExpr=datExpr, 
                            type=type, 
                            power = list_sft[[name]]$Power, 
                            corFnc = corFnc, 
                            corOptions = corOptions)
      consTomDS = TOMsimilarity(adjMat=adjacency,
                                TOMType=TOMType,
                                TOMDenom=TOMDenom,
                                verbose=verbose,
                                indent = indent)
      colnames(consTomDS) <- rownames(consTomDS) <- colnames(datExpr)
      save(consTomDS, file=sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", dirTmp, prefixData, prefixRun, name)) # Save TOM the way consensusTOM would have done
   }
    list_iterable = list("datExpr"=list_datExpr, "name"=sNames_1)
    invisible(safeParallel(fun=fun, 
                           list_iterable=list_iterable, 
               outfile=outfile#, 
               #type=type, 
               #list_sft=list_sft,
               #corFnc=corFnc, 
               #corOptions=corOptions, 
               #TOMType=TOMType, 
               #TOMDenom=TOMDenom, 
               #verbose=verbose, 
               #indent=indent, 
               #prefixData=prefixData, 
               #prefixRun=prefixRun, 
               #dirTmp=dirTmp
               ))
    list_consensus <- lapply(list_datExpr, function(datExpr)list("goodSamplesAndGenes"=list("goodGenes"=rep(TRUE,ncol(datExpr)))))
  }
  
  #rm(list_datExpr)
  
  #if (nRepTOM > 0) rm(list_multiExpr) #TODO: we'll need it for plotting permuted
  
  # Load TOMs into memory and convert to distance matrices
  fun = function(name) {
    load_obj(sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", dirTmp, prefixData, prefixRun, name))

  }
  list_iterable = list("name"=sNames_1)
  list_consTOM = safeParallel(fun=fun, 
                              list_iterable=list_iterable 
                              #dirTmp=dirTmp, 
                              #prefixData=prefixData, 
                              #prefixRun=prefixRun, 
                              #load_obj=load_obj
                              )
  
  # Filter datExpr to retain genes that were kept by goodSamplesGenesMS, called by consensusTOM 
  # Check whether this step works
  # fun = function(datExpr, consTOM) {
  #   datExpr[,match(colnames(consTOM), colnames(datExpr))]
  # }
  fun = function(datExpr, consensus, consTOM) {
    if (!is.null(consTOM)) {
      if (ncol(as.matrix(consTOM))==ncol(datExpr)) datExpr else datExpr[,consensus$goodSamplesAndGenes$goodGenes]
    } else {
      NULL
    }
  }
  
  list_iterable=list(datExpr=list_datExpr, consensus=list_consensus, consTOM=list_consTOM)
  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_datExpr_gg.txt")
  list_datExpr_gg <- safeParallel(fun=fun, list_iterable=list_iterable)
  
  rm(list_datExpr)
  
  # Convert proximity TOMs to distance matrices 
  list_iterable= list("X"=list_consTOM)
  fun = function(consTOM) { 
    1-as.dist(consTOM) 
    }
  list_dissTOM <- safeParallel(fun=fun, list_iterable=list_iterable)  
  rm(list_consTOM)
  # Save list_dissTOM to harddisk
  names(list_dissTOM) <- sNames_1
  
  saveRDS(list_dissTOM, file=sprintf("%s%s_%s_list_dissTOM.rds.gz", dirTmp, prefixData, prefixRun), compress = "gzip")
  
  ######################################################################
  ######################### CLEAR UP TOMS IN MEMORY AND HD #############
  ######################################################################
  
  rm(list_dissTOM)
  
  for (subsetName in sNames_1) {
    if (file.exists(sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", dirTmp, prefixData, prefixRun, subsetName))) {
      try(file.remove(sprintf("%s%s_%s_%s_consensusTOM-block.1.RData", dirTmp, prefixData, prefixRun, subsetName)))
    }
    # Delete individual TOMs from disk
    indivTOM_paths <- dir(path=dirTmp, pattern=paste0(prefixData, "_", prefixRun, "_", subsetName, "_individualTOM"), full.names = T)
    for (indivTOM in indivTOM_paths) {
      try(file.remove(indivTOM))
    }
  }
  
  ######################################################################
  ############################ CHECKPOINT ##############################
  ######################################################################
  
  resume="checkpoint_2"
  if (autosave) save.image( file=sprintf("%s%s_%s_checkpoint_2_image.RData.gz", dirTmp, prefixData, prefixRun), compress = "gzip")
  
} 

if (resume == "checkpoint_2") {
  
  ######################################################################
  #################### (UN) LOAD PACKAGES #############################
  ######################################################################
  
  # Free up DLLs
  invisible(R.utils::gcDLLs())
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("WGCNA"))
  
  #pkgs <- c("dplyr", "Biobase", "Matrix", "parallel", "reshape", "reshape2", "WGCNA")
  
  #ipak(pkgs)
  
  ######################################################################
  ##################### LOAD AND FILTER DISSTOMS #######################
  ######################################################################
  # TODO: Could filter earlier
  
  list_dissTOM_path <- dir(path = dirTmp, pattern = paste0(prefixData, "_", prefixRun, "_list_dissTOM"), full.names = T)
  list_dissTOM <- load_obj(list_dissTOM_path)
  
  # Remove NULL
  vec_logicalCellClustOK <- !sapply(list_dissTOM, is.null)

  if (!all(vec_logicalCellClustOK)) {
    log_entry <-paste0(sNames_1[!vec_logicalCellClustOK], ": Topological Overlap Matrix step failed, dropped from analysis")
    warning(log_entry)
    cat(log_entry, file = pathLogCellClustersDropped, append=T, sep = "\n")
  }
  # Edit objects
  list_dissTOM <- list_dissTOM[vec_logicalCellClustOK]
  list_datExpr_gg <- list_datExpr_gg[vec_logicalCellClustOK]
  sNames_2 <- sNames_1[vec_logicalCellClustOK]
  
  ######################################################################
  ####################### CLUSTER ON THE TOM ###########################
  ######################################################################
  
  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_parHclust.txt")
  fun <- function(dissTOM) {hclust(d=dissTOM, method=hclustMethod)}
  list_iterable = list("X"=list_dissTOM)
  
  list_geneTree <- safeParallel(fun=fun, 
                                list_iterable=list_iterable, 
                                outfile = outfile, 
                                #hclustMethod=hclustMethod
                                )
  names(list_geneTree) = sNames_2 # used for PlotDendro
  
  ######################################################################
  ######################### COMPUTE MODULES ############################
  ######################################################################
  
  message("Computing modules on the Topological Overlap Matrix")
  
  dims = sapply(list(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight), length) # list of number of values per parameter
  n_combs <- prod(dims) # Total number of combinations of parameter values
  list_comb = vector("list", length = n_combs)
  
  # Make a vector of different parameter combinations
  k=1
  for (size in 1:dims[1]) {
    for (ds in 1:dims[2]) { # Gandal et al 2018 use c(50,100, 200)
      for (pam in 1:dims[3]) {
        for (dthresh in 1:dims[4]) {
          list_comb[[k]] <- c(minClusterSize[size],deepSplit[ds], pamStage[pam], moduleMergeCutHeight[dthresh])
          k = k+1
        }
      }
    }
  }
  
  list_list_comb <- lapply(1:length(sNames_2), function(i) return(list_comb)) # just multiply..
  
  attr(x=list_list_comb, which = "names") <- sNames_2 # for some reason, just setting names(list_list_comb) <- sNames_2 filled the list with NULLs...
  
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_cutreeHybrid_for_vec.txt")
  
  fun <- function(geneTree,dissTOM) {
    lapply(list_comb, function(comb) {
      tryCatch({
        cutreeHybrid_for_vec(comb=comb, 
                                    geneTree=geneTree, 
                                    dissTOM = dissTOM,
                                    maxPamDist = maxPamDist, 
                                    useMedoids = useMedoids)
      }, error= function(err) {
      message(paste0(comb), ": cutreeHybrid failed")
      return(NULL)
     })
      }
    )}
  
  list_iterable = list("geneTree"=list_geneTree, "dissTOM"=list_dissTOM)
  
  list_list_cutree <- safeParallel(fun=fun, 
                                   list_iterable=list_iterable, 
                                   outfile=outfile#, 
                                   #list_comb=list_comb, 
                                   #maxPamDist=maxPamDist, 
                                   #useMedoids=useMedoids
                                   )
  
  # free up memory
  if (fuzzyModMembership=="kME") rm(list_dissTOM)
  
  # Produce labels for plotting the modules found by different parameters
  list_plot_label <- lapply(list_comb, function(comb) plotLabel_for_vec(comb=comb)) # list of labels for plot
  # Make a copy of the labels for each subset
  list_list_plot_label = list()
  list_list_plot_label <- lapply(sNames_2, function(name) list_list_plot_label[[name]] = list_plot_label)
  
  # Remove NULLs in nested lists where cutreeHybrid failed
  list_vec_idx.combcutreeNULL <- lapply(list_list_cutree, function(list_cutree) {
    sapply(list_cutree, is.null)
  })
  
  names(list_vec_idx.combcutreeNULL) = sNames_2
  
  # filter
  if (any(sapply(list_vec_idx.combcutreeNULL, function(idx) idx))) {
    
    list_list_cutree <- mapply(function(list_cutree, vec_idx.combcutreeNULL) {
      list_cutree[!vec_idx.combcutreeNULL]
    }, list_cutree=list_list_cutree, vec_idx.combcutreeNULL = list_vec_idx.combcutreeNULL, SIMPLIFY=F)
    
    list_list_comb <- mapply(function(list_comb, vec_idx.combcutreeNULL) {
      list_comb[!vec_idx.combcutreeNULL]
    }, list_comb=list_list_comb, vec_idx.combcutreeNULL = list_vec_idx.combcutreeNULL, SIMPLIFY=F)
    
    list_list_plot_label <- mapply(function(list_plot_label, vec_idx.combcutreeNULL) {
      list_plot_label[!vec_idx.combcutreeNULL]
    }, list_plot_label = list_list_plot_label, vec_idx.combcutreeNULL = list_vec_idx.combcutreeNULL, SIMPLIFY=F)
  }
  
  names(list_list_cutree) <- names(list_list_comb) <- names(list_list_plot_label) <- sNames_2
  
  # Filter at the cell cluster level
  vec_logicalCellClustOK <- sapply(list_list_cutree, function(x) length(x)>0, simplify=T) %>% unlist
  
  if (!all(vec_logicalCellClustOK)) {
    log_entry <-paste0(sNames_2[!vec_logicalCellClustOK], " : cutreeHybrid failed, removed from analysis")
    warning(log_entry)
    cat(log_entry, file = pathLogCellClustersDropped, append=T, sep = "\n")
  }
  
  list_list_cutree <- list_list_cutree[vec_logicalCellClustOK]# Filter(f=length, x= list_list_cutree)
  list_list_comb <-list_list_comb[vec_logicalCellClustOK] #Filter(f=length, x=list_list_comb)
  list_list_plot_label <-list_list_plot_label[vec_logicalCellClustOK]# Filter(f=length, x=list_list_plot_label)
  list_datExpr_gg <- list_datExpr_gg[vec_logicalCellClustOK]#list_datExpr_gg[names(list_datExpr_gg) %in% sNames_3]
  
  sNames_3 <- sNames_2[vec_logicalCellClustOK]# names(list_list_cutree) 
  
  ######################################################################
  ######################### MERGE CLOSE MODULES ########################
  ######################################################################
  
  message(paste0("Merging close modules"))
  
  if (fuzzyModMembership=="kME") {
    
    fun = function(list_comb,list_cutree, datExpr, cellType) {
      mapply(function(comb, cutree) {
        colors <- rep("grey", times=ncol(datExpr))
        MEs <- NULL 
        out <- list("colors"= colors, "MEs" = MEs)
        if (length(unique(cutree$labels))>1) {
          if (length(unique(cutree$labels))>2) {
            
            merged <- try({mergeCloseModules(exprData=as.matrix(datExpr), 
                                             colors = cutree$labels, 
                                             impute =T,
                                             corFnc = corFnc,
                                             corOptions = corOptions,
                                             cutHeight = comb[4],
                                             iterate = T,
                                             getNewMEs = F,
                                             getNewUnassdME = F)
            })
            
            if (class(merged)=="try-error") {
              colors = rep("grey", times=ncol(datExpr))
            } else {
              colors = labels2colors(merged$colors)
            }
            
          } else {
            warning(paste0(cellType, ": Only one proper module found, nothing to merge"))
            colors = labels2colors(cutree$labels) 
          }
          
          MEs <- try(moduleEigengenes_uv(expr = as.matrix(datExpr),
                                                colors=colors,
                                                excludeGrey = T))
          if (class(MEs)=="try-error") MEs <- NULL
          
        } else {
          warning(paste0(cellType, ": No modules found"))
        }
        out <- list("colors" = colors, "MEs"= MEs)
        out
      },
      cutree = list_cutree,
      comb = list_comb,
      SIMPLIFY=F)}
    
    list_iterable = list(list_comb=list_list_comb,
                list_cutree = list_list_cutree,
                datExpr = list_datExpr_gg,
                cellType = names(list_datExpr_gg))
    
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_mergeCloseModules.txt")
    
    list_list_merged <- safeParallel(fun=fun, 
                                     list_iterable=list_iterable, 
                                     outfile=outfile#, 
                                     #corFnc = corFnc,
                                     #corOptions=corOptions
                                     )
    
    names(list_list_merged) = sNames_3
    
    # get colors
    fun = function(list_merged,datExpr,cellType) {
      list_colors = lapply(list_merged, function(merged) {
        colors <- merged$colors
        names(colors) = colnames(datExpr)
        colors
      }) # list of merged colors
    } 
    
    list_iterable=list(list_merged=list_list_merged, 
              datExpr=list_datExpr_gg, 
              cellType=sNames_3)
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "getColors.txt")
    
    # Extract the colors from the list returned by mergeCloseModules
    list_list_colors <- safeParallel(fun=fun, 
                                     list_iterable=list_iterable, 
                                     outfile=outfile)
    names(list_list_colors) <- sNames_3
    
    # Extract the Module Eigengenes from the list returned by mergeCloseModules
    fun = function(list_merged) {lapply(list_merged, function(merged) {
      if (!is.null(merged$MEs)) merged$MEs$eigengenes else NULL})}
    list_iterable = list("X"=list_list_merged)
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "getMEs.txt")
    list_list_MEs <- safeParallel(fun=fun, 
                                  list_iterable=list_iterable, 
                                  outfile=outfile)
    names(list_list_MEs) <- sNames_3
    
    # Compute kMEs 
    fun = function(list_MEs, datExpr) {
      lapply(list_MEs, function(MEs) {
        kMEs <- NULL
        if (!is.null(MEs)) {
          kMEs <- signedKME(datExpr=as.matrix(datExpr),
                            datME=MEs,
                            outputColumnName="",
                            corFnc = corFnc )
          kMEs
        } else NULL
      })}
    list_iterable = list(list_MEs=list_list_MEs,  
                datExpr=list_datExpr_gg)
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "compute_kMEs.txt")
    
    list_list_kMs <- safeParallel(fun=fun, 
                                  list_iterable=list_iterable, 
                                  outfile=outfile)
    
  } else if (fuzzyModMembership=="kIM") { 
    
    fun = function(x) lapply(x, function(y) labels2colors(y$labels))
    list_iterable = list("X"=list_list_cutree)
    list_list_colors <- safeParallel(fun=fun, list_iterable=list_iterable)

    fun =  function(list_colors,datExpr) lapply(list_colors, function(colors) name_for_vec(to_be_named=colors, given_names = colnames(datExpr), dimension=NULL))
    list_iterable = list("list_colors"=list_list_colors, "datExpr"=list_datExpr_gg)
    list_list_colors <- safeParallel(fun=fun, list_iterable=list_iterable)

    
    fun = function(list_colors,list_plot_label) name_for_vec(to_be_named=list_colors, given_names = as.character(list_plot_label), dimension=NULL)
    list_iterable = list(list_colors = list_list_colors,
                list_plot_label = list_list_plot_label)
    list_list_colors <- safeParallel(fun=fun, list_iterable=list_iterable)

    
    names(list_list_colors) <- sNames_3
    
    # get IM embeddings in order to merge close modules using embeddings correlation
    fun = function(dissTOM, list_colors) lapply(list_colors, function(colors) kIM_eachMod_norm(dissTOM = dissTOM, 
                                                                                               colors = colors,
                                                                                               verbose=verbose,
                                                                                               excludeGrey=F,
                                                                                               do.par=F))
    list_iterable = list(dissTOM = list_dissTOM,
                list_colors = list_list_colors)
    list_list_kMs <- safeParallel(fun = fun, list_iterable=list_iterable)
    
    
    fun = function(list_colors,list_kMs,datExpr,dissTOM, cellType) {
      mapply(function(colors,kMs) {
        mergeCloseModskIM(datExpr=datExpr,
                          colors = colors,
                          kIMs  = kMs,
                          dissTOM = dissTOM,
                          moduleMergeCutHeight=moduleMergeCutHeight,
                          verbose=verbose,
                          cellType = cellType)},
        colors = list_colors,
        kMs = list_kMs,
        SIMPLIFY=F)}

    list_iterable = list(list_colors = list_list_colors,
                list_kMs = list_list_kMs,
                datExpr = list_datExpr_gg,
                dissTOM = list_dissTOM,
                cellType = names(list_list_colors))
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_mergeCloseModskIM.txt")
    list_list_merged <- safeParallel(fun=fun, list_iterable=list_iterable, outfile = outfile)
    
   
    fun = function(cellType) {tryCatch({lapply(list_list_merged[[cellType]], function(y) y$colors)}, 
                                       error = function(c) {
                                         warning(paste0(cellType, " failed"))
                                       })}
    list_iterable = list("X"=names(list_list_merged))
    list_list_colors <- safeParallel(fun=fun, list_iterable=list_iterable)

    fun = function(x) lapply(x, function(y) y[['colors']])
    list_iterable= list("X"=list_list_merged)
    list_list_colors <- safeParallel(fun=fun, list_iterable=list_iterable)

    fun = function(x) lapply(x, function(y) y[['kIMs']])
    list_iterable= list("X"=list_list_merged)
    list_list_kMs <- safeParallel(fun=fun, list_iterable=list_iterable)
    
    list_list_MEs <- NULL
    
    
  }
  names(list_list_colors) <- names(list_list_kMs) <- sNames_3
  
  ######################################################################
  ###################### REASSIGN GENES BASED ON kM ####################
  ######################################################################
  
  # see https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0420-6
  
  if (kMReassign) {
    message(paste0("Reassigning genes based on ", fuzzyModMembership))
    
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_par_kMReassign.txt")
    
    if (fuzzyModMembership=="kIM") {
      fun = function(list_colors,
                     dissTOM,
                     datExpr,
                     cellType) lapply(list_colors, function(cols) fnc_kMReassign(modcolors = cols,
                                                                                    fuzzyModMembership = fuzzyModMembership,
                                                                                    dissTOM = if (fuzzyModMembership=="kIM") dissTOM else NULL,
                                                                                    datExpr = if (fuzzyModMembership=="kME") datExpr else NULL,
                                                                                    corFnc = if (fuzzyModMembership=="kME") corFnc else NULL,
                                                                                    verbose=verbose,
                                                                                    max_iter = 3,
                                                                                    cellType = cellType)) 
      list_iterable = list(list_colors = list_list_colors,
                  dissTOM = list_dissTOM,
                  datExpr = list_datExpr_gg,
                  cellType = names(list_list_colors))
      
    } else if (fuzzyModMembership=="kME"){ # to avoid storing list_dissTOM
      fun = function(list_colors,
                     datExpr,
                     cellType) lapply(list_colors, function(cols) fnc_kMReassign(modcolors = cols,
                                                                                    fuzzyModMembership = fuzzyModMembership,
                                                                                    #dissTOM = if (fuzzyModMembership=="kIM") dissTOM else NULL,
                                                                                    datExpr = if (fuzzyModMembership=="kME") datExpr else NULL,
                                                                                    corFnc = if (fuzzyModMembership=="kME") corFnc else NULL,
                                                                                    verbose=verbose,
                                                                                    max_iter = 3,
                                                                                    cellType = cellType))
      list_iterable = list(list_colors = list_list_colors,
                  datExpr = list_datExpr_gg,
                  cellType = names(list_list_colors))
      
    }
    
    
    list_list_reassign = safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile#, 
                                      #fuzzyModMembership=fuzzyModMembership, corFnc=corFnc, verbose=verbose
                                      )
    
    
    # Extract new colors and kMs from the list of lists of lists returned by the vectorised kMReassign 
    list_list_colors_reassign <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$colors))
    list_list_kMs_reassign <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$kMs))
    list_list_reassign_log <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$log))
    
    #invisible(gc()); invisible(R.utils::gcDLLs())
  } else {
    list_list_colors_reassign <- list_list_colors
    list_list_kMs_reassign <- list_list_kMs
    list_list_reassign_log <- NULL
  }  
  
  ######################################################################
  #################### EXTRACT PRIMARY kMEs / kIMs #####################
  ######################################################################
  # Extract the primary kMs - i.e. those of each gene w.r.t. the module to which it belongs
  # use these for only submitting relatively 'core' genes for PPI validation
  
  message(paste0("Computing primary "), fuzzyModMembership, "s")
  
  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_par_primary_kMs_for_PPI.txt")
  fun = function(list_kMs,
                 list_colors) mapply(function(kMs, colors) {
                   if (length(unique(colors))>1) {
                     pkMs_fnc(kMs=kMs, colors=colors)  
                   } else {
                     pkMs<-numeric(length=length(colors))
                     names(pkMs) <- names(colors)
                     pkMs}
                 }, 
                 kMs = list_kMs, 
                 colors=list_colors, SIMPLIFY=F)
  
  list_iterable=list(list_kMs = list_list_kMs_reassign, 
            list_colors = list_list_colors_reassign)
  
  list_list_pkMs <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
  
  # The pkM vectors already have gene names. Now name each vector, and each list of vectors
  fun = function(list_pkMs,list_plot_label) name_for_vec(to_be_named = list_pkMs, 
                                                         given_names = as.character(list_plot_label), 
                                                         dimension = NULL)
  list_iterable = list(list_pkMs=list_list_pkMs, 
              list_plot_label=list_list_plot_label)
  list_list_pkMs <- safeParallel(fun=fun, list_iterable=list_iterable)
  
  names(list_list_pkMs) <- sNames_3

  ######################################################################
  ##### FILTER OUT GENES WITH INSIGNIFICANT GENE-MOD CORRELATION #######
  ######################################################################
  # Module expression as defined by ME - Module Eigengene or IM - IntraModular connectivity embedding
  
  if (kMSignifFilter) {
    
    message("Computing gene-module correlation t-test p-values")
    
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_par_t.test.txt")
    
    if (fuzzyModMembership=="kME") {
      
      fun = function(list_pkMs, list_colors, list_comb, cellType)  {
        mapply(function(pkMs, colors, comb) {
          tryCatch({
            if (length(unique(colors))>1){
              vec_t <- ifelse(colors!="grey", (pkMs*sqrt(length(pkMs)-2))/sqrt(1-pkMs^2), 0) # compute t.stats, for a vector of rho
              vec_p.val <- ifelse(colors!="grey", stats::pt(q=vec_t,  
                                                            lower.tail = F, 
                                                            df = length(pkMs)-2), 1) # compute p.values. We are genuinely only interested in upper tail 
              vec_q.val <- p.adjust(p = vec_p.val, method = "fdr")
              vec_idx_signif <- vec_q.val < pvalThreshold # compute boolean significance vector. Set "grey" to TRUE so we don't count them
              vec_idx_signif[colors=="grey"] <- TRUE
              out <- data.frame("p.val" = vec_p.val, "q.val" = vec_q.val, "signif" = vec_idx_signif)
              out } 
            else {
              out <- NULL
            }
          }, 
          error = function(err) {
            message(paste0(cellType, ", ", comb, ": ", "Failed to compute gene p-values with error: ", err))  
            out <- data.frame("p.val" = rep(0, length(colors)), "q.val" = rep(0, length(colors)), "signif" = rep(TRUE, length(colors)))
            out
          })
        },
        pkMs=list_pkMs, 
        colors=list_colors, 
        comb = list_comb,
        SIMPLIFY=F)}
      
      list_iterable = list(list_pkMs = list_list_pkMs, 
                  list_colors = list_list_colors_reassign, 
                  list_comb = list_list_comb,
                  cellType = sNames_3)
      
      list_list_geneMod_t.test <- safeParallel(fun=fun, 
                                               list_iterable=list_iterable#, 
                                               #pvalThreshold=pvalThreshold
                                               ) 
      
      
    } else if (fuzzyModMembership=="kIM") {
      
      fun = function(datExpr, list_colors, cellType, list_kMs, dissTOM, list_comb) { 
        mapply(function(colors, kMs, comb) {
          tryCatch({
            if (length(unique(colors))>1) {
              message(paste0(cellType, ": Computing module kIM embeddings"))
              cellModEmbed(datExpr=datExpr,
                           modcolors=colors,
                           latentGeneType="IM",
                           cellType=NULL, # to avoid having it in the embedding matrix column names
                           kMs = kMs,
                           dissTOM = dissTOM)
            } else {
              NULL
            }
          }, error = function(err) {
            message(paste0(cellType, ", ", comb, ": ", "failed to compute embedding matrix with error: ", err))
            NULL
          })
        },
        colors=list_colors,
        kMs = list_kMs,
        comb= list_comb,
        SIMPLIFY=F)
      }
      
      list_iterable  = list(datExpr = list_datExpr_gg,
                   list_colors = list_list_colors_reassign,
                   cellType = sNames_3,
                   list_kMs = list_list_kMs_reassign,
                   list_comb=list_list_comb,
                   dissTOM = list_dissTOM)
      
      list_list_embed_mat <- safeParallel(fun=fun, list_iterable=list_iterable)
      
      
      fun = function(datExpr, list_embed_mat, list_colors, list_comb, cellType) {
        mapply(function(embed_mat, colors, comb){
          tryCatch({
            if (length(unique(colors))>1) {
              message(paste0(cellType, ": Computing gene correlations with module kIM embeddings"))
              vec_p.val <- numeric(length=length(colors))
              names(vec_p.val) <- names(colors)
              vec_p.val[colors=="grey"] <- 0
              
              for (color in names(table(colors))[!grepl("^grey$",names(table(colors)))]) {
                vec_p.val[colors==color] <- apply(X = datExpr[,colors==color], MARGIN = 2, FUN = function(datExpr_col) cor.test(datExpr_col, embed_mat[, color,drop=F], 
                                                                                                                                alternative = "greater", 
                                                                                                                                method = "pearson", 
                                                                                                                                conf.level = 1-pvalThreshold)$p.value)
              }
              vec_q.val <- p.adjust(p = vec_p.val, method = "fdr")
              vec_idx_signif <- vec_q.val < pvalThreshold
              vec_idx_signif[colors=="grey"] <- TRUE# compute boolean significance vector. Set "grey" to TRUE so we don't count them
              out <- data.frame("p.val" = vec_p.val, "q.val" = vec_q.val, "signif" = vec_idx_signif)
              out } else out <- NULL
          }, error= function(err) {
            message(paste0(cellType, ", ", comb, ": ", "failed to compute t.test with error: ", err))
            out <- data.frame("p.val" = rep(0, length(colors)), "q.val" = rep(0, length(colors)), "signif" = rep(TRUE, length(colors)))
            out
          })
        },
        embed_mat=list_embed_mat,
        colors=list_colors,
        comb = list_comb,
        SIMPLIFY=F)}
      
      list_iterable = list(datExpr = list_datExpr_gg,
                  list_embed_mat = list_list_embed_mat,
                  list_colors = list_list_colors_reassign,
                  list_comb = list_list_comb,
                  cellType = sNames_3)
      
      list_list_geneMod_t.test <- safeParallel(fun=fun, 
                                               list_iterable=list_iterable, 
                                               outfile=outfile#, 
                                               #pvalThreshold=pvalThreshold
                                               )
      
      
      rm(list_list_embed_mat)
    } 
    #stopCluster(cl)
  } else {
    list_list_geneMod_t.test <- NULL 
  }  
  
  if (fuzzyModMembership=="kIM") rm(list_dissTOM)
  
  ########################################################################################
  #################### REASSIGN NON-SIGNIFICANT GENES TO "GREY" ##########################
  ########################################################################################
  
  list_list_colors_t.test <- list_list_colors_reassign

  if(kMSignifFilter) if (!all(sapply(list_list_geneMod_t.test, 
                                       function(list_geneMod_t.test) {
                                         sapply(list_geneMod_t.test, function(geneMod_t.test) {
                                           if (!is.null(geneMod_t.test)) all(geneMod_t.test$signif) else T}
                                           , simplify=T)}, simplify=T))) {
    message("Filtering out genes with insignificant gene-module correlation p-values in expression profile")

    # filter colors
    fun = function(list_colors_reassign, list_geneMod_t.test) {
      mapply(function(colors_reassign, geneMod_t.test) {
        
        if (!is.null(geneMod_t.test)) {
          colors_t.test <- colors_reassign
          colors_t.test[!geneMod_t.test$signif] <- rep("grey", times=sum(!geneMod_t.test$signif))
          return(colors_t.test)
        } else {
          colors_reassign
        }
      }, colors_reassign = list_colors_reassign, 
      geneMod_t.test = list_geneMod_t.test, 
      SIMPLIFY=F)
    }
    list_iterable = list(list_colors_reassign = list_list_colors_reassign, 
                list_geneMod_t.test = list_list_geneMod_t.test)
    list_list_colors_t.test <- safeParallel(fun=fun, list_iterable=list_iterable)

  } 
  
  ######################################################################
  ################### Filter out very small modules ####################
  ######################################################################
  
  if(kMSignifFilter) if (!all(sapply(list_list_geneMod_t.test, 
                                       function(list_geneMod_t.test) {
                                         sapply(list_geneMod_t.test, function(geneMod_t.test) {
                                           if (!is.null(geneMod_t.test)) all(geneMod_t.test$signif) else T}
                                           , simplify=T)}, simplify=T))) {
    
    message(paste0("Filtering out modules with fewer than ", minClusterSize, " genes left after t.testing"))
    

    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_par_minClusterSize_filter.txt")
    
    # filter colors last as we needed colors to filter pkMs
    fun = function(list_colors) {
      lapply(list_colors, function(colors) {
        mods_too_small <- names(table(colors))[table(colors) < minClusterSize]
        colors[colors %in% mods_too_small] <- rep("grey", times=sum(colors %in% mods_too_small))
        colors
      })
    }
    list_iterable= list("X"=list_list_colors_t.test)
    list_list_colors_t.test <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
   
  } 
  
  ######################################################################
  ######################### REMOVE ALL-GREY RUNS #######################
  ######################################################################
  
  message("Removing all-grey runs")
  
  # count the grey
  list_vec_n_grey <- lapply(X=list_list_colors_t.test, 
                            FUN = function(list_colors) sapply(X = list_colors, 
                                                               FUN=function(colors) sum(as.numeric(colors=="grey")), simplify = T))
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  list_logical_params_ok <- mapply(function(vec_n_grey,list_colors_t.test) {
    as.logical(mapply(function(n_grey,colors_t.test) {
      n_grey!=length(colors_t.test)
    }, 
    n_grey=vec_n_grey, 
    colors_t.test=list_colors_t.test, SIMPLIFY=F))
  }, 
  vec_n_grey=list_vec_n_grey, 
  list_colors_t.test=list_list_colors_t.test, 
  SIMPLIFY=F)
  
  vec_logicalCellClustOK <- sapply(list_logical_params_ok, any, simplify = T) %>% unlist
  
  # First, for each run, remove results of any parametrisation for which all the genes were assigned to grey
  # If for a subset all parameters gave only grey modules, take it out of the top level list.
  list_list_plot_label_ok <- mapply(function(list_plot_label,logical_params_ok) list_plot_label[logical_params_ok], 
                                    list_plot_label = list_list_plot_label, 
                                    logical_params_ok = list_logical_params_ok, 
                                    SIMPLIFY = F)
  
  list_list_cutree_ok <- mapply(function(list_cutree,logical_params_ok) list_cutree[logical_params_ok], 
                                list_cutree = list_list_cutree,
                                logical_params_ok = list_logical_params_ok,
                                SIMPLIFY = F)
  
  list_list_colors_t.test_ok <- mapply(function(list_colors_t.test,logical_params_ok) {
    list_colors_t.test[logical_params_ok]
  }, 
  list_colors_t.test=list_list_colors_t.test, 
  logical_params_ok=list_logical_params_ok, 
  SIMPLIFY = F)
  

  list_list_reassign_log_ok <- if (kMReassign) {mapply(function(x,y) x[y], x=list_list_reassign_log, y=list_logical_params_ok, SIMPLIFY = F)  
    } else NULL

  # Make a warning and report if any whole subsets produced no modules
  if (!all(vec_logicalCellClustOK)) {
    log_entry <- paste0(sNames_3[vec_logicalCellClustOK], " had no modules, dropped from the analysis")
    warning(log_entry)
    cat(log_entry, file = pathLogCellClustersDropped, append=T, sep = "\n")
  } 
  
  # filter
  list_list_plot_label_ok <- list_list_plot_label_ok[vec_logicalCellClustOK]
  list_list_cutree_ok <- list_list_cutree_ok[vec_logicalCellClustOK]
  list_datExpr_gg <- list_datExpr_gg[vec_logicalCellClustOK] # If for a subset all parameters gave only grey modules, take it out of the top level list.
  list_geneTree <- list_geneTree[vec_logicalCellClustOK]
  list_list_colors_t.test_ok <- list_list_colors_t.test_ok[vec_logicalCellClustOK]
  list_list_reassign_log_ok <- if (kMReassign) list_list_reassign_log_ok[vec_logicalCellClustOK] else NULL 
  
  # Get new subset names (after filtering)
  sNames_4 <- sNames_3[vec_logicalCellClustOK]
  
  # Assign gene names to each color vector
  list_list_colors_t.test_ok <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), 
                                                                                         dimension = NULL)), 
                                        x = list_list_colors_t.test_ok, 
                                        y = list_datExpr_gg, SIMPLIFY=F)
  
  rm(list_list_plot_label, list_list_cutree, list_list_colors_t.test, list_list_reassign_log, list_list_MEs)#, list_idx_mapped_gg)
  
  invisible(gc())
  
  ######################################################################
  ############################# CHECKPOINT #############################
  ######################################################################
  
  resume = "checkpoint_3"
  message("Reached checkpoint 3, saving session image")
  if (autosave) save.image( file=sprintf("%s%s_%s_checkpoint_3_image.RData.gz", dirTmp, prefixData, prefixRun), compress = "gzip")
} 

if (resume == "checkpoint_3") {
  
  ######################################################################
  ######################### LOAD PACKAGES ##############################
  ######################################################################
  
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("WGCNA"))
  suppressPackageStartupMessages(library("liger"))
  suppressPackageStartupMessages(library("boot"))

  ######################################################################
  ################ ORDER PARAMETER SETS BY NUMBER OF MODS ##############
  ######################################################################
  
  message("Selecting parameters with the highest number of modules")
  
  # Count how many modules were found under each parametrisation and order the parametrisations
  list_vec_n_mod <- lapply(list_list_colors_t.test_ok, function(list_colors) sapply(list_colors, function(vec_colors) length(unique(vec_colors)), simplify=T)) 
    
  # Order all the outputs by how many genes were assigned to a (non-grey) module
  if (kMReassign) {
    list_list_reassign_log_order <- mapply(function(list_reassign_log_ok,vec_n_mod) list_reassign_log_ok[order(vec_n_mod, decreasing=T)], 
                                           list_reassign_log_ok = list_list_reassign_log_ok, 
                                           vec_n_mod = list_vec_n_mod, SIMPLIFY=F) 
  } else {
    list_list_reassign_log_order <- NULL
  }
  
  list_list_plot_label_ok_order <- mapply(function(x,y) x[order(y, decreasing=T)], x = list_list_plot_label_ok, y = list_vec_n_mod, SIMPLIFY=F)
  list_list_colors_t.test_ok_order <- mapply(function(x,y) x[order(y, decreasing=T)], x  =  list_list_colors_t.test_ok , y = list_vec_n_mod, SIMPLIFY = F )
  list_list_geneMod_t.test_order <- if (kMSignifFilter) mapply(function(x,y) x[order(y, decreasing=T)], x = list_list_geneMod_t.test[names(list_list_geneMod_t.test) %in% sNames_4], y = list_vec_n_mod, SIMPLIFY=F) else NULL
  
  ######################################################################
  ######### FOR EACH CELLTYPE SELECT PARAMETERS WITH MOST MODULES  #####
  ######################################################################
  
  # Eliminate a layer of nesting by selecting only the 'best' parametrisation per celltype
  list_plot_label <- lapply(list_list_plot_label_ok_order, function(x) x[[1]])
  list_reassign_log <- if (kMReassign) lapply(list_list_reassign_log_order, function(x) x[[1]]) else NULL
  list_colors <- lapply(list_list_colors_t.test_ok_order, function(x) x[[1]])
  list_geneMod_t.test <- if (kMSignifFilter)  lapply(list_list_geneMod_t.test_order, function(x) x[[1]]) else NULL 

  # Name by cell clusters
  names(list_plot_label)  <-  names(list_colors) <- sNames_4

  if (kMReassign) names(list_reassign_log) <- sNames_4
  if (kMSignifFilter)  names(list_geneMod_t.test) <- sNames_4
  
  # Make list of list of final parameters
  param_names = c("minClusterSize", "deepSplit","pamStage", "moduleMergeCutHeight")
  list_list_cutree_params <- lapply(list_vec_n_mod, function(x) list_comb[order(x,decreasing = T)][[1]])
  list_list_cutree_params <- lapply(list_list_cutree_params, function(x) name_for_vec(to_be_named = x, given_names = param_names, dimension = NULL)) 
  
  ######################################################################
  ##################### MAKE MODULE COLORS UNIQUE ######################
  ######################################################################
  
  message("Making module colors unique across cell clusters")
  
  # Get nested list of modules
  list_mods <- lapply(list_colors, function(x) names(table(x)))
  
  list_mods <- lapply(list_mods, function(mods) {
    if (any(grepl("^grey$", mods))) {
      mods <- mods[mods!="grey"] 
    } else {
      mods
    }
  })
  
  mods <- unlist(list_mods)
  names(mods) <- NULL
  
  all_cols_nogrey_uniq <- unique(gsub("\\d+", "", colors()[-grep("grey|gray", colors())]))
  all_cols_nogrey <- colors()[-grep("grey|gray", colors())]
  
  # Replace colors with new unique colors
  if (length(mods) <= length(all_cols_nogrey_uniq) ) { # if there are enough unique colors without adding numbers
    mods_uniq <- all_cols_nogrey_uniq[sample(x=1:length(all_cols_nogrey_uniq), size=length(mods), replace=F)]
  } else if (length(mods) > length(all_cols_nogrey_uniq) & length(mods) < length(all_cols_nogrey) ) { # if there aren't enough unique colors unless they have numbers added
    mods_uniq <- all_cols_nogrey[sample(x=1:length(all_cols_nogrey), size=length(mods), replace=F)]
  } else if (length(mods) > length(all_cols_nogrey)) { # if there aren't enough unique colors in R
    mods_uniq <- paste0(all_cols_nogrey_uniq, "_", 1:(length(mods)))
  }
  
  list_mods_uniq <- vector(mode="list", length=length(list_mods))
  
  k=1
  for (j in 1:length(list_mods)) {
    for (i in 1:length(list_mods[[j]])) {
      list_mods_uniq[[j]][i] <- mods_uniq[k]
      k = k+1
    }
  }
  
  names(list_mods_uniq) <- names(list_colors)
  
  # Update colors
  list_colors_uniq <- mapply(function(cols,mods,mods_uniq) {
    newCols <- mods_uniq[match(cols, mods)]
    names(newCols) <- names(cols)
    newCols
    },
   cols = list_colors, 
   mod = list_mods, 
   mods_uniq = list_mods_uniq, SIMPLIFY = F)

  # Replace NAs (i.e. no match in non-grey modules) with "grey"
  list_colors_uniq <- lapply(list_colors_uniq, function(x) ifelse(is.na(x), yes="grey", no = x))
  
  ######################################################################
  ############# RECOMPUTE MEs, kMs, pkMs AFTER PPI FILTER ##############
  ######################################################################

  message("Computing kMs and pkMs")

  if (fuzzyModMembership=="kIM"){
    list_dissTOM_path <- dir(path = dirTmp, pattern = paste0(prefixData, "_", prefixRun, "_list_dissTOM"), full.names = T)
    list_dissTOM <- load_obj(list_dissTOM_path)
    list_dissTOM <- list_dissTOM[match(sNames_4,names(list_dissTOM))]
  }

  if (fuzzyModMembership == "kME") {

    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_kMEs.txt")
    fun = function(x,y,cellType) {
      out <- try({
        moduleEigengenes_uv(expr = as.matrix(x),
                                   colors=y,
                                   excludeGrey=T)#,

      })
      if (class(out) == "try-error"){
        warning(paste0(cellType, ": moduleEigengenes failed with the error: ", out))
        out <- list(eigengenes=NULL, u = NULL)
      } else {
        out
      }
    }
    list_iterable = list(x = list_datExpr_gg,
                y = list_colors_uniq,
                cellType = names(list_datExpr_gg))
    
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_MEs.txt")
    
    list_ModuleEigengenes_out <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

    list_MEs <- lapply(list_ModuleEigengenes_out, function(x) x$eigengenes)

    list_list_u <- lapply(list_ModuleEigengenes_out, function(x) x$u) # 

    names(list_MEs) <- names(list_list_u) <- sNames_4

    fun = function(x,y) {
      if(!is.null(y)) {
        signedKME(datExpr = as.matrix(x),
                  datME = y,
                  outputColumnName = "",
                  corFnc = corFnc)
      } else {
        NULL
      }
    }
    list_iterable = list(x=list_datExpr_gg,
                y=list_MEs)
    list_kMs <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

    # Remove 'ME' from eigengenes
    list_MEs <- lapply(list_MEs, function(x) {
      if (!is.null(x)) name_for_vec(to_be_named = x, given_names = gsub(pattern="ME" , replacement="", x = colnames(x), ignore.case = F), dimension = 2) else NULL
    })

    # add cell row names
    list_MEs <- mapply(function(x,y) {
      if (!is.null(x)) name_for_vec(to_be_named = x, given_names = rownames(y), dimension=1) else NULL
    },
    x = list_MEs,
    y = list_datExpr_gg,
    SIMPLIFY=F)

  } else if (fuzzyModMembership == "kIM"){

    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "log_kIMs.txt")
    fun = function(x,y) kIM_eachMod_norm(dissTOM = x,
                                         colors = y,
                                         excludeGrey=T,
                                         do.par=F)
    list_iterable = list(x = list_dissTOM,
                y = list_colors_uniq)
    list_kMs <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)



    names(list_kMs) <- sNames_4

    list_MEs <- NULL

    # Output list of list of kM gene weights, one vector per module
    fun = function(kMs, colors) {

      list_u <- lapply(colnames(kMs), function(module) {
        out <- kMs[match(names(colors)[colors==module], rownames(kMs)),module]
        names(out) = rownames(kMs)[colors==module]
        return(out)
      })
      names(list_u) <- colnames(kMs)
      return(list_u)
    }
    list_iterable = list(kMs = list_kMs,
                colors = list_colors_uniq)
    list_list_u <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

  }

  if (fuzzyModMembership == "kIM") rm(list_dissTOM)

  #####################################################################
  ########################### COMPUTE PRIMARY KMs #####################
  #####################################################################
  # Extract the primary kMs - i.e. those of each gene w.r.t. the module to which it belongs

  message(paste0("Computing primary "), fuzzyModMembership, "s")

  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "list_pkMs_PPI.txt")
  fun = function(kMs, colors) pkMs_fnc(kMs=kMs, colors=colors)
  list_iterable = list("kMs" = list_kMs, "colors"=list_colors_uniq)
  list_pkMs <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

  names(list_pkMs) <- sNames_4
  
  ######################################################################
  ##################### COMPUTE CELL x EIGENGENE MATRIX ################
  ######################################################################

  message("Computing all cell embeddings on all modules, across celltypes")
  
  path_datExpr <- dir(path = dirTmp, pattern = paste0(prefixData, "_", prefixRun, "_datExprNormFull.RDS.gz"), full.names = T)
  load_obj(path_datExpr[1]) %>% t -> datExpr 
  
  invisible(gc())
  
  latentGeneType <- if (fuzzyModMembership == "kME") "ME" else if (fuzzyModMembership=="kIM") "IM"
  
  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "make_cell_embed_mat.txt")
  fun = function(vec_colors,name,kMs) cellModEmbed(datExpr = datExpr, 
                                         modcolors = vec_colors,
                                         latentGeneType = latentGeneType,
                                         cellType = name,
                                         kMs = if (latentGeneType== "IM") kMs else NULL)#,

  list_iterable = list(vec_colors = list_colors_uniq, 
              name = names(list_colors_uniq),
              kMs = list_kMs)#,

  list_cellModEmbed_mat <- safeParallel(fun=fun, 
                                        list_iterable=list_iterable, 
                                        outfile=outfile)#, 
                                        #datExpr = datExpr, 
                                        #latentGeneType = latentGeneType)
  
  
  list_cellModEmbed_mat %>% Reduce(function(mat1, mat2) cbind(mat1, mat2), .) -> cellModEmbed_mat
  
  ### Z-score the embedding matrix 
  # cellModEmbed_mat <- apply(X=cellModEmbed_mat, 
  #                           MARGIN = 2, 
  #                           FUN = function(module_vec) { 
  #                             module_vec %>% '-'(mean(.)) %>% '/'(stats::sd(.)) 
  #                             })
  
  # TODO: does it make sense to Z-score it? Probably not..
  
  rownames(cellModEmbed_mat) <- NULL
  
  df_annot <- if (!is.null(mat_addAnnot)) { 
    data.frame("data"=rep(x = prefixData,  
                          times=nrow(cellModEmbed_mat)), 
                 mat_addAnnot,
                "cell_cluster"= vec_idents, 
                "cell_id" = rownames(datExpr), row.names = NULL) 
  } else {
    data.frame("data"=rep(x = prefixData,  
                          times=nrow(cellModEmbed_mat)), 
               "cell_cluster"= vec_idents, 
               "cell_id" = rownames(datExpr), row.names = NULL) 
  }
  
  cellModEmbed_mat_annot <- cbind(df_annot, cellModEmbed_mat)
                                 
  rm(datExpr)
  
  ######################################################################
  ############################# CHECKPOINT #############################
  ######################################################################
  
  resume = "checkpoint_4"
  message("Reached checkpoint 4, saving session image")
  if (autosave) save.image( file=sprintf("%s%s_%s_checkpoint_4_image.RData.gz", dirTmp, prefixData, prefixRun), compress = "gzip")
} 

if (resume == "checkpoint_4") {
  
  ##########################################################################
  ############################ (UN)LOAD PACKAGES ############################
  ##########################################################################
  
  invisible(R.utils::gcDLLs())
  
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("Biobase"))
  suppressPackageStartupMessages(library("Matrix"))
  suppressPackageStartupMessages(library("parallel"))
  suppressPackageStartupMessages(library("reshape"))
  suppressPackageStartupMessages(library("reshape2"))
  suppressPackageStartupMessages(library("readr"))
  
  #pkgs <- c("dplyr", "Biobase", "Matrix", "parallel", "reshape", "reshape2", "readr")
  
  #ipak(pkgs)
  
  ##########################################################################
  ######### PREPARE GENES LISTS AND DATAFRAME WITH MODULES, GENES ##########
  ##########################################################################

  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "write_outputs.txt")
  
  message("Preparing outputs")
  
  # prepare nested lists of module genes
  fun = function(a,b) lapply(b, function(x) names(a)[a==x])
  list_iterable = list(a=list_colors_uniq, b=list_mods_uniq)
  list_list_module_genes = safeParallel(fun=fun,list_iterable=list_iterable, outfile=outfile)
  
  fun = function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = NULL)
  list_iterable=  list(x=list_list_module_genes, y=list_mods_uniq)
  list_list_module_genes <- safeParallel(fun=fun,list_iterable=list_iterable,outfile=outfile)

  # make a copy for pkMs
  list_list_module_pkMs <- list_list_module_genes # just a template
  list_list_module_gene_loadings <- if (fuzzyModMembership == "kME") list_list_module_genes else NULL # just a template
  # order the gene lists by pkM and get pkM value

  # iterate over celltypes
  for (celltype in names(list_list_module_genes)) {
    # iterate over modules
    for (module in names(list_list_module_genes[[celltype]])) {
      pkMs = list_pkMs[[celltype]]
      if (fuzzyModMembership == "kME") mod_u <- list_list_u[[celltype]][[module]]
      mod_genes <- list_list_module_genes[[celltype]][[module]]
      mod_pkMs_sorted <- sort(pkMs[match(mod_genes,names(pkMs)), drop=F], decreasing=T)
      list_list_module_genes[[celltype]][[module]] <- names(mod_pkMs_sorted)  # sort the genes
      list_list_module_pkMs[[celltype]][[module]] <- mod_pkMs_sorted
      if (fuzzyModMembership == "kME") list_list_module_gene_loadings[[celltype]][[module]] <- mod_u[match(names(mod_pkMs_sorted), names(mod_u)), drop=F]
    }
  }
  
  message("Preparing module genes dataframe")
  
  cell_cluster <- rep(sNames_4, times=unlist(lapply(list_list_module_genes, FUN=function(x) sum(sapply(x, function(y) length(y), simplify=T)))))
  module <- unlist(lapply(list_list_module_genes, function(x) rep(names(x), sapply(x, function(y) length(y), simplify = T))), use.names = F)
  #ensembl <- if (mapToEnsembl) unlist(list_list_module_genes, recursive = T, use.names = F) else NULL 
  genes <- unlist(list_list_module_genes, recursive = T, use.names = F) 
  pkMs <- unlist(list_list_module_pkMs, recursive = T, use.names = F)
  if (fuzzyModMembership == "kME") gene_loadings <- unlist(list_list_module_gene_loadings, recursive = T, use.names = F)
  data <- rep(prefixData, times= length(pkMs))
  run <- rep(prefixRun, times= length(pkMs))    
  
  df_cell_cluster_module_genes <- data.frame(data, run, cell_cluster, module, genes, pkMs, row.names = NULL)

  if (fuzzyModMembership=="kME") df_cell_cluster_module_genes[["gene_loadings"]] <- gene_loadings
  
  ##########################################################################
  ############################# OUTPUT TABLES ##############################
  ##########################################################################
  
  message("Writing outputs to disk")
  
  ##################### WRITE OUT TABLES OF MODULE GENES ####################
  ###########################################################################
  
  write.csv(df_cell_cluster_module_genes, file = gzfile(sprintf("%s%s_%s_cell_cluster_module_genes.csv.gz",dirTables, prefixData, prefixRun)), quote = F, row.names = F)
  
  ################## SAVE KMs FOR ENRICHED MODULES #########################
  ##########################################################################
  
  # Prepare dfs with a gene column followed by kMEs / kIMs 
  list_kMs_out <- lapply(list_kMs, function(kMs) {
    kMs<-cbind(genes=as.character(rownames(kMs)), kMs)
    rownames(kMs) <- NULL
    kMs
    })

  # save kMs as a single outer join 
  kMs_out <- Reduce(f = function(df1, df2) {
    dplyr::full_join(df1, df2, by="genes")}, 
    x=list_kMs_out)
  
  write.csv(kMs_out, file=gzfile(sprintf("%s%s_%s_kMs_full_join.csv.gz", dirTables, prefixData, prefixRun)), 
            row.names=F, quote = F)
  

  ################## OUTPUT CELL MODULE EMBEDDINGS MATRIX ###################
  ###########################################################################
  
  invisible(write.csv(x=cellModEmbed_mat_annot, 
                      file=gzfile(sprintf("%s%s_%s_%s_cellModEmbed.csv.gz", dirTables, prefixData, prefixRun, fuzzyModMembership)), row.names=F, quote = F))
  
  ##################### WRITE PARAMETERS AND STATS TO FILES ################
  ##########################################################################
  
  t_finish <- as.character(Sys.time())
  
  sumstats_celltype_df <- matrix(data = NA_real_, nrow=length(sNames_0), ncol=13) %>% data.frame 
  colnames(sumstats_celltype_df) = c("prefixRun", 
                                     "subset", 
                                     "n_cells", 
                                     "n_genes", 
                                     "n_genes_used", 
                                     "softPower", 
                                     "SFT.R.sq", 
                                     "median.k.", 
                                     "plot_label_final", 
                                     "n_genes_reassign", 
                                     "prop_genes_t.test_fail", 
                                     "prop_genes_assign",
                                     "n_modules")#,

  
  sumstats_celltype_df[["prefixRun"]] = rep(prefixRun, times=length(sNames_0))
  sumstats_celltype_df[["subset"]] = sNames_0
  sumstats_celltype_df[["n_cells"]] = n_cells_subsets
  sumstats_celltype_df[["n_genes"]] = n_genes_subsets
  sumstats_celltype_df[["n_genes_used"]] = n_genes_subsets_used
  sumstats_celltype_df[["softPower"]][sNames_0 %in% sNames_1] = sapply(list_sft, function(x) x$Power)
  sumstats_celltype_df[["SFT.R.sq"]][sNames_0 %in% sNames_1] = sapply(list_sft, function(x) x$SFT.R.sq)
  sumstats_celltype_df[["median.k."]][sNames_0 %in% sNames_1] = sapply(list_sft, function(x) x$median.k.)
  sumstats_celltype_df[["plot_label_final"]][sNames_0 %in% sNames_4 ] = unlist(x=list_plot_label, use.names = F)
  if (kMReassign) sumstats_celltype_df[["n_genes_reassign"]][sNames_0 %in% sNames_4] = sapply(list_reassign_log, function(rlog) if(!is.null(rlog)) nrow(rlog) else 0, simplify=T)
  if (kMSignifFilter) sumstats_celltype_df[["prop_genes_t.test_fail"]][sNames_0 %in% sNames_4] = sapply(list_geneMod_t.test, function(x) {if (!is.null(x)) {sum(as.numeric(!x$signif))/length(x$signif)} else {NA_real_}}, simplify =T)
  sumstats_celltype_df[["prop_genes_assign"]][sNames_0 %in% sNames_4] = sapply(list_colors_uniq, function(x) round(sum(x!="grey")/length(x),2), simplify=T)
  sumstats_celltype_df[["n_modules"]][sNames_0 %in% sNames_4] = sapply(list_colors_uniq, function(x) length(unique(as.character(x)))-1, simplify=T)

  
  sumstats_run <- c(prefixRun = prefixRun,
                    tStart = tStart,
                    t_finish = t_finish,
                    mean_percent_mito = if (!is.null(riboMitoMeta)) mean(riboMitoMeta$percent.mito, na.rm=T) else NA,
                    mean_percent_ribo =  if (!is.null(riboMitoMeta)) mean(riboMitoMeta$percent.ribo, na.rm=T) else NA,
                    n_celltypes = length(sNames_0),
                    n_celltypes_used = length(sNames_4),
                    prop_genes_assign_mean = round(mean(sumstats_celltype_df$prop_genes_assign, na.rm=T),2),
                    prop_genes_t.test_fail_mean = if (kMSignifFilter) tryCatch({sum((sapply(list_geneMod_t.test, function(x)  {if (!is.null(x)) {sum(as.numeric(!x$signif))/length(x$signif)} else {NA_real_}}, simplify=T)))/ sum(sapply(list_geneMod_t.test, function(x) length(x$signif), simplify =T))}, error = function(err) {NA_real_}) else NA,
                    n_modules_total = sum(sumstats_celltype_df$n_modules))#,

  #params_run <- opt[-length(opt)] # all the user options/defaults
  
  # Convert sumstats and params from list to one-row data.frame
  matrix(unlist(sumstats_run), nrow=1, dimnames=list(NULL, c(names(sumstats_run)))) %>% as.data.frame(row.names=NULL, stringsAsFactors=F) -> sumstats_run_df
 
  # for (idx_null in which(sapply(params_run, is.null))) {
  #   params_run[idx_null] <- "NULL"
  # }
  
  #matrix(unlist(params_run,recursive = F), nrow=1, dimnames=list(NULL, c(names(params_run)))) %>% as.data.frame(row.names=NULL, stringsAsFactors=F) -> params_run_df
  
  # make path strings
  sumstats_celltype_path = sprintf("%s%s_%s_%s_sumstats_celltype.tab", dirLog, prefixData, prefixRun, flagDate)
  sumstats_run_path = sprintf("%s%s_%s_%s_sumstats_run.tab", dirLog, prefixData, prefixRun, flagDate)
  #params_run_path = sprintf("%s%s_%s_%s_params_run.tab", dirLog, prefixData, prefixRun, flagDate)
  
  # set append param values
  append_sumstats_celltype = F
  append_sumstats_run = F
  #append_params_run = F
  
  # write to file
  write.table(sumstats_celltype_df, file=sumstats_celltype_path, quote = F, sep = "\t", row.names=F, append = append_sumstats_celltype, col.names = !append_sumstats_celltype)
  write.table(sumstats_run_df, file=sumstats_run_path, quote = F, sep = "\t", row.names=F, append = append_sumstats_run, col.names = !append_sumstats_run)
  #write.table(params_run_df, file=params_run_path, quote = F, sep = "\t", row.names=F, append = append_params_run, col.names = !append_params_run)

  # save environment and parameters
  setwd(dirLog)
  saveMeta(doPrint=T)
  
  ######################################################################
  ################# SAVE SESSION IMAGE AND FINISH ######################
  ######################################################################
  
  message("Saving final session image")
  
  save.image(file=sprintf("%s%s_%s_final_session_image.RData.gz", dirRObjects, prefixData, prefixRun), compress = "gzip")
  
  tmp_file_paths <- dir(path = dirTmp, pattern = paste0(prefixData,"_", prefixRun), full.names = T)
  if (length(tmp_file_paths)>0) for (ch in tmp_file_paths) try(file.remove(ch))
  
  message("Script DONE!")
  
}

