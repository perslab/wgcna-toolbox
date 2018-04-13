# Title: Script to find robust WGCNA modules

######################################################################
############################## USAGE #################################
######################################################################

# e.g.


# RUNNING 
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_dev.R --data_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_neurons_sub.RData --dir_project /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-neurons-sub-19/ --data_prefix campbell-neurons-sub-19 --compare_params FALSE --do.center FALSE --genes_use PCA_5000 --corFnc bicor --networkType signed --anti_cor_action NULL --minClusterSize 15 --deepSplit 2 --moduleMergeCutHeight 0.20 --nPermutations 25 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --save_plots TRUE --plot_permuted F --n_cores 15

# TODO
# 
######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

# data_path = "/projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData"
# dir_project = "/projects/jonatan/tmp-rwgcna-tests/tmp-5/"
# data_prefix = "tmp-5"
# compare_params = F
# do.center = T
# genes_use = "PCA_5000"
# corFnc = "cor"
# networkType = "signed"
# anti_cor_action = NULL
# minClusterSize = 20
# deepSplit = 2
# moduleMergeCutHeight = 0.2
# replace = F
# nPermutations = 0
# STRINGdb_species = 10090
# ensembl_dataset = "mmusculus_gene_ensembl"
# save_plots = T
# plot_permuted = F
# n_cores = 4

######################################################################
########################### OptParse #################################
######################################################################

suppressMessages(library(optparse))

# specify our desired options in a list

option_list <- list(
  
  make_option("--data_path", type="character",
              help = "Provide full path to Rdata input file with Seurat object"),
  make_option("--dir_project", type="character", default=NULL,
              help = "Optional. Provide project directory. Must have subdirs RObjects, plots, tables. If not provided, assumed to be dir one level up from input data dir."),
  make_option("--data_prefix", type="character", default=paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_rWGCNA_run"),
              help = "Dataset prefix for output files"),
  make_option("--compare_params", type="logical", default=FALSE, 
              help="Compare many different cutreeHybrid parameters and select best. If FALSE, uses recommended parameters."), 
  make_option("--do.center", type="logical", default=T,
              help="Use centered data? In either case data is scaled and nUMI and mitochrondrial genes are regressed out."),
  make_option("--genes_use", type="character", default="PCA_5000",
              help="One of 'all', 'var.genes', 'hvg_<number of highly variable genes>', 'PCA_<number of high loading genes>'"), 
  make_option("--corFnc", type="character", default="cor",
              help="Use 'cor' for Pearson or 'bicor' for midweighted bicorrelation function."), 
  make_option("--networkType", type="character", default = "signed",
              help="'signed' scales correlations to [0:1]; 'unsigned' takes the absolute value (but the TOM can still be 'signed'); 'signed hybrid' sets negative correlations to zero."),
  make_option("--anti_cor_action", type="character", default=NULL, 
              help = "Optional. 'kME_reassign' reassigns genes with a negative kME more than 1.5 the kME w.r.t. their own (primary) module. Should be used only with networkType 'signed hybrid'."),
  make_option("--minClusterSize", type="integer", default=15L,
              help = "Minimum genes needed to form a module, recommended 5-25"),
  make_option("--deepSplit", type="integer", default=3L,
              help = "Controls the sensitivity of the cutreeDynamic/cutreeHybrid algorithm. Takes integer values 0-4, defaults to 3. Only applicable if test_params == F"),
  make_option("--moduleMergeCutHeight", type="double", default=0.2,
              help = "Cut-off level for the variable (1-correlation) for merging eigengenes. Recommended value range 0.1-0.2."),
  make_option("--nPermutations", type="integer", default=100L,
              help = "Number of times to permute the dataset, defaults to 100"),
  make_option("--replace", type="logical", default=TRUE,
              help = "Sample with replacement? Defaults to TRUE. If TRUE, uses all samples, if FALSE, uses 66% each time."),
  make_option("--STRINGdb_species", type="integer", default="10090",
              help = "Optional: species for which to retrieve protein data from STRINGdb to validate clusters. 10090 is mus musculus"),
  make_option("--ensembl_dataset", type="character", default="mmusculus_gene_ensembl",
              help = "Optional. Dataset for outputting colors with ensembl gene IDs for LD score regression."),
  make_option("--save_plots", type="logical", default=TRUE,
              help="Save plots?"),
  make_option("--plot_permuted", type="logical", default=FALSE,
              help="Compute and plot the modules on each resampled dataset? Good for visually inspecting how robust modules are to resampling."),
  make_option("--n_cores", type="integer", default=5,
              help = "Number of cores to use for parallelization [default %default]")
)

######################################################################
############################## LIBRARIES #############################
######################################################################

#suppressMessages(library(profvis)) #  for processing time and memory profiling
#suppressMessages(library(pryr)) # for processing time and memory profiling
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(WGCNA))
suppressMessages(library(STRINGdb))
suppressMessages(library(reshape))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(parallel))
#suppressMessages(library(biomaRt))

message("Libraries loaded")

######################################################################
################### GET COMMAND LINE OPTIONS #########################
######################################################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 

opt <- parse_args(OptionParser(option_list=option_list))

data_path <- opt$data_path 

dir_project <- opt$dir_project

data_prefix <- opt$data_prefix 

compare_params <- opt$compare_params
#if (compare_params==T) message("Sorry, compare_params is disabled in the current version of the script")
#compare_params <- F # Temporarily disable for bug testing TODO

do.center <- opt$do.center

genes_use <- opt$genes_use

corFnc <- opt$corFnc 

networkType <- opt$networkType

anti_cor_action <- opt$anti_cor_action

minClusterSize <- opt$minClusterSize 

deepSplit <- opt$deepSplit

moduleMergeCutHeight <- opt$moduleMergeCutHeight 

nPermutations <- opt$nPermutations

replace <- opt$replace

STRINGdb_species <- opt$STRINGdb_species

ensembl_dataset <- opt$ensembl_dataset

save_plots <- opt$save_plots

plot_permuted <- opt$plot_permuted

n_cores <- opt$n_cores  

######################################################################
############################## CONSTANTS #############################
######################################################################

if (!file.exists(data_path)) stop("Input data path not found")

# if no output directory provided, use the dir one level above that of the input file.

if (is.null(dir_project)) {
  pos <- tail(gregexpr("/", data_path)[[1]],2)[1]
  dir_project = substr(data_path, 1, pos)
}

# if specified output directory doesn't exist, create it 
if (!file.exists(dir_project)) {
  dir.create(dir_project) 
  message("Project directory not found, new one created")
}

saveIndividualTOMs = plot_permuted

dir_plots = paste0(dir_project,"plots/")
if (!file.exists(dir_plots)) dir.create(dir_plots) 

dir_tables = paste0(dir_project,"tables/")
if (!file.exists(dir_tables)) dir.create(dir_tables)

dir_RObjects = paste0(dir_project,"RObjects/")
if (!file.exists(dir_RObjects)) dir.create(dir_RObjects)

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

if (!(sapply(c("all", "var.genes", "hvg", "PCA"), function(x) grepl(x, genes_use, ignore.case=T)) %>% any())) {
  stop("genes_use must be one of 'all', 'var.genes', 'hvg_<number of highly variable genes>', 'PCA_<number of high loading genes>'")
}

if (!corFnc %in% c("cor", "bicor")) stop("corFnc must be one of 'cor' for Pearson's or 'bicor' for biweighted midcorrelation")

if (corFnc == "bicor" & do.center == T) warning("Using bicor with centered data will ignore non-positive values when detecting correlations")

if (!networkType %in% c('signed', 'unsigned', 'signed hybrid')) stop("networkType must be one of 'signed', 'unsigned' or 'signed hybrid'")

if (!minClusterSize %in% seq.int(from=5, to=100)) stop("minClusterSize must be an integer between 5 and 100")

if (!deepSplit %in% c(0L, 1L, 2L, 3L, 4L)) stop("deepSplit must be an integer between 0 and 4")

if (! (moduleMergeCutHeight > 0 & moduleMergeCutHeight < 0.5)) stop("moduleMergeCutHeight must be between 0 and 1, recommended range is between 0.1 and 0.2")
  
if (! (nPermutations >= 0 & nPermutations <= 100)) stop("nPermutations must be in the range 0-100")

if (! (n_cores >= 0 & n_cores <= 50)) stop("n_cores must be in the range 0-50")

# Load parameter values and utility functions
source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
source(file = "/projects/jonatan/functions-src/functions.R")


if (!is.null(ensembl_dataset)) {
  perslabMusEnsembl <- read.delim("/projects/jonatan/wgcna-src/rwgcna-pipeline/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
} 

if (nPermutations < 1) plot_permuted <- F 

# Set MAGMA variable values 
study_label = data_prefix
file_suffix = if (!is.null(STRINGdb_species)) "kMEs_PPI" else "kMEs"
output_label = flag_date

######################################################################
################# LOAD AND SUBSET SEURAT OBJECT ######################
######################################################################

message("Loading and subsetting seurat object..")

seurat_obj <- load_obj(data_path)

sNames <- names(table(seurat_obj@ident))

subsets <- lapply(sNames, function(x) SubsetData(seurat_obj,
                                                 ident.use = x, 
                                                 do.scale = F, 
                                                 do.center = F,
                                                 subset.raw = T))

names(subsets) <- sNames 

for (subsetName in sNames) {
  seurat_obj_sub <- subsets[[subsetName]]
  save(seurat_obj_sub, file=sprintf("%s%s_%s_seurat_subset_tmp_%s.RData", dir_RObjects, data_prefix, subsetName, flag_date))
}

rm(subsets)

message("Seurat object loaded and subsetted")

######################################################################
##################### SAVE PARAMETERS TO FILE ########################
######################################################################

subsets_n_cells <- data.frame(subset = sNames, n_cells = as.numeric(table(seurat_obj@ident)))

# Free up memory
rm(seurat_obj)

userParams <- cbind(list("data_path" = data_path,
                         "dir_project" = dir_project,
                         "data_prefix" = data_prefix,
                         "compare_params" = compare_params,
                         "do.center" = do.center,
                         "genes_use" = genes_use,
                         "corFnc" = corFnc,
                         "networkType" = networkType,
                         "anti_cor_action" = anti_cor_action,
                         "minClusterSize" = minClusterSize,
                         "deepSplit" = deepSplit,
                         "moduleMergeCutHeight" = moduleMergeCutHeight,
                         "nPermutations" = nPermutations,
                         "replace" = replace,
                         "STRINGdb_species" = STRINGdb_species,
                         "ensembl_dataset" = ensembl_dataset,
                         "save_plots" = save_plots,
                         "plot_permuted" = plot_permuted,
                         "n_cores" = n_cores))


builtInParams <- cbind(list("minCoreKME" = minCoreKME,
                            "minKMEtoStay" = minKMEtoStay,
                            "consensusQuantile" = consensusQuantile,
                            "replace" = replace, 
                            "fraction" = fraction,
                            "pamStage" =  if (compare_params==T) "compare_params" else pamStage,
                            "pamRespectsDendro" = if (compare_params==T) "compare_params" else pamRespectsDendro,
                            "hclustMethod" = hclustMethod,
                            "impute" = impute,
                            "nPC_seurat" = nPC_seurat,
                            "TOMType" = TOMType,
                            "PPI_pval_threshold" = PPI_pval_threshold))


# Save run params to file
write.csv(userParams, file=sprintf("%s%s_INFO_user_parameters_%s.csv", dir_tables, data_prefix, flag_date))
write.csv(builtInParams, file=sprintf("%s%s_INFO_built_in_parameters_%s.csv", dir_tables, data_prefix, flag_date))
write.csv(subsets_n_cells, file=sprintf("%s%s_INFO_subsets_n_cells_%s.csv", dir_tables, data_prefix, flag_date), row.names=F)

# clear up
rm(userParams)
rm(builtInParams)
rm(subsets_n_cells)

######################################################################
############# DEFINE THE FUNCTION FOR PARALLEL EXECUTION #############
######################################################################

parRWGCNA = function(sNames) { 

  disableWGCNAThreads() # to avoid nested parallelisation
  
  # retrieve subset from disk and delete it afterwards
  seurat_obj_sub <- load_obj(sprintf("%s%s_%s_seurat_subset_tmp_%s.RData", dir_RObjects, data_prefix, sNames, flag_date))
  file.remove(sprintf("%s%s_%s_seurat_subset_tmp_%s.RData", dir_RObjects, data_prefix, sNames, flag_date))
              
  ######################################################################
  ################# PROCESS SEURAT SUBSET OBJECT #######################
  ######################################################################
  
  # Filter out genes expressed in few cells
  if (min.cells > 0) {
    num.cells <- rowSums(seurat_obj_sub@data > 0)
    genes.use <- names(x = num.cells[which(x = num.cells >= min.cells)])
    seurat_obj_sub@raw.data <- seurat_obj_sub@raw.data[genes.use, ]
    seurat_obj_sub@data <- seurat_obj_sub@data[genes.use, ]
  }
  
  if (any(c("nUMI", "percent.mito") %in% names(seurat_obj_sub@meta.data))) {
    vars.to.filter_regress = c("nUMI", "percent.mito")[c("nUMI", "percent.mito") %in% names(seurat_obj_sub@meta.data)]
    seurat_obj_sub <- FilterCells(object = seurat_obj_sub, subset.names = vars.to.filter_regress, low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.2)) 
  } else {
    warning("nUMI and percent.mito not found in Seurat object meta data")
  }
  
  
  if (genes_use != "all"){ #if we do use all genes, we'll rune ScaleData below
    
    # Need to run ScaleData to get variable genes for PCA / hvg / var.genes
    seurat_obj_sub <- ScaleData(object = seurat_obj_sub, 
                                    #genes.use = NULL, # default: genes.use = all genes in @data
                                  vars.to.regress = if (any(c("nUMI", "percent.mito") %in% names(seurat_obj_sub@meta.data))) vars.to.filter_regress else NULL, 
                                  model.use="linear",
                                  do.par=F, #T,
                                  #num.cores = n_cores,
                                  do.scale=T,
                                  do.center=T)
        
    seurat_obj_sub <- FindVariableGenes(object = seurat_obj_sub,
                                        #mean.function = ExpMean,
                                        #dispersion.function = LogVMR,
                                        x.low.cutoff = 0.0125,
                                        x.high.cutoff = 3,
                                        y.cutoff = 0.5, 
                                        do.plot=F)
    
    
    if (grepl("PCA", genes_use, ignore.case = T)) {
      seurat_obj_sub <- RunPCA(object = seurat_obj_sub,
                               pcs.compute = min(nPC_seurat, 
                                                 length(seurat_obj_sub@var.genes)-1, 
                                                 ncol(seurat_obj_sub@data)-1),
                               do.print = F,
                               seed.use = randomSeed,
                               maxit = maxit,
                               fastpath=fastpath)
      
    }
 
  }

  # Store metadata
  metaData <- seurat_obj_sub@meta.data
  
  # Select genes

  if (do.center == F | genes_use == "all") {
    seurat_obj_sub <- ScaleData(object = seurat_obj_sub, 
                                #genes.use = NULL, # default: genes.use = all genes in @data
                                vars.to.regress = if (any(c("nUMI", "percent.mito") %in% names(seurat_obj_sub@meta.data))) vars.to.filter_regress else NULL, 
                                model.use="linear",
                                do.par=F,#T,
                                #num.cores = n_cores,
                                do.scale=T,
                                do.center=do.center) 
  }

  tryCatch({
    # Select a subset of the genes 
    if (genes_use == "all") {
      input_data <- seurat_obj_sub@scale.data
    } else if (genes_use == "var.genes") {
      input_data <- seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% seurat_obj_sub@var.genes,] 
    } else if (grepl("PCA", genes_use, ignore.case = T)) {
      n_genes <- as.numeric(gsub("[^0-9]", "", genes_use))
      seurat_obj_sub <- ProjectPCA(seurat_obj_sub, do.print = F, pcs.print = NULL, pcs.store = nPC_seurat, genes.print = NULL, replace.pc=F, do.center=T)
      max_loads <- apply(abs(seurat_obj_sub@dr$pca@gene.loadings.full), MARGIN=1, max, T)
      names_genes_use = names(sort(max_loads, decreasing=T)[1:n_genes])
      input_data <- seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,] 
    } else if (grepl("hvg", genes_use, ignore.case = T)) {
      n_genes <- as.numeric(gsub("[^0-9]", "", genes_use))
      names_genes_use = rownames(seurat_obj_sub@hvg.info)[order(seurat_obj_sub@hvg.info$gene.dispersion.scaled, decreasing=T)][1:n_genes]
      input_data <- seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,] 
    }
  }, error = function(c) {write.csv("ERROR", file=sprintf("%s%s_%s_input_data_selection_ERROR_%s.csv", dir_project, data_prefix, sNames, flag_date), row.names = F)
    stop("Input data error")})
  
    # Format the data: WGCNA needs a dataframe with the genes in columns
    datExpr <- t(input_data) 
    colnames(datExpr) <- rownames(input_data)

    # clear memory
    rm(seurat_obj_sub) # clear memory
    rm(input_data) #clear memory
    gc()


  ######################################################################
  ####### PICK SOFT THRESHOLD POWER FOR ADJACENCY MATRIX ###############
  ######################################################################
  
  softPower <- 8 # Set a default value as fall back
  
  tryCatch({
  
    disableWGCNAThreads()
    #powers = c(c(1:10), seq(from = 12, to=20, by=2))
    powers = c(1:20)

    sft = pickSoftThreshold(data=datExpr,
                            powerVector = powers,
                            blockSize = min(2500, ncol(datExpr)), #try to prevent crashing
                            corFnc = corFnc,
                            corOptions =  corOptions,
                            networkType = networkType,
                            verbose = verbose)
    
    if (save_plots==TRUE) {
      
      pdf(sprintf("%s%s_%s_pickSoftThresholdSFTFit_%s.pdf", dir_plots, data_prefix, sNames, flag_date),width=10,height=5)
      #par(mfrow = c(1,2));
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
      pdf(sprintf("%s%s_%s_pickSoftThresholdMeanCon_%s.pdf", dir_plots, data_prefix, sNames, flag_date),width=10,height=5)
      plot(sft$fitIndices[,1], sft$fitIndices[,5],
           xlab="Soft Threshold (power)",
           ylab="Mean Connectivity",
           type="n",
           main = paste("Mean connectivity"))
      text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
      dev.off()
    }
    
    # Select softPower: of lower .95 percentile connectivity, if several softPowers achieve 0.9 R.sq, 
    # take smallest softPower; else take best fitting one
    
    fitIndices <- as.data.frame(sft$fitIndices)
    fitIndices_filter <- fitIndices %>% dplyr::filter(median.k. <= quantile(median.k.,0.95, na.rm=T))
    if (sum(fitIndices_filter$SFT.R.sq >= 0.9) > 1) {
      softPower = min(fitIndices_filter$Power[fitIndices_filter$SFT.R.sq>=0.9])
    } else {
      softPower <- (fitIndices_filter %>% dplyr::arrange(desc(SFT.R.sq)) %>% dplyr::select(Power))[1,]
    }
  }, error = function(c) {
    write.csv("ERROR", file=sprintf("%s%s_%s_pickSoftThreshold_ERROR_set_to_8_as_default_%s.csv", dir_project, data_prefix, sNames, flag_date), row.names = F)})
  
  ######################################################################
  ############### RESAMPLE THE DATA FOR ROBUSTNESS #####################
  ######################################################################

  if (nPermutations > 0) {
    multiExpr <- bootstrap(datExpr=datExpr,
                           nPermutations = nPermutations,
                           replace = replace,
                           fraction = fraction,
                           randomSeed = randomSeed)
    
    ######################################################################
    ####################### SET WORKING DIRECTORY ########################
    ######################################################################
    
    # TOM files will be saved here by consensusTOM / blockwiseConsensusModules
    
    setwd(dir_RObjects)
    
    ######################################################################
    ####################### FIND CONSENSUS TOM ###########################
    ######################################################################
    
    consensus0 <- consensusTOM(
      # ... information needed to calculate individual TOMs
      multiExpr = multiExpr,
      
      # Data checking options
      checkMissingData = checkMissingData,
      
      # Blocking options
      #blocks = NULL,            
      maxBlockSize = maxBlockSize,
      blockSizePenaltyPower = blockSizePenaltyPower,
      #nPrecluste              ringCenters = NULL,
      randomSeed = randomSeed,
      
      # Network construction arguments: correlation options
      
      corType = corType,
      maxPOutliers = maxPOutliers,
      quickCor = quickCor,
      pearsonFallback = pearsonFallback,
      cosineCorrelation = cosineCorrelation,
      replaceMissingAdjacencies = replaceMissingAdjacencies,
      
      # Adjacency function options
      
      power = softPower,
      networkType = networkType,
      #checkPower = checkPower,
      
      # Topological overlap options
      
      TOMType = TOMType,
      TOMDenom = TOMDenom,
      
      # Save individual TOMs?
      
      saveIndividualTOMs = saveIndividualTOMs,
      individualTOMFileNames = paste0(sNames, "_individualTOM-Set%s-Block%b.RData"),
      
      # ... or individual TOM information
      
      #individualTOMInfo = NULL,
      #useIndivTOMSubset = NULL,
      
      ##### Consensus calculation options
      
      #useBlocks = NULL,
      
      # Network calibration
      networkCalibration = networkCalibration,
      #saveCalibratedIndividualTOMs = FALSE,
      #calibratedIndividualTOMFilePattern = "calibratedIndividualTOM-Set%s-Block%b.RData",
      
      # Simple quantile calibration options
      #calibrationQuantile = calibrationQuantile, # Only used if networkCalibration = "single quantile"
      sampleForCalibration = sampleForCalibration,
      sampleForCalibrationFactor = sampleForCalibrationFactor,
      getNetworkCalibrationSamples = getNetworkCalibrationSamples,
      
      # Consensus definition
      consensusQuantile = consensusQuantile,
      useMean = useMean,
      #setWeights = setWeights, # only if using a
      # Return options
      saveConsensusTOMs = saveConsensusTOMs,
      consensusTOMFilePattern = paste0(sNames,"_consensusTOM-block.%b.RData"),
      returnTOMs = F,
      # Internal handling of TOMs
      useDiskCache = T,
      #chunkSize = NULL, # automatic
      cacheDir = dir_RObjects,
      cacheBase = ".blockConsModsCache",
      #nThreads = n_cores,
      verbose = verbose,
      indent = indent)
        
    # blockwiseConsensusModules has filtered out some genes when creating the consensus TOM.
    # We need to only work with these from here on
    goodGenesTOM_idx <- consensus0$goodSamplesAndGenes$goodGenes
    datExpr_filter <- datExpr[,goodGenesTOM_idx] # keep only the genes kept by consensusTOM
    
    # Load the consensus TOM matrix 
    suppressMessages(load(sprintf("%s%s_consensusTOM-block.1.RData", dir_RObjects, sNames))) # Load the consensus TOM generated by blockwiseConsensusModules
    
    dissTOM <- 1-as.dist(consTomDS) # Convert proximity to distance
    
    # Clean up
    rm(datExpr)
    rm(consTomDS)
    gc()
    
  } else {
    
    adjacency = adjacency(datExpr, 
                          type=type, 
                          power = softPower, 
                          corFnc = corFnc, 
                          corOptions = corOptions)
    
    TOM = TOMsimilarity(adjMat=adjacency,
                        TOMType=TOMType,
                        TOMDenom=TOMDenom,
                        verbose=verbose,
                        indent = indent)
    
    dissTOM <- 1-as.dist(TOM) # Convert proximity to distance
    
    goodGenesTOM_idx <- rep("TRUE", ncol(datExpr))
    datExpr_filter <- datExpr #rename for consistency with script later

    rm(datExpr)
    rm(adjacency)
    rm(TOM)
    gc()
    
  }
  
  
  ######################################################################
  ####################### CLUSTERING ON THE TOM ########################
  ######################################################################
  
  # We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes.
  # We use Dynamic Tree Cut from the package dynamicTreeCut.    
   
  geneTree <- hclust(d=dissTOM,
                     method=hclustMethod) # Set up the dendrogram ready for cutting
  
  ######################################################################
  ########################### PLOT_PERMUTED ############################
  ######################################################################
  

  tryCatch({
    
    if (plot_permuted == T) {
      
      # NOTE: Alternatively use sampledBlockwiseModules() - see /projects/jonatan/wgcna-src/180228_robust_wgcna_module_detection.Rmd

      indiv_TOMs <- lapply(1:nPermutations, function(x) load_obj(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", dir_RObjects, sNames, x, 1))) # retrieve permuted dataset TOMs from disk
      dissTOMs <- lapply(indiv_TOMs, function(x) 1-as.dist(x)) # Convert permuted TOMs to distance matrices
      geneTrees <- lapply(dissTOMs, function(x) hclust(d=x, method=hclustMethod)) # Do hierarchical clustering 
      # cut out modules in each gene tree
      cutrees <- mapply(function(x,y) cutreeHybrid(dendro=x, cutHeight = NULL, minClusterSize = minClusterSize, distM=as.matrix(y), deepSplit=deepSplit,pamStage=pamStage, pamRespectsDendro=pamRespectsDendro), geneTrees, dissTOMs, SIMPLIFY=F, .scheduling="dynamic") 
  
      # merge modules with correlated eigengenes
      colors_MEs <- mapply(function(x,y) mergeCloseModules(exprData=as.matrix(x$data)[,goodGenesTOM_idx], colors = y$labels, cutHeight=moduleMergeCutHeight), multiExpr, cutrees, .scheduling="dynamic") # Merge modules whose eigengenes are highly correlated
      
      permutedColors = lapply(colors_MEs, function(x) labels2colors(x$colors)) 
 
      permutedLabels <- paste0("Permutation ", as.character(1:nPermutations)) 
 
      # Free up memory
      rm(indiv_TOMs)
      rm(dissTOMs)
      rm(colors_MEs)
      # plot permuted dataset modules later when we have the consensus modules
    }
  }, error = function(c) { write.csv("ERROR", file=sprintf("%s%s_%s_plot_permuted_ERROR_%s.csv", dir_project, data_prefix, flag_date, sNames, flag_date), row.names = F)
    stop("Error while computing modules on permuted datasets")})
  
  if (nPermutations>0) rm(multiExpr)
  
  # Delete individual TOMs produced by consensusTOMs / blockwiseConsensusModules from disk
  for (i in 1:nPermutations) {
    if (file.exists(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", dir_RObjects, sNames, i, 1))) {
      file.remove(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", dir_RObjects, sNames, i, 1))
    }
  }
  
  ######################################################################
  ########################### CUTREEHYBRID #############################
  ######################################################################
  
  # Call cutreeHybrid to cut out modules
  
  if (compare_params == T) {
    # Iterate over different combinations of parameters
    # Parallelised version of the for loop used in Gandal et al 2018 to test different cluster params
    
    minModSizeVals = c(10,20) # No reason to go higher as we can later merge close modules
    deepSplitVals =  c(2:4)
    pamVals = c(TRUE,FALSE) # sets both pamState and pamRespectsDendro
    dthreshVals = c(0.1, 0.2) # What distance of eigengene correlation to merge modules
    
    dims = sapply(list(minModSizeVals, deepSplitVals, pamVals, dthreshVals), length) # list of number of values per parameter
    n_combs <- prod(dims) # Total number of combinations of parameter values
    comb_list = vector("list", length = n_combs)
    
    k=1
    for (minModSize in 1:dims[1]) {
      for (ds in 1:dims[2]) { # Gandal et al 2018 use c(50,100, 200)
        for (pam in 1:dims[3]) {
          for (dthresh in 1:dims[4]) {
            comb_list[[k]] <- c(minModSizeVals[minModSize],deepSplitVals[ds], pamVals[pam], dthreshVals[dthresh])
            k = k+1
          }
        }
      }
    }
    
    cutrees <- lapply(comb_list, parCutreeHybrid) # Cut out different color assignments in the dendrogran
    
    colors_MEs <- mapply(parMergeCloseModules, cutrees, comb_list, SIMPLIFY=FALSE) # Merge modules whose eigengenes are highly correlated
    listColors = lapply(colors_MEs, function(x) x[[1]]) # Merged colors
    MEs = lapply(colors_MEs, function(x) x[[2]]) # Merged color eigengenes

    #if (save_plots == TRUE) {
    labels <- lapply(comb_list, parLabel) # Labels for plot
    
    ######################################################################
    ######### EVALUATE AND SELECT FINAL CLUSTERING PARAMETERS ############
    ######################################################################
    
    # Use module preservation statistics to identify the best set of clustering parameters
    # See sections 4-6 in Langfelder,...,Horvath_2011_Is My Network Module Preserved and Reproducible_S1
    
    n_combs <- length(listColors)
    n_modules <- sapply(listColors, FUN = function(x) length(unique(x)), simplify=T) # get the number of modules for each set of params
    min_modules <- quantile(n_modules, 0.66, na.rm=T) # Choose quantile as cut-off to continue only with parameters that found a higher number of modules
    idx_quantile <- n_modules >= min_modules
    
    if (save_plots == TRUE) {
      pdf(sprintf("%s%s_%s_%s_diffParams_n_modules.pdf", dir_plots, data_prefix, sNames, flag_date), width=10, height=5)
      plot(1:n_combs, n_modules, xlab = "Parameter setting", main="Number of modules found for different parameters")
      abline(h=min_modules)
      dev.off()
    }
    
    listDatExpr <- list("data"=datExpr_filter)
    names(listColors[idx_quantile]) <- labels[idx_quantile]
    
    tryCatch({
      
    clust_qual_stats <- wrapModulePreservation(listDatExpr = listDatExpr,
                                               listColors = listColors[idx_quantile],
                                               #labels = labels, # not needed since listColors has names
                                               dataIsExpr = dataIsExpr,
                                               networkType = networkType,
                                               corFnc = corFnc,
                                               corOptions = corOptions,
                                               # We actually just want to get quality stats on the reference network,
                                               # but the function requires testNetworks..
                                               nPermutations = nPermutations,
                                               includekMEallInSummary = includekMEallInSummary,
                                               restrictSummaryForGeneralNetworks = restrictSummaryForGeneralNetworks,
                                               calculateQvalue = calculateQvalue,
                                               randomSeed = randomSeed,
                                               maxGoldModuleSize = maxGoldModuleSize,
                                               maxModuleSize = maxModuleSize,
                                               quickCor = quickCor,
                                               ccTupletSize = ccTupletSize,
                                               #calculateCor.kIMall = calculateCor.kIMall,
                                               calculateClusterCoeff = calculateClusterCoeff,
                                               useInterpolation = useInterpolation,
                                               checkData = checkData,
                                               #greyName = greyName,
                                               savePermutedStatistics = savePermutedStatistics,
                                               loadPermutedStatistics = loadPermutedStatistics,
                                               permutedStatisticsFile = permutedStatisticsFile,
                                               plotInterpolation = plotInterpolation,
                                               interpolationPlotFile = interpolationPlotFile,
                                               discardInvalidOutput = discardInvalidOutput,
                                               parallelCalculation = parallelCalculation,
                                               verbose = verbose,
                                               indent = indent)
    
    # A statistic for evaluating the quality of a cluster/module assignment is the
    # proportion of variance in the expression of the genes that the cluster's first
    # principal component or 'eigengene' accounts for.
    # Here, for each set of parameters, we find the median module statistic value for each set of
    # parameters to rank the sets and pick a single set of parameters
    
    medianPropVarExpl <- numeric(length(listColors[idx_quantile]))
    for (i in 1:length(clust_qual_stats$quality$observed)[idx_quantile]) {
      medianPropVarExpl[i] <- median(clust_qual_stats$quality$observed[[i]]$'inColumnsAlsoPresentIn.test'$propVarExplained.qual, na.rm = T)
    }
    
    # Order the module assignments by the statistic above
    listColors_order <- listColors[idx_quantile][order(medianPropVarExpl, decreasing=T)]
    labels_order <- labels[idx_quantile][order(medianPropVarExpl, decreasing=T)]
    MEs_order <- MEs[idx_quantile][order(medianPropVarExpl, decreasing=T)]
    # plot the module assignments but only those in the selected quantile of number of modules,
    # and ranked by the mean proportion of variance within a module that the module eigengene explains
    
    
    diffParamsColors_order = matrix(0, ncol(datExpr_filter), n_combs)
    diffParamsColors_order[, 1] = listColors_order[[1]] # use the best colors as 'reference'
    
    for (r in 2:(n_combs))
    {
      diffParamsColors_order[, r] = matchLabels(source = listColors_order[[r-1]], # traverse gene colour assignments
                                              reference = diffParamsColors_order[, 1],
                                              pThreshold = 5e-2,
                                              na.rm = TRUE,
                                              extraLabels = paste0("Extra_", 0:100)) # labels for modules in source that cannot be matched to any modules in reference
    }
    
    pdf(sprintf("%s%s_%s_%s_diffParamsColors_ordered.pdf", dir_plots, data_prefix, sNames, flag_date),width=8,height=12)
    plotDendroAndColors(geneTree, 
                        labels2colors(diffParamsColors_order),
                        #matrix(unlist(permutedColors_final), ncol = length(permutedLabels_final), byrow = F), 
                        groupLabels=labels_order, 
                        addGuide= TRUE, 
                        dendroLabels=FALSE, 
                        main="Modules with diff. params ranked by the ME prop. var. expl.", 
                        cex.colorLabels=0.5)
    dev.off()
  
    # Select a set of preferred colors and add them as the first column to be shown at the top of the figure
    
    colors <- listColors_order[[1]]
    labels <- labels_order[[1]]
    MEs <- MEs_order[[1]]
    cutree_params_final = as.list(comb_list[idx_quantile][order(medianPropVarExpl,decreasing = T)][[1]]) # TODO CHECK THIS
    names(cutree_params_final) = c("minClusterSize", "deepSplit","pamStage", "cutHeight")
    
    }, error = function(c) { write.csv("ERROR", file=sprintf("%s%s_%s_clust_qual_stats_ERROR_%s.csv", dir_project, data_prefix, flag_date, sNames, flag_date), row.names = F)
      stop("Error while computing clustering quality statistics")})
    
    
    # (compare_params == T) chunk finished
    
  } else if (compare_params == F) { # Use given set of parameters
    
    cutree <- cutreeHybrid(dendro = geneTree,
                           cutHeight = NULL,
                           minClusterSize = minClusterSize,
                           distM = as.matrix(dissTOM),
                           deepSplit = deepSplit,
                           pamStage = pamStage,
                           pamRespectsDendro = pamRespectsDendro) # Cut out different color assignments in the dendrogran
    
    if (length(unique(cutree$labels)) == 1) {
      
      write.csv("ERROR", file=sprintf("%s%s_%s_stop_no_modules_found_%s.csv", dir_project, data_prefix, sNames, flag_date))
      stop("No modules found")
      
    } else {  
      
      merged = mergeCloseModules(exprData = as.matrix(datExpr_filter),
                                 colors = cutree$labels,
                                 cutHeight = moduleMergeCutHeight)
      
      colors = labels2colors(merged$colors)
      
      #MEs = merged$newMEs
      MEs = moduleEigengenes(expr = as.matrix(datExpr_filter),
                             colors,
                             excludeGrey = F)

    }
  }
  
  ######################################################################
  ######################################################################
  ######################################################################
  
  names(colors) <- colnames(datExpr_filter)
  
  # Clear up
  rm(dissTOM)
  gc()
  
  
  ######################################################################
  ######################## COMPUTE MODULE KMEs #########################
  ######################################################################
  
  tryCatch({
    # All module kMEs
    kMEs = signedKME(as.matrix(datExpr_filter),
                     MEs$eigengenes,
                     corFnc = corFnc)
    
    # Get primary module kMEs
    pkMEs <- vector(length=length(colors))
    
    for (i in 1:length(colors)) {
      pkMEs[i] <- as.numeric(unlist(kMEs %>% dplyr::select(matches(colors[[i]]))))[[i]]
    }
    
    # Delete grey kMEs 
    if (any(grepl("kMEgrey", colnames(kMEs)))) kMEs[['kMEgrey']] <- NULL

  }, error = function(c) {  write.csv("ERROR", file=sprintf("%s%s_%s_compute_kMEs_ERROR_%s.csv", dir_project, data_prefix, sNames, flag_date), row.names = F)
    stop("Error while computing kMEs")})

  ######################################################################
  #################### REASSIGN GENES BASED ON KMEs ####################
  ######################################################################
  
  if (anti_cor_action == "kME_reassign") {
    tryCatch({
      col_most_negkME = max.col(-1*kMEs) # indices of most negative kMEs
      rows_to_reassign <- logical(length=nrow(kMEs))
      
      for (i in nrow(kMEs)) {
        rows_to_reassign[i] <- kMEs[i, col_most_negkME[i]] < -kME_reassign_threshold*(pkMEs[i])
      }
      
      if (sum(rows_to_reassign) > 0) {
        
        colors[rows_to_reassign] <- colors[col_most_negkME][rows_to_reassign]
        
        # Recompute Module Eigengenes, kME, pkME
        MEs = moduleEigengenes(expr = as.matrix(datExpr_filter),
                               colors,
                               excludeGrey = F)
                               #softPower = softPower)
        
        kMEs = signedKME(as.matrix(datExpr_filter),
                         MEs$eigengenes,
                         corFnc = corFnc)
        
        pkMEs <- vector(length=length(colors))
        
        for (i in 1:length(colors)) {
          pkMEs[i] <- as.numeric(unlist(kMEs %>% dplyr::select(matches(colors[[i]]))))[[i]]
        }
      }
      message(paste0("Reassigned ", sum(rows_to_reassign), " genes to new modules based on a higher negative kME."))
    }, error = function(c) {
      write.csv("ERROR", file=sprintf("%s%s_%s_kMEs_reassign_ERROR_%s.csv", dir_project, data_prefix, sNames, flag_date), row.names = F)})
    } 
  ######################################################################
  ######################### NAME PKME VECTOR ###########################
  ######################################################################
  
  names(pkMEs) <- names(colors) 
  
  ######################################################################
  ######################### RUN PPI ENRICHMENT #########################
  ######################################################################
  
  message(paste0("Checking modules for Protein-Protein Interaction (PPI) enrichment via STRINGdb, adj p-value cut-off: ", PPI_pval_threshold))
  
  # TODO : try ensembleIDs instead?
  # get indices for genes for which primary kMEs are over a set threshold to validate
  # them with Protein-Protein Interaction database. Set it low since kME depends only on
  # the first PC
  
  module_PPI <- NULL
  unique_colors <- NULL
  unique_colors_PPI <- NULL
  colors_PPI <- NULL
  idx_pkME_ok <- NULL
  MEs_PPI <- NULL
  kMEs_PPI <- NULL
  pkMEs_PPI <- NULL 
  
  if (!is.null(STRINGdb_species)) {
      
    idx_pkME_ok <- abs(pkMEs) > 0.2
    unique_colors <- unique(colors[idx_pkME_ok])
    string_db <- STRINGdb$new(version="10", species = STRINGdb_species, score_threshold=0, input_directory="") # Default: 10090 (Mus musculus)
    module_PPI <- sapply(unique_colors, function(x) PPI_for_par(x, unique_colors = unique_colors, idx_pkME_ok = idx_pkME_ok, colors = colors, string_db = string_db), simplify=T) %>% t()
    module_PPI <- as.data.frame(module_PPI, row.names = unique_colors)
    
    ######################################################################
    ############### FILTER MODULES ON PPI ENRICHMENT  ####################
    ######################################################################
    
    if (!is.null(module_PPI)) {
     
      unique_colors_PPI = unique_colors[module_PPI$'p-value' < PPI_pval_threshold]
      genes_PPI_idx <- colors %in% unique_colors_PPI
      colors_PPI <- colors
      colors_PPI[!genes_PPI_idx] <- "grey"
      
      ######################################################################
      ############# RECOMPUTE MEs, kMEs, pkMEs AFTER PPI FILTER ############
      ######################################################################
      
      # TODO better to filter
      # recompute Module Eigengenes
      
      if (sum(genes_PPI_idx)!=0) {
        
        MEs_PPI = moduleEigengenes(expr = as.matrix(datExpr_filter),
                                   colors_PPI,
                                   excludeGrey = F) # T leads to errors when finding pkMEs
                                   #softPower = softPower)
        # recompute kMEs
        kMEs_PPI = signedKME(as.matrix(datExpr_filter),
                             MEs_PPI$eigengenes,
                             corFnc = corFnc)
        
        # recompute primary kMEs
        pkMEs_PPI <- vector(length=length(colors_PPI))
        for (i in 1:length(colors_PPI)) {
          pkMEs_PPI[i] <- as.numeric(unlist(kMEs_PPI %>% dplyr::select(matches(colors_PPI[[i]]))))[[i]]
        }
        
        names(pkMEs_PPI) <- names(colors_PPI)
        
        # Delete grey kME (needed it for pkME step though)
        if (any(grepl("kMEgrey", colnames(kMEs_PPI)))) kMEs_PPI[['kMEgrey']] <- NULL
        
      } else {
        write.csv("WARNING" , file=sprintf("%s%s_%s_WARNING_sum_genes_PPI_0_%s.csv", dir_project, data_prefix, sNames, flag_date))
      }
    } # end of if (!is.null(module_PPI)) 
  } # end of if (!is.null(STRINGdb_species)) 
  

  ######################################################################
  ######################### PLOT FINAL COLORS #########################
  ######################################################################
  
  tryCatch({
    all_colors = if (!is.null(STRINGdb_species) & !is.null(module_PPI)) cbind(colors, colors_PPI) else colors
    
    if (save_plots==T) {
      pdf(sprintf("%s%s_%s_%s_colors.pdf", dir_plots, data_prefix, sNames, flag_date),width=10,height=5)
      plotDendroAndColors(geneTree,
                          all_colors,
                          groupLabels = if (compare_params==T) labels else NULL,
                          addGuide=T,
                          dendroLabels=F,
                          main="Final colors",
                          cex.colorLabels=0.3)
      dev.off()
      
      if (plot_permuted==T) {  
        
        # Align the module coloring schemes 
        
        permutedColors_final = matrix(0, ncol(datExpr_filter), nPermutations + 1)
        permutedColors_final[, 1] = colors # use the consensus colors as 'reference'
        
        for (r in 2:(nPermutations+1))
        {
          permutedColors_final[, r] = matchLabels(source = permutedColors[[r-1]], # traverse gene colour assignments
                                    reference = permutedColors_final[, 1],
                                    pThreshold = 5e-2,
                                    na.rm = TRUE,
                                    extraLabels = paste0("Extra_", 0:100)) # labels for modules in source that cannot be matched to any modules in reference
        }

        permutedLabels_final <- c("Final colors", permutedLabels)

        pdf(sprintf("%s%s_%s_%s_diffpermutedColors.pdf", dir_plots, data_prefix, sNames, flag_date),width=8,height=8)
        plotDendroAndColors(geneTree, 
                            labels2colors(permutedColors_final),
                            #matrix(unlist(permutedColors_final), ncol = length(permutedLabels_final), byrow = F), 
                            groupLabels=permutedLabels_final, 
                            addGuide= TRUE, 
                            dendroLabels=FALSE, 
                            main="Modules: final and on resampled data", 
                            cex.colorLabels=0.5)
        dev.off()
      }
    }
  }, error = function(c) write.csv("ERROR", file=sprintf("%s%s_%s_plot_colors_ERROR_%s.csv", dir_project, data_prefix, sNames, flag_date), row.names = F))
  
  #rm(datExpr_filter)
  ######################################################################
  ############################ GET ENSEMBL IDs #########################
  ######################################################################
  # https://support.bioconductor.org/p/45386/
  # https://www.bioconductor.org/packages/3.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
  
  ensembl_out <- NULL
  ensembl_PPI_out <- NULL
  
  tryCatch({
    if (!is.null(ensembl_dataset)) {
      ensembl_IDs <- perslabMusEnsembl$ensembl_gene_id[match(names(colors), perslabMusEnsembl$gene_name_optimal)] 
      ensembl_out <- data.frame(module_name = colors, ensembl_ID = ensembl_IDs, row.names = NULL)
      prop_not_map <- paste0("EnsemblID_proportion_not_mapped:_", signif(as.double(sum(is.na(ensembl_IDs))) / as.double(1e-10+length(ensembl_IDs)),2))
      write.csv("", file=sprintf("%s%s_%s_%s_%s.csv", dir_tables, data_prefix, sNames, prop_not_map, flag_date))
      if (!is.null(colors_PPI)) {
        colors_PPI_noGrey <- colors_PPI[colors_PPI != grey]
        ensembl_PPI_IDs <- perslabMusEnsembl$ensembl_gene_id[match(names(colors_PPI_noGrey), perslabMusEnsembl$gene_name_optimal)] 
        ensembl_PPI_out <- data.frame(module_name_PPI = colors_PPI_noGrey, ensembl_ID = ensembl_PPI_IDs, row.names = NULL) 
        prop_not_map_PPI <- paste0("EnsemblID_PPI_proportion_not_mapped:_", signif(as.double(sum(is.na(ensembl_PPI_IDs))) / as.double(1e-10+length(ensembl_PPI_IDs)),2))
        write.csv("", file=sprintf("%s%s_%s_%s_%s.csv", dir_tables, data_prefix, sNames, prop_not_map_PPI, flag_date))
      }
    }
    # ensembl <- useMart("ensembl", dataset=ensembl_dataset) 
    # ensembl_query <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = ensembl)
    # ensembl_IDs <- ensembl_query$ensembl_gene_id[match(names(colors), ensembl_query$external_gene_name)]
    # ensembl_out <- data.frame(module_name = colors, ensembl_ID = ensembl_IDs, row.names = NULL)
    # if (!is.null(colors_PPI)) {
    # ensembl_PPI_IDs <- ensembl_query$ensembl_gene_id[match(names(colors_PPI), ensembl_query$external_gene_name)]
    # ensembl_PPI_out <- data.frame(module_name_PPI = colors_PPI, ensembl_ID = ensembl_PPI_IDs, row.names = NULL)
    # }
    
    
  }, error = function(c) write.csv("ERROR", file=sprintf("%s%s_%s_ensemblID_ERROR_%s.csv", dir_project, data_prefix, sNames, flag_date), row.names = F))

  ######################################################################
  ################### SAVE RUN OBJECTS, RESULTS AND PARAMS #############
  ######################################################################

  
  tryCatch({
    
    kMEs <- cbind(genes=names(colors), kMEs)
    pkMEs <- cbind(colors=colors, genes=names(colors), pkMEs = pkMEs)
    
    # save kMEs and kMEs primary to file
    write.csv(kMEs, file=sprintf("%s%s_%s_kMEs_%s.csv", dir_tables, data_prefix, sNames, flag_date), row.names=F)
    write.csv(pkMEs, file=sprintf("%s%s_%s_pkMEs_%s.csv", dir_tables, data_prefix, sNames, flag_date), row.names = F)
    
    if (!is.null(kMEs_PPI)) {
      kMEs_PPI <- cbind(genes=names(colors_PPI), kMEs_PPI) 
      pkMEs_PPI <- data.frame(colors_PPI = colors_PPI, genes = names(colors_PPI), pkMEs_PPI = pkMEs_PPI) 
      write.csv(kMEs_PPI, file=sprintf("%s%s_%s_kMEs_PPI_%s.csv", dir_tables, data_prefix, sNames, flag_date), row.names=F)
      write.csv(pkMEs_PPI, file=sprintf("%s%s_%s_pkMEs_PPI_%s.csv", dir_tables, data_prefix, sNames, flag_date), row.names=F)
    }
    
    # ensembl IDs for LD score regression
    if (!is.null(ensembl_out)) {
      write.csv(ensembl_out, file=sprintf("%s%s_%s_ensembl_out_%s.csv", dir_tables, data_prefix, sNames, flag_date), row.names=F)
      if (!is.null(ensembl_PPI_out)) write.csv(ensembl_PPI_out, file=sprintf("%s%s_%s_ensembl_PPI_out_%s.csv", dir_tables, data_prefix, sNames, flag_date), row.names=F)
    }   

    save.image(file=sprintf("%s%s_%s_rsession_%s.RData", dir_RObjects, data_prefix, sNames, flag_date))
    
  }, error = function(c) write.csv("ERROR", file=sprintf("%s%s_%s_save_results_ERROR_%s.csv", dir_project, data_prefix, sNames, flag_date), row.names = F))

  ##########################################################################
  ################################## MAGMA #################################
  ##########################################################################
  
  file_sep = ','

  project_path = "/projects/jonatan/tmp-bmi-brain/"
  magma_gwas_path = paste(project_path,"data/magma/",sep="/")
  figs_path = paste(project_path,"figs/",sep="/")
  mapping_hs_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz"
  mapping_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
  mapping_mm_synonyms_filepath = "/data/genetic-mapping/ncbi/Mus_musculus.gene_info_symbol2ensembl.gz"
  mapping_hs_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
  
  # Load WGCNA results
  log_not_mapped_filepath = paste0(project_path,"log/",output_label,"_magma_wgcna_not_mapped_",study_label,"_",sNames,".tab")
  modulekME = read.csv(file=sprintf("%s%s_%s_%s_%s.csv", dir_tables, study_label, sNames, file_suffix, output_label), row.names=1)#, check.names = FALSE, sep = file_sep)
  
  # Load MAGMA genes and remap to Ensembl gene IDs
  d = dir(path=magma_gwas_path, pattern="[.]genes.out", recursive = T)
  gwas = vector(mode="list")
  for(i in 1:length(d)) {
    gwas[[i]] = read.table(paste(magma_gwas_path, d[[i]],sep=""),head=T, check.names = FALSE)
  }
  names(gwas) = gsub(".genes.out", "", d)
  
  # Match and replace with ENSG 
  #genes = union(gwas[[1]]$GENE, gwas[[2]]$GENE)
  #for(i in 3:length(gwas)) genes = union(genes, gwas[[i]]$GENE)
  
  # Remapping from human Entrez to human Ensembl gene IDs
  mapping_hs_entrez2ensembl = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T)
  for(i in 1:length(gwas)) {
    idx = match(gwas[[i]]$GENE, mapping_hs_entrez2ensembl$entrezgene)
    mapping = data.frame(entrez=gwas[[i]]$GENE, ensembl=mapping_hs_entrez2ensembl$ensembl_gene_id[idx])
    gwas[[i]]$gene_name = mapping$ensembl
  }
  
  # Remapping the WGCNA data to human ensembl IDs (using synonyms)
  # Step 1: direct mapping
  mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T)
  mapping = data.frame(symbol=row.names(modulekME), ensembl.mouse=mapping_direct$ensembl_gene_id[ match(row.names(modulekME), mapping_direct$gene_name_optimal) ])
  
  # Step 2: map remaing using synonyms
  mapping_synonyms = read.csv(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T)
  mapping$ensembl.mouse[ which(is.na(mapping$ensembl.mouse)) ] = mapping_synonyms$ensembl[ match( mapping$symbol[which(is.na(mapping$ensembl.mouse)) ] ,mapping_synonyms$symbol) ]
  
  # Step 3: orthology mapping
  #mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #mapping_mm_orthologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), mart=mart)
  mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T)
  mapping$ensembl.human = mapping_orthology$ensembl_gene_id[ match(mapping$ensembl.mouse,mapping_orthology$mmusculus_homolog_ensembl_gene) ]
  #mapping$ensembl.human[mapping$ensembl.human == ""] = NA
  df_not_mapped = mapping[is.na(mapping$ensembl.human),]
  write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
  modulekME$symbol = mapping$symbol
  modulekME$ensembl = mapping$ensembl.human
  modulekME = na.omit(modulekME)
  tmp = within(modulekME, rm("symbol","ensembl"))
  
  # Average duplicated gene IDs
  modulekME_ens <-aggregate(tmp, by=list(modulekME$ensembl),FUN=mean, na.rm=TRUE)
  rownames(modulekME_ens) = modulekME_ens$Group.1
  modulekME_ens = within(modulekME_ens, rm("Group.1"))
  
  # Calculate spearman's correlation between gene module membership and GWAS gene significance
  colors = colnames(modulekME_ens)
  table.kme.cor.p = table.kme.cor.r<- matrix(NA,nrow=length(unique(colors)),ncol=length(gwas))
  rownames(table.kme.cor.r) = rownames(table.kme.cor.p) = unique(colors)
  colnames(table.kme.cor.r) = colnames(table.kme.cor.p) = names(gwas) 
  
  for(m in unique(colors)) {
    for(i in 1:length(gwas)) {
      #col = paste("kME", m, sep="")
      col = m
      genes = intersect(rownames(modulekME_ens),gwas[[i]]$gene_name)
      x = -log10(gwas[[i]]$P[match(genes, gwas[[i]]$gene_name)])
      y = modulekME_ens[match(genes,rownames(modulekME_ens)), col]
      cor = cor.test(x,y,method="spearman", exact=F)
      table.kme.cor.r[m,i] = cor$estimate
      table.kme.cor.p[m,i] = cor$p.value
    }
  }
  
  table.kme.cor.p.fdr = p.adjust(table.kme.cor.p, method="fdr")
  dim(table.kme.cor.p.fdr) = dim(table.kme.cor.p);  dimnames(table.kme.cor.p.fdr) = dimnames(table.kme.cor.p)
  
  d = -log10(table.kme.cor.p.fdr) * sign(table.kme.cor.r) 
  #pdf("SampleGraph.pdf",width=7,height=5)
  #sizeGrWindow(9,7)
  #par(mfrow = c(2,2))
  #par(mar = c(4, 5, 4, 6));
  labeledHeatmap(d,textMatrix = signif(table.kme.cor.r,1), xLabels = colnames(d), yLabels = rownames(d),invertColors = T, colors = blueWhiteRed(1000), main="GWAS - kME correlation", cex.text = 0.6)
  #dev.off()
  
  #dat = as.data.frame(table.kme.cor.p.fdr*sign(table.kme.cor.r))[c("tan","blue","yellow","purple","turquoise","green","greenyellow","salmon"),]
  dat = as.data.frame(table.kme.cor.p.fdr*sign(table.kme.cor.r))
  dat[dat<0]=1 #Only look for positive enrichment
  dat = -log10(dat)
  dat$module = gsub("kME","",rownames(dat))
  #dat$module = gsub("ME","",rownames(dat))
  dat2 = melt(dat)
  dat2$variable=as.character(dat2$variable)
  
  #p=ggplot(melt(dat),aes(x=variable,y=value,fill=colors)) + 
  p=ggplot(melt(dat),aes(x=variable,y=value,fill="blue")) + 
    geom_bar(stat="identity",position=position_dodge(),color="black") +
    scale_fill_manual(values=sort(unique(dat2$module)))+ theme_classic() +
    geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
    theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
  
  p=ggplot(melt(dat),aes(x=variable,y=value,fill=module)) + 
    geom_bar(stat="identity",position=position_dodge(),color="black") +
    scale_fill_manual(values=sort(unique(dat2$module)))+ theme_classic() +
    geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
    theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
  
  ggsave(p, filename = paste0(figs_path,"/",output_label,"_wgcna_magma_", file_suffix,"_", study_label,"_",sNames,".pdf") ,width=45,height=12)
  
}


##########################################################################
######################## RUN THE PARALLEL LOOP ###########################
##########################################################################

disableWGCNAThreads()

cl <- makeCluster(n_cores, type = "FORK")

# Check the libraries are installed on all cores in cluster
clusterEvalQ(cl, library(Matrix))
clusterEvalQ(cl, library(WGCNA))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(STRINGdb))
clusterEvalQ(cl, library(reshape))
clusterEvalQ(cl, library(reshape2))
clusterEvalQ(cl, library(ggplot2))

message("Running WGCNA and MAGMA..")

parLapply(cl, sNames, parRWGCNA)
#lapply(sNames, function(x) parMagmaWGCNA(x, study_label=study_label, file_suffix=file_suffix, output_label=output_label))

stopCluster(cl)

##########################################################################
################################ FINISH ##################################
##########################################################################

message("Script DONE!")
