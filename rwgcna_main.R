# Title: Script to find robust WGCNA modules

######################################################################
############################## USAGE #################################
######################################################################

# e.g.

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-14/ --magma_gwas_path /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix tmp-15 --meta.data_ID NULL --min.cells 5 --do.center TRUE --genes_use PCA_5000 --corFnc cor --networkType signed --anti_cor_action NULL --minClusterSize "c(15,20)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.20,0.25)" --nPermutations 10 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted T --n_cores 5

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_cell_ontology_id.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-9/ --magma_gwas_path "/projects/jonatan/tmp-bmi-brain/data/magma/" --data_prefix tmp-9 --meta.data_ID NULL --min.cells 5 --do.center TRUE --genes_use PCA_5000 --corFnc cor --networkType signed --anti_cor_action NULL --minClusterSize "c(10,20)" --deepSplit "c(1,2)" --pamStage "c(TRUE,FALSE)" --moduleMergeCutHeight "c(0.20)" --nPermutations 10 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted T --n_cores 10

######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

# seurat_path = "/projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData"
# project_dir = "/projects/jonatan/tmp-rwgcna-tests/tmp-16/"
# magma_gwas_path = "/projects/jonatan/tmp-bmi-brain/data/magma/"
# data_prefix = "tmp-16"
# meta.data_ID = NULL
# do.center = T
# genes_use = "PCA_5000"
# corFnc = "cor"
# networkType = "signed"
# anti_cor_action = "kME_reassign"
# minClusterSize = c(20)
# deepSplit = c(1,2)
# pamStage = c(TRUE)
# moduleMergeCutHeight = c(0.2)
# replace = T
# nPermutations = 2
# STRINGdb_species = 10090
# ensembl_dataset = "mmusculus_gene_ensembl"
# plot_permuted = T
# n_cores = 5

######################################################################
########################### OptParse #################################
######################################################################

suppressMessages(library(optparse))

option_list <- list(
  
  make_option("--seurat_path", type="character",
              help = "Provide full path to Rdata input file with Seurat object"),
  make_option("--project_dir", type="character", default=NULL,
              help = "Optional. Provide project directory. Must have subdirs RObjects, plots, tables. If not provided, assumed to be dir one level up from input data dir."),
  make_option("--magma_gwas_dir", type="character", default = '/projects/jonatan/tmp-bmi-brain/data/magma/',
              help = "MAGMA input GWAS data directory as a character, defaults to '/projects/jonatan/tmp-bmi-brain/data/magma/'"),
  make_option("--meta.data_ID", type="character", default=NULL,
              help = "Specify the name of a seurat@meta.data$... column to use for subsetting the Seurat object. If NULL (default) uses the @ident slot."),
  make_option("--data_prefix", type="character", default=paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_rWGCNA_run"),
              help = "Dataset prefix for output files"),
  make_option("--min.cells", type="integer", default=5L,
              help="What is the minimum number of cells in each subset in the data in which a gene should be detected to not be filtered out? Integer, defaults to 5."),
  make_option("--do.center", type="logical", default=T,
              help="Use centered data? In either case data is scaled and nUMI and mitochrondrial genes are regressed out."),
  make_option("--genes_use", type="character", default="PCA_5000",
              help="One of 'all', 'var.genes', 'hvg_<number of highly variable genes>', 'PCA_<number of genes with high loading on at least one PC>'"), 
  make_option("--corFnc", type="character", default="cor",
              help="Use 'cor' for Pearson or 'bicor' for midweighted bicorrelation function."), 
  make_option("--networkType", type="character", default = "signed",
              help="'signed' scales correlations to [0:1]; 'unsigned' takes the absolute value (but the TOM can still be 'signed'); 'signed hybrid' sets negative correlations to zero."),
  make_option("--anti_cor_action", type="character", default=NULL, 
              help = "Optional. 'kME_reassign' reassigns genes with a negative kME more than 1.5 the kME w.r.t. their own (primary) module. Should be used only with networkType 'signed hybrid'."),
  make_option("--minClusterSize", type="character", default="c(15L)",
              help = "Minimum genes needed to form a module, recommended 5-25. Takes a character with a vector of integer values to try, defaults to 'c(15L)'."),
  make_option("--deepSplit", type="character", default="c(3L)",
              help = "Controls the sensitivity of the cutreeDynamic/cutreeHybrid algorithm. Takes a character with a vector of integer values between 0-4 to try, defaults to 'c(3L)'"),
  make_option("--moduleMergeCutHeight", type="character", default="c(0.2)",
              help = "Cut-off level for the variable (1-correlation) for merging eigengenes. Takes a character with a vector of double values to try, e.g. 'c(0.1, 0.2)'. Defaults to 'c(0.2)'"),
  make_option("--pamStage", type="character", default="c(TRUE)",
              help = "cutreeHybrid. Perform additional Partition Around Medroids step? Takes a character with a vector of logicals to try, e.g. 'c(TRUE,FALSE)', default 'c(TRUE)'"),
  make_option("--nPermutations", type="integer", default=100L,
              help = "Number of times to permute the dataset, defaults to 100"),
  make_option("--replace", type="logical", default=TRUE,
              help = "Sample with replacement? Defaults to TRUE. If TRUE, uses all samples, if FALSE, uses 66% each time."),
  make_option("--STRINGdb_species", type="integer", default="10090",
              help = "Optional: species for which to retrieve protein data from STRINGdb to validate clusters. 10090 is mouse(mus musculus), 9606 is human"),
  make_option("--ensembl_dataset", type="character", default="mmusculus_gene_ensembl",
              help = "Optional. Dataset for outputting colors with ensembl gene IDs for LD score regression."),
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
suppressMessages(library(parallel))
suppressMessages(library(WGCNA))
suppressMessages(library(STRINGdb))
suppressMessages(library(reshape))
suppressMessages(library(reshape2))
# suppressMessages(library(ggplot2))

#suppressMessages(library(biomaRt))

message("Libraries loaded")

######################################################################
################### GET COMMAND LINE OPTIONS #########################
######################################################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 

opt <- parse_args(OptionParser(option_list=option_list))

seurat_path <- opt$seurat_path 

project_dir <- opt$project_dir

magma_gwas_dir <- opt$magma_gwas_dir 

data_prefix <- opt$data_prefix 

meta.data_ID <- opt$meta.data_ID

min.cells <- opt$min.cells

do.center <- opt$do.center

genes_use <- opt$genes_use

corFnc <- opt$corFnc 

networkType <- opt$networkType

anti_cor_action <- opt$anti_cor_action

minClusterSize <- eval(parse(text=opt$minClusterSize))

deepSplit <- eval(parse(text=opt$deepSplit))

moduleMergeCutHeight <- eval(parse(text=opt$moduleMergeCutHeight))

pamStage <- eval(parse(text=opt$pamStage))

nPermutations <- opt$nPermutations

replace <- opt$replace

STRINGdb_species <- opt$STRINGdb_species

ensembl_dataset <- opt$ensembl_dataset

plot_permuted <- opt$plot_permuted

n_cores <- opt$n_cores  

######################################################################
############################## CONSTANTS #############################
######################################################################

if (!file.exists(seurat_path)) stop("Input data path not found")

# if no output directory provided, use the dir one level above that of the input file.

if (is.null(project_dir)) {
  pos <- tail(gregexpr("/", seurat_path)[[1]],2)[1]
  project_dir = substr(seurat_path, 1, pos)
}

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

saveIndividualTOMs = plot_permuted

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

# Load parameter values and utility functions
source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
source(file = "/projects/jonatan/functions-src/functions.R")

if (!is.null(ensembl_dataset)) {
  perslabMusEnsembl <- read.delim("/projects/jonatan/wgcna-src/rwgcna-pipeline/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
} 

if (nPermutations < 1) plot_permuted <- F 

# Set MAGMA variable values
file_suffix = if (!is.null(STRINGdb_species)) "kMEs_PPI" else "kMEs"

######################################################################
############################ CHECK INPUT #############################
######################################################################

if (min.cells < 0) stop("min.cells must be a non-negative integer") 
if (min.cells > 10 | min.cells < 3) warning("Recommended values for min.cells are between 3 and 10")

if (!(sapply(c("all", "var.genes", "hvg", "PCA"), function(x) grepl(x, genes_use, ignore.case=T)) %>% any())) {
  stop("genes_use must be one of 'all', 'var.genes', 'hvg_<number of highly variable genes>', 'PCA_<number of high loading genes>'")
}

if (!corFnc %in% c("cor", "bicor")) stop("corFnc must be one of 'cor' for Pearson's or 'bicor' for biweighted midcorrelation")

if (corFnc == "bicor" & do.center == T) warning("Using bicor with centered data will ignore non-positive values when detecting correlations")

if (!networkType %in% c('signed', 'unsigned', 'signed hybrid')) stop("networkType must be one of 'signed', 'unsigned' or 'signed hybrid' (not 'signed_hybrid')")

if (length(c(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight))>4) warning("Comparing different parameters increases processing time")

if (min(minClusterSize) < 5) stop("minClusterSize must be a vector of integers over 5")

if (max(deepSplit) > 4 | min(deepSplit) < 0) stop("deepSplit must be a vector of integers between 0 and 4")

if (min(moduleMergeCutHeight) < 0 | max(moduleMergeCutHeight) > 1) stop("moduleMergeCutHeight must be a vector of doubles between 0 and 1, recommended range is between 0.1 and 0.2")
  
if (! (nPermutations >= 0 & nPermutations <= 100)) stop("nPermutations must be in the range 0-100")

if (plot_permuted == T) warning("plotting modules on permuted datasets increases processing time")

if (! (n_cores >= 0 & n_cores <= 100)) stop("n_cores must be in the range 0-100")

######################################################################
################# LOAD AND SUBSET SEURAT OBJECT ######################
######################################################################

message("Loading and subsetting seurat object..")

seurat_obj <- load_obj(seurat_path)

if (!is.null(meta.data_ID)) {
  if (!any(grepl(meta.data_ID, names(seurat_obj@meta.data), ignore.case = T))) { 
    stop(sprintf("meta.data_ID not found in Seurat object at %s", seurat_path)) } else {
      seurat_obj <- SetAllIdent(seurat_obj, id = meta.data_ID)
    }
  }

sNames <- names(table(seurat_obj@ident))

subsets <- lapply(sNames, function(x) SubsetData(seurat_obj,
                                                 ident.use = x, 
                                                 do.scale = F, 
                                                 do.center = F,
                                                 subset.raw = T))
names(subsets) <- sNames 

for (subsetName in sNames) {
  seurat_obj_sub <- subsets[[subsetName]]
  save(seurat_obj_sub, file=sprintf("%s%s_%s_seurat_subset_tmp_%s.RData", RObjects_dir, data_prefix, subsetName, flag_date))
}

rm(subsets)

message("Seurat object loaded and subsetted")

######################################################################
##################### SAVE PARAMETERS TO FILE ########################
######################################################################

subsets_n_cells <- data.frame(subset = sNames, n_cells = as.numeric(table(seurat_obj@ident)))

# Free up memory
rm(seurat_obj)

userParams <- cbind(list("seurat_path" = seurat_path,
                         "project_dir" = project_dir,
                         "magma_gwas_dir" = magma_gwas_dir,
                         "data_prefix" = data_prefix,
                         "meta.data_ID" = meta.data_ID,
                         "min.cells" = min.cells,
                         "do.center" = do.center,
                         "genes_use" = genes_use,
                         "corFnc" = corFnc,
                         "networkType" = networkType,
                         "anti_cor_action" = anti_cor_action,
                         "minClusterSize" = minClusterSize,
                         "deepSplit" = deepSplit,
                         "pamStage" = pamStage,
                         "moduleMergeCutHeight" = moduleMergeCutHeight,
                         "nPermutations" = nPermutations,
                         "replace" = replace,
                         "STRINGdb_species" = STRINGdb_species,
                         "ensembl_dataset" = ensembl_dataset,
                         "plot_permuted" = plot_permuted,
                         "n_cores" = n_cores))

builtInParams <- cbind(list("minCoreKME" = minCoreKME,
                            "minKMEtoStay" = minKMEtoStay,
                            "consensusQuantile" = consensusQuantile,
                            "fraction" = fraction,
                            "hclustMethod" = hclustMethod,
                            "impute" = impute,
                            "nPC_seurat" = nPC_seurat,
                            "TOMType" = TOMType,
                            "PPI_pval_threshold" = PPI_pval_threshold))

# Save run params to file
write.csv(userParams, file=sprintf("%s%s_INFO_user_parameters_%s.csv", tables_dir, data_prefix, flag_date))
write.csv(builtInParams, file=sprintf("%s%s_INFO_built_in_parameters_%s.csv", tables_dir, data_prefix, flag_date))
write.csv(subsets_n_cells, file=sprintf("%s%s_INFO_subsets_n_cells_%s.csv", tables_dir, data_prefix, flag_date), row.names=F)

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
  seurat_obj_sub <- load_obj(sprintf("%s%s_%s_seurat_subset_tmp_%s.RData", RObjects_dir, data_prefix, sNames, flag_date))
  
  if (file.exists(sprintf("%s%s_%s_seurat_subset_tmp_%s.RData", RObjects_dir, data_prefix, sNames, flag_date))) {
    invisible(file.remove(sprintf("%s%s_%s_seurat_subset_tmp_%s.RData", RObjects_dir, data_prefix, sNames, flag_date)))
  }
  
  ######################################################################
  ################# PROCESS SEURAT SUBSET OBJECT #######################
  ######################################################################
  
  # Store metadata
  metaData <- seurat_obj_sub@meta.data
  
  # Filter out genes expressed in few cells
  if (min.cells > 0) {
    num.cells <- rowSums(seurat_obj_sub@data > 0)
    genes.use <- names(x = num.cells[which(x = num.cells >= min.cells)])
    seurat_obj_sub@raw.data <- seurat_obj_sub@raw.data[genes.use, ]
    seurat_obj_sub@data <- seurat_obj_sub@data[genes.use, ]
  }
  
  if (any(c("nUMI", "percent.mito") %in% names(seurat_obj_sub@meta.data))) {
    vars.to.filter_regress = c("nUMI", "percent.mito")[c("nUMI", "percent.mito") %in% names(seurat_obj_sub@meta.data)]
    #seurat_obj_sub <- FilterCells(object = seurat_obj_sub, subset.names = vars.to.filter_regress, low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.2)) 
  } else {
    warning("nUMI and percent.mito not found in Seurat object meta data")
  }
  
  
  if (genes_use != "all"){ #if we do use all genes, we'll rune ScaleData below
    
    # Need to run ScaleData to get variable genes for PCA / hvg / var.genes
    seurat_obj_sub <- ScaleData(object = seurat_obj_sub, 
                                    #genes.use = NULL, # default: genes.use = all genes in @data
                                  vars.to.regress = if (any(c("nUMI", "percent.mito") %in% names(seurat_obj_sub@meta.data))) vars.to.filter_regress else NULL, 
                                  model.use="linear",
                                  do.par=F,#F,#T, #T,
                                  #num.cores = n_cores,
                                  do.scale=T,
                                  do.center=T)
        
    seurat_obj_sub <- FindVariableGenes(object = seurat_obj_sub,
                                        x.low.cutoff = 0.0125,
                                        x.high.cutoff = 3,
                                        y.cutoff = 0.5, 
                                        do.plot=F)
    
    if (grepl("PCA", genes_use, ignore.case = T)) {
      seurat_obj_sub <- RunPCA(object = seurat_obj_sub,
                               pcs.compute = min(nPC_seurat, 
                                                 length(seurat_obj_sub@var.genes)-1, 
                                                 ncol(seurat_obj_sub@data)-1),
                               weight.by.var = F,#T,#F,
                               do.print = F,
                               seed.use = randomSeed,
                               maxit = maxit,
                               fastpath=fastpath)
      
      # Select significant PCs using empirical p-value based on JackStraw resampling
      # Source: https://rdrr.io/cran/Seurat/src/R/plotting.R
      # score the PCs using JackStraw resampling to get an empirical null distribution to get p-values for the PCs based on the p-values of gene loadings
      
      # seurat_obj_sub <- JackStraw(object = seurat_obj_sub,
      #                             num.pc = nPC_seurat,
      #                             num.replicate = 100,
      #                             display.progress = FALSE,
      #                             do.par = F)#T,
      #                             #num.cores = 10)
      # 
      # score.thresh = (5e-2)/nPC_seurat # TODO: Is this a suitable multiple testing adjustment 
      # 
      # pAll <- GetDimReduction(seurat_obj_sub, reduction.type = "pca", slot = "jackstraw")@emperical.p.value
      # pAll <- pAll[, 1:nPC_seurat, drop = FALSE]
      # pAll <- as.data.frame(pAll)
      # pAll$Contig <- rownames(x = pAll)
      # pAll.l <- melt(data = pAll, id.vars = "Contig")
      # colnames(x = pAll.l) <- c("Contig", "PC", "Value")
      # 
      # score.df <- NULL
      # 
      # for (i in (1:nPC_seurat)) {
      #   pc.score <- suppressWarnings(prop.test( # test if there ais
      #     x = c(length(x = which(x = pAll[, i] <= score.thresh)), floor(x = nrow(x = pAll) * score.thresh)),
      #     n = c(nrow(pAll), nrow(pAll)))$p.val)
      #   if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
      #     pc.score <- 1
      #   }
      #   if (is.null(x = score.df)) {
      #     score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
      #   } else {
      #     score.df <- rbind(score.df, data.frame(PC = paste0("PC",i), Score = pc.score))
      #   }
      # }
      # #PC_select <- score.df %>% filter(Score < score.thresh) %>% select(PC)
      # PC_select_idx <- which(score.df$Score < score.thresh)

      # Only keep gene loadings of significant PCs
      
    }

  }

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
      #seurat_obj_sub@dr$pca@gene.loadings.full <- seurat_obj_sub@dr$pca@gene.loadings.full[,PC_select_idx]
      max_loads <- apply(abs(seurat_obj_sub@dr$pca@gene.loadings.full), MARGIN=1, max, T)
      names_genes_use = names(sort(max_loads, decreasing=T)[1:n_genes])
      input_data <- seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,] 
    } else if (grepl("hvg", genes_use, ignore.case = T)) {
      n_genes <- as.numeric(gsub("[^0-9]", "", genes_use))
      names_genes_use = rownames(seurat_obj_sub@hvg.info)[order(seurat_obj_sub@hvg.info$gene.dispersion.scaled, decreasing=T)][1:n_genes]
      input_data <- seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,] 
    }
  }, error = function(c) {write.csv("ERROR", file=sprintf("%s%s_%s_input_data_selection_ERROR_%s.csv", project_dir, data_prefix, sNames, flag_date), row.names = F)
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
    
    powers = c(1:30)

    sft = pickSoftThreshold(data=datExpr,
                            powerVector = powers,
                            blockSize = min(2500, ncol(datExpr)), #try to prevent crashing
                            corFnc = corFnc,
                            corOptions =  corOptions,
                            networkType = networkType,
                            verbose = verbose)
    
    pdf(sprintf("%s%s_%s_pickSoftThresholdSFTFit_%s.pdf", plots_dir, data_prefix, sNames, flag_date),width=10,height=5)
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
    pdf(sprintf("%s%s_%s_pickSoftThresholdMeanCon_%s.pdf", plots_dir, data_prefix, sNames, flag_date),width=10,height=5)
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
    fitIndices_filter <- fitIndices %>% dplyr::filter(median.k. <= quantile(median.k.,0.95, na.rm=T))
    if (sum(fitIndices_filter$SFT.R.sq >= 0.9) > 1) {
      softPower = min(fitIndices_filter$Power[fitIndices_filter$SFT.R.sq>=0.9])
    } else {
      softPower <- (fitIndices_filter %>% dplyr::arrange(desc(SFT.R.sq)) %>% dplyr::select(Power))[1,]
    }
  }, error = function(c) {
    write.csv("ERROR", file=sprintf("%s%s_%s_pickSoftThreshold_ERROR_set_to_8_as_default_%s.csv", project_dir, data_prefix, sNames, flag_date), row.names = F)})
  
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
    
    setwd(RObjects_dir)
    
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
      cacheDir = RObjects_dir,
      cacheBase = ".blockConsModsCache",
      #nThreads = n_cores,
      verbose = verbose,
      indent = indent)
        
    # blockwiseConsensusModules has filtered out some genes when creating the consensus TOM.
    # We need to only work with these from here on
    goodGenesTOM_idx <- as.logical(consensus0$goodSamplesAndGenes$goodGenes)
    datExpr_filter <- datExpr[,goodGenesTOM_idx] # keep only the genes kept by consensusTOM
    
    # Load the consensus TOM matrix 
    suppressMessages(load(sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, sNames))) # Load the consensus TOM generated by blockwiseConsensusModules
    
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
  ########################### CUTREEHYBRID #############################
  ######################################################################

  # Iterate over different combinations of parameters
  # Parallelised version of the for loop used in Gandal et al 2018 to test different cluster params
  dims = sapply(list(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight), length) # list of number of values per parameter
  n_combs <- prod(dims) # Total number of combinations of parameter values
  comb_list = vector("list", length = n_combs)
  
  k=1
  for (minModSize in 1:dims[1]) {
    for (ds in 1:dims[2]) { # Gandal et al 2018 use c(50,100, 200)
      for (pam in 1:dims[3]) {
        for (dthresh in 1:dims[4]) {
          comb_list[[k]] <- c(minClusterSize[minModSize],deepSplit[ds], pamStage[pam], moduleMergeCutHeight[dthresh])
          k = k+1
        }
      }
    }
  }
  
  list_cutrees <- lapply(comb_list, function(x) parCutreeHybrid(comb=x, geneTree=geneTree, dissTOM = dissTOM)) # Cut out different color assignments in the dendrogran
  
  list_merged <- mapply(function(x,y) parMergeCloseModules(cutree=x,comb=y, datExpr_filter=datExpr_filter), list_cutrees, comb_list, SIMPLIFY=FALSE) # Merge modules whose eigengenes are highly correlated
  
  list_MEs = lapply(list_merged, function(x) x$MEGs) # list of merged color eigengenes
  
  labels <- lapply(comb_list, function(x) parLabel(comb=x)) # list of labels for plot
  
  # Get the colors and put them in a list
  list_colors = lapply(list_merged, function(x) extract_and_name_colors(merged=x, datExpr_filter=datExpr_filter)) # list of merged colors
  
  if (length(list_colors)>1) {
    # Match colors between iterations by reference to the set of colors with the highest number of unique colors
    max_colors_idx = which.max(lapply(list_MEs, function(x) length(x$eigengenes)))
  
    diffParams_colors = matrix(0, nrow=ncol(datExpr_filter), ncol=n_combs)
    diffParams_colors[, max_colors_idx] = list_colors[[max_colors_idx]] # use the set of colors with the highest number of modules as 'reference'
  
    for (j in (1:n_combs)[-max_colors_idx]) {
      diffParams_colors[,j] <- matchLabels(source = list_colors[[j]], # traverse gene colour assignments
                              reference = diffParams_colors[,max_colors_idx],
                              pThreshold = 5e-2,
                              na.rm = TRUE,
                              extraLabels = standardColors())#paste0("Extra_", 0:500)) # labels for modules in source that cannot be matched to any modules in reference
  
    }
    
    # Split the matrix of matched colors back into a list 
    list_colors_matched <- split(diffParams_colors, rep(1:ncol(diffParams_colors), each = nrow(diffParams_colors))) 
    
  } else {
    diffParams_colors <- NULL
    list_colors_matched <- list_colors
  }  
    
  names(list_colors_matched) <- labels
  
  for (j in 1:length(list_colors_matched)) {
    names(list_colors_matched[[j]]) <- colnames(datExpr_filter)
  }
  
  # Compute kMEs and pkMEs for each set of colors corresponding to distinct parameters
  list_kMEs <- vector(mode="list", length=length(list_MEs))
  list_pkMEs <- vector(mode="list", length=length(list_MEs))
  
  for (j in 1:length(list_MEs)) {
    # Rename eigengene columns according to new color names adopted above 
    names(list_MEs[[j]]$eigengenes) <- paste0("ME",sort(names(table(list_colors_matched[j])), decreasing = F))
    
    # get kMEs for each set of eigengenes
    list_kMEs[[j]] = signedKME(as.matrix(datExpr_filter),
                     list_MEs[[j]]$eigengenes,
                     outputColumnName="",
                     corFnc = corFnc)
    
    # Get the 'principal kMEs', i.e. kME of each gene to the module to which it is allocated
    pkMEs <- vector(mode="numeric",length=nrow(diffParams_colors))
    
    # Loop over every gene and get its pkME
    for (i in 1:length(pkMEs)) {
      #pkMEs[i] <- list_kMEs[[j]][diffParams_colors[i,j]][i,]
      pkMEs[i] <- list_kMEs[[j]][list_colors_matched[[j]][i]][i,]
    }
    
    list_pkMEs[[j]] <- pkMEs 
    
    # Delete grey kMEs
    for (k in 1:length(list_kMEs)) {
      if (any(grepl("grey", colnames(list_kMEs[[k]])))) list_kMEs[[k]][['grey']] <- NULL
    }
  }
  
  ######################################################################
  #################### REASSIGN GENES BASED ON KMEs ####################
  ######################################################################
  
  if (anti_cor_action == "kME_reassign") {
    
    list_list_kME_reassign_in <- list(list_MEs, list_kMEs, list_pkMEs, list_colors)
    
    tryCatch({
      list_kME_reassign_out <- lapply(1:length(list_kMEs), 
                                      FUN= function(i) kME_reassign_fnc(MEs = list_list_kME_reassign_in[[1]][[i]],
                                                                        kMEs = list_list_kME_reassign_in[[2]][[i]], 
                                                                        pkMEs = list_list_kME_reassign_in[[3]][[i]], 
                                                                        kME_reassign_threshold = kME_reassign_threshold, 
                                                                        colors = list_list_kME_reassign_in[[4]][[i]], 
                                                                        datExpr_filter = datExpr_filter, 
                                                                        corFnc = corFnc))
      
      list_colors = lapply(list_kME_reassign_out, function(x) x$colors)
      list_MEs = lapply(list_kME_reassign_out, function(x) x$MEs)
      list_kMEs = lapply(list_kME_reassign_out, function(x) x$kMEs)
      list_pkMEs = lapply(list_kME_reassign_out, function(x) x$pkMEs)
      list_reassigned = lapply(list_kME_reassign_out, function(x) x$reassigned)
      
    }, error = function(c) {
      write.csv("ERROR", file=sprintf("%s%s_%s_kME_reassign_ERROR_%s.csv", project_dir, data_prefix, sNames, flag_date), row.names = F)
    })
    
    # Free up memory
    rm(list_kME_reassign_out)
    
  }
  
  names(list_pkMEs) <- labels
  
  ######################################################################
  #################### CHECK MODULES FOR PPI ENRICHMENT ################
  ######################################################################
  
  list_colors_PPI <- mapply(function(x,y) PPI_for_par_outer(pkMEs = x,
                                                            colors = y,
                                                            STRINGdb_species = STRINGdb_species,
                                                            PPI_pval_threshold = PPI_pval_threshold,
                                                            project_dir = project_dir, 
                                                            data_prefix = data_prefix, 
                                                            subsetName = sNames, 
                                                            flag_date = flag_date),
                            list_pkMEs,
                            list_colors_matched,
                            SIMPLIFY = F)
  
  # If we have more than one set of colourings, count the number of assigned genes per set and rank the sets by it
  
  if (length(list_colors_PPI) > 1) {
    
    n_grey <- sapply(list_colors_PPI, function(x) sum(x=="grey"), simplify=T)
    
    list_colors_PPI_order <- list_colors_PPI[order(n_grey, decreasing=F)]
    labels_order <- labels[order(n_grey, decreasing=F)]
    MEs_order <- list_MEs[order(n_grey, decreasing=F)]
    
    pdf(sprintf("%s%s_%s_diffParams_colors_PPI_ordered_%s.pdf", plots_dir, data_prefix, sNames, flag_date),width=8,height=6+length(list_colors_PPI)/2)
    plotDendroAndColors(geneTree, 
                        #labels2colors(list_colors_PPI_order),
                        #matrix(unlist(list_colors_PPI_order), ncol = length(list_colors_PPI_order), byrow = F), 
                        diffParams_colors,
                        groupLabels=labels_order, 
                        addGuide= TRUE, 
                        dendroLabels=FALSE, 
                        main="Diff. params modules, PPI enriched, ranked by n assigned genes.", 
                        cex.colorLabels=0.5)
    dev.off()
    
    # Select a set of preferred colors and add them as the first column to be shown at the top of the figure
    colors <- list_colors_matched[order(n_grey, decreasing=F)][[1]]
    colors_PPI <- list_colors_PPI_order[[1]]
    labels <- labels_order[[1]]
    MEs <- MEs_order[[1]]
    cutree_params_final = as.list(comb_list[order(n_grey,decreasing = F)][[1]]) 
 
  } else {
    colors <- list_colors[[1]]
    colors_PPI <- list_colors_PPI[[1]]
    labels = labels[[1]]
    MEs = list_MEs[[1]]
    cutree_params_final = as.list(comb_list[[1]])
  }
  
  names(cutree_params_final) = c("minClusterSize", "deepSplit","pamStage", "moduleMergeCutHeight")
  
  ######################################################################
  ############# RECOMPUTE MEs, kMEs, pkMEs AFTER PPI FILTER ############
  ######################################################################
  
  # TODO better to filter
  # recompute Module Eigengenes
  
  if (unique(colors_PPI)!=1) {
    
    MEs_PPI = moduleEigengenes(expr = as.matrix(datExpr_filter),
                               colors_PPI,
                               excludeGrey = F) # T leads to errors when finding pkMEs
    #softPower = softPower)
    # recompute kMEs
    kMEs_PPI = signedKME(as.matrix(datExpr_filter),
                         MEs_PPI$eigengenes,
                         outputColumnName = "",
                         corFnc = corFnc)
    
    # recompute primary kMEs
    pkMEs_PPI <- vector(length=length(colors_PPI))
    
    for (i in 1:length(colors_PPI)) {
      pkMEs_PPI[i] <- kMEs_PPI[colors_PPI[i]][i,]#as.numeric(unlist(kMEs_PPI %>% dplyr::select(matches(colors_PPI[[i]]))))[[i]]
    }
    
    #names(pkMEs_PPI) <- names(colors_PPI) # don't use them or write then to csv..
    
    # Delete grey kME (needed it for pkME step though)
    if (any(grepl("grey", colnames(kMEs_PPI)))) kMEs_PPI[['grey']] <- NULL 
    
  } else {
    write.csv("WARNING" , file=sprintf("%s%s_%s_WARNING_sum_genes_PPI_0_%s.csv", project_dir, data_prefix, sNames, flag_date))
  }
  
  
  ######################################################################
  ######################### PLOT FINAL COLORS #########################
  ######################################################################
  
  tryCatch({
    
    all_colors = if (!is.null(STRINGdb_species) & !is.null(colors_PPI)) cbind(colors, colors_PPI) else colors
    labels[[2]] <- "PPI"
    
    pdf(sprintf("%s%s_%s_colors_%s.pdf", plots_dir, data_prefix, sNames, flag_date),width=10,height=5)
    plotDendroAndColors(geneTree,
                          all_colors,
                          groupLabels = if (!is.null(compare_params)) labels else NULL,
                          addGuide=T,
                          dendroLabels=F,
                          main="Final colors",
                          cex.colorLabels=0.3)
    dev.off()
      
    if (plot_permuted==T) {  
      
      # NOTE: Alternatively use sampledBlockwiseModules() - see /projects/jonatan/wgcna-src/180228_robust_wgcna_module_detection.Rmd
      bootstrap_seq <- if (replace==T) 1:(nPermutations+1) else 1:nPermutations
      
      bootstrap_indiv_TOMs <- lapply(bootstrap_seq, function(x) load_obj(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", RObjects_dir, sNames, x, 1))) # retrieve permuted dataset TOMs from disk
      
      # Delete individual TOMs from disk
      for (i in bootstrap_seq) {
        if (file.exists(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", RObjects_dir, sNames, i, 1))) {
          file.remove(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", RObjects_dir, sNames, i, 1))
        }
      }
      
      bootstrap_dissTOMs <- lapply(bootstrap_indiv_TOMs, function(x) 1-as.dist(x)) # Convert permuted TOMs to distance matrices
      
      bootstrap_geneTrees <- lapply(bootstrap_dissTOMs, function(x) hclust(d=x, method=hclustMethod)) # Do hierarchical clustering 
      
      # cut out modules in each gene tree
      bootstrap_cutrees <- mapply(function(x,y) cutreeHybrid(dendro = x, distM = as.matrix(y), cutHeight = NULL, minClusterSize = cutree_params_final$minClusterSize, deepSplit = cutree_params_final$deepSplit, pamStage = cutree_params_final$pamStage, pamRespectsDendro = cutree_params_final$pamStage), bootstrap_geneTrees, bootstrap_dissTOMs, SIMPLIFY=F) 
      
      # merge modules with correlated eigengenes
      bootstrap_merged <- mapply(function(x,y) mergeCloseModules(exprData = as.matrix(x$data)[,goodGenesTOM_idx], colors = y$labels, cutHeight=cutree_params_final$moduleMergeCutHeight), multiExpr, bootstrap_cutrees, SIMPLIFY=FALSE) # Merge modules whose eigengenes are highly correlated
      
      bootstrap_colors = lapply(bootstrap_merged, function(x) labels2colors(x$colors)) 
      
      bootstrap_labels <- c("Original data", paste0("Permutation ", as.character(bootstrap_seq)[-max(bootstrap_seq)]))
      
      # Free up memory
      rm(bootstrap_indiv_TOMs)
      rm(bootstrap_dissTOMs)
      
      # Align the module coloring schemes 
      bootstrap_colors_final = matrix(0, ncol(datExpr_filter), nPermutations + 1)
      bootstrap_colors_final[, 1] = colors # use the consensus colors as 'reference'
      
      for (r in 2:(nPermutations+1))
      {
        bootstrap_colors_final[, r] = matchLabels(source = bootstrap_colors[[r-1]], # traverse gene colour assignments
                                  reference = bootstrap_colors_final[, 1],
                                  pThreshold = 5e-2,
                                  na.rm = TRUE,
                                  extraLabels = standardColors())# paste0("Extra_", 0:100)) # labels for modules in source that cannot be matched to any modules in reference
      }

      bootstrap_labels_final <- c("Consensus", bootstrap_labels)

      pdf(sprintf("%s%s_%s_diffbootstrap_colors_%s.pdf", plots_dir, data_prefix, sNames, flag_date),width=8,height=8)
      plotDendroAndColors(geneTree, 
                          labels2colors(bootstrap_colors_final),
                          #matrix(unlist(bootstrap_colors_final), ncol = length(bootstrap_labels_final), byrow = F), 
                          groupLabels=bootstrap_labels_final, 
                          addGuide= TRUE, 
                          dendroLabels=FALSE, 
                          main="Modules: consensus and on each resampling", 
                          cex.colorLabels=0.5)
      dev.off()
    }
  }, error = function(c) write.csv("ERROR", file=sprintf("%s%s_%s_plot_colors_ERROR_%s.csv", project_dir, data_prefix, sNames, flag_date), row.names = F))
  
  # clear up
  if (nPermutations>0) rm(multiExpr)
  
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
      write.csv("", file=sprintf("%s%s_%s_%s_%s.csv", tables_dir, data_prefix, sNames, prop_not_map, flag_date))
      
      if (!is.null(colors_PPI)) {
        colors_PPI_noGrey <- colors_PPI[colors_PPI != grey]
        ensembl_PPI_IDs <- perslabMusEnsembl$ensembl_gene_id[match(names(colors_PPI_noGrey), perslabMusEnsembl$gene_name_optimal)] 
        ensembl_PPI_out <- data.frame(module_name_PPI = colors_PPI_noGrey, ensembl_ID = ensembl_PPI_IDs, row.names = NULL) 
        prop_not_map_PPI <- paste0("EnsemblID_PPI_proportion_not_mapped:_", signif(as.double(sum(is.na(ensembl_PPI_IDs))) / as.double(1e-10+length(ensembl_PPI_IDs)),2))
        write.csv("", file=sprintf("%s%s_%s_%s_%s.csv", tables_dir, data_prefix, sNames, prop_not_map_PPI, flag_date))
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
    
    
  }, error = function(c) write.csv("ERROR", file=sprintf("%s%s_%s_ensemblID_ERROR_%s.csv", project_dir, data_prefix, sNames, flag_date), row.names = F))

  ######################################################################
  ################### SAVE RUN OBJECTS, RESULTS AND PARAMS #############
  ######################################################################

  
  tryCatch({
    
    kMEs <- cbind(genes=names(colors), kMEs)
    pkMEs <- cbind(colors=colors, genes=names(colors), pkMEs = pkMEs)
    
    # save kMEs and kMEs primary to file
    write.csv(kMEs, file=sprintf("%s%s_%s_kMEs_%s.csv", tables_dir, data_prefix, sNames, flag_date), row.names=F)
    write.csv(pkMEs, file=sprintf("%s%s_%s_pkMEs_%s.csv", tables_dir, data_prefix, sNames, flag_date), row.names = F)
    
    if (!is.null(kMEs_PPI)) {
      kMEs_PPI <- cbind(genes=names(colors_PPI), kMEs_PPI) 
      pkMEs_PPI <- data.frame(colors_PPI = colors_PPI, genes = names(colors_PPI), pkMEs_PPI = pkMEs_PPI) 
      write.csv(kMEs_PPI, file=sprintf("%s%s_%s_kMEs_PPI_%s.csv", tables_dir, data_prefix, sNames, flag_date), row.names=F)
      write.csv(pkMEs_PPI, file=sprintf("%s%s_%s_pkMEs_PPI_%s.csv", tables_dir, data_prefix, sNames, flag_date), row.names=F)
    }
    
    # ensembl IDs for LD score regression
    if (!is.null(ensembl_out)) {
      write.csv(ensembl_out, file=sprintf("%s%s_%s_ensembl_out_%s.csv", tables_dir, data_prefix, sNames, flag_date), row.names=F)
      if (!is.null(ensembl_PPI_out)) write.csv(ensembl_PPI_out, file=sprintf("%s%s_%s_ensembl_PPI_out_%s.csv", tables_dir, data_prefix, sNames, flag_date), row.names=F)
    }   

    save.image(file=sprintf("%s%s_%s_rsession_%s.RData", RObjects_dir, data_prefix, sNames, flag_date))
    
  }, error = function(c) write.csv("ERROR", file=sprintf("%s%s_%s_save_results_ERROR_%s.csv", project_dir, data_prefix, sNames, flag_date), row.names = F))
}
  
##########################################################################
######################## RUN THE WGCNA LOOP ##############################
##########################################################################

cl <- makeCluster(n_cores, type = "FORK")

# Check the libraries are installed on all cores in cluster
clusterEvalQ(cl, library(Matrix))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(WGCNA))
clusterEvalQ(cl, library(STRINGdb))

message("Running WGCNA")
parLapply(cl, sNames, parRWGCNA)

stopCluster(cl)

##########################################################################
######################## RUN THE MAGMA LOOP ##############################
##########################################################################

message("WGCNA all done, running MAGMA..")

#cl <- makeCluster(n_cores, type = "FORK")

#clusterEvalQ(cl, library(Matrix))
#clusterEvalQ(cl, library(plyr))
#clusterEvalQ(cl, library(WGCNA))
#clusterEvalQ(cl, library(reshape))
#clusterEvalQ(cl, library(reshape2))
#clusterEvalQ(cl, library(ggplot2))

#parLapply(cl, sNames, function(x) parMagma(subsetName = x, 
#                                   project_dir = project_dir, 
#                                   plots_dir = plots_dir, 
#                                   log_dir = log_dir,
#                                   tables_dir = tables_dir,
#                                   magma_gwas_dir = magma_gwas_dir, 
#                                   data_prefix = data_prefix, 
#                                   file_suffix = file_suffix, 
#                                   flag_date = flag_date))

lapply(sNames, function(x) parMagma(subsetName = x, 
                                   project_dir = project_dir, 
                                   plots_dir = plots_dir, 
                                   log_dir = log_dir,
                                   tables_dir = tables_dir,
                                   magma_gwas_dir = magma_gwas_dir, 
                                   data_prefix = data_prefix, 
                                   file_suffix = file_suffix, 
                                   flag_date = flag_date))

#stopCluster(cl)

##########################################################################
################################ FINISH ##################################
##########################################################################

message("Script DONE!")
