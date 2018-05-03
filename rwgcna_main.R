# Title: Script to find robust WGCNA modules

######################################################################
############################## USAGE #################################
######################################################################

# e.g.

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_n05_to_n10.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell_n05_to_n10-9/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-n05_to_n10-9 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.20)" --num.replicate 10 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 6
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_neurons_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-neurons-12/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-neurons-12 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.1,0.20)" --num.replicate 250 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 14
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_s_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-glia-12/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-glia-12 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.1,0.20)" --num.replicate 250 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 11
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-8/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-AgRP-8 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-9/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-AgRP-9 --min.cells 3 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 4

# ===================RUNNING=========================
#  time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_pancreas.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-pancreas-6/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix maca-pancreas-6 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 7
#  - time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-17/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-AgRP-17 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 4
# tmux 0 - time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-19/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-AgRP-19 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 350 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 4
# tmux 0:3 - time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-20/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-AgRP-20 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 4
# ===================================================

# Campbell AgRP 6 - var.genes, num.replicate 0
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-6/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-AgRP-6 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(1,2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 4

# tmux 8 - RUNNING - neurons 16, var.genes, num.replicate 0
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_neurons_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-neurons-16/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-neurons-16 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(1,2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 7

# tmux 9 - RUNNING - glia 16, var.genes, num.replicate 0
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_s_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-glia-16/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-glia-16 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(1,2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 6
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_s_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-glia-16/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-glia-16 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(1,2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --num.replicate 0 --nPermutations 0 --replace T --STRINGdb_species 10090 --ensembl_dataset mmusculus_gene_ensembl --plot_permuted F --n_cores 6


######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

if (FALSE) { 
  seurat_path = "/projects/jonatan/tmp-holst-hsl/RObjects/campbell_n05_to_n10.RData"
  project_dir = "/projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-12/"
  magma_gwas_dir = "/projects/jonatan/tmp-bmi-brain/data/magma/"
  data_prefix = "campbell-AgRP-12"
  resume = NULL 
  meta.data_subset_col = NULL
  meta.data_corr_col = NULL #c("nUMI", "nGene", "percent.mito", "percent.ribo")
  use.imputed = F
  min.cells = 5
  do.center = T
  genes_use = "PCA"
  pca_genes = 'all'
  corFnc = "cor"
  networkType = "signed"
  anti_cor_action = NULL
  minClusterSize = c(20)
  deepSplit = c(2)
  pamStage = c(TRUE)
  moduleMergeCutHeight = c(0.2)
  replace = T
  num.replicate = 0
  nPermutations = 0
  STRINGdb_species = 10090
  ensembl_dataset = "mmusculus_gene_ensembl"
  plot_permuted = F
  n_cores = 6
}

######################################################################
########################### OptParse #################################
######################################################################

suppressMessages(library(optparse))

option_list <- list(
  
  make_option("--seurat_path", type="character",
              help = "Provide full path to Rdata input file with Seurat object with lognormalized expression data in the @data slot"),
  make_option("--project_dir", type="character", default=NULL,
              help = "Optional. Provide project directory. Must have subdirs RObjects, plots, tables. If not provided, assumed to be dir one level up from input data dir."),
  make_option("--magma_gwas_dir", type="character", default = '/projects/jonatan/tmp-bmi-brain/data/magma/',
              help = "MAGMA input GWAS data directory as a character, defaults to '/projects/jonatan/tmp-bmi-brain/data/magma/'"),
  make_option("--data_prefix", type="character", default=paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_rWGCNA_run"),
              help = "Dataset prefix for output files"),
  make_option("--resume", type="character", default=NULL,
              help = "Resume from a previous session image? Must have same path and data_prefix. Options are 'checkpoint_1' - 'checkpoint_7'"),
  make_option("--meta.data_subset_col", type="character", default=NULL,
              help = "Specify a seurat@meta.data$... column to use for subsetting the Seurat object. If NULL (default) uses the @ident slot."),
  make_option("--meta.data_corr_col", type="character", default=NULL,
              help = "Specify seurat@meta.data$... column(s) for which to compute correlations with gene modules. Takes a character with a vector of meta.data column names e.g. 'nUMI' or 'c('nUMI', 'AGE')'. For factor or character metadata, each levels is analysed as a dummy variable, so exercise caution. Defaults to NULL."),
  make_option("--use.imputed", type="logical", default=F,
              help="Use data in the obj@imputed slot for the computations to replace the @data slot? If the @imputed slot is empty, will revert to the default (FALSE)"),
  make_option("--min.cells", type="integer", default=5L,
              help="What is the minimum number of cells in each subset in the data in which a gene should be detected to not be filtered out? Integer, defaults to 5."),
  make_option("--do.center", type="logical", default=T,
              help="Use centered data? In either case data is scaled and nUMI and mitochrondrial genes are regressed out."),
  make_option("--genes_use", type="character", default="PCA",
              help="One of 'all', 'var.genes' or 'PCA' for genes with significant loading on at least one significant PC. 'All' is not recommended. Defaults to 'PCA'"), 
  make_option("--pca_genes", type="character", default="var.genes",
              help="'all' or 'var.genes'. 'All' is computationally expensive but allows for selecting genes based on PC loading p-values rather than magnitudes as with the default (var.genes)"), 
  make_option("--corFnc", type="character", default="cor",
              help="Use 'cor' for Pearson or 'bicor' for midweighted bicorrelation function."), 
  make_option("--networkType", type="character", default = "signed",
              help="'signed' scales correlations to [0:1]; 'unsigned' takes the absolute value (but the TOM can still be 'signed'); 'signed hybrid' sets negative correlations to zero."),
  make_option("--anti_cor_action", type="character", default=NULL, 
              help = "Optional. 'kME_reassign' reassigns genes with a negative kME more than 1.5 the kME w.r.t. their own (primary) module. Should be used only with networkType 'signed hybrid'."),
  make_option("--minClusterSize", type="character", default="15L",
              help = "Minimum genes needed to form a module, recommended 5-25. Takes a character with a vector of integer values to try 'c(15,20,25)', defaults to 'c(15L)'."),
  make_option("--deepSplit", type="character", default="3L",
              help = "Controls the sensitivity of the cutreeDynamic/cutreeHybrid algorithm. Takes a character with a vector of integer values between 0-4 to try, e.g. 'c(1,2,3)', defaults to '3L'"),
  make_option("--moduleMergeCutHeight", type="character", default="c(0.2)",
              help = "Cut-off level for the variable (1-correlation) for merging eigengenes. Takes a character with a vector of double values to try, e.g. 'c(0.1, 0.2)'. Defaults to '0.2'"),
  make_option("--pamStage", type="character", default="c(TRUE)",
              help = "cutreeHybrid. Perform additional Partition Around Medroids step? Takes a character with a vector of logicals to try, e.g. 'c(TRUE,FALSE)', default 'c(TRUE)'"),
  make_option("--num.replicate", type="integer", default=300L,
              help = "Number of times to resample to make null distributions for empirical significance tests in JackStraw and other functions. Non-negative integer, 0 to skip Jackstraw, defaults to 300."),
  make_option("--nPermutations", type="integer", default=100L,
              help = "Number of times to resample the dataset, defaults to 100"),
  make_option("--replace", type="logical", default=TRUE,
              help = "Sample with replacement? Defaults to TRUE. If TRUE, uses all samples, if FALSE, uses 66% each time."),
  make_option("--STRINGdb_species", type="integer", default="10090",
              help = "Optional. Species for which to retrieve protein data from STRINGdb to validate clusters. Mouse (mus musculus) is 10090, homo sapiens is 9606"),
  make_option("--ensembl_dataset", type="character", default=NULL,
              help = "Optional. Dataset for mapping gene symbols to ensembl_ID. Mouse is mmusculus_gene_ensembl, homo sapiens is hsapiens_gene_ensembl"),
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
suppressMessages(library(gProfileR))
suppressMessages(library(reshape))
suppressMessages(library(reshape2))
suppressMessages(library(corrplot))
suppressMessages(library(pheatmap))

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

resume <- opt$resume

meta.data_subset_col <- opt$meta.data_subset_col

meta.data_corr_col <- opt$meta.data_corr_col

if (!is.null(meta.data_corr_col)) {
  meta.data_corr_col <- eval(parse(text=meta.data_corr_col))
}

use.imputed <- opt$use.imputed

min.cells <- opt$min.cells

do.center <- opt$do.center

genes_use <- opt$genes_use

pca_genes <- opt$pca_genes

corFnc <- opt$corFnc 

networkType <- opt$networkType

anti_cor_action <- opt$anti_cor_action

minClusterSize <- eval(parse(text=opt$minClusterSize))

deepSplit <- eval(parse(text=opt$deepSplit))

moduleMergeCutHeight <- eval(parse(text=opt$moduleMergeCutHeight))

pamStage <- eval(parse(text=opt$pamStage))

num.replicate <- opt$num.replicate

nPermutations <- opt$nPermutations

replace <- opt$replace

STRINGdb_species <- opt$STRINGdb_species

ensembl_dataset <- opt$ensembl_dataset

plot_permuted <- opt$plot_permuted

n_cores <- opt$n_cores  

######################################################################
############################## CONSTANTS #############################
######################################################################

if (!file.exists(seurat_path) & is.null(resume)) stop("Input data path not found and no previous session to resume given")
    
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

if (!is.null(ensembl_dataset)) if (ensembl_dataset == "hsapiens_gene_ensembl") {
  perslabEnsembl <- read.delim("/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz")
} else if (ensembl_dataset == "mmusculus_gene_ensembl") {
  perslabEnsembl <- read.delim("/projects/jonatan/wgcna-src/rwgcna-pipeline/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz") 
}

# Gprofiler
if (!is.null(ensembl_dataset)) {
    if (ensembl_dataset == "hsapiens_gene_ensembl") {
    organism = "hsapiens" 
  } else if (ensembl_dataset == "mmusculus_gene_ensembl") {
    organism = "mmusculus"
  } else {
    organism = NULL
  }
} else if (is.null(ensembl_dataset)) {
  organism <- NULL
}

if (nPermutations < 1) plot_permuted <- F 

######################################################################
############################ CHECK INPUT #############################
######################################################################

if (min.cells < 0) stop("min.cells must be a non-negative integer") 
if (min.cells > 25 | min.cells < 5) warning("Recommended values for min.cells between 5 and 10")

if (!(sapply(c("all", "var.genes", "PCA"), function(x) grepl(x, genes_use, ignore.case=T)) %>% any())) {
  stop("genes_use must be one of 'all', 'var.genes' or 'PCA'")
}

if (!pca_genes %in% c('all', 'var.genes')) {
  warning("Invalid pca_genes argument, reverted to var.genes")
  pca_genes <- var.genes
}

if (!corFnc %in% c("cor", "bicor")) stop("corFnc must be one of 'cor' for Pearson's or 'bicor' for biweighted midcorrelation")

if (corFnc == "bicor" & do.center == T) warning("Using bicor with centered data will ignore non-positive values when detecting correlations")

if (!networkType %in% c('signed', 'unsigned', 'signed hybrid')) stop("networkType must be one of 'signed', 'unsigned' or 'signed hybrid' (not 'signed_hybrid')")

if (length(c(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight))>4) warning("Comparing different parameters increases processing time")

if (min(minClusterSize) < 5) stop("minClusterSize must be a vector of integers over 5")

if (max(deepSplit) > 4 | min(deepSplit) < 0) stop("deepSplit must be a vector of integers between 0 and 4")

if (min(moduleMergeCutHeight) < 0 | max(moduleMergeCutHeight) > 1) stop("moduleMergeCutHeight must be a vector of doubles between 0 and 1, recommended range is between 0.1 and 0.2")

if (! (num.replicate >= 0)) stop("num.replicate must be 0 or higher")

if (! (nPermutations >= 0)) stop("nPermutations must be non-negative")

if (plot_permuted == T) warning("plotting modules on permuted datasets increases processing time")

if (! (n_cores >= 0 & n_cores <= 100)) stop("n_cores must be in the range 0-100")

######################################################################
################# LOAD AND SUBSET SEURAT OBJECT ######################
######################################################################

if (is.null(resume)) {
  message("Loading and subsetting seurat object..")
  
  seurat_obj <- load_obj(f=seurat_path)
  
  # if (!is.null(meta.data_subset_col)) {
  #   if (!any(grepl(meta.data_subset_col, names(seurat_obj@meta.data), ignore.case = T))) {
  #     stop(sprintf("meta.data_subset_col not found in Seurat object at %s", seurat_path)) } else {
  #       seurat_obj <- SetAllIdent(seurat_obj, id = meta.data_subset_col)
  #     }
  #   }
  
  if (is.null(seurat_obj@imputed) | any(dim(seurat_obj@imputed)) == FALSE) use.imputed <- F
  if (use.imputed==T) seurat_obj@data <- seurat_obj@imputed
  
  sNames <- names(table(seurat_obj@ident))
  
  all_genes <- rownames(seurat_obj@data)
  
  # Convert any character or factor meta.data to numeric dummy variables each level with its own numeric column
  
  if (!is.null(meta.data_corr_col)) {
    #meta.data_corr_col_numeric <- character()
    #for (column in meta.data_corr_col) {
    if (all(sapply(meta.data_corr_col, function(x) !is.null(seurat_obj@meta.data[[x]])))) {
      meta.data <- SeuratFactorToIndicator(obj = seurat_obj,
                                            meta.colnames = meta.data_corr_col)
    } else {
      message('Did not find all columns in meta.data_corr_col in the Seurat object meta.data slot')
      meta.data <- NULL
    }
  }
    #if (!is.null(seurat_obj@meta.data[[column]])) if (class(seurat_obj@meta.data[[column]]) %in% c("character", "factor")) {
  
      
        #meta.data_corr_col_numeric <- c(meta.data_corr_col_numeric,levels(seurat_obj@meta.data[[column]]))
    # } else if (class(seurat_obj@meta.data[[column]]) %in% c("integer", "numeric")) {
    #     meta.data_corr_col_numeric <- c(meta.data_corr_col_numeric,column)
    #   }
    # }
    # Take the updated meta data columns from the seurat object (so we can remove it)
    #meta.data <- sapply(meta.data_corr_col_numeric, function(x) seurat_obj@meta.data[[x]], simplify = T)
  
  #}
  
  
  subsets <- lapply(sNames, function(x) SubsetData(seurat_obj,
                                                   ident.use = x, 
                                                   do.scale = F, 
                                                   do.center = F,
                                                   subset.raw = if (!is.null(seurat_obj@raw.data)) T else F))
  
  names(subsets) <- sNames 
  
  message("Seurat object loaded and subsetted")
  
  ######################################################################
  ##################### SAVE PARAMETERS TO FILE ########################
  ######################################################################
  
  subsets_n_cells <- data.frame(subset = sNames, n_cells = as.numeric(table(seurat_obj@ident)))
  
  userParams <- cbind(list("seurat_path" = seurat_path,
                           "project_dir" = project_dir,
                           "magma_gwas_dir" = magma_gwas_dir,
                           "data_prefix" = data_prefix,
                           "meta.data_subset_col" = meta.data_subset_col,
                           "meta.data_corr_col" = meta.data_corr_col,
                           "use.imputed" = use.imputed,
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
                           "num.replicate" = num.replicate,
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
  
  ######################################################################
  ##################### SEURAT PROCESSING ##############################
  ######################################################################
  
  ### 180502_v1.7_1
  subsets <- lapply(subsets, function(x) FilterGenes(x, min.cells = min.cells))
  
  ### 180504_v1.7_3
  subsets <- lapply(subsets, function(x) FilterCells(object = x,
                                                   subset.names = c("nGene", "percent.mito"),
                                                   low.thresholds = c(200,-Inf),
                                                   high.thresholds = c(7500, 0.2)))
  
  vars.to.regress = c("nUMI", "percent.mito", "percent.ribo")[c("nUMI", "percent.mito", "percent.ribo") %in% names(seurat_obj@meta.data)]

  subsets <- lapply(subsets, function(x) ScaleData(object = x,
                                                  vars.to.regress = if (length(vars.to.regress) > 0) vars.to.regress else NULL,
                                                  model.use="linear",
                                                  do.par=T,
                                                  num.cores = min(n_cores, detectCores()-1),
                                                  do.scale=T,
                                                  do.center=do.center))
  ###
  
  
  ### 180502_v1.7_2
  # subsets <- lapply(subsets, wrapFilterCellsScaleData)
  # subsets <- lapply(subsets, function(x) wrapFilterCellsScaleData(x, 
  #                                                                 n_cores = n_cores, 
  #                                                                 do.center=do.center)) 
  
  ### 180504_v2.7_4
  # subsets <- lapply(subsets, function(x) FindVariableGenes(object = x,
  #                                                          # x.low.cutoff = 0.0125,
  #                                                          # x.high.cutoff = 3,
  #                                                          # y.cutoff = 0.5,
  #                                                          do.plot=F))
  subsets <- lapply(subsets, function(x) FindVariableGenes(object = x,
                                                           x.low.cutoff = 0.0125,
                                                           x.high.cutoff = 8,
                                                           y.cutoff = 1,
                                                           do.plot=F,
                                                           display.progress=F))
  ###
  invisible(gc())

  if (grepl("PCA", genes_use, ignore.case = T)) {

    message("Performing Principal Component Analysis..")
    
    start_time <- Sys.time()
  
    cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "log_parPCA.txt"))
    ### EDIT_180501_1
    #subsets <- parLapply(cl, subsets, function(x) parPCA(x, nPC_seurat = nPC_seurat)) # if use.imputed=T the @imputed slot has been copied to @data
    #subsets <- parLapply(cl, subsets, function(x) parPCAfull(x, nPC_seurat = nPC_seurat)) # if use.imputed=T the @imputed slot has been copied to @data
    
    tryCatch({
      subsets <- parLapply(cl, subsets, function(x) RunPCA(object = x,
                                                           pc.genes = if (pca_genes == 'all') rownames(x@data) else x@var.genes,
                                                           pcs.compute = min(nPC_seurat, (if (pca_genes == 'all') nrow(x@data) else length(x@var.genes)) %/% 2, ncol(x@data) %/% 2),
                                                           use.imputed = F, # if use.imputed=T the @imputed slot has been copied to @data
                                                           weight.by.var = F,
                                                           do.print = F,
                                                           seed.use = randomSeed,
                                                           maxit = maxit, # set to 500 as default
                                                           fastpath = fastpath)) 
      
    }, warning = function(c) {
      subsets <- parLapply(cl, subsets, function(x) RunPCA(object = x,
                                                           pc.genes = if (pca_genes == 'all') rownames(x@data) else x@var.genes,
                                                           pcs.compute = min(nPC_seurat, (if (pca_genes == 'all') nrow(x@data) else length(x@var.genes)) %/% 2, ncol(x@data) %/% 2),
                                                           use.imputed = F, # if use.imputed=T the @imputed slot has been copied to @data
                                                           weight.by.var = F,
                                                           do.print = F,
                                                           seed.use = randomSeed,
                                                           maxit = maxit*2, # set to 500 as default
                                                           fastpath = F))
    })
    ### 
    stopCluster(cl)
    invisible(gc())
    end_time <- Sys.time()
   
    message(sprintf("PCA done, time elapsed: %s seconds", round(end_time - start_time,2)))
    
    # Select significant PCs using empirical p-value based on JackStraw resampling
    # Source: https://rdrr.io/cran/Seurat/src/R/plotting.R
    # score the PCs using JackStraw resampling to get an empirical null distribution to get p-values for the PCs based on the p-values of gene loadings
    
    #message(sprintf("Performing JackStraw with %s replications to select genes that load on significant PCs", num.replicate))
    list_datExpr <- lapply(subsets, function(x) wrapJackStraw(seurat_obj_sub = x, n_cores = n_cores, num.replicate = num.replicate))
    
    invisible(gc())
    
  } else if (genes_use == "all") {
    list_datExpr <- lapply(subsets, function(x) seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rep(TRUE, nrow(x@scale.data))))
  } else if (genes_use == "var.genes") {
    list_datExpr <- lapply(subsets, function(x) seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rownames(x@scale.data) %in% x@var.genes))
  } 
  
  # just in case, since JackStraw is slooow
  save.image(file=sprintf("%s%s_checkpoint_1_image.RData", RObjects_dir, data_prefix))
} else if (!is.null(resume)) {
  if (resume == "checkpoint_1") {
    load((file=sprintf("%s%s_checkpoint_1_image.RData", RObjects_dir, data_prefix)))
    source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
    source(file = "/projects/jonatan/functions-src/functions_v1.7.R")
  }
}

### 180503_v1.7
# Better to do this later so we can run MAGMA on the original
# if (!is.null(organism)) {
#   ######################################################################
#   ######################### MAP to ENSEMBL IDs #########################
#   ######################################################################
#   
#   # The whole dataset
#   all_genes_ensembl_IDs <- perslabEnsembl$ensembl_gene_id[match(all_genes, perslabEnsembl$gene_name_optimal)]
#   # Where the mapping failed, revert to the previous name
#   all_genes_IDs_mix <- ifelse(is.na(all_genes_ensembl_IDs), yes = all_genes, no = all_genes_ensembl_IDs)
#   # Each subset
#   list_ensembl_IDs <- lapply(list_datExpr, function(x) perslabEnsembl$ensembl_gene_id[match(colnames(x), perslabEnsembl$gene_name_optimal)]) 
#   list_idx_mapped <- lapply(list_ensembl_IDs, function(x) !is.na(x))
#   list_prop_not_map <- lapply(list_idx_mapped, function(x) signif(as.double(sum(!x)) / as.double(1e-10+length(x)),2))
#   # Where the mapping failed, revert to the previous name
#   list_ensembl_IDs_mix <- mapply(FUN = function(x,y) ifelse(is.na(x), yes = colnames(y), no = x), x=list_ensembl_IDs, y=list_datExpr, SIMPLIFY=F)
# 
#   # Rename the expression matrix colmns
#   list_datExpr <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = 2),  x=list_datExpr,  y=list_ensembl_IDs_mix, SIMPLIFY=F)
# }
    
###
# Clear up
### 180503 v1.7
#rm(perslabEnsembl, subsets, userParams, builtInParams, subsets_n_cells, seurat_obj)
rm(subsets, userParams, builtInParams, subsets_n_cells, seurat_obj)
###
invisible(gc())
  
  ######################################################################
  ####### PICK SOFT THRESHOLD POWER FOR ADJACENCY MATRIX ###############
  ######################################################################
  message("Computing soft threshold powers to maximise the fit of a scale free topology to the adjacency matrix")
  
### 180503 v1.7
  # invisible(gc())
  # cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_sft_for_par.txt"))
  # list_softPower <- clusterMap(cl, function(x,y) sft_for_par(datExpr=x, subsetName=y), list_datExpr, sNames, SIMPLIFY=F)
  # stopCluster(cl)
  # invisible(gc())

  list_softPower <- mapply(function(x,y) sft_for_par(datExpr=x, subsetName=y), list_datExpr, sNames, SIMPLIFY=F)

  ###
  
  names(list_softPower) <- sNames
  
  # TOM files will be saved here by consensusTOM / blockwiseConsensusModules
  setwd(RObjects_dir)
  
  if (nPermutations > 0) {
    
    ######################################################################
    ############### RESAMPLE THE DATA FOR ROBUSTNESS #####################
    ######################################################################
    
    message(sprintf("Resampling the data %s times with replace = %s", nPermutations, replace))
    
    invisible(gc())
    
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_bootStrap.txt"))
    
    list_multiExpr <- parLapply(cl, list_datExpr, function(x) bootstrap(datExpr=x, 
                                                                        nPermutations = nPermutations,
                                                                        replace = replace,
                                                                        fraction = fraction,
                                                                        randomSeed = randomSeed))
    stopCluster(cl)
    
    invisible(gc())
    
    names(list_multiExpr) <- sNames 
  
    ######################################################################
    ######### RUN CONSENSUSTOM ON THE RESAMPLED DATASETS #################
    ######################################################################
    
    message(sprintf("Computing the consensus Topological Overlap Matrix with %s permutations", nPermutations))
  
    invisible(gc())
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_consensusTOM.txt"))
    list_consensus <- clusterMap(cl, function(x,y)  consensusTOM(multiExpr = x, 
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
                                                                power = list_softPower[[y]],
                                                                networkType = networkType,
                                                                TOMDenom = TOMDenom,
                                                                saveIndividualTOMs = saveIndividualTOMs,
                                                                individualTOMFileNames = paste0(y, "_individualTOM-Set%s-Block%b.RData"),
                                                                networkCalibration = networkCalibration,
                                                                sampleForCalibration = sampleForCalibration,
                                                                sampleForCalibrationFactor = sampleForCalibrationFactor,
                                                                getNetworkCalibrationSamples = getNetworkCalibrationSamples,
                                                                consensusQuantile = consensusQuantile,
                                                                useMean = useMean,
                                                                saveConsensusTOMs = saveConsensusTOMs,
                                                                consensusTOMFilePattern = paste0(y,"_consensusTOM-block.%b.RData"),
                                                                returnTOMs = F,
                                                                useDiskCache = T,
                                                                cacheDir = RObjects_dir,
                                                                cacheBase = ".blockConsModsCache",
                                                                verbose = verbose,
                                                                indent = indent), 
                                 x=list_multiExpr, 
                                 y=sNames, 
                                 SIMPLIFY=F)
    stopCluster(cl)
    invisible(gc())
    
    # For each consensusTOM, get a logical vector where 'good_genes' that were used are TRUE
    list_goodGenesTOM_idx <- lapply(list_consensus, function(x) as.logical(x$goodSamplesAndGenes$goodGenes))
    
    # Use the goodGenesTOM_idx to filter the datExpr matrices
    #list_datExpr_filter <- mapply(FUN=parDatExpr_filter, list_datExpr, list_goodGenesTOM_idx, SIMPLIFY=F)
    list_datExpr_gg <- mapply(FUN=function(x,y) x[,y], x=list_datExpr, y=list_goodGenesTOM_idx, SIMPLIFY=F)
    
    ### 180503 v1.7
    #list_idx_mapped_gg <- mapply(FUN=function(x,y) x[y], x=list_idx_mapped, y=list_goodGenesTOM_idx, SIMPLIFY=F)
    ###
    
      #   resume="checkpoint_2"
      #   
      # }, error = function(c) {
      #   save.image(file=sprintf("%s%s_checkpoint_1_image.RData", RObjects_dir, data_prefix))
      #   stop(sprintf("ConsensusTOM computation failed for %s",sNames))
        
        # warning(sprintf("ConsensusTOM computation failed for %s, computing the TOM on the unpermuted data instead",sNames))
        # 
        # try({
        #   rm(list_multiExpr)
        #   stopCluster(cl)
        #   invisible(gc())
        # }, silent = T)
        # 
        # invisible(gc())
        # cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_TOM_for_par.txt"))
        # list_goodGenesTOM_idx <- clusterMap(cl, function(x,y) TOM_for_par(datExpr=x, subsetName=y, softPower=list_softPower[[y]]), x=list_datExpr, y=sNames, SIMPLIFY=F)
        # stopCluster(cl)
        # 
        # invisible(gc())
        # 
        # list_datExpr_gg <- list_datExpr
        # list_idx_mapped_gg <- list_idx_mapped
      # })
    
    # just in case
    save.image(file=sprintf("%s%s_checkpoint_2_image.RData", RObjects_dir, data_prefix))
    
  } else if (nPermutations==0) {
    
    message("Computing the Topological Overlap Matrix")
    
    ######################################################################
    ##### COMPUTE THE ADJACENCY MATRIX AND TOM WITHOUT RESAMPLING ########
    ######################################################################
    
    invisible(gc())
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_TOM_for_par.txt"))
    list_goodGenesTOM_idx <- clusterMap(cl, function(x,y) TOM_for_par(datExpr=x, subsetName=y, softPower=list_softPower[[y]]), x=list_datExpr, y=sNames, SIMPLIFY=F)
    stopCluster(cl)
    invisible(gc())
    
    list_datExpr_gg <- list_datExpr
    ### 180503 v1.7
    #list_idx_mapped_gg <- list_idx_mapped
    ###
}

message("Done computing the consensus Topological Overlap Matrix")
  

# resume = "checkpoint_2"
# message("Reached checkpoint 2, saving session image")
# save.image(file=sprintf("%s%s_checkpoint_2_image.RData", RObjects_dir, data_prefix))

# } else if (resume == "checkpoint_2") {
#   load(file=sprintf("%s%s_checkpoint_2_image.RData", RObjects_dir, data_prefix))
#   # Load parameter values and utility functions anew (in case a bug was fixed)
#   source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
#   source(file = "/projects/jonatan/functions-src/functions.R")
# }


if (nPermutations > 0) rm(list_multiExpr) #TODO: we'll need it for plotting permuted
invisible(gc())


######################################################################
####################### CLUSTER ON THE TOM ###########################
######################################################################

# Convert TOM to distance matrix
list_dissTOM <- lapply(sNames, dissTOM_for_par)

# Cluster
invisible(gc())
cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parHclust.txt"))
list_geneTree <- parLapply(cl, list_dissTOM, function(x) hclust(d=x, method=hclustMethod))
stopCluster(cl)

names(list_geneTree) = sNames # used for PlotDendro

######################################################################
######################### CUT OUT MODULES ############################
######################################################################
# Iterate over different combinations of parameters
# NB: the following section uses vectorised computation nested within a parallel loop
# Credit: Gandal et al 2018 

message("Computing modules on the Topological Overlap Matrix")

dims = sapply(list(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight), length) # list of number of values per parameter
n_combs <- prod(dims) # Total number of combinations of parameter values
comb_list = vector("list", length = n_combs)
  
# Make a vector of different parameter combinations
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

invisible(gc())
cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_cutreeHybrid_for_vec.txt"))
# nb: comb_list is in the environment - no need to pass it through clusterMap
list_list_cutree <- clusterMap(cl, function(x,y) lapply(comb_list, function(z) cutreeHybrid_for_vec(comb=z, geneTree=x, dissTOM = y)), x=list_geneTree, y=list_dissTOM, SIMPLIFY=F)
stopCluster(cl)

names_list_list_cutree = sNames

# Merge close modules
invisible(gc())
cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parvecMergeCloseModules.txt"))
list_list_merged <- parLapply(cl, sNames, function(x) lapply(1:n_combs, function(y) mergeCloseModules_for_vec(cutree=list_list_cutree[[x]][[y]], comb=comb_list[[y]], datExpr=list_datExpr_gg[[x]], excludeGrey = excludeGrey)))
stopCluster(cl)
invisible(gc())

names(list_list_merged) = sNames

# Produce labels for plotting the modules found by different parameters
list_plot_label <- lapply(comb_list, function(x) plotLabel_for_vec(comb=x)) # list of labels for plot
# Make a copy for each subset
list_list_plot_label <- list()
list_list_plot_label <- lapply(sNames, function(x) list_list_plot_label$x = list_plot_label)

names(list_list_plot_label) = sNames

# Extract the Module Eigengenes from the list returned by mergeCloseModules
list_list_MEs <- lapply(list_list_merged, parGetMEs) # nested lapply

names(list_list_MEs) = sNames

# Extract the colors from the list returned by mergeCloseModules
list_list_colors <- mapply(parGetColors, list_list_merged, list_datExpr_gg, SIMPLIFY=F)
names(list_list_colors) <- sNames

# Test if for a given vector all the genes are grey
# count the grey
list_vec_n_grey <- lapply(list_list_colors, count_grey_in_list_of_vec)
# For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, we're ok.

#list_idx_ok <- mapply(function(x,y) as.logical(mapply(function(a,b) which(a!=length(b)), a=x, b=y, SIMPLIFY=F)), x=list_vec_n_grey, y=list_list_colors_matched, SIMPLIFY=F)
list_idx_ok <- mapply(function(x,y) as.logical(mapply(function(a,b) a!=length(b), a=x, b=y, SIMPLIFY=F)), x=list_vec_n_grey, y=list_list_colors, SIMPLIFY=F)
subsets_ok <- sapply(list_idx_ok, any, simplify = T)

# First, for each run, remove results of any parametrisation for which all the genes were assigned to grey
# If for a subset all parameters gave only grey modules, take it out of the top level list.
#if (!all(unlist(list_idx_ok), recursive = T)) {
list_list_plot_label_filter <- mapply(function(x,y) x[y], x = list_list_plot_label, y = list_idx_ok, SIMPLIFY = F) %>% Filter(f=length)
list_list_cutree_filter <- mapply(function(x,y) x[y], x = list_list_cutree, y = list_idx_ok, SIMPLIFY = F) %>% Filter(f=length)
list_list_MEs_filter <- mapply(function(x,y) x[y], x = list_list_MEs, y = list_idx_ok, SIMPLIFY = F) %>% Filter(f=length)
list_list_colors_filter <- mapply(function(x,y) x[y], x=list_list_colors, y=list_idx_ok, SIMPLIFY = F) %>% Filter(f=length)

sNames_filter <- sNames[subsets_ok]
# list_datExpr_filter <- list_datExpr_gg[names(list_datExpr_gg) %in% names(list_list_colors_filter)] # If for a subset all parameters gave only grey modules, take it out of the top level list.
# list_geneTree_filter <- list_geneTree[names(list_geneTree) %in% names(list_datExpr_gg)]
# list_list_plot_label_filter <- list_list_plot_label[names(list_list_plot_label) %in% names(list_list_colors_filter)]
list_datExpr_filter <- list_datExpr_gg[subsets_ok] # If for a subset all parameters gave only grey modules, take it out of the top level list.
list_geneTree_filter <- list_geneTree[subsets_ok]
#list_list_plot_label_filter <- list_list_plot_label[subsets_ok] # we do this above already

### 180503 v1.7
#list_idx_mapped_gg_filter <- list_idx_mapped_gg[subsets_ok]
###

list_dissTOM_filter <- list_dissTOM[subsets_ok]
#}

# Assign gene names to each color vector
list_list_colors_filter <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), dimension = NULL)), x = list_list_colors_filter, y = list_datExpr_filter, SIMPLIFY=F)

# Make a warning and report if any whole subsets produced no modules
discarded <- setdiff(sNames, sNames_filter)
if (length(discarded > 0)) {
  warning(sprintf("The %s were dropped because the script didn't find any modules", discarded))
  fileConn<-file(sprintf("%s%s_%s_ERROR_report_%s.txt", log_dir,  data_prefix, discarded, flag_date ))
  writeLines(sprintf("The %s were dropped because the script didn't find any modules", discarded), fileConn)
  close(fileConn)
}

# clear up
### 180503 v1.7
rm(list_datExpr, list_list_cutree, list_list_MEs, list_datExpr_gg, list_geneTree, list_list_plot_label, list_dissTOM)#, list_idx_mapped_gg)
###
invisible(gc())
  
######################################################################
###################### compute kMEs & pkMEs ##########################
######################################################################

message("Computing kMEs: connection (correlation) between each gene and each module 'eigengene', the primary principal component of a module's gene expression")
# Compute kMEs
invisible(gc())
cl <- makeCluster(spec=n_cores, type = "FORK", outfile = paste0(log_dir, "log_parkMEs_parPkMEs.txt"))
list_list_kMEs_filter <- clusterMap(cl, function(x,y) parkMEs(list_MEs=x, datExpr=y), x=list_list_MEs_filter, y=list_datExpr_filter, SIMPLIFY = F)
#list_list_kMEs_filter <- mapply(function(x,y) parkMEs(list_MEs=x, datExpr=y), x=list_list_MEs_filter, y=list_datExpr_filter, SIMPLIFY = F)
# Extract the primary kMEs - i.e. those of each gene w.r.t. the module to which it belongs
list_list_pkMEs_filter <- clusterMap(cl, function(x,y) parPkMEs(list_kMEs=x, list_colors=y), x=list_list_kMEs_filter, y=list_list_colors_filter, SIMPLIFY = F)
#list_list_pkMEs_filter <- mapply(function(x,y) parPkMEs(list_kMEs=x, list_colors=y), x=list_list_kMEs_filter, y=list_list_colors_filter, SIMPLIFY = F)
stopCluster(cl)
invisible(gc())

######################################################################
########################## compute kIMs ##############################
######################################################################

message("Computing normalised Intra-module connectivity, kIM: each gene's degree with respect to each module")

invisible(gc())

cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_kIMs.txt"))

list_list_kIMs_filter <- clusterMap(cl, function(x,y) lapply(y, function(z) kIM_eachMod_norm(dissTOM = x, 
                                                                              colors = z,
                                                                              genes = names(z))),
                     x = list_dissTOM_filter,
                     y = list_list_colors_filter, 
                     SIMPLIFY=F)

stopCluster(cl)
invisible(gc())


######################################################################
############## Match colors between parameter settings ###############
######################################################################

message("Matching module color labels between parameter settings")

# Align / match the colors so similar modules found with different sets of parameters have the same name
invisible(gc())
cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parMatchColors.txt"))
list_list_colors_filter_matched <- parLapply(cl, list_list_colors_filter, parMatchColors)
stopCluster(cl)

names(list_list_colors_filter_matched) = sNames_filter
# Name each entry of the vectors of color assignments with the corresponding genes
list_list_colors_filter_matched <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), dimension = NULL)), x=list_list_colors_filter_matched, y=list_datExpr_filter, SIMPLIFY=F)

# Rename the eigengenes according to the matched color scheme - nope, not necessary since we won't be using them
#list_list_MEs_matched_filter <- mapply(function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x$eigengenes, given_names= sort(names(table(y)))), x=a, y=b, SIMPLIFY=F), a= list_list_MEs_filter, b = list_list_colors_matched_filter, SIMPLIFY=F)

# Rename the columns of the kMEs and kIMs according to the matched color scheme
list_list_kMEs_filter_matched <- mapply(function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x, given_names = sort(names(table(y))), dimension = 2), x = a, y = b, SIMPLIFY=F), a = list_list_kMEs_filter, b = list_list_colors_filter_matched, SIMPLIFY=F)
list_list_kIMs_filter_matched <- mapply(function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x, given_names = sort(names(table(y))), dimension = 2), x = a, y = b, SIMPLIFY=F), a = list_list_kIMs_filter, b = list_list_colors_filter_matched, SIMPLIFY=F)


### EDIT_180430_2
# Get rid of any 'fake' colors

# list_list_colors_filter_matched <- lapply(list_list_colors_filter_matched, function(x) lapply(x, function(y) replace_unplottable_colors(y)))
# 
# list_list_kMEs_filter_matched <- mapply(function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(table(y)), dimension = 2), 
#                                                              x = a, y = b, SIMPLIFY = F), 
#                                         a = list_list_kMEs_filter_matched, b = list_list_colors_filter_matched, SIMPLIFY = F)
# 
# list_list_kIMs_filter_matched <- mapply(function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(table(y)), dimension = 2), 
#                                                              x = a, y = b, SIMPLIFY = F), 
#                                         a = list_list_kIMs_filter_matched, b = list_list_colors_filter_matched, SIMPLIFY = F)
# ###

### EDIT_180428_1
# # Delete grey kMEs and grey kIMs
# list_list_kMEs_filter_matched_noGrey <- lapply(list_list_kMEs_filter_matched, deleteGrey)
# list_list_kIMs_filter_matched_noGrey <- lapply(list_list_kIMs_filter_matched, deleteGrey)
###

### EDIT_180428_1
# Name
# names(list_list_kMEs_filter_matched_noGrey) <- sNames_filter
# list_list_kMEs_filter_matched_noGrey <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_kMEs_filter_matched_noGrey, y=list_list_plot_label_filter, SIMPLIFY=F)
# names(list_list_kIMs_filter_matched_noGrey) <- sNames_filter
# list_list_kIMs_filter_matched_noGrey <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_kIMs_filter_matched_noGrey, y=list_list_plot_label_filter, SIMPLIFY=F)
names(list_list_kMEs_filter_matched) <- sNames_filter
list_list_kMEs_filter_matched <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_kMEs_filter_matched, y=list_list_plot_label_filter, SIMPLIFY=F)
names(list_list_kIMs_filter_matched) <- sNames_filter
list_list_kIMs_filter_matched <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_kIMs_filter_matched, y=list_list_plot_label_filter, SIMPLIFY=F)
###

# Name primary (p-) kMEs
names(list_list_pkMEs_filter) <- sNames_filter
list_list_pkMEs_filter <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_pkMEs_filter, y=list_list_plot_label_filter, SIMPLIFY=F)

# clear up 
rm(list_list_colors_filter)
invisible(gc())

if (!is.null(STRINGdb_species)) {
  ######################################################################
  #################### CHECK MODULES FOR PPI ENRICHMENT ################
  ######################################################################
  message("Checking modules for significant Protein-Protein Interactions through STRINGdb")

  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_PPI_outer_for_vec.txt"))
  list_list_colors_PPI <- clusterMap(cl, function(a,b) mapply(function(x,y) PPI_outer_for_vec(colors = x,
                                                                                              pkMEs = y,
                                                                                              STRINGdb_species = STRINGdb_species, 
                                                                                              PPI_pval_threshold = PPI_pval_threshold, 
                                                                                              project_dir = project_dir, 
                                                                                              data_prefix = data_prefix, 
                                                                                              flag_date = flag_date), 
                                                              x = a, y = b, SIMPLIFY=F), 
                                     a = list_list_colors_filter_matched, 
                                     b = list_list_pkMEs_filter, 
                                     SIMPLIFY=F)
  stopCluster(cl)
  invisible(gc())







  
  ######################################################################
  ########### ORDER PARAMETER SETS BY PPI ENRICHMENT ###################
  ######################################################################
  
  message("Selecting parameters with the highest number of genes assigned to modules significantly enriched for Protein-Protein Interactions")
    
  # Count how many genes were assigned to grey under each parametrisation and order the parametrisations
  list_vec_n_grey <- lapply(list_list_colors_PPI, count_grey_in_list_of_vec)
  
  # Order all the outputs by how many genes were assigned to a (non-grey) module
  list_list_plot_label_filter_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_plot_label_filter, y = list_vec_n_grey, SIMPLIFY=F)
  names(list_list_plot_label_filter_order) = sNames_filter
  
  list_list_colors_filter_matched_order <- mapply(function(x,y) x[order(y, decreasing=F)], x  =  list_list_colors_filter_matched , y = list_vec_n_grey, SIMPLIFY = F )
  names(list_list_colors_filter_matched_order) = sNames_filter
  
  list_list_MEs_filter_order <- mapply(function(x,y) x[order(y, decreasing=F)], x =  list_list_MEs_filter, y = list_vec_n_grey, SIMPLIFY=F)
  names(list_list_MEs_filter_order) = sNames_filter
  
  ### EDIT_180428_1
  # list_list_kMEs_filter_matched_noGrey_order <- mapply(function(x,y) x[order(y, decreasing=F)],  x = list_list_kMEs_filter_matched_noGrey, y = list_vec_n_grey, SIMPLIFY=F)
  # names(list_list_kMEs_filter_matched_noGrey_order) = sNames_filter
  # list_list_kIMs_filter_matched_noGrey_order <- mapply(function(x,y) x[order(y, decreasing=F)],  x = list_list_kIMs_filter_matched_noGrey, y = list_vec_n_grey, SIMPLIFY=F)
  # names(list_list_kIMs_filter_matched_noGrey_order) = sNames_filter
  list_list_kMEs_filter_matched_order <- mapply(function(x,y) x[order(y, decreasing=F)],  x = list_list_kMEs_filter_matched, y = list_vec_n_grey, SIMPLIFY=F)
  names(list_list_kMEs_filter_matched_order) = sNames_filter
  list_list_kIMs_filter_matched_order <- mapply(function(x,y) x[order(y, decreasing=F)],  x = list_list_kIMs_filter_matched, y = list_vec_n_grey, SIMPLIFY=F)
  names(list_list_kIMs_filter_matched_order) = sNames_filter
  
  ###
  list_list_pkMEs_filter_order <- mapply(function(x,y) x[order(y, decreasing=F)],  x  = list_list_pkMEs_filter, y = list_vec_n_grey, SIMPLIFY=F)
  names(list_list_pkMEs_filter_order) = sNames_filter
  
  list_list_colors_PPI_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_colors_PPI, y = list_vec_n_grey, SIMPLIFY=F)
  names(list_list_colors_PPI_order) = sNames_filter

  ### EDIT_180430 remove any '.' from color names to avoid plotting error:
  # Error in rect(xl, yb, xr, yt, col = as.character(C[, j]), border = as.character(C[,  :  invalid color name 'darkorange.1' 
  # list_list_colors_filter_matched_order <- lapply(list_list_colors_filter_matched_order, function(x) lapply(x, function(y) gsub(pattern="\\.", replacement="", y)))
  # list_list_colors_PPI_order <- lapply(list_list_colors_PPI_order, function(x) lapply(x, function(y) gsub(pattern="\\.", replacement="", y)))
  ###
  
  # For each subset, plot the ranked colorings
  # Pre-PPI test
  invisible(mapply(function(x,y) parPlotDiffParams(list_colors=x, subsetName=y), x=list_list_colors_filter_matched_order, y=sNames_filter, SIMPLIFY = F))
  # Post-PPI test
  invisible(mapply(function(x,y) parPlotDiffParams_PPI(list_colors=x, subsetName=y), x=list_list_colors_PPI_order, y=sNames_filter, SIMPLIFY = F))

  ######################################################################
  ##### FOR EACH SUBSET SELECT PARAMETERS WITH BEST PPI ENRICHMENT #####
  ######################################################################

  # In each subset, select PPI enriched colors for the best parametrisation. 
  # Eliminate a layer of nesting.
  list_plot_label_final <- lapply(list_list_plot_label_filter_order, function(x) x[[1]])
  names(list_plot_label_final) <- sNames_filter
  
  list_colors_PPI_final <- lapply(list_list_colors_PPI_order, function(x) x[[1]])
  list_colors_PPI_final <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = colnames(y), dimension = NULL), 
                                  x = list_colors_PPI_final, 
                                  y = list_datExpr_filter, 
                                  SIMPLIFY=F)

  names(list_colors_PPI_final) <- sNames_filter
  
  # Also select corresponding colors, MEs, kMEs, kIMs, pkMEs for module assignment 
  # prior to PPI enrichment test
  # !! Note that we keep the grey color assignments in each of the list_colors_final vectors
  # but not the corresponding values for list_kMEs_final (noGrey)!!
  
  list_colors_final <- lapply(list_list_colors_filter_matched_order, function(x) x[[1]])
  names(list_colors_final) <- sNames_filter
  
  list_MEs_final <- lapply(list_list_MEs_filter_order, function(x) x[[1]])
  names(list_MEs_final) <- sNames_filter
  
  ### EDIT 180428_1
  # list_kMEs_final <- lapply(list_list_kMEs_filter_matched_noGrey_order, function(x) x[[1]])
  # names(list_kMEs_final) <- sNames_filter
  # list_kIMs_final <- lapply(list_list_kIMs_filter_matched_noGrey_order, function(x) x[[1]])
  # names(list_kIMs_final) <- sNames_filter  
  
  list_kMEs <- lapply(list_list_kMEs_filter_matched_order, function(x) x[[1]])
  names(list_kMEs) <- sNames_filter
  list_kIMs <- lapply(list_list_kIMs_filter_matched_order, function(x) x[[1]])
  names(list_kIMs) <- sNames_filter  

  # Name the rows
  # list_kMEs_final <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kMEs_final, y=list_colors_final, SIMPLIFY=F)
  # list_kIMs_final <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kIMs_final, y=list_colors_final, SIMPLIFY=F)
  list_kMEs <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kMEs, y=list_colors_final, SIMPLIFY=F)
  list_kIMs <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kIMs, y=list_colors_final, SIMPLIFY=F)
  
  #list_pkMEs_final <- lapply(list_list_pkMEs_filter_order, function(x) x[[1]])
  #names(list_pkMEs_final) <- sNames_filter
  
  ###
  
  param_names = c("minClusterSize", "deepSplit","pamStage", "moduleMergeCutHeight")
  list_list_cutree_params_final <- lapply(list_vec_n_grey, function(x) comb_list[order(x,decreasing = F)][[1]])
  list_list_cutree_params_final <- lapply(list_list_cutree_params_final, function(x) name_for_vec(to_be_named = x, given_names = param_names, dimension = NULL)) 

  ######################################################################
  ############# RECOMPUTE MEs, kMEs, kIMs, pkMEs AFTER PPI FILTER ############
  ######################################################################
  message("Re-computing kME, primary kMEs and kIMs after filtering modules for significant Protein-Protein Interactions")
  
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_moduleEigengenes_kMEs_pkMEs_PPI.txt"))
  
  list_MEs_PPI_final <- clusterMap(cl, function(x,y) moduleEigengenes(expr=as.data.frame(x,col.names=col.names(x)),
                                                                      colors=y,
                                                                      excludeGrey=excludeGrey), 
                                   x= list_datExpr_filter, 
                                   y= list_colors_PPI_final, 
                                   SIMPLIFY = F)
  
  #names(list_MEs_PPI_final) <- sNames_filter
  list_kMEs_PPI <- clusterMap(cl, function(x,y) signedKME(as.matrix(x),
                      y$eigengenes,
                      outputColumnName = "",
                      corFnc = corFnc), x=list_datExpr_filter, y=list_MEs_PPI_final)
 
  list_kIMs_PPI <- clusterMap(cl, function(x,y) kIM_eachMod_norm(dissTOM=x, colors=y, genes=names(y)),
                              x = list_dissTOM_filter,
                              y = list_colors_PPI_final,
                              SIMPLIFY=F)
  
  stopCluster(cl)
  invisible(gc())
  ### EDIT_180428_1
  # Give row names to kME and kIM vectors
  # list_kMEs_PPI_final <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kMEs_PPI, y=list_colors_PPI_final, SIMPLIFY=F)
  # list_kIMs_PPI_final <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kIMs_PPI, y=list_colors_PPI_final, SIMPLIFY=F)
  # 
  # list_kIMs_PPI_final <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = colnames(y), dimension = 2),
  #                               x = list_kIMs_PPI_final,
  #                               y = list_kMEs_PPI_final, SIMPLIFY=F)
  list_kMEs_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kMEs_PPI, y=list_colors_PPI_final, SIMPLIFY=F)
  list_kIMs_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kIMs_PPI, y=list_colors_PPI_final, SIMPLIFY=F)
  # Give kIM columns names of kMEs
  list_kIMs_PPI <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = colnames(y), dimension = 2),
                                x = list_kIMs_PPI,
                                y = list_kMEs_PPI, SIMPLIFY=F)
  ###
  
  # None of these two function calls work at present - but then we don't need to output pkMEs?
  #list_pkMEs_PPI_final <- parPkMEs_2(list_kMEs=list_kMEs_PPI_final, list_colors=list_colors_PPI_final)
  #clusterMap(cl, function(x,y) getPkMEs(colors=x, kMEs = y), 
  #                                   x = list_colors_PPI_final, 
  #                                   y = list_kMEs_PPI_final, 
  #                                   SIMPLIFY = F)
  #
  #names(list_pkMEs_PPI_final) <- sNames_filter
  
  # Check and remove grey kMEs/kIMs
  #list_kMEs_PPI_final <- lapply(list_kMEs_PPI, checkGrey)
  #list_kIMs_PPI_final <- lapply(list_kIMs_PPI, checkGrey)
  
  #list_kMEs_PPI_final_filter <- list_kMEs_PPI_final[sapply(list_kMEs_PPI_final, function(x) length(colnames(x))>0, simplify=T)]
  #list_kIMs_PPI_final_filter <- list_kIMs_PPI_final[sapply(list_kIMs_PPI_final, function(x) length(colnames(x))>0, simplify=T)]
  
  #sNames_PPI_filter <- sNames_filter[names(list_kMEs_PPI_final) %in% names(list_kMEs_PPI_final_filter)]
  
  #list_colors_PPI_final_filter <- list_colors_PPI_final[sNames_filter %in% sNames_PPI_filter]
  
  # Also update the pre-PPI lists
  #list_colors_final <- list_colors_final[sNames_filter %in% sNames_PPI_filter]
  #list_MEs_final <- list_MEs_final[sNames_filter %in% sNames_PPI_filter]
  #list_kMEs_final <- list_kMEs_final[sNames_filter %in% sNames_PPI_filter]
  #list_kIMs_final <- list_kIMs_final[sNames_filter %in% sNames_PPI_filter]
  #list_geneTree_filter <- list_geneTree_filter[sNames_filter %in% sNames_PPI_filter]

} else if (is.null(STRINGdb_species)) {
  # Just take the first of the colours in each case (since number assigned may be zero..)
  # Can be improved
  list_colors_final <- lapply(list_list_colors_filter_matched, function(x) x[[1]])
  names(list_colors_final) <- sNames_filter
  list_colors_final <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = colnames(y), dimension = NULL), x=list_colors_final, y=list_datExpr_filter, SIMPLIFY=F)
  
  list_MEs_final <- lapply(list_list_MEs_filter, function(x) x[[1]])
  names(list_MEs_final) <- sNames_filter
  
  ### EDIT 180428_1
  # list_kMEs_final <- lapply(list_list_kMEs_filter_matched_noGrey, function(x) x[[1]])
  # names(list_kMEs_final) <- sNames_filter
  # list_kMEs_final <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = colnames(y), dimension = 1), x=list_kMEs_final, y=list_datExpr_filter, SIMPLIFY=F)
  # 
  # list_kIMs_final <- lapply(list_list_kIMs_filter_matched_noGrey, function(x) x[[1]])
  # names(list_kIMs_final) <- sNames_filter
  # list_kIMs_final <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = colnames(y), dimension = 1), x=list_kIMs_final, y=list_datExpr_filter, SIMPLIFY=F)
   
  list_kMEs <- lapply(list_list_kMEs_filter_matched, function(x) x[[1]])
  names(list_kMEs) <- sNames_filter
  list_kMEs <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = colnames(y), dimension = 1), x=list_kMEs, y=list_datExpr_filter, SIMPLIFY=F)
  
  list_kIMs <- lapply(list_list_kIMs_filter_matched, function(x) x[[1]])
  names(list_kIMs) <- sNames_filter
  list_kIMs <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = colnames(y), dimension = 1), x=list_kIMs, y=list_datExpr_filter, SIMPLIFY=F)
}  

# clear memory
rm(list_dissTOM_filter, list_list_plot_label_filter_order, list_list_colors_PPI_order, list_list_colors_filter_matched_order,
   list_list_MEs_filter_order, list_list_kMEs_filter_matched_order, list_list_kIMs_filter_matched_order, list_list_pkMEs_filter_order)
invisible(gc()) 

    #list_pkMEs_final <- lapply(list_list_pkMEs_filter, function(x) x[[1]])
    #names(list_pkMEs_final) <- sNames_filter
  
  #}  
  
  # resume = "checkpoint_6"
  # message("Reached checkpoint 6, saving session image")
  # save.image(file=sprintf("%s%s_checkpoint_5_image.RData", RObjects_dir, data_prefix))

# } else if (resume == "checkpoint_5") {
#   load(file=sprintf("%s%s_checkpoint_5_image.RData", RObjects_dir, data_prefix))
#   # Load parameter values and utility functions anew (in case a bug was fixed)
#   source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
#   source(file = "/projects/jonatan/functions-src/functions.R")
# } 

# if (plot_permuted==T) {  
#   
#   # NOTE: Alternatively use sampledBlockwiseModules() - see /projects/jonatan/wgcna-src/180228_robust_wgcna_module_detection.Rmd
#   bootstrap_seq <- if (replace==T) 1:(nPermutations+1) else 1:nPermutations
#   
#   bootstrap_indiv_TOMs <- lapply(bootstrap_seq, function(x) load_obj(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", RObjects_dir, sNames, x, 1))) # retrieve permuted dataset TOMs from disk
#   
#   # Delete individual TOMs from disk
#   for (i in bootstrap_seq) {
#     if (file.exists(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", RObjects_dir, sNames, i, 1))) {
#       file.remove(sprintf("%s%s_individualTOM-Set%s-Block%s.RData", RObjects_dir, sNames, i, 1))
#     }
#   }
#   
#   bootstrap_dissTOMs <- lapply(bootstrap_indiv_TOMs, function(x) 1-as.dist(x)) # Convert permuted TOMs to distance matrices
#   
#   bootstrap_geneTrees <- lapply(bootstrap_dissTOMs, function(x) hclust(d=x, method=hclustMethod)) # Do hierarchical clustering 
#   
#   # cut out modules in each gene tree
#   bootstrap_cutrees <- mapply(function(x,y) cutreeHybrid(dendro = x, distM = as.matrix(y), cutHeight = NULL, minClusterSize = cutree_params_final$minClusterSize, deepSplit = cutree_params_final$deepSplit, pamStage = cutree_params_final$pamStage, pamRespectsDendro = cutree_params_final$pamStage), bootstrap_geneTrees, bootstrap_dissTOMs, SIMPLIFY=F) 
#   
#   # merge modules with correlated eigengenes
#   bootstrap_merged <- mapply(function(x,y) mergeCloseModules(exprData = as.matrix(x$data)[,goodGenesTOM_idx], colors = y$labels, cutHeight=cutree_params_final$moduleMergeCutHeight), multiExpr, bootstrap_cutrees, SIMPLIFY=FALSE) # Merge modules whose eigengenes are highly correlated
#   
#   bootstrap_colors = lapply(bootstrap_merged, function(x) labels2colors(x$colors)) 
#   
#   bootstrap_plot_labels <- c("Original data", paste0("Permutation ", as.character(bootstrap_seq)[-max(bootstrap_seq)]))
#   
#   # Free up memory
#   rm(bootstrap_indiv_TOMs)
#   rm(bootstrap_dissTOMs)
#   
#   # Align the module coloring schemes 
#   bootstrap_colors_final = matrix(0, ncol(datExpr_filter), nPermutations + 1)
#   bootstrap_colors_final[, 1] = colors # use the consensus colors as 'reference'
#   
#   for (r in 2:(nPermutations+1))
#   {
#     bootstrap_colors_final[, r] = matchLabels(source = bootstrap_colors[[r-1]], # traverse gene colour assignments
#                               reference = bootstrap_colors_final[, 1],
#                               pThreshold = 5e-2,
#                               na.rm = TRUE,
#                               extraLabels = standardColors())# paste0("Extra_", 0:100)) # labels for modules in source that cannot be matched to any modules in reference
#   }
# 
#   bootstrap_plot_labels_final <- c("Consensus", bootstrap_plot_labels)
# 
#   pdf(sprintf("%s%s_%s_diffbootstrap_colors_%s.pdf", plots_dir, data_prefix, sNames, flag_date),width=9,height=8)
#   plotDendroAndColors(geneTree, 
#                       labels2colors(bootstrap_colors_final),
#                       #matrix(unlist(bootstrap_colors_final), ncol = length(bootstrap_plot_labels_final), byrow = F), 
#                       groupLabels=bootstrap_plot_labels_final, 
#                       addGuide= TRUE, 
#                       dendroLabels=FALSE, 
#                       main="Modules: consensus and on each resampling", 
#                       cex.colorLabels=0.5)
#   dev.off()


######################################################################
############ ADD GENE ONTOLOGY TERMS FOR MODULES #####################
######################################################################

list_list_gprofiles <- NULL
list_list_gprofiles_PPI <- NULL


if (FALSE) {#(!is.null(organism)) {
  message("Searching the Gene Ontology database for terms corresponding to the gene modules")
  
  list_list_module_genes <- lapply(list_colors_final, function(x) lapply(names(table(x)), function(y) names(x)[x==y])) 
  
  # name each list entry
  list_list_module_genes <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names= colnames(y), dimension = NULL),  x=list_list_module_genes,  y=list_kIMs, SIMPLIFY=F)
    
  # Order genes by kIM to their own module. This allows us to submit them as ranked queries to gprofiler for GSEA style p-values (https://cran.r-project.org/web/packages/gProfileR/gProfileR.pdf)
  list_list_module_genes_order <- mapply(function(a,b) lapply(names(a), 
                                                              function(x) a[[x]]  [ order(b[[x]] [rownames(b) %in% a[[x]]], decreasing = T) ]  ), a = list_list_module_genes,  b = list_kIMs, SIMPLIFY=F)

  ### EDIT 180429_1 Not necessary: https://biit.cs.ut.ee/gprofiler/help.cgi?help_id=0#help_id_3
  # Only include genes that have been mapped to ensembl ENSG earlier
  ### EDIT_180428_2
  #list_list_module_genes_order_ENSG <- mapply(function(x,y) lapply(x, function(z) z[z %in% y]), x=list_list_module_genes_order, y = list_ensembl_IDs, SIMPLIFY = F)
  ###
  ###
  
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", 
                     outfile = paste0(log_dir, "log_gprofiler.txt"))
  list_list_gprofiles <- parLapply(cl, list_list_module_genes_order, function(x) lapply(x, function(y) gprofiler(query=y, 
                              organism = organism, 
                              sort_by_structure = T,
                              ordered_query = T,
                              significant = T,
                              exclude_iea = F, # TODO check with Dylan
                              underrep = F,
                              evcodes = F,
                              region_query = F,
                              max_p_value = 1,
                              min_set_size = 0,
                              max_set_size = 0,
                              min_isect_size = 0,
                              correction_method = "analytical",
                              hier_filtering = "none",
                              domain_size = "annotated",
                              custom_bg = "",
                              numeric_ns = "",
                              png_fn = NULL,
                              include_graph = F,
                              src_filter = NULL)))
  
  stopCluster(cl)
  invisible(gc())
  
  # Assign names to each list with gprofiler dataframes
  list_list_gprofiles <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names= colnames(y), dimension = NULL),  x=list_list_gprofiles,  y=list_kIMs, SIMPLIFY=F)
  

  if (!is.null(STRINGdb_species)) {
    
    list_list_module_PPI_genes <- lapply(list_colors_PPI_final, function(x) lapply(names(table(x)), function(y) names(x)[x==y])) 
    list_list_module_PPI_genes <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names= colnames(y), dimension = NULL),  x=list_list_module_PPI_genes,  y=list_kIMs_PPI_final, SIMPLIFY=F)
  
    # Order genes by kIM to their own module. This allows us to submit them as ranked queries to gprofiler for GSEA style p-values (https://cran.r-project.org/web/packages/gProfileR/gProfileR.pdf)
    list_list_module_PPI_genes_order <- mapply(function(a,b) lapply(names(a), 
                                                                function(x) a[[x]]  [ order(b[[x]] [rownames(b) %in% a[[x]]], decreasing = T) ]  ), a = list_list_module_PPI_genes,  b = list_kIMs_PPI, SIMPLIFY=F)
    
    ### EDIT 180429_1 Not necessary: https://biit.cs.ut.ee/gprofiler/help.cgi?help_id=0#help_id_3
    # Only include genes that have been mapped to ensembl ENSG earlier
    ### EDIT_180428_2
    #list_list_module_PPI_genes_order_ENSG <- mapply(function(x,y) lapply(x, function(z) z[z %in% y]), x=list_list_module_PPI_genes_order, y = list_ensembl_IDs, SIMPLIFY = F)
    ### 
    ###
    invisible(gc())
    cl <- makeCluster(n_cores, type = "FORK", 
                      outfile = paste0(log_dir, "log_gprofiler_PPI.txt"))
    
    list_list_gprofiles_PPI <- parLapply(cl, list_list_module_PPI_genes_order, function(x) lapply(x, function(y) gprofiler(query=y, 
                                                                                 organism = organism, 
                                                                                 sort_by_structure = T,
                                                                                 ordered_query = T, 
                                                                                 significant = T, 
                                                                                 exclude_iea = F, # TODO check with Dylan
                                                                                 underrep = F,
                                                                                 evcodes = F, 
                                                                                 region_query = F, 
                                                                                 max_p_value = 1, 
                                                                                 min_set_size = 0,
                                                                                 max_set_size = 0, 
                                                                                 min_isect_size = 0, 
                                                                                 correction_method = "analytical",
                                                                                 hier_filtering = "none", 
                                                                                 domain_size = "annotated", 
                                                                                 custom_bg = "",
                                                                                 numeric_ns = "", 
                                                                                 png_fn = NULL, 
                                                                                 include_graph = F, 
                                                                                 src_filter = NULL)))
    
  stopCluster(cl)
  invisible(gc())
  # Assign names to each list with gprofiler dataframes
  list_list_gprofiles_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names= colnames(y), dimension = NULL),  x=list_list_gprofiles_PPI,  y=list_kIMs_PPI, SIMPLIFY=F)
  
  } # if (!is.null(STRINGdb_species)) {} 
} #if (!is.null(organism)) 

######################################################################
######################### FORMAT OUTPUTS #############################
######################################################################


list_kMEs_noGrey <- deleteGrey(list_kMEs)
list_kIMs_noGrey <- deleteGrey(list_kIMs)

# Reformat results for writing to csv
### REV_180428: Changed suffix from _out to _final
list_kMEs_out <- lapply(list_kMEs_noGrey, function(x) cbind(genes=rownames(x), x))
list_kIMs_out <- lapply(list_kIMs_noGrey, function(x) cbind(genes=rownames(x), x))

if (!is.null(STRINGdb_species)) {
  ### REV_180428: Changed suffix from _out to _final
  list_kMEs_PPI_noGrey <- deleteGrey(list_kMEs_PPI)
  list_kIMs_PPI_noGrey <- deleteGrey(list_kIMs_PPI)
  
  list_kMEs_PPI_out <- lapply(list_kMEs_PPI_noGrey, function(x) cbind(genes=rownames(x), x))
  list_kIMs_PPI_out <- lapply(list_kIMs_PPI_noGrey, function(x) cbind(genes=rownames(x), x))
}  

### 180503 v1.7
if (!is.null(organism)) {
  # The whole dataset
  #all_genes_ensembl_IDs <- perslabEnsembl$ensembl_gene_id[match(all_genes, perslabEnsembl$gene_name_optimal)]
  # Where the mapping failed, revert to the previous name
  #all_genes_IDs_mix <- ifelse(is.na(all_genes_ensembl_IDs), yes = all_genes, no = all_genes_ensembl_IDs)
  # Each subset
  list_ensembl_IDs <- lapply(list_datExpr_filter, function(x) perslabEnsembl$ensembl_gene_id[match(colnames(x), perslabEnsembl$gene_name_optimal)])
  list_idx_mapped <- lapply(list_ensembl_IDs, function(x) !is.na(x))
  list_prop_not_map <- lapply(list_idx_mapped, function(x) signif(as.double(sum(!x)) / as.double(1e-10+length(x)),2)) # TODO currently unused
  # Where the mapping failed, revert to the previous name
  #list_ensembl_IDs_mix <- mapply(FUN = function(x,y) ifelse(is.na(x), yes = colnames(y), no = x), x=list_ensembl_IDs, y=list_datExpr, SIMPLIFY=F)

  # Rename the expression matrix colmns
  #list_datExpr <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = 2),  x=list_datExpr,  y=list_ensembl_IDs_mix, SIMPLIFY=F)
  list_ensembl_out <- mapply(function(x,y) data.frame(module_name = x, symbol = names(x), ensembl_ID = y, row.names = NULL), x=list_colors_final, y=list_idx_mapped, SIMPLIFY=F)
  if(!is.null(STRINGdb_species)) {
    list_ensembl_PPI_out <- mapply(function(x,y) data.frame(module_name = x, symbol = names(x), ensembl_ID = y, row.names = NULL), list_colors_PPI_final, list_idx_mapped, SIMPLIFY=F)
  }
}

# if (!is.null(ensembl_dataset)) {
#   list_ensembl_out <- mapply(function(x,y) data.frame(module_name = x, ensembl_ID = names(x), mapped_ENSG_ok = y, row.names = NULL), list_colors_final, list_idx_mapped_gg_filter, SIMPLIFY=F)
#   if(!is.null(STRINGdb_species)) {
#     list_ensembl_PPI_out <- mapply(function(x,y) data.frame(module_name = x, ensembl_ID = names(x), mapped_ENSG_ok = y, row.names = NULL), list_colors_PPI_final, list_idx_mapped_gg_filter, SIMPLIFY=F)
#   }
# }

# clear memory
rm(list_kMEs_noGrey, list_kIMs_noGrey, list_kMEs_PPI_noGrey, list_kIMs_PPI_noGrey)
invisible(gc())

######################################################################
########################## SAVE TABLES ###############################
######################################################################
#list_pkMEs_out <- mapply(function(x,y) data.frame(colors_final = x, genes = names(x), pkMEs = y), list_colors_final, list_pkMEs_final, SIMPLIFY = F)
# save kMEs and kMIs primary to file
invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kMEs_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_kMEs_out, sNames_filter, SIMPLIFY = F))
invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kIMs_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_kIMs_out, sNames_filter, SIMPLIFY = F))
#invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_pkMEs_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_pkMEs_final, sNames_filter, SIMPLIFY = F))

if (!is.null(STRINGdb_species)) { 
  #list_pkMEs_PPI_final <- mapply(function(x,y) data.frame(colors_PPI_final = x, genes = names(x), pkMEs_PPI = y), list_colors_PPI_final, list_pkMEs_PPI_final, SIMPLIFY = F)
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kMEs_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_kMEs_PPI_out, sNames_filter, SIMPLIFY = F))
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kIMs_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_kIMs_PPI_out, sNames_filter, SIMPLIFY = F))
  #invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_pkMEs_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_pkMEs_PPI_final, sNames_filter, SIMPLIFY = F))
}

# ensembl IDs for LD score regression
if (!is.null(ensembl_dataset)) {
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_ensembl_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_ensembl_out, sNames_filter, SIMPLIFY = F))
  if(!is.null(ensembl_dataset)) {
    invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_ensembl_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_ensembl_PPI_out, sNames_filter, SIMPLIFY = F))
  }  
}

# gprofiler results 
if (!is.null(list_list_gprofiles)) {
  # invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_gprofiles_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_list_gprofiles, sNames_filter, SIMPLIFY = F))
  # if(!is.null(STRINGdb_species) & !is.null(list_list_gprofiles_PPI)) {
  #   invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_gprofiles_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_list_gprofiles_PPI, sNames_filter, SIMPLIFY = F))
  # } 
  invisible(mapply(function(x,y) mapply(function(a,b) write.csv(a, file=sprintf("%s%s_%s_gprofiles_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), a=x, b=names(x)), list_list_gprofiles, sNames_filter, SIMPLIFY = F))
  if(!is.null(STRINGdb_species) & !is.null(list_list_gprofiles_PPI)) {
    invisible(mapply(function(x,y) mapply(function(a,b) write.csv(a, file=sprintf("%s%s_%s_gprofiles_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), a=x, b=names(x)), list_list_gprofiles_PPI, sNames_filter, SIMPLIFY = F))
  } 
}

invisible(gc())
    
##########################################################################
#################### COMPUTE MAGMA GWAS ENRICHMENT #######################
##########################################################################
message("Scoring modules for enrichment with genes linked by GWAS to phenotypes of interest")

file_suffix = "kMEs"
cl <- makeCluster(n_cores, type = "FORK", 
                  outfile = paste0(log_dir, "log_MAGMA.txt"))
invisible(parLapply(cl, sNames_filter, function(x) magma_for_par(subsetName = x,
                                                            project_dir = project_dir,
                                                            plots_dir = plots_dir,
                                                            log_dir = log_dir,
                                                            tables_dir = tables_dir,
                                                            magma_gwas_dir = magma_gwas_dir,
                                                            data_prefix = data_prefix,
                                                            file_suffix = file_suffix,
                                                            flag_date = flag_date)))

stopCluster(cl)
invisible(gc())

if (!is.null(STRINGdb_species)) {
  
  file_suffix = "kMEs_PPI"
  cl <- makeCluster(n_cores, type = "FORK", 
                    outfile = paste0(log_dir, "log_MAGMA_PPI.txt"))
  
  invisible(parLapply(cl, sNames_filter, function(x) magma_for_par(subsetName = x,
                                                              project_dir = project_dir,
                                                              plots_dir = plots_dir,
                                                              log_dir = log_dir,
                                                              tables_dir = tables_dir,
                                                              magma_gwas_dir = magma_gwas_dir,
                                                              data_prefix = data_prefix,
                                                              file_suffix = file_suffix,
                                                              flag_date = flag_date)))
  stopCluster(cl)
  invisible(gc())
}

    
######################################################################
############################# SAVE PLOTS #############################
######################################################################

message("Plotting")

if (!is.null(STRINGdb_species)) {
  list_colors_final_both = mapply(FUN=cbind, list_colors_final, list_colors_PPI_final, SIMPLIFY = F)
  names(list_colors_final_both) = sNames_filter
  invisible(lapply(sNames_filter, plotFinalColors))
  
} else if (is.null(STRINGdb_species)) {
  list_colors_final_both = list_colors_final
  names(list_colors_final_both) = sNames_filter
  invisible(lapply(sNames_filter, plotFinalColors))
}

if (FALSE) {  # TODO
  
  ######################################################################
  ############# PLOT kME - kIM CORRELATIONS WITHIN CELLTYPES ###########
  ######################################################################
  # This is a diagnostic plot to see whether kIMs are highly correlated with kMEs as expected
  
  list_corr_kMEs_kIMs <- mapply(function(x,y) cor(x,y), x=list_kMEs, y=list_kIMs, SIMPLIFY=F)
  invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
                                                  vectorNames = "kME-kIM", 
                                                  subsetName=y, 
                                                  diag=T,
                                                  is.corr = T), 
                   x=list_corr_kMEs_kIMs, y=sNames_filter, SIMPLIFY=F))
  
  if (!is.null(STRINGdb_species)){
    list_corr_kMEs_kIMs_PPI <- mapply(function(x,y) cor(x,y), x=list_kMEs_PPI, y=list_kIMs_PPI, SIMPLIFY=F)
    invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
                                                    vectorNames = "kME-kIM_PPI", 
                                                    subsetName=y, 
                                                    diag=T,
                                                    is.corr = T), 
                     x=list_corr_kMEs_kIMs_PPI, y=sNames_filter, SIMPLIFY=F))
  }
  ######################################################################
  ######## PLOT kME-kME & kIM-kIM CORRELATIONS WITHIN CELLTYPES ########
  ######################################################################
  
  list_corr_kMEs_kMEs <- lapply(list_kMEs, function(x) WGCNA::cor(x,x, method="spearman", verbose=verbose))
  list_corr_kIMs_kIMs <- lapply(list_kIMs, function(x) WGCNA::cor(x,x, method="spearman", verbose=verbose))
  
  # EDIT_180429_4
  # The names screw up the plots
  # # name the correlation matrix columns
  # list_corr_kMEs_kMEs <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = 2), 
  #                               x = list_corr_kMEs_kMEs,
  #                               y = list_kMEs, SIMPLIFY=F)
  # list_corr_kIMs_kIMs <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = 2), 
  #                               x = list_corr_kIMs_kIMs,
  #                               y = list_kIMs, SIMPLIFY=F)
  # 
  # # name the correlation matrix rows
  # list_corr_kMEs_kMEs <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension=1), 
  #                               x = list_corr_kMEs_kMEs,
  #                               y = list_kMEs, SIMPLIFY=F)
  # list_corr_kIMs_kIMs <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension=1), 
  #                               x = list_corr_kIMs_kIMs,
  #                               y = list_kIMs, SIMPLIFY=F)
  
  # plot
  invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
                                                  vectorNames = "kME-kME", 
                                                  subsetName=y, 
                                                  diag=F,
                                                  is.corr=T,
                                                  order = "hclust",
                                                  hclust.method = hclustMethod), 
                   x=list_corr_kMEs_kMEs, y=sNames_filter, SIMPLIFY=F))
  
  invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
                                                  vectorNames = "kIM-kIM", 
                                                  subsetName=y, 
                                                  diag=F,
                                                  is.corr=T,
                                                  order = "hclust",
                                                  hclust.method = hclustMethod), 
                   x=list_corr_kIMs_kIMs, y=sNames_filter, SIMPLIFY=F))

  
  if (!is.null(STRINGdb_species)){
    list_corr_kMEs_kMEs_PPI <- lapply(list_kMEs_PPI, function(x) WGCNA::cor(x,x, method="spearman", verbose=verbose))
    list_corr_kIMs_kIMs_PPI <- lapply(list_kIMs_PPI, function(x) WGCNA::cor(x,x, method="spearman", verbose=verbose))
    
    # Name the correlation matrix columns
    # list_corr_kMEs_kMEs_PPI <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = 2), 
    #                                   x = list_corr_kMEs_kMEs_PPI,
    #                                   y = list_kMEs_PPI, SIMPLIFY=F)
    # list_corr_kIMs_kIMs_PPI <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = 2), 
    #                                   x = list_corr_kIMs_kIMs_PPI,
    #                                   y = list_kIMs_PPI, SIMPLIFY=F)
    # 
    # # Name the correlation matrix rows
    # list_corr_kMEs_kMEs_PPI <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension=1), 
    #                                   x = list_corr_kMEs_kMEs_PPI,
    #                                   y = list_kMEs_PPI, SIMPLIFY=F)
    # list_corr_kIMs_kIMs_PPI <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension=1), 
    #                                   x = list_corr_kIMs_kIMs_PPI,
    #                                   y = list_kIMs_PPI, SIMPLIFY=F)
    # 
    # plot
    invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
                                                    vectorNames = "kME-kME_PPI", 
                                                    subsetName=y, 
                                                    diag=F,
                                                    is.corr=T,
                                                    order = "hclust",
                                                    hclust.method = hclustMethod), 
                     x=list_corr_kMEs_kMEs_PPI, y=sNames_filter, SIMPLIFY=F))
    
    invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
                                                    vectorNames = "kIM-kIM_PPI", 
                                                    subsetName=y, 
                                                    diag=F,
                                                    is.corr=T,
                                                    order = "hclust",
                                                    hclust.method = hclustMethod), 
                     x=list_corr_kIMs_kIMs_PPI, y=sNames_filter, SIMPLIFY=F))
  }

  ######################################################################
  ######################### LOAD FULL DATA #############################
  ######################################################################
  
  #datExpr_full <- t(load_obj(seurat_path)@data)
  seurat_obj <- load_obj(seurat_path)
  datExpr_full <- t(seurat_obj@data)
  ######################################################################
  ####### PLOT kME-kME and kIM-kIM correlation across all subsets ######
  ######################################################################

  # Include all genes in kMEs / kIMs with value 0 if they were not included in the subset analysis
  ### 180503_v1.7
  # need to use original names
  # list_kMEs_all_genes <- lapply(list_kMEs, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes_IDs_mix))
  # list_kIMs_all_genes <- lapply(list_kIMs, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes_IDs_mix))
  list_kMEs_all_genes <- lapply(list_kMEs, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes))
  list_kIMs_all_genes <- lapply(list_kIMs, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes))
  ###
  
  # Bind all celltype module kMEs / kIMs in single dataframes
  
  for (i in 1:length(list_kMEs_all_genes)) {
    if (i == 1) {
      kMEs_across_subsets <- list_kMEs_all_genes[[i]] 
    } else { kMEs_across_subsets <- cbind(kMEs_across_subsets,list_kMEs_all_genes[[i]])
    }
  }
  
  for (i in 1:length(list_kIMs_all_genes)) {
    if (i == 1) {
      kIMs_across_subsets <- list_kIMs_all_genes[[i]] 
    } else { kIMs_across_subsets <- cbind(kIMs_across_subsets,list_kIMs_all_genes[[i]])
    }
  }
  #kMEs_across_subsets <- cbind(eval(parse(text=sprintf("c(%s)", paste0("list_kMEs_all_genes[[", 1:length(list_kMEs_all_genes), "]]")))))
  #kIMs_across_subsets <- cbind(eval(parse(text='c(paste0("list_kIMs_all_genes[[", 1:length(list_kIMs_all_genes), "]]"))')))

  # name the dataframes and matrices so as to distinguish modules from different subsets   
  prefixes = unlist(mapply(function(x,y) rep(x = x, times= ncol(y)), x=sNames_filter, y= list_kMEs_all_genes, SIMPLIFY=T)) # same for kME and kIM

  colnames(kMEs_across_subsets) = paste(prefixes, colnames(kMEs_across_subsets), sep="_")
  colnames(kIMs_across_subsets) = paste(prefixes, colnames(kIMs_across_subsets), sep="_")
  
  # correlation matrices
  corr_kMEs_across_subsets <- WGCNA::cor(as.matrix(kMEs_across_subsets), method="spearman", verbose = verbose)
  corr_kIMs_across_subsets <- WGCNA::cor(as.matrix(kIMs_across_subsets), method="spearman", verbose = verbose)
  
  # Name the matrices
  # rownames(corr_kMEs_across_subsets) = paste(prefixes, colnames(kMEs_across_subsets), sep="_")
  # rownames(corr_kIMs_across_subsets) = paste(prefixes, colnames(kIMs_across_subsets), sep="_")
  # 
  # colnames(corr_kMEs_across_subsets) = paste(prefixes, colnames(kMEs_across_subsets), sep="_")
  # colnames(corr_kIMs_across_subsets) = paste(prefixes, colnames(kIMs_across_subsets), sep="_")
  # 
  # Include all genes in kIMs with value 0 if they were not included in the subset analysis
  # Bind all celltype module kMEs in a single dataframe
  # name the dataframe and matrix so as to distinguish subsets 
  
  # plot
  pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "kME-kME", flag_date))#,width=ncol(kMEs_across_subsets) %/% 4,height=ncol(kMEs_across_subsets) %/% 4)
  corrplot(corr = corr_kMEs_across_subsets, 
           method = "color",
           is.corr=T,
           title=sprintf("%s %s corr across subsets", data_prefix, "kME-kME"),
           #addCoef.col = F,
           number.digits = NULL)
  
  invisible(dev.off())
  
 # plot
  pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "kIM-kIM", flag_date),width=ncol(kIMs_across_subsets) %/% 4,height=ncol(kIMs_across_subsets) %/% 4)
  corrplot(corr = corr_kIMs_across_subsets, 
           add = F,
           method = "color",
           is.corr=T,
           #order = "hclust",
           #hclust.method = hclustMethod,
           title=sprintf("%s %s corr across subsets", data_prefix, "kIM-kIM"),
           #addCoef.col = F,
           number.digits = NULL)
  
  invisible(dev.off())
  
  if (!is.null(STRINGdb_species)){
    #### 180503_v1.7
    # Include all genes in kMEs / kIMs with value 0 if they were not included in the subset analysis
    # list_kMEs_PPI_all_genes <- lapply(list_kMEs_PPI, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes_IDs_mix))
    # list_kIMs_PPI_all_genes <- lapply(list_kIMs_PPI, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes_IDs_mix))
    list_kMEs_PPI_all_genes <- lapply(list_kMEs_PPI, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes))
    list_kIMs_PPI_all_genes <- lapply(list_kIMs_PPI, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes))
    
    ###
    
    # Bind all celltype module kMEs / kIMs in single dataframes
    
    for (i in 1:length(list_kMEs_PPI_all_genes)) {
      if (i == 1) {
        kMEs_PPI_across_subsets <- list_kMEs_PPI_all_genes[[i]] 
      } else { kMEs_PPI_across_subsets <- cbind(kMEs_PPI_across_subsets,list_kMEs_PPI_all_genes[[i]])
      }
    }
    
    for (i in 1:length(list_kIMs_PPI_all_genes)) {
      if (i == 1) {
        kIMs_PPI_across_subsets <- list_kIMs_PPI_all_genes[[i]] 
      } else { kIMs_PPI_across_subsets <- cbind(kIMs_PPI_across_subsets,list_kIMs_PPI_all_genes[[i]])
      }
    }
    
    # Bind all celltype module kMEs / kIMs in single dataframes
    # kMEs_PPI_across_subsets <- cbind(eval(parse(text=c(paste0("list_kMEs_PPI_all_subsets[[", 1:length(list_kMEs_PPI_all_genes), "]]")))))
    # kMEs_PPI_across_subsets <- cbind(eval(parse(text=c(paste0("list_kMEs_PPI_all_subsets[[", 1:length(list_kMEs_PPI_all_genes), "]]")))))
    # 
    # name the dataframes and matrices so as to distinguish modules from different subsets   
    prefixes = unlist(mapply(function(x,y) rep(x = x, times= ncol(y)), x=sNames_filter, y= list_kMEs_PPI_all_genes, SIMPLIFY=T)) # same for kME and kIM
    
    colnames(kMEs_PPI_across_subsets) = paste(prefixes, colnames(kMEs_PPI_across_subsets), sep="_")
    colnames(kMEs_PPI_across_subsets) = paste(prefixes, colnames(kMEs_PPI_across_subsets), sep="_")
    
    # correlation matrices
    corr_kMEs_PPI_across_subsets <- WGCNA::cor(as.matrix(kMEs_PPI_across_subsets), method="spearman", verbose=verbose)
    corr_kMEs_PPI_across_subsets <- WGCNA::cor(as.matrix(kMEs_PPI_across_subsets), method="spearman", verbose=verbose)
    
    # Name the matrices
    # rownames(corr_kMEs_PPI_across_subsets) = paste(prefixes, colnames(kMEs_PPI_across_subsets), sep="_")
    # rownames(corr_kMEs_PPI_across_subsets) = paste(prefixes, colnames(kMEs_PPI_across_subsets), sep="_")
    # 
    # colnames(corr_kMEs_PPI_across_subsets) = paste(prefixes, colnames(kMEs_PPI_across_subsets), sep="_")
    # colnames(corr_kMEs_PPI_across_subsets) = paste(prefixes, colnames(kMEs_PPI_across_subsets), sep="_")
    # 
    # Include all genes in kIMs with value 0 if they were not included in the subset analysis
    # Bind all celltype module kMEs in a single dataframe
    # name the dataframe and matrix so as to distinguish subsets 
    
    # plot
    pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "kME-kME_PPI", flag_date),width=ncol(kMEs_PPI_across_subsets) %/% 4,height=ncol(kMEs_PPI_across_subsets) %/% 4)
    corrplot(corr = corr_kMEs_PPI_across_subsets, 
             method = "color",
             is.corr=T,
             title=sprintf("%s %s corr across subsets", data_prefix, "kME-kME (PPI)"))#,
             #order = "hclust",
             #hclust.method = hclustMethod)
    
    invisible(dev.off())
    
    # plot
    pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "kIM-kIM_PPI", flag_date),width=ncol(kMEs_PPI_across_subsets) %/% 4,height=ncol(kMEs_PPI_across_subsets) %/% 4)
    corrplot(corr = corr_kMEs_PPI_across_subsets, 
             method = "color",
             is.corr=T,
             title=sprintf("%s %s corr across subsets", data_prefix, "kIM-kIM (PPI)"))#,
             #order = "hclust",
             #hclust.method = hclustMethod)
    
    invisible(dev.off())
  }
  ######################################################################
  ########### PLOT KME - & KIM - METADATA CORRELATIONS #################
  ######################################################################
  
  if (!is.null(meta.data)) {
    
    # Compute cell embedding
    cell_embeddings_kME <- datExpr_full %*% as.matrix(kMEs_across_subsets) # cell x gene %*% gene x module -> cell x module 
    cell_embeddings_kIM <- datExpr_full %*% as.matrix(kIMs_across_subsets) # cell x gene %*% gene x module -> cell x module 
    
    # correlation matrices 
    corr_kME_meta.data <- WGCNA::cor(x=cell_embeddings_kME, y=meta.data, method="spearman", verbose=verbose)
    corr_kIM_meta.data <- WGCNA::cor(x=cell_embeddings_kME, y=meta.data, method="spearman", verbose=verbose)
    
    # name the matrices
    # rownames(corr_kME_meta.data) = meta.data_corr_col_numeric
    # rownames(corr_kIM_meta.data) = meta.data_corr_col_numeric
    # 
    # colnames(corr_kME_meta.data) = paste(prefixes, colnames(kMEs_PPI_across_subsets), sep="_")
    # colnames(corr_kIM_meta.data) = paste(prefixes, colnames(kIMs_PPI_across_subsets), sep="_")
    
    # plot
    pdf(sprintf("%s%s_corr_%s_with_metadata_%s.pdf", plots_dir, data_prefix, "kME_embeddings", flag_date),width=ncol(meta.data) %/% 4,height=ncol(kMEs_across_subsets) %/% 4)
    corrplot(corr = corr_kME_meta.data, 
             method = "color",
             is.corr=T,
             title=sprintf("%s corr: module %s with metadata", data_prefix, "kME embeddings"))#,
             #order = "hclust",
             #hclust.method = hclustMethod)
    
    dev.off()
    
    pdf(sprintf("%s%s_corr_%s_with_metadata_%s.pdf", plots_dir, data_prefix, "kIM_embeddings", flag_date),width=ncol(meta.data) %/% 4,height=ncol(kIMs_across_subsets) %/% 4)
    corrplot(corr = corr_kIM_meta.data, 
             method = "color",
             is.corr=T,
             title=sprintf("%s corr: module %s with metadata", data_prefix, "kIM embeddings"))#,
             #order = "hclust",
             #hclust.method = hclustMethod)
    
    dev.off()
    
    if (!is.null(STRINGdb_species)) {

      # Compute cell embedding
      cell_embeddings_kMEs_PPI <- datExpr_full %*% as.matrix(kMEs_PPIs_across_subsets) # cell x gene %*% gene x module -> cell x module 
      cell_embeddings_kIMs_PPI <- datExpr_full %*% as.matrix(kIMs_PPIs_across_subsets) # cell x gene %*% gene x module -> cell x module 
      
      # correlation matrices 
      corr_kMEs_PPI_meta.data <- WGCNA::cor(x=cell_embeddings_kMEs_PPI, y=meta.data, method="spearman", verbose=verbose)
      corr_kIMs_PPI_meta.data <- WGCNA::cor(x=cell_embeddings_kMEs_PPI, y=meta.data, method="spearman", verbose=verbose)
      
      # name the matrices
      # rownames(corr_kMEs_PPI_meta.data) = meta.data_corr_col_numeric
      # rownames(corr_kIMs_PPI_meta.data) = meta.data_corr_col_numeric
      # 
      # colnames(corr_kMEs_PPI_meta.data) = paste(prefixes, colnames(kMEs_PPIs_PPI_across_subsets), sep="_")
      # colnames(corr_kIMs_PPI_meta.data) = paste(prefixes, colnames(kIMs_PPI_across_subsets), sep="_")
      # 
      # plot
      pdf(sprintf("%s%s_corr_%s_with_metadata_%s.pdf", plots_dir, data_prefix, "kMEs_PPI_embeddings", flag_date),width = max(4, ncol(meta.data) %/% 4), height = ncol(kMEs_PPIs_across_subsets) %/% 4)
      corrplot(corr = corr_kMEs_PPI_meta.data, 
               method = "color",
               is.corr=T,
               title=sprintf("%s corr: module %s with metadata", data_prefix, "kMEs_PPI embeddings"))#,
               #order = "hclust",
               #hclust.method = hclustMethod)
      
      dev.off()
      
      pdf(sprintf("%s%s_corr_%s_with_metadata_%s.pdf", plots_dir, data_prefix, "kIMs_PPI_embeddings", flag_date),width = max(4, ncol(meta.data) %/% 4), height = ncol(kIMs_PPIs_across_subsets) %/% 4)
      corrplot(corr = corr_kIMs_PPI_meta.data, 
               method = "color",
               is.corr=T,
               title=sprintf("%s corr: module %s with metadata", data_prefix, "kIMs_PPI embeddings"))#,
               #order = "hclust",
               #hclust.method = hclustMethod)
      
      dev.off()
    }
  }
  
  ######################################################################
  ##### PLOT CELL EXPRESSION ON MODULES WITH HIGH GWAS SIGNIFICANCE ####
  ######################################################################
  
  gwas_sub_dirs = list.dirs(path = magma_gwas_dir, full.names = FALSE, recursive = FALSE)
  
  file_suffix = "kMEs"
  # upper level list of subsets, nested list with tables for each GWAS repo
  list_list_magma_tables <- lapply(sNames_filter, function(x) lapply(gwas_sub_dirs, function(y) load_obj(f=paste0(tables_dir, data_prefix,"_", x, "_", file_suffix, "_magma_GWAS_", y, "_", flag_date, ".csv"))))
  list_list_idx_signif_modules <- lapply(list_list_magma_tables, function(x) lapply(x, function(y) x[x[,max.col(x[,-which(colnames(x)=="module")])] > -log10(5^-2),][["module"]]))
  list_MAGMA_signif_modules <- sapply(list_list_idx_signif_modules, function(x) unique(unlist(x)))
  
  # Compute cell embeddings
  list_cell_embeddings_kME_MAGMA <- lapply(list_MAGMA_signif_modules, function(x) as.matrix(kMEs_across_subsets[,colnames(kMEs_across_subsets) %in% x]) %*% t(datExpr_full)) # module x gene %*% gene x cell -> module x cell
  
  # name the matrix
  #colnames
  # list_cell_embeddings_kME_MAGMA <- lapply(list_cell_embeddings_kME_MAGMA, function(x) name_for_vec(to_be_named = x, given_names = rownames(datExpr_full), dimension = 2)) # cell names
  # #rownames
  # list_cell_embeddings_kME_MAGMA <- lapply(list_cell_embeddings_kME_MAGMA, function(x) name_for_vec(to_be_named = x, given_names = paste(prefixes, colnames(x), sep="_"), dimension = 1))
  
  # group the columns of the embeddings matrix by cell type (@ident)
  celltypes_idx <- integer(0)
  for (celltype in names(table(seurat_obj@ident))) {
    celltypes_idx <- c(celltypes_idx, which(seurat_obj@ident == celltype)) 
  }
  
  list_cell_embeddings_kME_MAGMA <- lapply(list_cell_embeddings_kME_MAGMA, function(x) x[,celltypes_idx]) 
  
  # plot
  for (i in 1:length(gwas_sub_dirs)) {
    pheatmap(mat = list_cell_embeddings_kME_MAGMA[[i]], 
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))[0:100],
             kmeans_k = NA,
             breaks = NA, 
             cluster_rows = F,
             cluster_cols = F,
             #clustering_method = hclustMethod,
             show_rownames = T,
             show_colnames = F,
             main = sprintf("%s_%s_%s_%s", data_prefix, "module_embeddings_MAGMA", gwas_sub_dirs[i], flag_date),
             filename=sprintf("%s%s_%s_%s_%s.pdf", plots_dir, data_prefix, "module_embeddings_MAGMA", gwas_sub_dirs[i], flag_date),
             width = 10,
             height = length(list_MAGMA_signif_modules) %/% 4,
             silent = T)
  }
  
  if (!is.null(STRINGdb_species)) {
    
    file_suffix = "kMEs_PPI"
    # upper level list of subsets, nested list with tables for each GWAS repo
    list_list_PPI_magma_tables <- lapply(sNames_filter, function(x) lapply(gwas_sub_dirs, function(y) load_obj(file=paste0(tables_dir, data_prefix,"_", x, "_", file_suffix, "_magma_GWAS_", y, "_", flag_date, ".csv"))))
    list_list_idx_PPI_signif_modules <- lapply(list_list_PPI_magma_tables, function(x) lapply(x, function(y) x[x[,max.col(x[,-which(colnames(x)=="module")])] > -log10(5^-2),][["module"]]))
    list_PPI_MAGMA_signif_modules <- sapply(list_list_idx_PPI_signif_modules, function(x) unique(unlist(x)))
    
    # Compute cell embeddings
    list_cell_embeddings_kME_PPI_MAGMA <- lapply(list_PPI_MAGMA_signif_modules, function(x) as.matrix(kMEs_PPI_across_subsets[,colnames(kMEs_PPI_across_subsets) %in% x]) %*% t(datExpr_full)) # module x gene %*% gene x cell -> module x cell
    
    # name the matrix
    #colnames
    # list_cell_embeddings_kME_PPI_MAGMA <- lapply(list_cell_embeddings_kME_PPI_MAGMA, function(x) name_for_vec(to_be_named = x, given_names = rownames(datExpr_full), dimension = 2)) #cell names
    # #rownames
    # list_cell_embeddings_kME_PPI_MAGMA <- lapply(list_cell_embeddings_kME_PPI_MAGMA, function(x) name_for_vec(to_be_named = x, given_names = paste(prefixes, colnames(x), sep="_"), dimension = 1))
     
    # group the columns of the embeddings matrix by cell type (@ident)
    
    list_cell_embeddings_kME_PPI_MAGMA <- lapply(list_cell_embeddings_kME_PPI_MAGMA, function(x) x[,celltypes_idx]) 
    
    for (i in 1:length(gwas_sub_dirs)) {
      # plot
      pheatmap(mat = list_cell_embeddings_kME_PPI_MAGMA[[i]], 
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))[0:100],
               kmeans_k = NA,
               breaks = NA, 
               cluster_rows = F,
               cluster_cols = F,
               #clustering_method = hclustMethod,
               show_rownames = T,
               show_colnames = F,
               main = sprintf("%s_%s_%s_%s", data_prefix, "module_embeddings_PPI_MAGMA", gwas_sub_dirs[i], flag_date),
               filename=sprintf("%s%s_%s_%s_%s.pdf", plots_dir, data_prefix, "module_embeddings_PPI_MAGMA", gwas_sub_dirs[i], flag_date),
               width = 10,
               height = length(list_MAGMA_signif_modules) %/% 4,
               silent = T)
    }
  }
  ######################################################################
  ######## PLOT MODULE PRESERVATION IN DIFFERENT TISSUES ###############
  ######################################################################
  
  # Check preservation of each set of colors in all celltypes (including itself)
  module_preservation_celltypes <- sapply(list_colors_final, function(x) wrapModulePreservation(listDatExpr = list_datExpr_filter,
                                                                                                  listColors = list(x))[['preservation']][['log.pBonf']],
                                                            simplify = T)
  
  colnames(module_preservation_celltypes) = sNames_filter
  rownames(module_preservation_celltypes) = paste(prefixes, names(unlist(list_colors_final, recursive = F, use.names = T)), sep="_")
  
  # plot
  pheatmap(mat = module_preservation_celltypes, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))[0:100],
           kmeans_k = NA,
           breaks = NA, 
           cluster_rows = F,
           cluster_cols = F,
           #clustering_method = hclustMethod,
           show_rownames = T,
           show_colnames = T,
           main = sprintf("%s_%s_%s", data_prefix, "module preservation across celltypes", flag_date),
           filename=sprintf("%s%s_%s_%s.pdf", plots_dir, data_prefix, "module_preservation", flag_date),
           width = length(list_datExpr_filter) %/% 2,
           height = length(unlist(list_colors_final, recursive = F)) %/% 4,
           silent = T)
  
  if (!is.null(STRINGdb_species)) {
    module_PPI_preservation_celltypes <- sapply(list_colors_PPI_final, function(x) wrapModulePreservation(listDatExpr = list_datExpr_filter,
                                                                                                  listColors = list(x))[['preservation']][['log.pBonf']],
                                            simplify = T)
    
    colnames(module_PPI_preservation_celltypes) = sNames_filter
    rownames(module_PPI_preservation_celltypes) = paste(prefixes, names(unlist(list_colors_final_PPI, recursive = F, use.names = T)), sep="_")
    
    # plot
    pheatmap(mat = module_PPI_preservation_celltypes, 
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))[0:100],
             kmeans_k = NA,
             breaks = NA, 
             cluster_rows = F,
             cluster_cols = F,
             #clustering_method = hclustMethod,
             show_rownames = T,
             show_colnames = T,
             main = sprintf("%s_%s_%s", data_prefix, "PPI-enriched modules: preservation across celltypes", flag_date),
             filename=sprintf("%s%s_%s_%s.pdf", plots_dir, data_prefix, "PPI_module_preservation", flag_date),
             width = length(list_datExpr_filter) %/% 2,
             height = length(unlist(list_colors_PPI_final, recursive = F)) %/% 4,
             silent = T)
  }
} # end of if (FALSE) {make all the plots} TODO: remove

######################################################################
########################### CLEAR UP #################################
######################################################################

save.image(file=sprintf("%s%s_final_session_image.RData", RObjects_dir, data_prefix, flag_date))

# Delete checkpoint session images from disk

for (subsetName in sNames) {
  if (file.exists(sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, subsetName))) {
    file.remove(sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, subsetName))
  }
}

# for (i in 1:11) {
#   if (file.exists(sprintf("%s%s_checkpoint_%s_image.RData", RObjects_dir, data_prefix, i))) {
#     file.remove(sprintf("%s%s_checkpoint_%s_image.RData", RObjects_dir, data_prefix, i))
#   }
# }

rm(list=ls())
invisible(gc())

##########################################################################
################################ FINISH ##################################
##########################################################################

message("Script DONE!")
