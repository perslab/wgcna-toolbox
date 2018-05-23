# Title: Script to find robust WGCNA modules

######################################################################
############################## USAGE #################################
######################################################################

# e.g.
# CURRENT 180515
# TEST v1.8_dev3
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_n01_to_n04.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell_n01_to_n04-11/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-n01_to_n04-11 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.20)" --jackstraw.num.replicate 0 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_clust_all_sub_n.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-clust-all-sub-n-1/ --meta.data_corr_col 'c("nGene", "nUMI", "percent.mito", "percent.ribo", "X2.group", "X3.batches", "X4.sex", "X5.Diet", "X6.FvF")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-clust-all-sub-n-1 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 400 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 12 --resume checkpoint_2

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-21/ --meta.data_corr_col 'c("nGene", "nUMI", "percent.mito", "percent.ribo", "X2.group", "X3.batches", "X4.sex", "X5.Diet", "X6.FvF")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-AgRP-21 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 0 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 4 --resume checkpoint_3


# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_cell_ontology_class_filterRibo0.05_min20cell.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-5/ --meta.data_corr_col 'c("nGene", "nUMI")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix maca-5 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 500 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 20 --resume checkpoint_3
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_pancreas.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-pancreas-10/ --meta.data_corr_col 'c("percent.mito", "percent.ribo")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix maca-pancreas-10 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 500 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 14 --resume checkpoint_2
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/EX_so2.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-5/ --meta.data_corr_col 'c("dataset", "percent.mito", "nUMI")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix epilepsy-5 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 500 --TOM.num.replicate 0 --replace T --organism hsapiens --plot_permuted F --n_cores 10


# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_clust_all.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-clust-all-2/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-clust-all-1 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 400 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 21 --resume checkpoint_4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_clust_all_sub_n.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-clust-all-sub-n-2/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-clust-all-sub-n-2 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 300 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 12
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_pancreas.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-pancreas-9/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix maca-pancreas-9 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 400 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 14 --resume checkpoint_4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-4/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix epilepsy-4 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(3,4)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 500 --TOM.num.replicate 0 --replace T --organism hsapiens --plot_permuted F --n_cores 10
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_n01_to_n04.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell_n01_to_n04-10/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-n01_to_n04-10 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.20)" --jackstraw.num.replicate 0 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_n05_to_n10.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell_n05_to_n10-9/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-n05_to_n10-9 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.20)" --jackstraw.num.replicate 0 --TOM.num.replicate 0 --replace T --organism mmusculus --plot_permuted F --n_cores 6
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-6/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-AgRP-6 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(1,2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 0 --TOM.num.replicate 0 --replace T --organism mmusculus--plot_permuted F --n_cores 4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_neurons_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-neurons-17/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-neurons-17 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 500 --TOM.num.replicate 0 --replace T --organism mmusculus--plot_permuted F --n_cores 14 --resume checkpoint_4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v1.8_dev3.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_s_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-glia-17/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/ --data_prefix campbell-glia-17 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstraw.num.replicate 500 --TOM.num.replicate 0 --replace T --organism mmusculus--plot_permuted F --n_cores 11 --resume checkpoint_4
  
######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

if (FALSE) { 
  seurat_path = "/projects/jonatan/tmp-holst-hsl/RObjects/campbell_n05_to_n10.RData"
  project_dir = "/projects/jonatan/tmp-rwgcna-tests/tmp-campbell-n05_to_n10-10/"
  magma_gwas_dir = "/projects/jonatan/tmp-bmi-brain/data/magma/"
  data_prefix = "campbell-n10_to_n10-10"
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
  hminClusterSize = c(20)
  deepSplit = c(2)
  pamStage = c(TRUE)
  moduleMergeCutHeight = c(0.2)
  replace = T
  jackstraw.num.replicate = 0
  TOM.num.replicate = 0
  ### 180507 v1.8_dev2
  # STRINGdb_species = 10090
  # ensembl_dataset = "mmusculus_gene_ensembl"
  ##
  organism = "mmusculus"
  plot_permuted = F
  n_cores = 6
}

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

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
              help = "Resume from a previous session image? Must have same path and data_prefix. Options are 'checkpoint_1' - 'checkpoint_4'"),
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
  make_option("--jackstraw.num.replicate", type="integer", default=300L,
              help = "Number of times to resample to make null distributions for empirical significance tests in JackStraw and other functions. Non-negative integer, 0 to skip Jackstraw, defaults to 300."),
  make_option("--TOM.num.replicate", type="integer", default=100L,
              help = "Number of times to resample the dataset when finding the consensus TOM, defaults to 100"),
  make_option("--replace", type="logical", default=TRUE,
              help = "Sample with replacement? Defaults to TRUE. If TRUE, uses all samples, if FALSE, uses 66% each time."),
  # 180507 v1.8_dev2
  # make_option("--STRINGdb_species", type="integer", default="10090",
  #             help = "Optional. Species for which to retrieve protein data from STRINGdb to validate clusters. Mouse (mus musculus) is 10090, homo sapiens is 9606"),
  # make_option("--ensembl_dataset", type="character", default=NULL,
  #             help = "Optional. Dataset for mapping gene symbols to ensembl_ID. Mouse is mmusculus_gene_ensembl, homo sapiens is hsapiens_gene_ensembl"),
  make_option("--organism", type="character", default="mmusculus",
              help = "'hsapiens' or 'mmusculus'"),
  ###  
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
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(WGCNA)) 
#suppressPackageStartupMessages(library(STRINGdb)) 
suppressPackageStartupMessages(library(liger))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(qusage))

#suppressPackageStartupMessages(library(data.table))

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

jackstraw.num.replicate <- opt$jackstraw.num.replicate

TOM.num.replicate <- opt$TOM.num.replicate

replace <- opt$replace

### 180507_v1.8_dev2
# STRINGdb_species <- opt$STRINGdb_species
# 
# ensembl_dataset <- opt$ensembl_dataset
organism <- opt$organism
###
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
source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params_v1.8_dev3.R")
source(file = "/projects/jonatan/functions-src/functions_v1.8_dev3.R")

if (organism == "mmusculus") {
  # source: South Dakota State University 
  # http://ge-lab.org/gskb/
  mm_geneset_path <- "/projects/jonatan/genesets/mGSKB_Ensembl.gmt"
} 

if (TOM.num.replicate < 1) plot_permuted <- F 
  
# Checkpoint 0
if (is.null(resume)) {
  message("Loading and subsetting seurat object..")
  
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
  
  if (! (jackstraw.num.replicate >= 0)) stop("jackstraw.num.replicate must be 0 or higher")
  
  if (! (TOM.num.replicate >= 0)) stop("TOM.num.replicate must be non-negative")
  
  if (plot_permuted == T) warning("plotting modules on permuted datasets increases processing time")
  
  if (! (n_cores >= 0 & n_cores <= 100)) stop("n_cores must be in the range 0-100")
  
  ######################################################################
  ################# LOAD AND SUBSET SEURAT OBJECT ######################
  ######################################################################

  seurat_obj <- load_obj(f=seurat_path)
  
  # 180507 v1.8_dev2
  
  if (!any(c("ENSG", "ENSMUSG") %in% rownames(seurat_obj@data))) { 
    
    if (is.null(organism)) {
      stop("Error: seurat object gene names are not in ensembl ID format. To remap from hgcn to ensembl, please provide an organism 'mmusculus' or 'hsapiens'")
    
      } else if (!is.null(organism)) {
      # Map to ensembl_id
      if (organism == "mmusculus") { 
        
        STRINGdb_species <- 10090
        
        mapping_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
        mapping_mm_synonyms_filepath = "/data/genetic-mapping/ncbi/Mus_musculus.gene_info_symbol2ensembl.gz"
        # Step 1: direct mapping
        mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T)
        mapping = data.frame(symbol=row.names(seurat_obj@data), ensembl = mapping_direct$ensembl_gene_id[ match(row.names(seurat_obj@data), mapping_direct$gene_name_optimal) ])
        # Step 2: map remaing using synonyms
        
        ### 180508_v1.8_dev2
        #mapping_synonyms = read.csv(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T)
        mapping_synonyms = read.delim(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T)
        ###
        mapping$ensembl[ which(is.na(mapping$ensembl)) ] = mapping_synonyms$ensembl[ match( mapping$symbol[which(is.na(mapping$ensembl)) ] ,mapping_synonyms$symbol) ]
        
      } else if (organism == "hsapiens") {
        
        STRINGdb_species <- 9606

        mapping_hs_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz" # columns: ensembl_gene_id, entrezgene, hgcn_symbol
        mapping_direct = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T) # columns: ensembl_gene_id, entrezgene, hgcn_symbol
        # Step 1: direct mapping
        mapping = data.frame(symbol=row.names(seurat_obj@data), ensembl = mapping_direct$ensembl_gene_id[ match(row.names(seurat_obj@data), mapping_direct$hgcn_symbol) ])
   
      }

      # Make a log of unmapped genes 
      log_not_mapped_filepath = sprintf("%s%s_not_mapped_to_ensembl_%s", log_dir, data_prefix, flag_date,".tab")
      log_duplicate_filepath = sprintf("%s%s_duplicate_ensembl_id_genenames_%s", log_dir, data_prefix, flag_date,".tab")
      
      df_not_mapped = mapping[is.na(mapping$ensembl),]
      write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
      
      # Make a log of duplicate genes
      idx_duplicate_genes <- duplicated(mapping$ensembl)
      df_duplicate <- mapping[idx_duplicate_genes,]
      write.table(df_duplicate,log_duplicate_filepath,quote=F,sep="\t",row.names=F)
      
      # Filter out unmapped and duplicate genes from Seurat object
      seurat_obj@data <- seurat_obj@data[!is.na(mapping$ensembl) & !idx_duplicate_genes,]
      seurat_obj@raw.data <- seurat_obj@raw.data[!is.na(mapping$ensembl) & !idx_duplicate_genes,]

      # rename Seurat object rows where mapping was successful to ensembl ID
      rownames(seurat_obj@data) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
      rownames(seurat_obj@raw.data) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
      
      # Save ensembl seurat object for later use
      save(seurat_obj, file=sprintf("%sseurat_obj_ensembl.RData",RObjects_dir))
      
    }
  } else if (any(grepl("ENSG", rownames(seurat_obj@data)))) {
    message("Homo sapiens ensembl id gene names detected, skipping remapping from hgnc to ensembl")
    organism <- "hsapiens"
    
  } else if (any(grepl("ENSMUSG", rownames(seurat_obj@data)))) {
    message("mus musculus ensembl id gene names detected, skipping remapping from hgnc to ensembl")
    organism <- "mmusculus"
  }
  
  # Use imputed data if user opted to do so 
  if (is.null(seurat_obj@imputed) | any(dim(seurat_obj@imputed)) == FALSE) use.imputed <- F
  if (use.imputed==T) seurat_obj@data <- seurat_obj@imputed
  
  # Subset the seurat object
  if (!is.null(meta.data_subset_col)) {
    seurat_obj <- SetAllIdent(object = seurat_obj,id = meta.data_subset_col)
  }
  sNames <- names(table(seurat_obj@ident))
  
  all_genes <- rownames(seurat_obj@data)
  
  # Convert any character or factor meta.data to numeric dummy variables each level with its own numeric column
  if (!is.null(meta.data_corr_col)) {
    if (all(sapply(meta.data_corr_col, function(x) !is.null(seurat_obj@meta.data[[x]])))) {
      meta.data <- SeuratFactorToIndicator(obj = seurat_obj,
                                            meta.colnames = meta.data_corr_col)
    } else {
      message('Did not find all columns in meta.data_corr_col in the Seurat object meta.data slot')
      meta.data <- NULL
    }
  }

  # Subset the seurat object
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
                           "jackstraw.num.replicate" = jackstraw.num.replicate,
                           "TOM.num.replicate" = TOM.num.replicate,
                           "replace" = replace,
                           
                           ### 180507 v1.8_dev2
                           #"STRINGdb_species" = STRINGdb_species,
                           #"ensembl_dataset" = ensembl_dataset,
                           "organism" = organism,
                           ###
                           
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

  vars.to.regress = c("nUMI", "percent.mito", "percent.ribo")[c("nUMI", "percent.mito", "percent.ribo") %in% names(seurat_obj@meta.data)]

  subsets <- lapply(subsets, function(x) ScaleData(object = x,
                                                  vars.to.regress = if (length(vars.to.regress) > 0) vars.to.regress else NULL,
                                                  model.use="linear",
                                                  do.par=T,
                                                  num.cores = min(n_cores, detectCores()-1),
                                                  do.scale=T,
                                                  do.center=do.center))

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
    if (jackstraw.num.replicate > 0) message(sprintf("Performing JackStraw with %s replications to select genes that load on significant PCs", jackstraw.num.replicate))
    list_datExpr <- lapply(subsets, function(x) wrapJackStraw(seurat_obj_sub = x, n_cores = n_cores, jackstraw.num.replicate = jackstraw.num.replicate))
    
    invisible(gc())
    
  } else if (genes_use == "all") {
    list_datExpr <- lapply(subsets, function(x) seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rep(TRUE, nrow(x@scale.data))))
  } else if (genes_use == "var.genes") {
    list_datExpr <- lapply(subsets, function(x) seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rownames(x@scale.data) %in% x@var.genes))
  } 
  
  # Clear up
  rm(subsets, userParams, builtInParams, subsets_n_cells, seurat_obj)
  invisible(gc())
  
  # Save or load session image 
  resume = "checkpoint_1"
  # just in case, since JackStraw is slow
  save.image(file=sprintf("%s%s_checkpoint_1_image.RData", RObjects_dir, data_prefix))
  
  # TODO: unload packages
  
  
} else if (!is.null(resume)) {
  
  if (resume == "checkpoint_1") {
    
    load((file=sprintf("%s%s_checkpoint_1_image.RData", RObjects_dir, data_prefix)))
    source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params_v1.8_dev3.R")
    source(file = "/projects/jonatan/functions-src/functions_v1.8_dev3.R")
    TOM.num.replicate=0
  }
  
}

if (resume == "checkpoint_1") {

  ######################################################################
  ####### PICK SOFT THRESHOLD POWER FOR ADJACENCY MATRIX ###############
  ######################################################################
  message("Computing soft threshold powers to maximise the fit of a scale free topology to the adjacency matrix")
  
    list_softPower <- mapply(function(x,y) sft_for_par(datExpr=x, subsetName=y), list_datExpr, sNames, SIMPLIFY=F)
  
    ###
    invisible(gc())
    names(list_softPower) <- sNames
    
    # TOM files will be saved here by consensusTOM / blockwiseConsensusModules
    setwd(RObjects_dir)
    
    if (TOM.num.replicate > 0) {
      
      ######################################################################
      ############### RESAMPLE THE DATA FOR ROBUSTNESS #####################
      ######################################################################
      
      message(sprintf("Resampling the data %s times with replace = %s", TOM.num.replicate, replace))
      
      invisible(gc())
      
      cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_bootStrap.txt"))
      
      list_multiExpr <- parLapply(cl, list_datExpr, function(x) bootstrap(datExpr=x, 
                                                                          ### 180508_v1.8_dev2
                                                                          nPermutations = TOM.num.replicate,
                                                                          ###
                                                                          replace = replace,
                                                                          fraction = fraction,
                                                                          randomSeed = randomSeed))
      stopCluster(cl)
      
      invisible(gc())
      
      names(list_multiExpr) <- sNames 
    
      ######################################################################
      ######### RUN CONSENSUSTOM ON THE RESAMPLED DATASETS #################
      ######################################################################
      
      message(sprintf("Computing the consensus Topological Overlap Matrix with %s permutations", TOM.num.replicate))
    
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
      list_datExpr_gg <- mapply(FUN=function(x,y) x[,y], x=list_datExpr, y=list_goodGenesTOM_idx, SIMPLIFY=F)
      
    } else if (TOM.num.replicate==0) {
      
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
  }
  
  message("Done computing the consensus Topological Overlap Matrix")
  
  resume="checkpoint_2"
    
  if (TOM.num.replicate > 0) rm(list_multiExpr) #TODO: we'll need it for plotting permuted
  
  save.image(file=sprintf("%s%s_checkpoint_2_image.RData", RObjects_dir, data_prefix))

} else if (!is.null(resume)) {
  if (resume == "checkpoint_2") {
    load((file=sprintf("%s%s_checkpoint_2_image.RData", RObjects_dir, data_prefix)))
    source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params_v1.8_dev3.R")
    source(file = "/projects/jonatan/functions-src/functions_v1.8_dev3.R")
  }
}

if (resume == "checkpoint_2") {
 
  ######################################################################
  ####################### CLUSTER ON THE TOM ###########################
  ######################################################################
  # in:
  #   consTomDS (from disk!)
  #
  # out:
  #   list_geneTree
  
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
  # in:
  #   list_dissTOM
  #   list_geneTree
  #
  # out:
  #   list_list_plot_label_ok 
  #   list_list_MEs_ok 
  #   list_list_colors_ok 
  #   sNames_ok 
  #   list_datExpr_ok 
  #   list_geneTree_ok
  #   list_dissTOM_ok 

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
  
  names(list_list_cutree) = sNames
  
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
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  list_logical_params_ok <- mapply(function(x,y) as.logical(mapply(function(a,b) a!=length(b), a=x, b=y, SIMPLIFY=F)), x=list_vec_n_grey, y=list_list_colors, SIMPLIFY=F)
  logical_subsets_ok <- sapply(list_logical_params_ok, any, simplify = T)
  
  # First, for each run, remove results of any parametrisation for which all the genes were assigned to grey
  # If for a subset all parameters gave only grey modules, take it out of the top level list.
  ### 180510_v1.8_dev2_1
  list_list_plot_label_ok <- mapply(function(x,y) x[y], x = list_list_plot_label, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_cutree_ok <- mapply(function(x,y) x[y], x = list_list_cutree, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_MEs_ok <- mapply(function(x,y) x[y], x = list_list_MEs, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_colors_ok <- mapply(function(x,y) x[y], x=list_list_colors, y=list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]

  sNames_ok <- sNames[logical_subsets_ok]
  list_datExpr_ok <- list_datExpr_gg[logical_subsets_ok] # If for a subset all parameters gave only grey modules, take it out of the top level list.
  list_geneTree_ok <- list_geneTree[logical_subsets_ok]
  list_dissTOM_ok <- list_dissTOM[logical_subsets_ok]
  
  # Assign gene names to each color vector
  list_list_colors_ok <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), dimension = NULL)), x = list_list_colors_ok, y = list_datExpr_ok, SIMPLIFY=F)
  
  # Make a warning and report if any whole subsets produced no modules
  discarded <- setdiff(sNames, sNames_ok)
  if (length(discarded > 0)) {
    warning(sprintf("The %s were dropped because the script didn't find any modules", discarded))
    fileConn<-file(sprintf("%s%s_%s_ERROR_report_%s.txt", log_dir,  data_prefix, discarded, flag_date ))
    writeLines(sprintf("The %s were dropped because the script didn't find any modules", discarded), fileConn)
    close(fileConn)
  }
  
  # clear up
  rm(list_datExpr, list_list_cutree, list_list_MEs, list_datExpr_gg, list_geneTree, list_list_plot_label, list_dissTOM)#, list_idx_mapped_gg)
  invisible(gc())
    
  ######################################################################
  ###################### compute kMEs & pkMEs ##########################
  ######################################################################
  # in:
  #   list_list_MEs_ok
  #   list_datExpr_ok
  #
  # out:
  #   list_list_kMEs
  #   list_list_pkMEs
  
  message("Computing kMEs: correlation between each gene and each module 'eigengene'")
 
  # Compute kMEs
  invisible(gc())
  cl <- makeCluster(spec=n_cores, type = "FORK", outfile = paste0(log_dir, "log_parkMEs_parPkMEs.txt"))
  list_list_kMEs <- clusterMap(cl, function(x,y) parkMEs(list_MEs=x, datExpr=y), x=list_list_MEs_ok, y=list_datExpr_ok, SIMPLIFY = F)
  
  # Extract the primary kMEs - i.e. those of each gene w.r.t. the module to which it belongs
  # We use these to filter genes that we submit to STRINGdb 
  list_list_pkMEs <- clusterMap(cl, function(x,y) parPkMEs(list_kMEs=x, list_colors=y), x=list_list_kMEs, y=list_list_colors_ok, SIMPLIFY = F)
  stopCluster(cl)
  invisible(gc())
  
  list_list_pkMEs <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_pkMEs, y=list_list_plot_label_ok, SIMPLIFY=F)
  names(list_list_pkMEs) <- sNames_ok
  
  ######################################################################
  ########################## compute kIMs ##############################
  ######################################################################
  # in:
  #   list_dissTOM_ok
  #   list_list_colors_ok
  #
  # out:
  #   list_list_kIMs
  
  # message("Computing normalised Intra-module connectivity, kIM: each gene's degree with respect to each module")
  # 
  # invisible(gc())
  # 
  # cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_kIMs.txt"))
  # 
  # list_list_kIMs <- clusterMap(cl, function(x,y) lapply(y, function(z) kIM_eachMod_norm(dissTOM = x, 
  #                                                                               colors = z,
  #                                                                               genes = names(z))),
  #                      x = list_dissTOM_ok,
  #                      y = list_list_colors_ok, 
  #                      SIMPLIFY=F)
  # 
  # stopCluster(cl)
  # invisible(gc())

  ######################################################################
  ############## Match colors between parameter settings ###############
  ######################################################################
  # in:
  #   list_list_colors_ok
  #   list_list_MEs_ok
  #   list_list_kMEs
  #   list_list_kIMs
  #
  # out:
  #   list_list_colors_ok_matched
  #   list_list_MEs_ok_matched
  #   list_list_kMEs_matched
  #   list_list_kIMs_matched
  
  message("Matching module color labels between parameter settings")

  # Align / match the colors so similar modules found with different sets of parameters have the same name
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parMatchColors.txt"))
  list_list_colors_ok_matched <- parLapply(cl, list_list_colors_ok, parMatchColors)
  stopCluster(cl)
  invisible(gc())
  
  # Rename highest level list entries (cell clusters)
  names(list_list_colors_ok_matched) = sNames_ok
  
  # Name each entry of the vectors of color assignments with the corresponding genes
  list_list_colors_ok_matched <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), dimension = NULL)), x=list_list_colors_ok_matched, y=list_datExpr_ok, SIMPLIFY=F)
  
  ### no need to match MEs, kMEs and kIMs as within each cell cluster we just need pkMEs for PPI, then
  # we remove all but one parameter set anyway
  #
  # Rename the eigengenes according to the matched color scheme - nope, not necessary since we won't be using them
  #list_list_MEs_ok_matched <- mapply(function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x$eigengenes, given_names= sort(names(table(y))), dimension = 2), x=a, y=b, SIMPLIFY=F), a= list_list_MEs_ok, b = list_list_colors_ok_matched, SIMPLIFY=F)
  # Rename the columns of the kMEs and kIMs according to the matched color scheme
  #list_list_kMEs_matched <- mapply(function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x, given_names = sort(names(table(y))), dimension = 2), x = a, y = b, SIMPLIFY=F), a = list_list_kMEs, b = list_list_colors_ok_matched, SIMPLIFY=F)
  #list_list_kIMs_matched <- mapply(function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x, given_names = sort(names(table(y))), dimension = 2), x = a, y = b, SIMPLIFY=F), a = list_list_kIMs, b = list_list_colors_ok_matched, SIMPLIFY=F)
  # Rename the nested list entries
  #list_list_kMEs_matched <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_kMEs_matched, y=list_list_plot_label_ok, SIMPLIFY=F)
  #list_list_kIMs_matched <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_kIMs_matched, y=list_list_plot_label_ok, SIMPLIFY=F)

  # Rename the upper level list entries, i.e. the nested lists
  #names(list_list_kMEs_matched) <- sNames_ok
  #names(list_list_kIMs_matched) <- sNames_ok
  
  # No need to rename the pkME modules since we don't query or display them, the values are all that matters
  
  # clear up 
  rm(list_list_colors_ok)
  invisible(gc())

  ######################################################################
  #################### CHECK MODULES FOR PPI ENRICHMENT ################
  ######################################################################
  # in:
  #   list_list_colors_ok_matched
  #   list_list_pkMEs
  #
  # out:
  #   list_list_colors_ok_matched_PPI
  
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
                                     a = list_list_colors_ok_matched, 
                                     b = list_list_pkMEs, 
                                     SIMPLIFY=F)
  stopCluster(cl)
  invisible(gc())
  
  # assign the gene names to the color vectors
  # list_list_colors_PPI <- lapply(list_list_colors_PPI, function(a) mapply(function(x,y) name_for_vec(to_be_named = x, given_names = colnames(y), dimension = NULL), 
  #                           x = a, 
  #                           y = list_datExpr_ok, 
  #                           SIMPLIFY=F))
  list_list_colors_PPI <- mapply(function(x,y) lapply(x, function(a) name_for_vec(to_be_named = a, given_names = colnames(y), dimension = NULL)), 
                                                                          x = list_list_colors_PPI, 
                                                                          y = list_datExpr_ok, 
                                                                          SIMPLIFY=F)
  ######################################################################
  ########### ORDER PARAMETER SETS BY PPI ENRICHMENT AND PLOT ##########
  ######################################################################
  # in:
  #   list_list_colors_ok_matched
  #   list_list_colors_ok_matched_PPI
  #   list_list_MEs_ok_matched
  #   list_list_kMEs_matched
  #   list_list_pkMEs
  #   list_list_kIMs_matched
  #   list_list_plot_label_ok
  #   
  # out:
  #   list_list_colors_ok_matched_PPI_order
  #   list_list_plot_label_ok_order
  #   list_list_colors_ok_matched_order
  #   list_list_MEs_ok_matched_order
  #   list_list_kMEs_matched_order
  #   list_list_pkMEs_order
  #   list_list_kIMs_matched_order
  #   plots of pre-PPI ranked colorings
  #   plots of post-PPI ranked colorings

  message("Selecting parameters with the highest number of genes assigned to modules significantly enriched for Protein-Protein Interactions")
    
  # Count how many genes were assigned to grey under each parametrisation and order the parametrisations
  list_PPI_vec_n_grey <- lapply(list_list_colors_PPI, count_grey_in_list_of_vec)

  # Order all the outputs by how many genes were assigned to a (non-grey) module
  list_list_plot_label_ok_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_plot_label_ok, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_colors_ok_matched_order <- mapply(function(x,y) x[order(y, decreasing=F)], x  =  list_list_colors_ok_matched , y = list_PPI_vec_n_grey, SIMPLIFY = F )
  # list_list_MEs_ok_matched_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_MEs_ok_matched, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  # list_list_kMEs_matched_order <- mapply(function(x,y) x[order(y, decreasing=F)],  x = list_list_kMEs_matched, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  # list_list_kIMs_matched_order <- mapply(function(x,y) x[order(y, decreasing=F)],  x = list_list_kIMs_matched, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  # list_list_pkMEs_order <- mapply(function(x,y) x[order(y, decreasing=F)],  x  = list_list_pkMEs, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_colors_PPI_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_colors_PPI, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  
  # # Name upper level list entries (cell clusters)
  # for (obj in list(list_list_plot_label_ok_order, 
  #                  list_list_colors_ok_matched_order, 
  #                  list_list_MEs_ok_matched_order, 
  #                  list_list_kMEs_matched_order,
  #                  list_list_kIMs_matched_order,
  #                  list_list_pkMEs_order,
  #                  list_list_colors_ok_matched_PPI_order)) {
  #   names(obj) <- sNames_ok
  # }
  
  ######################################################################
  ############## PLOT THE RANKED MODULE ASSIGNMENTS ####################
  ######################################################################

  # Pre-PPI 
  invisible(mapply(function(x,y) parPlotDiffParams(list_colors=x, subsetName=y), x=list_list_colors_ok_matched_order, y=sNames_ok, SIMPLIFY = F))
  # Post-PPI 
  invisible(mapply(function(x,y) parPlotDiffParams_PPI(list_colors=x, subsetName=y), x=list_list_colors_PPI_order, y=sNames_ok, SIMPLIFY = F))
  
  ######################################################################
  ##### FOR EACH SUBSET SELECT PARAMETERS WITH BEST PPI ENRICHMENT #####
  ######################################################################
  
  # Eliminate a layer of nesting by selecting only the best parametrisation per celltype
  list_plot_label_final <- lapply(list_list_plot_label_ok_order, function(x) x[[1]])
  list_colors_PPI <- lapply(list_list_colors_PPI_order, function(x) x[[1]])
  list_colors <- lapply(list_list_colors_ok_matched_order, function(x) x[[1]])
  # list_MEs <- lapply(list_list_MEs_ok_order, function(x) x[[1]])
  # list_kMEs <- lapply(list_list_kMEs_matched_order, function(x) x[[1]])
  # list_kIMs <- lapply(list_list_kIMs_matched_order, function(x) x[[1]])
  
  # Name by cell clusters
  names(list_plot_label_final) <- sNames_ok
  names(list_colors_PPI) <- sNames_ok
  names(list_colors) <- sNames_ok
  # names(list_MEs) <- sNames_ok
  # names(list_kMEs) <- sNames_ok
  # names(list_kIMs) <- sNames_ok  

  # Name the kME and KIM rows 
  # list_kMEs <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kMEs, y=list_colors, SIMPLIFY=F)
  # list_kIMs <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kIMs, y=list_colors, SIMPLIFY=F)

  # Make list of list of final parameters
  param_names = c("minClusterSize", "deepSplit","pamStage", "moduleMergeCutHeight")
  list_list_cutree_params_final <- lapply(list_PPI_vec_n_grey, function(x) comb_list[order(x,decreasing = F)][[1]])
  list_list_cutree_params_final <- lapply(list_list_cutree_params_final, function(x) name_for_vec(to_be_named = x, given_names = param_names, dimension = NULL)) 
  
  ######################################################################
  ######### PLOT FINAL MODULE ASSIGNMENT PRE- & POST-PPI FILTER ########
  ######################################################################
  
  message("Plotting")

  list_colors_both = mapply(FUN=cbind, list_colors, list_colors_PPI, SIMPLIFY = F)
  names(list_colors_both) = sNames_ok
  invisible(lapply(sNames_ok, plotFinalColors))
  
  ######################################################################
  ################### DELETE PPI RUNS WITH ONLY GREY ###################
  ######################################################################
  # in:
  #   list_colors_PPI
  #   list_datExpr_ok
  #   sNames_ok
  #
  # out:
  #   list_colors_PPI_ok
  #   list_datExpr_PPI_ok
  #   sNames_PPI_ok
  
  PPI_vec_n_grey <- count_grey_in_list_of_vec(list_colors_PPI)
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  logical_subsets_PPI_ok <- as.logical(mapply(function(x,y) x != length(y), x = PPI_vec_n_grey, y = list_colors_PPI, SIMPLIFY = F))
  
  sNames_PPI_ok <- sNames_ok[logical_subsets_PPI_ok]
  list_colors_PPI_ok <- list_colors_PPI[logical_subsets_PPI_ok]
  list_datExpr_PPI_ok <- list_datExpr_ok[logical_subsets_PPI_ok]
  list_dissTOM_PPI_ok <- list_dissTOM_ok[logical_subsets_PPI_ok]
  
  ######################################################################
  ########## RECOMPUTE MEs, kMEs, kIMs, pkMEs AFTER PPI FILTER #########
  ######################################################################
  # in:
  #   list_colors_PPI_ok
  #   list_datExpr_PPI_ok
  #   list_dissTOM_PPI_ok
  #
  # out:
  #   list_MEs_PPI
  #   list_kMEs_PPI
  #   list_kIMs_PPI
  
  message("Re-computing kME, primary kMEs and kIMs after filtering modules for significant Protein-Protein Interactions")
  
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_moduleEigengenes_kMEs_pkMEs_PPI.txt"))

  list_MEs_PPI <- clusterMap(cl, function(x,y) moduleEigengenes(expr=as.data.frame(x,col.names=col.names(x)),
                                                                      colors=y,
                                                                      excludeGrey=T), 
                                   x = list_datExpr_PPI_ok, 
                                   y = list_colors_PPI_ok, 
                                   SIMPLIFY = F)
  
  #names(list_MEs_PPI_final) <- sNames_ok
  list_kMEs_PPI <- clusterMap(cl, function(x,y) signedKME(as.matrix(x),
                      y$eigengenes,
                      outputColumnName = "",
                      corFnc = corFnc), 
                      x=list_datExpr_PPI_ok, 
                      y=list_MEs_PPI)
 
  list_kIMs_PPI <- clusterMap(cl, function(x,y) kIM_eachMod_norm(dissTOM=x, colors=y, genes=names(y)),
                              x = list_dissTOM_PPI_ok,
                              y = list_colors_PPI_ok,
                              SIMPLIFY=F)
  
  #list_kIMs_PPI <- deleteGrey(list_kIMs_PPI)
  
  stopCluster(cl)
  invisible(gc())
  
  # Name 
  list_kMEs_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kMEs_PPI, y=list_colors_PPI_ok, SIMPLIFY=F)
  list_kIMs_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kIMs_PPI, y=list_colors_PPI_ok, SIMPLIFY=F)

    
  # clear memory
  rm(list_dissTOM_PPI_ok, 
     list_list_plot_label_ok_order, 
     list_list_colors_PPI_order)
  invisible(gc()) 
  
  ###180509_v1.8_dev2
  # There should be no need since we have excludeGrey = T. 
  # But it may be necessary because there are definitely still grey modules
  list_kMEs_PPI_ok <- deleteGrey(list_kMEs_PPI) %>% Filter(f=length)
  list_kIMs_PPI_ok <- deleteGrey(list_kIMs_PPI) %>% Filter(f=length)

  ### 180515 v1.8_dev3 
  ### only output BMI module kMEs
  #list_kMEs_PPI_ok <- lapply(list_kMEs_PPI_ok, function(x) cbind(genes=rownames(x), x))
  #list_kIMs_PPI_ok <- lapply(list_kIMs_PPI_ok, function(x) cbind(genes=rownames(x), x))
  #rownames(list_kMEs_PPI_out) <- NULL
  #rownames(list_kIMs_PPI_out) <- NULL
  
  #list_mapping_PPI_out <- lapply(list_colors_PPI_ok, function(x) data.frame(module = x, ensembl_ID = names(x), hgnc_symbol = mapping$symbol[match(names(x), mapping$ensembl)], row.names = NULL))

  #invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kMEs_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_kMEs_PPI_out, sNames_PPI_ok, SIMPLIFY = F))
  #invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kIMs_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_kIMs_PPI_out, sNames_PPI_ok, SIMPLIFY = F))
  
  #invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_mapping_PPI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_mapping_PPI_out, sNames_PPI_ok, SIMPLIFY = F))

  # checkpoint
  resume = "checkpoint_3"
  message("Reached checkpoint 3, saving session image")
  save.image(file=sprintf("%s%s_checkpoint_3_image.RData", RObjects_dir, data_prefix))
  
  # Delete consensus TOMs from disk
  for (subsetName in sNames) {
    if (file.exists(sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, subsetName))) {
      file.remove(sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, subsetName))
    }
  }
  
} else if (resume == "checkpoint_3") {
  load(file=sprintf("%s%s_checkpoint_3_image.RData", RObjects_dir, data_prefix))
  # Load parameter values and utility functions anew (in case a bug was fixed)
  source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params_v1.8_dev3.R")
  source(file = "/projects/jonatan/functions-src/functions_v1.8_dev3.R")
}  

if (resume == "checkpoint_3") {

  ##########################################################################
  #################### COMPUTE MAGMA GWAS ENRICHMENT #######################
  ##########################################################################
  list_kMEs_PPI_ok <- deleteGrey(list_kMEs_PPI) %>% Filter(f=length)
  list_kIMs_PPI_ok <- deleteGrey(list_kIMs_PPI) %>% Filter(f=length)
  
  message("Scoring PPI enriched modules for enrichment with genes linked by GWAS to phenotypes of interest")
    
  file_suffix = "kMEs_PPI"    
  
  #list_dat_BMI <- magma_par(subsetNames = sNames_PPI_ok,
  list_dat_BMI <- magma_par(subsetNames = names(list_kMEs_PPI_ok),                         
                            list_kMEs = list_kMEs_PPI_ok,
                            project_dir = project_dir,
                            plots_dir = plots_dir,
                            log_dir = log_dir,
                            tables_dir = tables_dir,
                            magma_gwas_dir = magma_gwas_dir,
                            data_prefix = data_prefix,
                            file_suffix = file_suffix,
                            flag_date = flag_date,
                            organism = organism)
  
  if (is.null(list_dat_BMI)) stop("No modules enriched for BMI")

  dat_BMI_fdr <-  list_dat_BMI[["fdr"]]
  dat_BMI_corrCoef <-  list_dat_BMI[["corrCoef"]]

  # get a new vector of names of cell clusters with BMI enrichment
  modRegExpr <- ""
  for (i in 1:length(unique(unlist(list_colors_PPI_ok)))) {
    if (i<length(unique(unlist(list_colors_PPI_ok)))) {
        suffix <- "|"
    } else if (i == length(unique(unlist(list_colors_PPI_ok)))) {
      suffix <- ""
    }
    modRegExpr <- paste0(modRegExpr, unique(unlist(list_colors_PPI_ok))[i], suffix)
  }
  
  sNames_BMI <- unique(gsub(paste0(modRegExpr), "", dat_BMI_fdr$module))
  sNames_BMI <- gsub("_$", "", sNames_BMI)
  
  ### 180515 v1.8_dev3
  # rather make a data frame
  #~~Prepare nested lists of BMI-enriched module genes in both ensembl and HGNC symbol~~
  # Prepare 'melted' dataframe with modules and genes in ensembl and symbol
  ### 
  
  list_modules_BMI <- vector(mode="list", length=length(sNames_BMI))
  names(list_modules_BMI) <- sNames_BMI

  for (name in sNames_BMI) {
    idx <- grep(pattern = paste0(name, "_.*"), x = dat_BMI_fdr$module)
    list_modules_BMI[[name]] <- sapply(idx, function(j) gsub(".*_", "", dat_BMI_fdr$module[j]), simplify = T)
    #sNames_BMI[j] <- sNames_PPI_ok[j]
  }

  list_colors_BMI <- list_colors_PPI_ok[names(list_colors_PPI_ok) %in% sNames_BMI]
  
  # Prepare nested lists of module genes
  list_list_module_BMI_genes <- mapply(function(a,b) lapply(b, function(x) names(a)[a==x]), a=list_colors_BMI, b=list_modules_BMI, SIMPLIFY=F)
  list_list_module_BMI_genes <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = NULL), x=list_list_module_BMI_genes, y=list_modules_BMI)

  ######################################################################
  ########## PLOT CELL EXPRESSION ON MODULES HEAT MAP ##################
  ######################################################################
  
  # Output heatmap pdf and return nested lists of genes of BMI-enriched modules
  eigen_mat <- eigenPlotAllSamples(RObjects_dir = RObjects_dir,
                                   list_list_module_BMI_genes = list_list_module_BMI_genes,
                                   list_colors_BMI = list_colors_BMI, 
                                   list_modules_BMI = list_modules_BMI,
                                   dat_BMI_fdr = dat_BMI_fdr,
                                   dat_BMI_corrCoef = dat_BMI_corrCoef,
                                   n_cores = n_cores,
                                   plots_dir = plots_dir, 
                                   data_prefix = data_prefix, 
                                   flag_date = flag_date)

  ##########################################################################
  ######### PREPARE GENES LISTS AND DATAFRAME WITH MODULES, GENES ##########
  ##########################################################################
  
  list_list_module_BMI_genes_hgnc <- lapply(list_list_module_BMI_genes, function(x) lapply(x, function(y) mapping$symbol[match(y, mapping$ensembl)]))
  list_list_module_BMI_genes_hgnc <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = NULL), x=list_list_module_BMI_genes_hgnc, y=list_modules_BMI)
  
  # Prepare module genes dataframe
  cell_cluster <- rep(sNames_BMI, unlist(sapply(list_list_module_BMI_genes, FUN=function(x) sum(sapply(x, function(y) length(y), simplify=T)), simplify=T)))
  module <- unlist(sapply(list_list_module_BMI_genes, function(x) rep(names(x), sapply(x, function(y) length(y), simplify = T)), simplify=T), use.names = F)
  ensembl <- unlist(list_list_module_BMI_genes, recursive = T, use.names = F)
  hgnc <- unlist(list_list_module_BMI_genes_hgnc, recursive=T, use.names=F)
  
  df_BMI_module_genes <- data.frame(cell_cluster, module, ensembl, hgnc, row.names = NULL)
  
  # Also filter kME lists
  # upper level
  list_kMEs_BMI <- list_kMEs_PPI_ok[names(list_kMEs_PPI_ok) %in% sNames_BMI]
  list_kIMs_BMI <- list_kIMs_PPI_ok[names(list_kIMs_PPI_ok) %in% sNames_BMI]
  # nested level
  list_kMEs_BMI <- mapply(function(x,y) x[names(x) %in% names(y)], x=list_kMEs_BMI, y=list_list_module_BMI_genes, SIMPLIFY=F)
  list_kIMs_BMI <- mapply(function(x,y) x[names(x) %in% names(y)], x=list_kIMs_BMI, y=list_list_module_BMI_genes, SIMPLIFY=F)
  
  
  ######################################################################
  ############# DO GENE SET ENRICHMENT ANALYSIS (GSEA) #################
  ######################################################################
  
  # read in gene set information
  
  if (organism == "mmusculus") {
    gene.sets <- read.gmt(mm_geneset_path)
  } else if (organism == "hsapiens") {
    data("org.Hs.GO2Symbol.list")
    gene.sets <- org.Hs.GO2Symbol.list
  }
  
  ### 180522 can't install dependencies for clusterProfiler  
  
  list_list_kMEs_BMI <- lapply(list_kMEs_BMI, function(x) as.list(x))
  list_list_kMEs_BMI <- mapply(function(a,b) lapply(a, function(x) name_for_vec(to_be_named = x, given_names = rownames(b), dimension = NULL)),
                               a = list_list_kMEs_BMI,
                               b = list_kMEs_BMI,
                               SIMPLIFY = F)
  
  # list_list_gsea <- list()
  # for (i in 1:length(list_list_kMEs_BMI)) {
  #   list_list_gsea
  #   for (j in 1:length(list_list_kMEs_BMI[[i]])) {
  #     
  #     list_list_gsea <- gsea(values = list_list_kMEs_BMI[[i]][[j]], mc.cores=1, plot=T)
  #   }
  # }
  
  gene.sets.select <- gene.sets[grepl("^GO_|REACTOME|KEGG",  names(gene.sets))]
  
  invisible(gc())
  cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "log_GSEA.txt"))
  list_list_gsea <- parLapply(cl, list_list_kMEs_BMI, function(x) lapply(x, function(y) iterative.bulk.gsea(values=y, 
                                                                                                            set.list=gene.sets.select)))
  stopCluster(cl)
  invisible(gc())
  
  names(list_list_gsea) <- names(list_list_kMEs_BMI)
  
  signif.threshold <- 5e-2 / length(gene.sets.select)
  
  list_list_gsea_signif <- lapply(list_list_gsea, function(x) lapply(x, function(y) y[y$p.val<signif.threshold,])) 
  
  # TODO: Output in a single table?
  save(list_list_gsea_signif, sprintf("%s%s_list_list_gsea_signif_%s.RData",dir_RObjects, data_prefix, flag_date))

  ##########################################################################
  ##################### COMPUTE MODULE-MODULE CORRELATIONS #################
  ##########################################################################
  
  if (length(list_modules_BMI) > 1) {
    
    corr_eigengenes <- WGCNA::cor(x=eigen_mat, method=c("pearson"), verbose=verbose)
    
    # plot
    pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "BMI-enriched module eigengenes", flag_date),width=ncol(eigen_mat) %/% 3,height=ncol(eigen_mat) %/% 3)
    corrplot(corr = corr_eigengenes, 
             method = "color",
             diag = F,
             is.corr=T,
             title=sprintf("%s %s corr across subsets", data_prefix, "BMI-enriched module eigengenes"),#,
             order = "hclust",
             hclust.method = hclustMethod)
    
    invisible(dev.off())
    
    if (FALSE) {
      
      ##########################################################################
      #####################  MERGE CORRELATED MODULES ##########################
      ##########################################################################
      
      # TODO: first verify that merging makes biological sense 
      
      # Identify highly correlated modules
      corr_eigengenes_lt <- corr_eigengenes
      corr_eigengenes_lt[!lower.tri(corr_eigengenes, diag = F)] <- NA
      idx_toMerge <- which(corr_eigengenes_lt>(1-moduleMergeCutHeight), arr.ind = T)
      
      
      # Merge highly correlated modules
      
      modules_merged <- dat_BMI_fdr$module
      list_modules_BMI_merged <- list_modules_BMI # this is a nested version of the above list
      list_list_module_BMI_genes_merged  <- list_list_module_BMI_genes
      list_colors_BMI_merged <- list_colors_BMI
      
      if (nrow(idx_toMerge) > 0) {
        
        for (i in 1:nrow(idx_toMerge)) {
          
          message(paste0("Merging ", dat_BMI_fdr$module[idx_toMerge[i,1]], " with ", dat_BMI_fdr$module[idx_toMerge[i,2]], " - eigengene correlation is ", round(corr_eigengenes[idx_toMerge[i,1],idx_toMerge[i,2]],2)))
          subset_from <- gsub("_.*", "", dat_BMI_fdr$module[idx_toMerge[i,1]])
          module_from <- gsub(".*_", "", dat_BMI_fdr$module[idx_toMerge[i,1]])
          subset_to <- gsub("_.*", "", dat_BMI_fdr$module[idx_toMerge[i,2]])
          module_to <- gsub(".*_", "", dat_BMI_fdr$module[idx_toMerge[i,2]])
          
          list_list_module_BMI_genes_merged[[subset_to]] [[module_to]] <- union(list_list_module_BMI_genes[[subset_to]] [[module_to]] , list_list_module_BMI_genes[[subset_from]][[module_from]])
          
          list_list_module_BMI_genes_merged[[subset_from]] [[module_from]] <- NULL
          
          
          list_colors_BMI_merged[[subset_from]] [ list_colors_BMI_merged[[subset_from]] == module_from] <- module_to 
          
          modules_merged[idx_toMerge[i,2]] <- paste0(dat_BMI_fdr$module[idx_toMerge[i,2]], "_", dat_BMI_fdr$module[idx_toMerge[i,1]])
          modules_merged[idx_toMerge[i,1]] <- NA
          
        }
        
        modules_merged <- modules_merged[!is.na(modules_merged)]
        
        # get a new vector of names of cell clusters with BMI enrichment after merging
        sNames_BMI <- unique(gsub("_.*", "", modules_merged))
        
        # Transform the vector of modules into a nested list, as before
        list_modules_BMI <- vector(mode="list", length=length(sNames_BMI_merged))
        names(list_modules_BMI) <- sNames_BMI
        
        for (name in sNames_BMI_merged) {
          idx <- grep(pattern = paste0(name, "_.*"), x = modules_merged)
          list_modules_BMI[[name]] <- sapply(idx, function(j) gsub(".*_", "", modules_merged[j]), simplify = T)
        }
      }
    }
  }
  
  
  ##########################################################################
  ############ SAVE GENE LISTS AND KMEs FOR BMI-ENRICHED MODULES ###########
  ##########################################################################
  
  # output module gene df 
  invisible(write.csv(df_BMI_module_genes, file = sprintf("%s%s_BMI_modules_genes_%s.csv",tables_dir, data_prefix, flag_date), row.names = F))
  
  # Prepare dfs with a gene column followed by kMEs / kIMs
  list_kMEs_BMI_out <- lapply(list_kMEs_BMI, function(x) cbind(genes=rownames(x), x))
  list_kIMs_BMI_out <- lapply(list_kIMs_BMI, function(x) cbind(genes=rownames(x), x))
  rownames(list_kMEs_BMI_out) <- NULL
  rownames(list_kIMs_BMI_out) <- NULL
  
  # output kMEs and kIMs for BMI enriched modules
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kMEs_BMI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_kMEs_BMI_out, sNames_BMI, SIMPLIFY = F))
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kIMs_BMI_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F), list_kIMs_BMI_out, sNames_BMI, SIMPLIFY = F))
  
  ### 180504_v1.8_dev1
  resume = "checkpoint_4"
  message("Reached checkpoint 4, saving session image")
  save.image(file=sprintf("%s%s_checkpoint_4_image.RData", RObjects_dir, data_prefix))
  
} else if (resume == "checkpoint_4") {
  load(file=sprintf("%s%s_checkpoint_4_image.RData", RObjects_dir, data_prefix))
  # Load parameter values and utility functions anew (in case a bug was fixed)
  source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params_v1.8_dev3.R")
  source(file = "/projects/jonatan/functions-src/functions_v1.8_dev3.R")
}  
###

######################################################################
############### PLOT EIGENGENE - METADATA CORRELATION ################
######################################################################

if (!is.null(meta.data)) {
  # compute metadata correlations with eigengene embeddings
  # arrange meta.data rows to match eigenvector embeddings
  meta.data <- meta.data[match(rownames(meta.data), rownames(eigen_mat)),]
  eigen_meta.data_corr <- WGCNA::cor(x=meta.data, y=eigen_mat, method = c("pearson"), verbose = 0)
  # plot
  pdf(sprintf("%s%s_corr_eigengenes_metadata_%s.pdf", plots_dir, data_prefix, flag_date))
  corrplot(corr = eigen_meta.data_corr,
           add = F,
           method = "color",
           is.corr=T,
           #order = "hclust",
           #hclust.method = hclustMethod,
           title=sprintf("%s eigenvector - meta.data correlations", data_prefix),
           #addCoef.col = F,
           number.digits = NULL)
  invisible(dev.off())
}

######################################################################
############# PLOT kME - kIM CORRELATIONS ############################
######################################################################
# This is a diagnostic plot to see whether kIMs are highly correlated with kMEs as expected

### 180516 v1.8_dev3 
# Not really so useful divided into celltypes. better to do across celltypes

# list_corr_kMEs_kIMs <- mapply(function(x,y) cor(x,y), x=list_kMEs, y=list_kIMs, SIMPLIFY=F)
# invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
#                                                 vectorNames = "kME-kIM", 
#                                                 subsetName=y, 
#                                                 diag=T,
#                                                 is.corr = T), 
#                  x=list_corr_kMEs_kIMs, y=sNames_ok, SIMPLIFY=F))
###

# list_corr_kMEs_kIMs_BMI <- mapply(function(x,y) cor(x,y), x=list_kMEs_BMI, y=list_kIMs_BMI, SIMPLIFY=F)
# invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
#                                                 vectorNames = "module kME-kIM (BMI enriched)", 
#                                                 subsetName=y, 
#                                                 diag=T,
#                                                 is.corr = T), 
#                  x=list_corr_kMEs_kIMs_PPI, y=sNames_BMI, SIMPLIFY=F))
# #}

######################################################################
######## PLOT kME-kME & kIM-kIM CORRELATIONS WITHIN CELLTYPES ########
######################################################################
  
# list_corr_kMEs_kMEs <- lapply(list_kMEs, function(x) WGCNA::cor(x,x, method="spearman", verbose=verbose))
# list_corr_kIMs_kIMs <- lapply(list_kIMs, function(x) WGCNA::cor(x,x, method="spearman", verbose=verbose))
# 
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
# invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
#                                                 vectorNames = "kME-kME", 
#                                                 subsetName=y, 
#                                                 diag=F,
#                                                 is.corr=T,
#                                                 order = "hclust",
#                                                 hclust.method = hclustMethod), 
#                  x=list_corr_kMEs_kMEs, y=sNames_ok, SIMPLIFY=F))
# 
# invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
#                                                 vectorNames = "kIM-kIM", 
#                                                 subsetName=y, 
#                                                 diag=F,
#                                                 is.corr=T,
#                                                 order = "hclust",
#                                                 hclust.method = hclustMethod), 
#                  x=list_corr_kIMs_kIMs, y=sNames_ok, SIMPLIFY=F))

### 180607_v1.8_dev2
#if (!is.null(STRINGdb_species)){
# list_corr_kMEs_kMEs_PPI <- lapply(list_kMEs_PPI, function(x) WGCNA::cor(x,x, method="spearman", verbose=verbose))
# list_corr_kIMs_kIMs_PPI <- lapply(list_kIMs_PPI, function(x) WGCNA::cor(x,x, method="spearman", verbose=verbose))
# 
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
# invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
#                                                 vectorNames = "kME-kME_PPI", 
#                                                 subsetName=y, 
#                                                 diag=F,
#                                                 is.corr=T,
#                                                 order = "hclust",
#                                                 hclust.method = hclustMethod), 
#                  x=list_corr_kMEs_kMEs_PPI, y=sNames_ok, SIMPLIFY=F))

# invisible(mapply(function(x,y) plotCorr_for_par(corrMatrix=x, 
#                                                 vectorNames = "kIM-kIM_PPI", 
#                                                 subsetName=y, 
#                                                 diag=F,
#                                                 is.corr=T,
#                                                 order = "hclust",
#                                                 hclust.method = hclustMethod), 
#                  x=list_corr_kIMs_kIMs_PPI, y=sNames_ok, SIMPLIFY=F))
#}

######################################################################
######################### LOAD FULL DATA #############################
######################################################################
  
#datExpr_full <- t(load_obj(seurat_path)@data)
#seurat_obj <- load_obj(f=sprintf("%sseurat_obj_ensembl.RData",RObjects_dir))

#datExpr_full <- t(seurat_obj@data)

######################################################################
####### PLOT kME-kME and kIM-kIM correlation across all subsets ######
######################################################################

# Include all genes in kMEs / kIMs with value 0 if they were not included in the subset analysis
### 180503_v1.7
# need to use original names
# list_kMEs_all_genes <- lapply(list_kMEs, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes))
# list_kIMs_all_genes <- lapply(list_kIMs, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes))
# ###

# Bind all celltype module kMEs / kIMs in single dataframes

# for (i in 1:length(list_kMEs_all_genes)) {
#   if (i == 1) {
#     kMEs_across_subsets <- list_kMEs_all_genes[[i]] 
#   } else { kMEs_across_subsets <- cbind(kMEs_across_subsets,list_kMEs_all_genes[[i]])
#   }
# }

# for (i in 1:length(list_kIMs_all_genes)) {
#   if (i == 1) {
#     kIMs_across_subsets <- list_kIMs_all_genes[[i]] 
#   } else { kIMs_across_subsets <- cbind(kIMs_across_subsets,list_kIMs_all_genes[[i]])
#   }
# }
#kMEs_across_subsets <- cbind(eval(parse(text=sprintf("c(%s)", paste0("list_kMEs_all_genes[[", 1:length(list_kMEs_all_genes), "]]")))))
#kIMs_across_subsets <- cbind(eval(parse(text='c(paste0("list_kIMs_all_genes[[", 1:length(list_kIMs_all_genes), "]]"))')))

# name the dataframes and matrices so as to distinguish modules from different subsets   
#prefixes = unlist(mapply(function(x,y) rep(x = x, times= ncol(y)), x=sNames_ok, y= list_kMEs_all_genes, SIMPLIFY=T)) # same for kME and kIM

#colnames(kMEs_across_subsets) = paste(prefixes, colnames(kMEs_across_subsets), sep="_")
#colnames(kIMs_across_subsets) = paste(prefixes, colnames(kIMs_across_subsets), sep="_")

# correlation matrices
#corr_kMEs_across_subsets <- WGCNA::cor(as.matrix(kMEs_across_subsets), method="spearman", verbose = verbose)
#corr_kIMs_across_subsets <- WGCNA::cor(as.matrix(kIMs_across_subsets), method="spearman", verbose = verbose)

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
#  pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "kME-kME", flag_date))#,width=ncol(kMEs_across_subsets) %/% 4,height=ncol(kMEs_across_subsets) %/% 4)
#  corrplot(corr = corr_kMEs_across_subsets, 
#           method = "color",
#           is.corr=T,
#           title=sprintf("%s %s corr across subsets", data_prefix, "kME-kME"),
#           #addCoef.col = F,
#           number.digits = NULL)
#  
#  invisible(dev.off())
#  
# # plot
#  pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "kIM-kIM", flag_date),width=ncol(kIMs_across_subsets) %/% 4,height=ncol(kIMs_across_subsets) %/% 4)
#  corrplot(corr = corr_kIMs_across_subsets, 
#           add = F,
#           method = "color",
#           is.corr=T,
#           #order = "hclust",
#           #hclust.method = hclustMethod,
#           title=sprintf("%s %s corr across subsets", data_prefix, "kIM-kIM"),
#           #addCoef.col = F,
#           number.digits = NULL)
#  
#  invisible(dev.off())

### 180607_v1.8_dev2
#if (!is.null(STRINGdb_species)){
#### 180503_v1.7
# Include all genes in kMEs / kIMs with value 0 if they were not included in the subset analysis
# list_kMEs_PPI_all_genes <- lapply(list_kMEs_PPI, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes_IDs_mix))
# list_kIMs_PPI_all_genes <- lapply(list_kIMs_PPI, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes_IDs_mix))
# list_kMEs_BMI_all_genes <- lapply(list_kMEs_BMI, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes))
# list_kIMs_BMI_all_genes <- lapply(list_kIMs_BMI, function(x) get_kMEs_all_genes(kMEs=x, all_genes=all_genes))

###

# Bind all celltype module kMEs / kIMs in single dataframes

# for (i in 1:length(list_kMEs_BMI_all_genes)) {
#   if (i == 1) {
#     kMEs_BMI_across_subsets <- list_kMEs_BMI_all_genes[[i]] 
#   } else { kMEs_BMI_across_subsets <- cbind(kMEs_BMI_across_subsets,list_kMEs_BMI_all_genes[[i]])
#   }
# }
  
# for (i in 1:length(list_kIMs_BMI_all_genes)) {
#   if (i == 1) {
#     kIMs_BMI_across_subsets <- list_kIMs_BMI_all_genes[[i]] 
#   } else { kIMs_BMI_across_subsets <- cbind(kIMs_BMI_across_subsets,list_kIMs_BMI_all_genes[[i]])
#   }
# }
  
# Bind all celltype module kMEs / kIMs in single dataframes
# kMEs_PPI_across_subsets <- cbind(eval(parse(text=c(paste0("list_kMEs_PPI_all_subsets[[", 1:length(list_kMEs_PPI_all_genes), "]]")))))
# kMEs_PPI_across_subsets <- cbind(eval(parse(text=c(paste0("list_kMEs_PPI_all_subsets[[", 1:length(list_kMEs_PPI_all_genes), "]]")))))
# 
# name the dataframes and matrices so as to distinguish modules from different subsets   
# prefixes = unlist(mapply(function(x,y) rep(x = x, times= ncol(y)), x=sNames_BMI, y= list_kMEs_BMI_all_genes, SIMPLIFY=T)) # same for kME and kIM
# 
# colnames(kMEs_BMI_across_subsets) = paste(prefixes, colnames(kMEs_BMI_across_subsets), sep="_")
# colnames(kMEs_BMI_across_subsets) = paste(prefixes, colnames(kMEs_BMI_across_subsets), sep="_")
# 
# correlation matrices
#corr_kMEs_BMI_across_subsets <- WGCNA::cor(as.matrix(kMEs_BMI_across_subsets), method="spearman", verbose=verbose)
#corr_kMEs_BMI_across_subsets <- WGCNA::cor(as.matrix(kMEs_BMI_across_subsets), method="spearman", verbose=verbose)


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
#   

### 180617 Moved up so we can merge correlated modules
# corr_eigengenes_BMI_across_subsets <- WGCNA::cor(x=eigen_mat, method=c("pearson"), verbose=verbose)
# 
# # plot
# pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "BMI-enriched module eigengenes", flag_date),width=ncol(eigen_mat) %/% 3,height=ncol(eigen_mat) %/% 3)
# corrplot(corr = corr_eigengenes_BMI_across_subsets, 
#          method = "color",
#          diag = F,
#          is.corr=T,
#          title=sprintf("%s %s corr across subsets", data_prefix, "BMI-enriched module eigengenes"),#,
#          order = "hclust",
#          hclust.method = hclustMethod)
# 
# invisible(dev.off())
###

# plot
# pdf(sprintf("%s%s_%s_corr_across_subsets_%s.pdf", plots_dir, data_prefix, "kIM-kIM_BMI", flag_date),width=ncol(kMEs_BMI_across_subsets) %/% 4,height=ncol(kMEs_BMI_across_subsets) %/% 4)
# corrplot(corr = corr_kMEs_BMI_across_subsets, 
#          method = "color",
#          is.corr=T,
#          title=sprintf("%s %s corr across subsets", data_prefix, "module kIM-kIM (BMI enriched)"))#,
#          #order = "hclust",
#          #hclust.method = hclustMethod)
# 
# invisible(dev.off())
#}

######################################################################
########### PLOT KME - & KIM - METADATA CORRELATIONS #################
######################################################################

# if (!is.null(meta.data)) {
#   
#   # Compute cell embedding
#   cell_embeddings_kME <- datExpr_full %*% as.matrix(kMEs_across_subsets) # cell x gene %*% gene x module -> cell x module 
#   cell_embeddings_kIM <- datExpr_full %*% as.matrix(kIMs_across_subsets) # cell x gene %*% gene x module -> cell x module 
#   
#   # correlation matrices 
#   corr_kME_meta.data <- WGCNA::cor(x=cell_embeddings_kME, y=meta.data, method="spearman", verbose=verbose)
#   corr_kIM_meta.data <- WGCNA::cor(x=cell_embeddings_kME, y=meta.data, method="spearman", verbose=verbose)
#   
#   # name the matrices
#   # rownames(corr_kME_meta.data) = meta.data_corr_col_numeric
#   # rownames(corr_kIM_meta.data) = meta.data_corr_col_numeric
#   # 
#   # colnames(corr_kME_meta.data) = paste(prefixes, colnames(kMEs_PPI_across_subsets), sep="_")
#   # colnames(corr_kIM_meta.data) = paste(prefixes, colnames(kIMs_PPI_across_subsets), sep="_")
#   
#   # plot
#   pdf(sprintf("%s%s_corr_%s_with_metadata_%s.pdf", plots_dir, data_prefix, "kME_embeddings", flag_date),width=ncol(meta.data) %/% 4,height=ncol(kMEs_across_subsets) %/% 4)
#   corrplot(corr = corr_kME_meta.data, 
#            method = "color",
#            is.corr=T,
#            title=sprintf("%s corr: module %s with metadata", data_prefix, "kME embeddings"))#,
#            #order = "hclust",
#            #hclust.method = hclustMethod)
#   
#   dev.off()
#   
#   pdf(sprintf("%s%s_corr_%s_with_metadata_%s.pdf", plots_dir, data_prefix, "kIM_embeddings", flag_date),width=ncol(meta.data) %/% 4,height=ncol(kIMs_across_subsets) %/% 4)
#   corrplot(corr = corr_kIM_meta.data, 
#            method = "color",
#            is.corr=T,
#            title=sprintf("%s corr: module %s with metadata", data_prefix, "kIM embeddings"))#,
#            #order = "hclust",
#            #hclust.method = hclustMethod)
#   
#   dev.off()
#   
#   ### 180607_v1.8_dev2
#   #if (!is.null(STRINGdb_species)) {
# 
#   # Compute cell embedding
#   cell_embeddings_kMEs_PPI <- datExpr_full %*% as.matrix(kMEs_PPIs_across_subsets) # cell x gene %*% gene x module -> cell x module 
#   cell_embeddings_kIMs_PPI <- datExpr_full %*% as.matrix(kIMs_PPIs_across_subsets) # cell x gene %*% gene x module -> cell x module 
#   
#   # correlation matrices 
#   corr_kMEs_PPI_meta.data <- WGCNA::cor(x=cell_embeddings_kMEs_PPI, y=meta.data, method="spearman", verbose=verbose)
#   corr_kIMs_PPI_meta.data <- WGCNA::cor(x=cell_embeddings_kMEs_PPI, y=meta.data, method="spearman", verbose=verbose)
#   
#   # name the matrices
#   # rownames(corr_kMEs_PPI_meta.data) = meta.data_corr_col_numeric
#   # rownames(corr_kIMs_PPI_meta.data) = meta.data_corr_col_numeric
#   # 
#   # colnames(corr_kMEs_PPI_meta.data) = paste(prefixes, colnames(kMEs_PPIs_PPI_across_subsets), sep="_")
#   # colnames(corr_kIMs_PPI_meta.data) = paste(prefixes, colnames(kIMs_PPI_across_subsets), sep="_")
#   # 
#   # plot
#   pdf(sprintf("%s%s_corr_%s_with_metadata_%s.pdf", plots_dir, data_prefix, "kMEs_PPI_embeddings", flag_date),width = max(4, ncol(meta.data) %/% 4), height = ncol(kMEs_PPIs_across_subsets) %/% 4)
#   corrplot(corr = corr_kMEs_PPI_meta.data, 
#            method = "color",
#            is.corr=T,
#            title=sprintf("%s corr: module %s with metadata", data_prefix, "kMEs_PPI embeddings"))#,
#            #order = "hclust",
#            #hclust.method = hclustMethod)
#   
#   dev.off()
#   
#   pdf(sprintf("%s%s_corr_%s_with_metadata_%s.pdf", plots_dir, data_prefix, "kIMs_PPI_embeddings", flag_date),width = max(4, ncol(meta.data) %/% 4), height = ncol(kIMs_PPIs_across_subsets) %/% 4)
#   corrplot(corr = corr_kIMs_PPI_meta.data, 
#            method = "color",
#            is.corr=T,
#            title=sprintf("%s corr: module %s with metadata", data_prefix, "kIMs_PPI embeddings"))#,
#            #order = "hclust",
#            #hclust.method = hclustMethod)
#   
#   dev.off()
  #}
#} # end of if (!is.null(meta.data))

######################################################################
##### PLOT CELL EXPRESSION ON MODULES WITH HIGH GWAS SIGNIFICANCE ####
######################################################################
# 
# gwas_sub_dirs = list.dirs(path = magma_gwas_dir, full.names = FALSE, recursive = FALSE)
# 
# file_suffix = "kMEs"
# # upper level list of subsets, nested list with tables for each GWAS repo
# list_list_magma_tables <- lapply(sNames_ok, function(x) lapply(gwas_sub_dirs, function(y) load_obj(f=paste0(tables_dir, data_prefix,"_", x, "_", file_suffix, "_magma_GWAS_", y, "_", flag_date, ".csv"))))
# list_list_idx_signif_modules <- lapply(list_list_magma_tables, function(x) lapply(x, function(y) x[x[,max.col(x[,-which(colnames(x)=="module")])] > -log10(5^-2),][["module"]]))
# list_MAGMA_signif_modules <- sapply(list_list_idx_signif_modules, function(x) unique(unlist(x)))
# 
# # Compute cell embeddings
# list_cell_embeddings_kME_MAGMA <- lapply(list_MAGMA_signif_modules, function(x) as.matrix(kMEs_across_subsets[,colnames(kMEs_across_subsets) %in% x]) %*% t(datExpr_full)) # module x gene %*% gene x cell -> module x cell
# 
# # name the matrix
# #colnames
# # list_cell_embeddings_kME_MAGMA <- lapply(list_cell_embeddings_kME_MAGMA, function(x) name_for_vec(to_be_named = x, given_names = rownames(datExpr_full), dimension = 2)) # cell names
# # #rownames
# # list_cell_embeddings_kME_MAGMA <- lapply(list_cell_embeddings_kME_MAGMA, function(x) name_for_vec(to_be_named = x, given_names = paste(prefixes, colnames(x), sep="_"), dimension = 1))
# 
# # group the columns of the embeddings matrix by cell type (@ident)
# celltypes_idx <- integer(0)
# for (celltype in names(table(seurat_obj@ident))) {
#   celltypes_idx <- c(celltypes_idx, which(seurat_obj@ident == celltype)) 
# }
# 
# list_cell_embeddings_kME_MAGMA <- lapply(list_cell_embeddings_kME_MAGMA, function(x) x[,celltypes_idx]) 
# 
# # plot
# for (i in 1:length(gwas_sub_dirs)) {
#   pheatmap(mat = list_cell_embeddings_kME_MAGMA[[i]], 
#            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))[0:100],
#            kmeans_k = NA,
#            breaks = NA, 
#            cluster_rows = F,
#            cluster_cols = F,
#            #clustering_method = hclustMethod,
#            show_rownames = T,
#            show_colnames = F,
#            main = sprintf("%s_%s_%s_%s", data_prefix, "module_embeddings_MAGMA", gwas_sub_dirs[i], flag_date),
#            filename=sprintf("%s%s_%s_%s_%s.pdf", plots_dir, data_prefix, "module_embeddings_MAGMA", gwas_sub_dirs[i], flag_date),
#            width = 10,
#            height = length(list_MAGMA_signif_modules) %/% 4,
#            silent = T)
# }
# 
# ### 180607_v1.8_dev2
# #if (!is.null(STRINGdb_species)) {
#   
# file_suffix = "kMEs_PPI"
# # upper level list of subsets, nested list with tables for each GWAS repo
# list_list_PPI_magma_tables <- lapply(sNames_ok, function(x) lapply(gwas_sub_dirs, function(y) load_obj(file=paste0(tables_dir, data_prefix,"_", x, "_", file_suffix, "_magma_GWAS_", y, "_", flag_date, ".csv"))))
# list_list_idx_PPI_signif_modules <- lapply(list_list_PPI_magma_tables, function(x) lapply(x, function(y) x[x[,max.col(x[,-which(colnames(x)=="module")])] > -log10(5^-2),][["module"]]))
# list_PPI_MAGMA_signif_modules <- sapply(list_list_idx_PPI_signif_modules, function(x) unique(unlist(x)))
# 
# # Compute cell embeddings
# list_cell_embeddings_kME_PPI_MAGMA <- lapply(list_PPI_MAGMA_signif_modules, function(x) as.matrix(kMEs_PPI_across_subsets[,colnames(kMEs_PPI_across_subsets) %in% x]) %*% t(datExpr_full)) # module x gene %*% gene x cell -> module x cell
# 
# # name the matrix
# #colnames
# # list_cell_embeddings_kME_PPI_MAGMA <- lapply(list_cell_embeddings_kME_PPI_MAGMA, function(x) name_for_vec(to_be_named = x, given_names = rownames(datExpr_full), dimension = 2)) #cell names
# # #rownames
# # list_cell_embeddings_kME_PPI_MAGMA <- lapply(list_cell_embeddings_kME_PPI_MAGMA, function(x) name_for_vec(to_be_named = x, given_names = paste(prefixes, colnames(x), sep="_"), dimension = 1))
#  
# # group the columns of the embeddings matrix by cell type (@ident)
# 
# list_cell_embeddings_kME_PPI_MAGMA <- lapply(list_cell_embeddings_kME_PPI_MAGMA, function(x) x[,celltypes_idx]) 
# 
# for (i in 1:length(gwas_sub_dirs)) {
#   # plot
#   pheatmap(mat = list_cell_embeddings_kME_PPI_MAGMA[[i]], 
#            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))[0:100],
#            kmeans_k = NA,
#            breaks = NA, 
#            cluster_rows = F,
#            cluster_cols = F,
#            #clustering_method = hclustMethod,
#            show_rownames = T,
#            show_colnames = F,
#            main = sprintf("%s_%s_%s_%s", data_prefix, "module_embeddings_PPI_MAGMA", gwas_sub_dirs[i], flag_date),
#            filename=sprintf("%s%s_%s_%s_%s.pdf", plots_dir, data_prefix, "module_embeddings_PPI_MAGMA", gwas_sub_dirs[i], flag_date),
#            width = 10,
#            height = length(list_MAGMA_signif_modules) %/% 4,
#            silent = T)
# }
# #}
###
######################################################################
######## PLOT MODULE PRESERVATION IN DIFFERENT TISSUES ###############
######################################################################

if (FALSE) {
  
  suppressPackageStartupMessages(library(pheatmap))
  
  # Check preservation of each set of colors in all celltypes (including itself)
  module_preservation_celltypes <- sapply(list_colors, function(x) wrapModulePreservation(listDatExpr = list_datExpr_ok,
                                                                                                  listColors = list(x))[['preservation']][['log.pBonf']],
                                                            simplify = T)
  
  colnames(module_preservation_celltypes) = sNames_ok
  rownames(module_preservation_celltypes) = paste(prefixes, names(unlist(list_colors, recursive = F, use.names = T)), sep="_")
  
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
           width = length(list_datExpr_ok) %/% 2,
           height = length(unlist(list_colors, recursive = F)) %/% 4,
           silent = T)
  
  ### 180607_v1.8_dev2
  #if (!is.null(STRINGdb_species)) {
  module_PPI_preservation_celltypes <- sapply(list_colors_PPI, function(x) wrapModulePreservation(listDatExpr = list_datExpr_ok,
                                                                                                listColors = list(x))[['preservation']][['log.pBonf']],
                                          simplify = T)
  
  colnames(module_PPI_preservation_celltypes) = sNames_ok
  rownames(module_PPI_preservation_celltypes) = paste(prefixes, names(unlist(list_colors_PPI, recursive = F, use.names = T)), sep="_")
  
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
           width = length(list_datExpr_ok) %/% 2,
           height = length(unlist(list_colors_PPI, recursive = F)) %/% 4,
           silent = T)
  #}
  ###
} # end of if (FALSE) {make all the plots} TODO: remove

######################################################################
########################### CLEAR UP #################################
######################################################################
### 180503_v1.8_dev1
save.image(file=sprintf("%s%s_final_session_image.RData", RObjects_dir, data_prefix, flag_date))

rm(list=ls())
invisible(gc())

##########################################################################
################################ FINISH ##################################
##########################################################################

message("Script DONE!")
