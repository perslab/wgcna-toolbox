# Title: Script to find robust WGCNA modules

######################################################################
############################## USAGE #################################
#######################o###############################################

# e.g.
# export R_MAX_NUM_DLLS=999
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R  --seurat_path /raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-wgcna_modules/maca.seurat_obj.facs.figshare_180126.top29_BMI.RData  --project_dir /projects/jonatan/tmp-rwgcna-tests/maca_top29_BMI_2/  --data_prefix maca_top29_BMI_2  --data_organism mmusculus  --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --min.cells 10  --genes_use PCA  --pca_genes var.genes  --corFnc cor  --networkType signed  --hclustMethod average  --minClusterSize "c(15)"  --deepSplit "c(3)"  --pamStage "c(TRUE)"  --moduleMergeCutHeight "c(0.1)"  --jackstrawnReplicate 500  --TOMnReplicate 0  --n_cores 13 --autosave T --fuzzyModMembership kME --checkPPI T --resume checkpoint_2

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R  --seurat_path /raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-wgcna_modules/maca.seurat_obj.facs.figshare_180126.top29_BMI.RData  --project_dir /projects/jonatan/tmp-rwgcna-tests/maca_top29_BMI_2/  --data_prefix maca_top29_BMI_2  --data_organism mmusculus  --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --min.cells 10  --genes_use PCA  --pca_genes var.genes  --corFnc cor  --networkType signed  --hclustMethod average  --minClusterSize "c(15)"  --deepSplit "c(3)"  --pamStage "c(TRUE)"  --moduleMergeCutHeight "c(0.1)"  --jackstrawnReplicate 500  --TOMnReplicate 0  --n_cores 13 --autosave T --fuzzyModMembership kME --checkPPI T --resume checkpoint_2
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R  --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_all_author_ids.RData   --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell_3/  --data_prefix campbell_3 --data_organism mmusculus  --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --min.cells 10  --genes_use PCA  --pca_genes var.genes  --corFnc cor  --networkType signed  --hclustMethod average  --minClusterSize "c(15)"  --deepSplit "c(2)"  --pamStage "c(TRUE)"  --moduleMergeCutHeight "c(0.1)"  --jackstrawnReplicate 500  --TOMnReplicate 0  --n_cores 7  --autosave T --fuzzyModMembership kME --checkPPI T  --resume checkpoint_2
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R  --seurat_path /projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-L5.top10_cts.seurat_obj.RData  --project_dir /projects/jonatan/tmp-rwgcna-tests/mousebrain_top10_4/ --metadata_subset_col ClusterName --data_prefix mousebrain_top10_5 --data_organism mmusculus  --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --min.cells 10 --genes_use PCA  --pca_genes var.genes  --corFnc cor  --networkType signed  --hclustMethod average  --minClusterSize "c(15)"  --deepSplit "c(3)"  --pamStage "c(TRUE)"  --moduleMergeCutHeight "c(0.1)"  --jackstrawnReplicate 500  --TOMnReplicate 0  --n_cores 5  --autosave T --fuzzyModMembership kME --checkPPI T
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main_v2.3_dev.R  --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter_ensembl_UPF.RData  --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-UPF-4/ --metadata_subset_col "c('res.0.6')" --data_prefix EP-UPF-4 --data_organism hsapiens  --magma_gwas_dir /projects/jonatan/tmp-epilepsy/data/magma/ilae-lancet-2014/ --min.cells 10 --genes_use PCA  --pca_genes all  --corFnc cor  --networkType signed  --hclustMethod average  --minClusterSize "c(15)"  --deepSplit "c(2)"  --pamStage "c(TRUE)"  --moduleMergeCutHeight "c(0.1)"  --jackstrawnReplicate 500  --TOMnReplicate 0  --n_cores 10  --autosave T --fuzzyModMembership kIM --checkPPI F

######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

TODO = T # controls flows around broken / unfinished code

if (FALSE) { 
  seurat_path = "/projects/jonatan/tmp-holst-hsl/RObjects/campbell_n05_to_n10.RData"
  project_dir = "/projects/jonatan/tmp-rwgcna-tests/tmp-campbell-n05_to_n10-10/"
  data_prefix = "campbell-n10_to_n10-10"
  resume = NULL 
  autosave = T
  metadata_subset_col = NULL
  metadata_corr_col = NULL #c("nUMI", "nGene", "percent.mito", "percent.ribo")
  metadata_corr_filter_vals = "UPF"
  use.imputed = F
  regress_out = c("nUMI", "percent.mito", "percent.ribo")
  min.cells = 5
  genes_use = "PCA"
  pca_genes = 'all'
  corFnc = "cor"
  networkType = "signed"
  minClusterSize = c(20)
  deepSplit = c(2)
  pamStage = c(TRUE)
  moduleMergeCutHeight = c(0.2)
  jackstrawnReplicate = 0
  TOMnReplicate = 0
  data_organism = "mmusculus"
  magma_gwas_dir = "/projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/"
  gwas_filter_traits <- "c('BMI', 't1d', 't2d')"
  fuzzyModMembership = "kIM"
  scale_MEs_by_kIMs = F
  n_cores = 6
}

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  
  make_option("--seurat_path", type="character",
              help = "Provide full path to Rdata input file with Seurat object with lognormalized expression data in the @data slot "),
  make_option("--project_dir", type="character", default=NULL,
              help = "Optional. Provide project directory. Must have subdirs RObjects, plots, tables. If not provided, assumed to be dir one level up from input data dir. [default %default]"),
  make_option("--data_prefix", type="character", default=paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_rWGCNA_run"),
              help = "Dataset prefix for output files"),
  make_option("--autosave", type="logical", default=T,
              help = "Autosave session images at regular intervals to resume later? [default %default]."),
  make_option("--resume", type="character", default=NULL,
              help = "Resume from a checkpoint? Must have same path and data_prefix as provided. Options are 'checkpoint_1' - 'checkpoint_4' [default %default]"),
  make_option("--quit_session", type="character", default=NULL,
              help = "Quit after saving session image at a checkpoint? Options are 'checkpoint_1' - 'checkpoint_4'"),
  make_option("--metadata_subset_col", type="character", default=NULL,
              help = "Specify a seurat@meta.data$... column to use for subsetting the Seurat object. If NULL (default) uses the @ident slot."),
  make_option("--metadata_corr_col", type="character", default='NULL',
              help = "Specify seurat_obj@meta.data$... column(s) for which to compute correlations with gene modules. Takes a character with a vector of seurat_obj@meta.data column names e.g. 'nUMI' or 'c('nUMI', 'Age')'. For factor or character metadata, each levels is analysed as a dummy variable, so exercise caution.  [default %default]"),
  make_option("--metadata_corr_filter_vals", type="character", default='NULL',
              help = "Specify one or more values within the seurat@meta.data$... column(s). Retain only modules which are significantly (anti-) correlated (at present, there is no threshold value for the correlation). Takes a character with a vector of meta.data column names e.g. 'Female' or 'c('fasted', 'HFD')'. Case-insensitive [default %default]."),
  make_option("--use.imputed", type="logical", default=F,
              help="Use data in the obj@imputed slot for the computations to replace the @data slot? If the @imputed slot is empty, will revert to [default %default]"),
  make_option("--regress_out", type="character", default='c("nUMI", "percent.mito", "percent.ribo")',
              help="Provide arguments to Seurat's ScaleData function in the form of a vector in quotes, defaults to c('nUMI', 'percent.mito', 'percent.ribo') (not recommended, however) [default %default]"),
  make_option("--min.cells", type="integer", default=5L,
              help="What is the minimum number of cells in each subset in the data in which a gene should be detected to not be filtered out? Integer, [default %default]."),
  make_option("--genes_remove_dir", type="character", default=NULL,
              help="Directory with genesets to remove at the beginning of the analysis, e.g. mitochrondrial or ribosomal genes, in ensembl format. Takes .gz, .tab, .csv and .txt (delim) files. [default %default]"), 
  make_option("--genes_use", type="character", default="PCA",
              help="One of 'all', 'var.genes' or 'PCA' for genes with significant loading on at least one significant PC. 'All' is not recommended. [default %default]"), 
  make_option("--pca_genes", type="character", default="var.genes",
              help="'all' or 'var.genes'. 'All' is computationally expensive but allows for selecting genes based on PC loading p-values rather than magnitudes [default %default]"), 
  make_option("--corFnc", type="character", default="bicor",
              help="Use 'cor' for Pearson or 'bicor' for midweighted bicorrelation function (https://en.wikipedia.org/wiki/Biweight_midcorrelation). [default %default]"), 
  make_option("--networkType", type="character", default = "signed",
              help="'signed' scales correlations to [0:1]; 'unsigned' takes the absolute value (but the TOM can still be 'signed'); 'signed hybrid' sets negative correlations to zero. [default %default]"),
  make_option("--hclustMethod", type="character", default="average",
              help = "Hierarchical clustering agglomeration method. One of 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA), 'median' (= WPGMC) or 'centroid' (= UPGMC). See hclust() documentation for further information. [default %default]"),
  make_option("--minClusterSize", type="character", default="15L",
              help = "Minimum genes needed to form a module, or an initial cluster before the Partitioning Around Medoids-like step. WGCNA authors recommend decreasing the minimum cluster size when using higher settings of deepSplit. Takes a character with a vector of integer values to try 'c(15,20,25)' [default %default]"),
  make_option("--deepSplit", type="character", default="3L",
              help = "Controls the sensitivity of the cutreeDynamic/cutreeHybrid algorithm. Takes a character with a vector of integer values between 0-4 to try, e.g. 'c(1,2,3)', [default %default]"),
  make_option("--moduleMergeCutHeight", type="character", default="c(0.1)",
              help = "Cut-off level for 1-cor(eigen-gene, eigen-gene) for merging modules. Takes a character with a vector of double values to try, e.g. 'c(0.1, 0.2)', [default %default]"),
  make_option("--pamStage", type="character", default="c(TRUE)",
              help = "cutreeHybrid. Perform additional Partition Around Medroids step? Users favoring specificity over sensitivity (that is, tight clusters with few mis-classifications) should set pamStage = FALSE, or at least specify the maximum allowable object-cluster dissimilarity to be zero (in which case objects are assigned to clusters only if the object-cluster dissimilarity is smaller than the radius of the cluster). Takes a character with a vector of logicals to try, e.g. 'c(TRUE,FALSE)', [default %default]"),
  make_option("--kM_reassign", type="logical", default="TRUE",
              help = "Following hierarchical clustering, do additional k-means clustering? [default %default]"),
  make_option("--jackstrawnReplicate", type="integer", default=500L,
              help = "Number of times to re-run PCA after permuting a small proportion of genes to perform empirical significance tests, i.e. the `JackStraw` procedure (see `pca_genes` above), [default %default]."),
  make_option("--TOMnReplicate", type="integer", default=100L,
              help = "Number of times to resample the dataset when finding the consensus TOM [default %default]"),
  make_option("--fuzzyModMembership", type="character", default="kIM",
              help="Which 'fuzzy' measure of gene membership in a module should be used? Options are 'kME' (correlation between gene expression and module PC1 expression) and 'kIM' (sum of edges between a gene and genes in the module, normalized by number of genes in the module; i.e. average distance in the TOM between a gene and genes in the module [default %default]."),  
  make_option("--scale_MEs_by_kIMs", type="logical", default="F",
              help="If fuzzyModMembership == kME, should the gene weight in eigengene computation be scaled by kIM? [default %default]"),  
  make_option("--checkPPI", type="logical", default=T,
              help="Valiate gene modules using Protein-Protein Interactions?"),
  make_option("--data_organism", type="character", default="mmusculus",
              help = "'hsapiens' or 'mmusculus', [default %default]"),
  make_option("--magma_gwas_dir", type="character", default = NULL,
              help = "MAGMA input GWAS data directory, e.g. '/projects/jonatan/tmp-epilepsy/data/magma/ilae-lancet-2014/', '/projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/'. NULL skips the magma GWAS step [default %default]"),
  make_option("--gwas_filter_traits", type="character", default = NULL,
              help = "Filter out modules not significantly correlated with matching gwas studies within the magma_gwas_dir. Takes a character with a vector of character names to match within the filename of the GWAS , e.g. 'body_BMI_Locke2015' or 'c('BMI', 'T1D', 'T2D')'. Case-insensitive. [default %default]"),
  make_option("--n_cores", type="integer", default=5,
              help = "Number of cores to use for parallelization [default %default]")
)

######################################################################
############################## PACKAGES #############################
######################################################################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(WGCNA)) 

message("Packages loaded")

######################################################################
################### GET COMMAND LINE OPTIONS #########################
######################################################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 

opt <- parse_args(OptionParser(option_list=option_list))

seurat_path <- opt$seurat_path 

project_dir <- opt$project_dir

data_prefix <- opt$data_prefix 

autosave <- opt$autosave

resume <- opt$resume

quit_session <- opt$quit_session

metadata_subset_col <- opt$metadata_subset_col

metadata_corr_col <- eval(parse(text=opt$metadata_corr_col))

metadata_corr_filter_vals <- eval(parse(text=opt$metadata_corr_filter_vals))

use.imputed <- opt$use.imputed

regress_out <- eval(parse(text=opt$regress_out))

min.cells <- opt$min.cells

genes_remove_dir <- opt$genes_remove_dir

genes_use <- opt$genes_use

pca_genes <- opt$pca_genes

corFnc <- opt$corFnc 

networkType <- opt$networkType

hclustMethod <- opt$hclustMethod

minClusterSize <- eval(parse(text=opt$minClusterSize))

deepSplit <- eval(parse(text=opt$deepSplit))

moduleMergeCutHeight <- eval(parse(text=opt$moduleMergeCutHeight))

pamStage <- eval(parse(text=opt$pamStage))

kM_reassign <- opt$kM_reassign

jackstrawnReplicate <- opt$jackstrawnReplicate

TOMnReplicate <- opt$TOMnReplicate

fuzzyModMembership <- opt$fuzzyModMembership

scale_MEs_by_kIMs <- opt$scale_MEs_by_kIMs
   
checkPPI <- opt$checkPPI

magma_gwas_dir <- opt$magma_gwas_dir 

gwas_filter_traits <- eval(parse(text=opt$gwas_filter_traits))

data_organism <- opt$data_organism

n_cores <- opt$n_cores  

######################################################################
############################## CONSTANTS #############################
######################################################################

if (!file.exists(seurat_path) & is.null(resume)) stop("Input data path not found and no previous session to resume")

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

scratch_dir = "/scratch/tmp-wgcna/"

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

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

# Load parameter values and utility functions
source(file = paste0(current.dir, "rwgcna_params.R"))
source(file = paste0(current.dir, "rwgcna_functions.R"))

# Entrezgene to ensemb_gene_id mapping file path for MAGMA
mapping_hs_entrez2ensembl_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz"

if (data_organism == "mmusculus") {
  # For mapping symbol to ensembl
  mapping_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
  # Synonyms
  mapping_mm_synonyms_filepath = "/data/genetic-mapping/ncbi/Mus_musculus.gene_info_symbol2ensembl.gz"
  # Mouse to human ortholog mapping  
  mapping_hs_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
  
  # mendelian and coding variants geneset path
  variants_path <- "/projects/jonatan/genesets/variants_ensembl_mm.RData"
  
} else if (data_organism == "hsapiens") {
  
  mapping_hs_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz" # columns: ensembl_gene_id, entrezgene, hgnc_symbol
  
  # mendelian and coding variants geneset path
  variants_path <- "/projects/jonatan/genesets/variants_ensembl_hs.RData"
  
}

try(disableWGCNAThreads())

if (is.null(resume)) {
  
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
  
  if (!networkType %in% c('signed', 'unsigned', 'signed hybrid')) stop("networkType must be one of 'signed', 'unsigned' or 'signed hybrid' (not 'signed_hybrid')")
  
  if (length(c(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight))>4) warning("Comparing different parameters increases processing time")
  
  if (min(minClusterSize) < 5) stop("minClusterSize must be a vector of integers over 5")
  
  if (max(deepSplit) > 4 | min(deepSplit) < 0) stop("deepSplit must be a vector of integers between 0 and 4")
  
  if (min(moduleMergeCutHeight) < 0 | max(moduleMergeCutHeight) > 1) stop("moduleMergeCutHeight must be a vector of doubles between 0 and 1, recommended range is between 0.1 and 0.2")
  
  if (!fuzzyModMembership %in% c("kME", "kIM")) {
    warning("Invalid fuzzeModuleMembership value -  reset to kME (default)")
    fuzzyModMembership <- "kME"
  }
  
  if (! (jackstrawnReplicate >= 0)) stop("jackstrawnReplicate must be 0 or higher")
  
  if (! (TOMnReplicate >= 0)) stop("TOMnReplicate must be non-negative")
  
  if (! (n_cores >= 0 & n_cores <= 100)) stop("n_cores must be in the range 0-100")
  
  if (!is.null(genes_remove_dir)) {
    if (!file.exists(genes_remove_dir)) {
      warning("genes_remove_dir not found")
      genes_remove_dir <- NULL
    } 
  }
  
  if (fuzzyModMembership=="kIM" & scale_MEs_by_kIMs==T) scale_MEs_by_kIMs <- F
  
  if (TODO) {
    scale_MEs_by_kIMs = F
  }
  ######################################################################
  ################# LOAD AND SUBSET SEURAT OBJECT ######################
  ######################################################################
  
  message("Loading and subsetting seurat object..")
  
  seurat_obj <- load_obj(f=seurat_path)
  
  ######################################################################
  ######################## USE IMPUTED DATA? ###########################
  ######################################################################
  
  # Use imputed data if user opted to do so 
  if (is.null(seurat_obj@imputed) | any(dim(seurat_obj@imputed)) == FALSE) use.imputed <- F
  if (use.imputed==T) seurat_obj@raw.data <- seurat_obj@imputed
  
  ######################################################################
  ################### COMPUTE PERCENT MITO, PERCENT RIBO ###############
  ######################################################################
  
  
  # Identify mitochrondrial and ribosomal genes
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = seurat_obj@raw.data), value = TRUE, ignore.case=T)
  ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = seurat_obj@raw.data), value = TRUE)
  percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito.genes, ])/Matrix::colSums(seurat_obj@raw.data)
  percent.ribo <- Matrix::colSums(seurat_obj@raw.data[ribo.genes, ])/Matrix::colSums(seurat_obj@raw.data)
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.ribo, col.name = "percent.ribo")
  
  ######################################################################
  #################### FILTER OUT GENES IN genes_remove_dir ############
  ######################################################################
  
  if (!is.null(genes_remove_dir)) {
    genes_remove_paths <- dir(genes_remove_dir, pattern = "\\.csv|\\.tab|\\.txt", full.names = T)
    
    if (length(genes_remove_paths) > 0) {
      
      remove_idx_all <- rep(FALSE, times = nrow(seurat_obj@raw.data))
      
      for (path in genes_remove_paths) {
        
        geneset <- NULL
        
        if (grepl(pattern = "\\.gz", x = path)) {
          if (grepl(pattern = "\\.csv", x = path)) {
            geneset <- read.csv(file=gzfile(geneset),sep="\t", stringsAsFactors = F) 
          } else if (grepl(pattern = "\\.tab", x = path)) {
            geneset <- read.table(file=gzfile(geneset),sep="\t", stringsAsFactors = F)
          } else if (grepl(pattern = "\\.txt", x = path)) {
            geneset <- read.delim(file=gzfile(geneset), sep = "\t", stringsAsFactors = F)
          }
        } else if (grepl(pattern = "\\.csv", x = path)) {
          geneset <- read.csv(file=geneset,sep="\t", stringsAsFactors = F) 
        } else if (grepl(pattern = "\\.tab", x = path)) {
          geneset <- read.table(file=geneset,sep="\t", stringsAsFactors = F)
        } else if (grepl(pattern = "\\.txt", x = path)) {
          geneset <- read.delim(file=geneset,sep="\t", stringsAsFactors = F)
        }
        
        if (!is.null(geneset)) {
          
          remove_idx <- rownames(seurat_obj@raw.data) %in% geneset
          remove_idx_all <- remove_idx | remove_idx_all 
          
        }
      }
    }
    seurat_obj@raw.data <- seurat_obj@raw.data[!remove_idx,]
  }
  
  ######################################################################
  #################### SET SEURAT IDENT AND SUBSET NAMES ###############
  ######################################################################
  
  # Set seurat object ident
  if (!is.null(metadata_subset_col)) {
    seurat_obj <- SetAllIdent(object = seurat_obj,id = metadata_subset_col)
  }
  
  sNames <- names(table(seurat_obj@ident))
  
  ######################################################################
  ################# IF SYMBOL, REMAP TO ENSEMBL ID #####################
  ######################################################################
  
  if (!any(grepl("ENSG|ENSMUSG",rownames(seurat_obj@raw.data)))) {
    
    if (is.null(data_organism)) {
      stop("Error: seurat object gene names are not in ensembl ID format. To remap from hgnc to ensembl, please provide a data_organism 'mmusculus' or 'hsapiens'")
      
    } else if (!is.null(data_organism)) {
      
      if (data_organism == "mmusculus") { 
        
        # Step 1: direct mapping
        mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
        mapping = data.frame(symbol=row.names(seurat_obj@raw.data), ensembl = mapping_direct$ensembl_gene_id[ match(toupper(row.names(seurat_obj@raw.data)), toupper(mapping_direct$gene_name_optimal)) ])
        
        # Step 2: map remaining using synonyms
        mapping_synonyms = read.delim(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T, stringsAsFactors = F)
        mapping$ensembl[ which(is.na(mapping$ensembl)) ] = mapping_synonyms$ensembl[ match( toupper(mapping$symbol[which(is.na(mapping$ensembl)) ]) , toupper(mapping_synonyms$symbol)) ]
        
        # save mapping file for reference and later use
        write.csv(mapping, file=sprintf("%s%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, data_organism), quote = F, row.names = F)
        
      } else if (data_organism == "hsapiens") {
        
        mapping_direct = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T, stringsAsFactors = F) # columns: ensembl_gene_id, entrezgene, hgnc_symbol
        # Step 1: direct mapping
        mapping = data.frame(symbol=row.names(seurat_obj@raw.data), ensembl = mapping_direct$ensembl_gene_id[ match(toupper(row.names(seurat_obj@raw.data)), toupper(mapping_direct$hgnc_symbol)) ])
        
        # save mapping file for reference and later use
        write.csv(mapping, file=sprintf("%s%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, data_organism), quote = F, row.names = F)
        
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
      #seurat_obj@data <- seurat_obj@data[!is.na(mapping$ensembl) & !idx_duplicate_genes,]
      seurat_obj@raw.data <- seurat_obj@raw.data[!is.na(mapping$ensembl) & !idx_duplicate_genes,]
      
      # rename Seurat object rows where mapping was successful to ensembl ID
      #rownames(seurat_obj@data) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
      rownames(seurat_obj@raw.data) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
      
    }
  } else if (any(grepl("ENSG", rownames(seurat_obj@raw.data)))) {
    
    message("Homo sapiens ensembl id gene names detected")
    data_organism <- "hsapiens"
    
  } else if (any(grepl("ENSMUSG", rownames(seurat_obj@raw.data)))) {
    
    message("mus musculus ensembl id gene names detected")
    data_organism <- "mmusculus"
    
  }
  
  ######################################################################
  ######## DO SEURAT PROCESSING ON FULL EXPRESSION MATRIX ##############
  ######################################################################
  
  # Normalise 
  seurat_obj <- NormalizeData(object = seurat_obj, display.progress = T)
  
  # Scale and regress out confounders
  vars.to.regress = if (!is.null(regress_out)) regress_out[regress_out %in% names(seurat_obj@meta.data)] else NULL
  seurat_obj <- ScaleData(object = seurat_obj, 
                          vars.to.regress = vars.to.regress,
                          model.use="linear",
                          do.par=T,
                          num.cores = min(n_cores, detectCores()-1),
                          do.scale=T,
                          do.center=T)
  
  scale_data <- seurat_obj@scale.data
  ident <- seurat_obj@ident
  
  # Save scale and regressed whole expression matrix with ensembl rownames for later use
  save(scale_data, file = sprintf("%stmp_%s_scale_regr_data_ensembl.RData", scratch_dir, data_prefix))
  save(ident, file = sprintf("%stmp_%s_ident_%s.RData", scratch_dir, data_prefix, flag_date))
  
  ######################################################################
  ######## EXTRACT METADATA AND CONVERT FACTORS TO MODEL MATRIX ########
  ######################################################################
  
  # Convert any character or factor meta.data to numeric dummy variables each level with its own numeric column
  if (!is.null(metadata_corr_col)) {
    if (any(colnames(seurat_obj@meta.data) %in% metadata_corr_col)) {
      metadata <- matrix(NA, nrow=nrow(seurat_obj@meta.data), ncol=1)
      include <- seurat_obj@meta.data[,colnames(seurat_obj@meta.data) %in% metadata_corr_col, drop=F]
      for (i in 1:ncol(include)) {
        if (class(include[,i]) %in% c("factor", "character")) {
          metadata <- cbind(metadata, factorToIndicator(include[,i, drop=T]))        
        } else {
          metadata <- cbind(metadata, include[,i, drop=F])
        }
      }
      #rownames(metadata) <- rownames(seurat_obj@meta.data)
      metadata <- metadata[,-1, drop=F]
      if (!is.null(metadata_corr_filter_vals)) metadata <- metadata[, toupper(colnames(metadata)) %in% toupper(metadata_corr_filter_vals), drop = F]
      metadata <- as.data.frame(metadata)
      
      # Filter out any metadata columns where all the values are identical
      metadata <- metadata[apply(metadata, MARGIN=2, FUN = function(x) length(unique(x))>1)]
      if (ncol(metadata) == 0) metadata <- NULL    
      rownames(metadata) = rownames(seurat_obj@meta.data)
    } else metadata <- NULL
  } else metadata <- NULL
  
  ######################################################################
  ######## DO SEURAT PROCESSING ON SUBSETTED EXPRESSION MATRICES #######
  ######################################################################
  
  # Subset the seurat object
  # NB: We avoid parallel subsetting due to the RAM considerations
  
  subsets <- lapply(sNames, function(x) SubsetData(seurat_obj,
                                                   ident.use = x, 
                                                   do.scale = F, 
                                                   do.center = F,
                                                   subset.raw = T,
                                                   do.clean = F))
  
  ### Filter out cell types that have too few cells (<20).
  # We do this do avoid downstream problems with Seurat or WGCNA. 
  # E.g. Seurat will ScaleData will fail if regressing out variables when there are only 2 cells in the data.
  
  subsets_ok_idx <- sapply(subsets, function(x) {
    if(ncol(x@data)<=20) {
      warning(sprintf("%s cell cluster filtered out because it contains <20 cells",unique(as.character(x@ident))[1] ))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }, simplify = T)
  
  subsets <- subsets[subsets_ok_idx]
  sNames <- sNames[subsets_ok_idx]
  names(subsets) <- sNames 
  
  # Filter genes expressed in fewer than min.cells in a subset as these will also lead to spurious associations and computational difficulties
  subsets <- lapply(subsets, function(x) FilterGenes(x, min.cells = min.cells))
  
  # Scale data and regress out confounders
  vars.to.regress = if (!is.null(regress_out)) regress_out[regress_out %in% names(seurat_obj@meta.data)] else NULL
  
  subsets <- lapply(subsets, function(x) ScaleData(object = x,
                                                   vars.to.regress = vars.to.regress,
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
      subsets <- parLapplyLB(cl, subsets, function(x) RunPCA(object = x,
                                                             pc.genes = if (pca_genes == 'all') rownames(x@data) else x@var.genes,
                                                             pcs.compute = min(nPC_seurat, (if (pca_genes == 'all') nrow(x@data) else length(x@var.genes)) %/% 2, ncol(x@data) %/% 2),
                                                             use.imputed = F, # if use.imputed=T the @imputed slot has been copied to @data
                                                             weight.by.var = F,
                                                             do.print = F,
                                                             seed.use = randomSeed,
                                                             maxit = maxit, # set to 500 as default
                                                             fastpath = fastpath)) 
      
    }, warning = function(c) {
      subsets <- parLapplyLB(cl, subsets, function(x) RunPCA(object = x,
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
    invisible(gc()); invisible(R.utils::gcDLLs())
    end_time <- Sys.time()
    
    message(sprintf("PCA done, time elapsed: %s seconds", round(end_time - start_time,2)))
    
    # Select significant PCs using empirical p-value based on JackStraw resampling
    # Source: https://rdrr.io/cran/Seurat/src/R/plotting.R
    # score the PCs using JackStraw resampling to get an empirical null distribution to get p-values for the PCs based on the p-values of gene loadings
    if (jackstrawnReplicate > 0 & genes_use == "PCA") message(sprintf("Performing JackStraw with %s replications to select genes that load on significant PCs", jackstrawnReplicate))
    list_datExpr <- lapply(subsets, function(x) wrapJackStraw(seurat_obj_sub = x, n_cores = n_cores, jackstrawnReplicate = jackstrawnReplicate, pvalThreshold = pvalThreshold))
    
    invisible(gc())
    
  } else if (genes_use == "all") {
    list_datExpr <- lapply(subsets, function(x) 
      seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rep(TRUE, nrow(x@scale.data))))
  } else if (genes_use == "var.genes") {
    list_datExpr <- lapply(subsets, function(x) seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rownames(x@scale.data) %in% x@var.genes))
  } 
  
  ######################################################################
  ##################### SAVE PARAMETERS TO FILE ########################
  ######################################################################
  
  diagnostic_stats <- data.frame(subset = sNames, 
                                 n_cells = sapply(subsets, function(x) ncol(x@data), simplify = T),
                                 n_genes = sapply(list_datExpr, function(x) ncol(x), simplify=T))
  
  
  userParams <- cbind(list("seurat_path" = seurat_path,
                           "project_dir" = project_dir,
                           "data_prefix" = data_prefix,
                           "autosave" = autosave,
                           "resume" = resume,
                           "quit_session" = quit_session,
                           "metadata_subset_col" = metadata_subset_col,
                           "metadata_corr_col" = metadata_corr_col,
                           "use.imputed" = use.imputed,
                           "regress_out" = regress_out,
                           "min.cells" = min.cells,
                           "genes_remove_dir" = genes_remove_dir,
                           "genes_use" = genes_use,
                           "pca_genes" = pca_genes,
                           "corFnc" = corFnc,
                           "networkType" = networkType,
                           "hclustMethod" = hclustMethod,
                           "minClusterSize" = minClusterSize,
                           "deepSplit" = deepSplit,
                           "pamStage" = pamStage,
                           "kM_reassign" = kM_reassign,
                           "moduleMergeCutHeight" = moduleMergeCutHeight,
                           "fuzzyModMembership" = fuzzyModMembership,
                           "scale_MEs_by_kIMs" = scale_MEs_by_kIMs,
                           "jackstrawnReplicate" = jackstrawnReplicate,
                           "TOMnReplicate" = TOMnReplicate,
                           "checkPPI" = checkPPI,
                           "data_organism" = data_organism,
                           "magma_gwas_dir" = magma_gwas_dir,
                           "n_cores" = n_cores))
  
  builtInParams <- cbind(list("minCoreKME" = minCoreKME,
                              "minKMEtoStay" = minKMEtoStay,
                              "consensusQuantile" = consensusQuantile,
                              "fraction" = fraction,
                              "do.center" = do.center,
                              "replace" = replace,
                              "impute" = impute,
                              "nPC_seurat" = nPC_seurat,
                              "TOMType" = TOMType,
                              "pvalThreshold" = pvalThreshold))
  
  # Save run params to file
  write.csv(userParams, file=sprintf("%s%s_INFO_user_parameters_%s.csv", log_dir, data_prefix, flag_date), quote = F, row.names = T)
  write.csv(builtInParams, file=sprintf("%s%s_INFO_built_in_parameters_%s.csv", log_dir, data_prefix, flag_date), quote = F, row.names = T)
  
  
  
  # Clear up
  rm(subsets, userParams, builtInParams, seurat_obj)
  invisible(gc())
  
  # Save or load session image 
  resume = "checkpoint_1"
  if (autosave==T) save.image(file=sprintf("%s%s_checkpoint_1_image.RData", scratch_dir, data_prefix))
  if (!is.null(quit_session)) if (quit_session=="checkpoint_1") quit(save="no")
  
} else if (!is.null(resume)) {
  
  if (resume == "checkpoint_1") {
    
    tryCatch({load((file=sprintf("%s%s_checkpoint_1_image.RData", scratch_dir, data_prefix)))},
             error = function(x) {stop(paste0(resume, " session image file not found in ", scratch_dir))})
    
    source(file = paste0(current.dir, "rwgcna_params.R"))
    source(file = paste0(current.dir, "rwgcna_functions.R"))
    options(stringsAsFactors = F)
    disableWGCNAThreads() 
  }
}

if (resume == "checkpoint_1") {
  
  ######################################################################
  ###################### (UN) LOAD PACKAGES ############################
  ######################################################################
  
  # Unload packages
  try(detach("package:Seurat", unload=TRUE))
  
  # Free up DLLs
  invisible(R.utils::gcDLLs())
  # avoid multithreading in WGCNA
  suppressPackageStartupMessages(library(WGCNA)) 
  disableWGCNAThreads()
  
  ######################################################################
  ####### PICK SOFT THRESHOLD POWER FOR ADJACENCY MATRIX ###############
  ######################################################################
  message("Computing soft threshold powers to maximise the fit of a scale free topology to the adjacency matrix")
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  softPower <- 8 # Set a default value as fall back
  
  cl <- makeCluster(n_cores, 
                    type = "FORK", 
                    outfile = paste0(log_dir, "log_compute_softPowers.txt"))
  
  list_softPower <- clusterMap(cl, function(datExpr,
                                            subsetName) { 
    
    powers = c(1:30)
    sft = pickSoftThreshold(data=datExpr,
                            powerVector = powers,
                            blockSize = min(maxBlockSize , ncol(datExpr)), #try to prevent crashing
                            corFnc = corFnc,
                            corOptions =  corOptions,
                            networkType = networkType,
                            verbose = verbose)
    pdf(sprintf("%s%s_%s_pickSoftThresholdSFTFit_%s.pdf", plots_dir, data_prefix, subsetName, flag_date),width=10,height=5)
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
    pdf(sprintf("%s%s_%s_pickSoftThresholdMeanCon_%s.pdf", plots_dir, data_prefix, subsetName, flag_date),width=10,height=5)
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
    if (sum(fitIndices_filter$SFT.R.sq >= 0.9) > 1) {
      softPower = min(fitIndices_filter$Power[fitIndices_filter$SFT.R.sq>=0.9])
    } else {
      (fitIndices_filter %>% dplyr::arrange(desc(SFT.R.sq)) %>% dplyr::select(Power))[1,] -> softPower 
    }
    return(softPower)
  }, 
  datExpr = list_datExpr, 
  subsetName = sNames, 
  SIMPLIFY=F,
  .scheduling = c("dynamic"))
  
  stopCluster(cl)
  
  ###
  invisible(gc()); invisible(R.utils::gcDLLs())
  names(list_softPower) <- sNames
  
  # TOM files will be saved here by consensusTOM / blockwiseConsensusModules
  setwd(scratch_dir)
  
  if (TOMnReplicate > 0) {
    
    ######################################################################
    ############### RESAMPLE THE DATA FOR ROBUSTNESS #####################
    ######################################################################
    
    message(sprintf("Resampling the samples %s times with replace = %s", TOMnReplicate, replace))
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_bootStrap.txt"))
    
    list_multiExpr <- parLapplyLB(cl, list_datExpr, function(x) bootstrap(datExpr=x, 
                                                                          ### 180508_v1.8_dev2
                                                                          nPermutations = TOMnReplicate,
                                                                          ###
                                                                          replace = replace,
                                                                          fraction = fraction,
                                                                          randomSeed = randomSeed))
    stopCluster(cl)
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    names(list_multiExpr) <- sNames 
    
    ######################################################################
    ######### RUN CONSENSUSTOM ON THE RESAMPLED DATASETS #################
    ######################################################################
    
    message(sprintf("Computing the consensus Topological Overlap Matrix with %s permutations", TOMnReplicate))
    
    invisible(gc())
    
    # Set a ceiling on n_cores to avoid depleting memory
    
    cl <- makeCluster(min(n_cores,10), type = "FORK", outfile = paste0(log_dir, "log_consensusTOM.txt"))
    list_consensus <- clusterMap(cl, function(x,y) consensusTOM(multiExpr = x, 
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
                                                                cacheDir = scratch_dir,
                                                                cacheBase = ".blockConsModsCache",
                                                                verbose = verbose,
                                                                indent = indent), 
    x=list_multiExpr, 
    y=sNames, 
    SIMPLIFY=F,
    .scheduling = c("dynamic"))
    stopCluster(cl)
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    # For each consensusTOM, get a logical vector where 'good_genes' that were used are TRUE
    list_goodGenesTOM_idx <- lapply(list_consensus, function(x) as.logical(x$goodSamplesAndGenes$goodGenes))
    
    # Use the goodGenesTOM_idx to filter the datExpr matrices
    list_datExpr_gg <- mapply(FUN=function(x,y) x[,y], x=list_datExpr, y=list_goodGenesTOM_idx, SIMPLIFY=F)
    
  } else if (TOMnReplicate==0) {
    
    message("Computing the Topological Overlap Matrix")
    
    ######################################################################
    ##### COMPUTE THE ADJACENCY MATRIX AND TOM WITHOUT RESAMPLING ########
    ######################################################################
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_TOM_for_par.txt"))
    list_goodGenesTOM_idx <- clusterMap(cl, function(x,y) TOM_for_par(datExpr=x, 
                                                                      subsetName=y, 
                                                                      softPower=list_softPower[[y]]), 
                                        x=list_datExpr, 
                                        y=sNames, 
                                        SIMPLIFY=F,
                                        .scheduling = c("dynamic"))
    stopCluster(cl)
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    list_datExpr_gg <- list_datExpr
    
  }
  
  if (TOMnReplicate > 0) rm(list_multiExpr) #TODO: we'll need it for plotting permuted
  
  ######################################################################
  ############################ CHECKPOINT ##############################
  ######################################################################
  
  resume="checkpoint_2"
  
  if (autosave==T) save.image( file=sprintf("%s%s_checkpoint_2_image.RData", scratch_dir, data_prefix))
  if (!is.null(quit_session)) if (quit_session=="checkpoint_2") quit(save="no")
  
} else if (!is.null(resume)) {
  if (resume == "checkpoint_2") {
    
    tryCatch({load((file=sprintf("%s%s_checkpoint_2_image.RData", scratch_dir, data_prefix)))},
             error = function(x) {stop(paste0(resume, " session image file not found in ", scratch_dir))})
    
    source(file = paste0(current.dir, "rwgcna_params.R"))
    source(file = paste0(current.dir, "rwgcna_functions.R"))
    options(stringsAsFactors = F)
  }
}

if (resume == "checkpoint_2") {
  
  ######################################################################
  #################### (UN) LOAD PACKAGES #############################
  ######################################################################
  
  # Free up DLLs
  invisible(R.utils::gcDLLs())
  
  ######################################################################
  ####################### CLUSTER ON THE TOM ###########################
  ######################################################################
  # in:
  #   consTomDS (from disk!)
  #
  # out:
  #   list_geneTree
  
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parHclust.txt"))
  # Convert TOM to distance matrix
  list_dissTOM <- parLapplyLB(cl, sNames, dissTOM_for_par)
  # Cluster
  list_geneTree <- parLapplyLB(cl, list_dissTOM, function(x) hclust(d=x, method=hclustMethod))
  stopCluster(cl)
  invisible(gc()); invisible(R.utils::gcDLLs())
  names(list_geneTree) = sNames # used for PlotDendro
  
  ######################################################################
  ######################### COMPUTE MODULES ############################
  ######################################################################
  
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
  
  invisible(gc()); invisible(R.utils::gcDLLs())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_cutreeHybrid_for_vec.txt"))
  # nb: comb_list is in the environment - no need to pass it through clusterMap
  list_list_cutree <- clusterMap(cl, function(x,y) lapply(comb_list, function(z) cutreeHybrid_for_vec(comb=z, 
                                                                                                      geneTree=x, 
                                                                                                      dissTOM = y,
                                                                                                      maxPamDist = maxPamDist, 
                                                                                                      useMedoids = useMedoids)), 
                                 x=list_geneTree,
                                 y=list_dissTOM, 
                                 SIMPLIFY=F,
                                 .scheduling = c("dynamic"))
  stopCluster(cl)
  invisible(gc()); invisible(R.utils::gcDLLs())
  names(list_list_cutree) = sNames
  
  # Produce labels for plotting the modules found by different parameters
  list_plot_label <- lapply(comb_list, function(x) plotLabel_for_vec(comb=x)) # list of labels for plot
  
  # Make a copy of the labels for each subset
  list_list_plot_label <- list()
  list_list_plot_label <- lapply(sNames, function(x) list_list_plot_label$x = list_plot_label)
  names(list_list_plot_label) = sNames
  
  ######################################################################
  ######################### MERGE CLOSE MODULES ########################
  ######################################################################
  
  message(paste0("Merging modules with a distance below ", moduleMergeCutHeight))
    
  if (fuzzyModMembership=="kME") {
    
    # Merge close modules
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parvecMergeCloseModules.txt"))
    
    
    
    list_list_merged <- parLapplyLB(cl, sNames, function(x) lapply(1:n_combs, function(y) mergeCloseModules_for_vec(cutree=list_list_cutree[[x]][[y]],
                                                                                                                    comb=comb_list[[y]],
                                                                                                                    datExpr=list_datExpr_gg[[x]],
                                                                                                                    excludeGrey = excludeGrey,
                                                                                                                    scale_MEs_by_kIMs = scale_MEs_by_kIMs,
                                                                                                                    dissTOM = if (scale_MEs_by_kIMs==T) list_dissTOM[[x]] else NULL)))

    names(list_list_merged) = sNames
    
    # Extract the colors from the list returned by mergeCloseModules
    list_list_colors <- mapply(parGetColors, list_list_merged, list_datExpr_gg, SIMPLIFY=F)
    
    # Extract the Module Eigengenes from the list returned by mergeCloseModules
    list_list_MEs <- lapply(list_list_merged, function(y) lapply(y, function(x) x$MEs)) 
    
    names(list_list_MEs) <- sNames
    
    # Compute kMEs 
    
    list_list_kMs <- clusterMap(cl, 
                                function(x,y) parkMEs(list_MEs=x, 
                                                      datExpr=y), 
                                x=list_list_MEs, 
                                y=list_datExpr_gg, 
                                SIMPLIFY = F,
                                .scheduling = c("dynamic"))
    
    stopCluster(cl)
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    
    
  } else if (fuzzyModMembership=="kIM") { 
    
    list_list_colors <- lapply(list_list_cutree, function(x) lapply(x, function(y) labels2colors(y$labels)))
    
    list_list_colors <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names = colnames(y), dimension=NULL)),
                               x = list_list_colors,
                               y = list_datExpr_gg,
                               SIMPLIFY=F)
    
    list_list_colors <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = y, dimension=NULL),
                               x = list_list_colors,
                               y = list_plot_label,
                               SIMPLIFY=F)
    
    names(list_list_colors) <- sNames
    
    # Merge close modules using correlation between IM embeddings
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_par_mergeCloseModskIM.txt"))
    
    list_list_kMs <- clusterMap(cl, function(a,b) lapply(b, function(x) kIM_eachMod_norm(dissTOM = a, 
                                                                                         colors = x,
                                                                                         verbose=verbose,
                                                                                         excludeGrey=T)),
                                a = list_dissTOM,
                                b = list_list_colors,
                                SIMPLIFY=F,
                                .scheduling = c("dynamic"))
    


    list_list_merged <- mapply(function(a,b,c,d,e) clusterMap(cl, function(x,y) mergeCloseModskIM(datExpr=c,
                                                                                        colors = x,
                                                                                        kIMs  = y,
                                                                                        dissTOM = d,
                                                                                        moduleMergeCutHeight=moduleMergeCutHeight,
                                                                                        verbose=verbose,
                                                                                        celltype = e),
                                                        x = a,
                                                        y = b,
                                                        SIMPLIFY=F,
                                                        .scheduling = c("dynamic")),
                               a = list_list_colors,
                               b = list_list_kMs,
                               c = list_datExpr_gg,
                               d = list_dissTOM,
                               e = names(list_list_colors),
                               SIMPLIFY=F)

    stopCluster(cl)
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    list_list_colors <- lapply(list_list_merged, function(x) lapply(x, function(y) y$colors))
    list_list_kMs <- lapply(list_list_merged, function(x) lapply(x, function(y) y$kIMs))
    list_list_MEs <- NULL
  }
  
  names(list_list_colors) <- names(list_list_kMs) <- sNames
  
  ######################################################################
  ###################### REASSIGN GENES BASED ON kM ####################
  ######################################################################
  # see https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0420-6
  if (kM_reassign ==T ) {
    message(paste0("Reassigning genes based on ", fuzzyModMembership))
    
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_par_kM_reassign.txt"))
    
    list_list_reassign = mapply(function(a,b,c,d) parLapplyLB(cl, a, function(x) kM_reassign_fnc(colors = x,
                                                                                             fuzzyModMembership = fuzzyModMembership,
                                                                                             dissTOM = if (fuzzyModMembership=="kIM" | scale_MEs_by_kIMs==T) b else NULL,
                                                                                             datExpr = if (fuzzyModMembership=="kME") c else NULL,
                                                                                             corFnc = if (fuzzyModMembership=="kME") corFnc else NULL,
                                                                                             verbose=verbose,
                                                                                             max_iter = 10,
                                                                                             celltype = d)),
                                         a = list_list_colors,
                                         b = list_dissTOM,
                                         c = list_datExpr_gg,
                                         d = names(list_list_colors),
                                         SIMPLIFY=F)
    
    stopCluster(cl)
    
    # Extract new colors and kMs from the list of lists of lists returned by the vectorised kM_reassign 
    list_list_colors_reassign <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$colors))
    list_list_kMs_reassign <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$kMs))
    list_list_reassign_log <- lapply(list_list_reassign, function(x) lapply(x, function(y) y$log))
  
    invisible(gc()); invisible(R.utils::gcDLLs())
  } else {
    list_list_colors_reassign <- list_list_colors
    list_list_kMs_reassign <- list_list_kMs
    list_list_reassign_log <- NULL
  }  
  ######################################################################
  ################## extract primary kMEs / kIMs #######################
  ######################################################################
  
  # Extract the primary kMs - i.e. those of each gene w.r.t. the module to which it belongs
  # We use these to filter genes that we submit to STRINGdb

  message(paste0("Computing primary "), fuzzyModMembership, "s")
  
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_par_primary_kMs.txt"))

  list_list_pkMs <- clusterMap(cl, function(x,y) parPkMs(list_kMs=x, 
                                                         list_colors=y), 
                               x=list_list_kMs_reassign, 
                               y=list_list_colors_reassign, 
                               SIMPLIFY = F)
  
  stopCluster(cl)
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  # The pkM vectors already have gene names. Now name each vector, and each list of vectors
  list_list_pkMs <- mapply(function(x,y) name_for_vec(to_be_named = x, 
                                                      given_names = y[[1]], 
                                                      dimension = NULL), 
                           x=list_list_pkMs, 
                           y=list_list_plot_label, 
                           SIMPLIFY=F)
  
  names(list_list_pkMs) <- sNames
  
  ######################################################################
  ############## Match colors between parameter settings ###############
  ######################################################################
  # This step is purely for plotting. 
  # For each celltype, we will anyway take only modules found using one set of parameters
  
  message("Matching module color labels between parameter settings")
  
  # Align / match the colors so similar modules found with different sets of parameters have the same name
  invisible(gc()); invisible(R.utils::gcDLLs())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parMatchColors.txt"))
  list_list_colors_matched <- parLapplyLB(cl, list_list_colors_reassign, function(x) parMatchColors(list_colors = x))
  stopCluster(cl)
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  # Rename highest level list entries (cell clusters)
  names(list_list_colors_matched) = sNames
  
  # Name each entry of the vectors of color assignments with the corresponding genes
  list_list_colors_matched <- mapply(function(x,y) lapply(x, 
                                                          function(z) name_for_vec(to_be_named=z, 
                                                                                   given_names=colnames(y), 
                                                                                   dimension = NULL)), 
                                     x=list_list_colors_matched, 
                                     y=list_datExpr_gg, 
                                     SIMPLIFY=F)
  
  # clear up 
  rm(list_list_colors)
  invisible(gc())
  
  ######################################################################
  ######################### REMOVE ALL-GREY RUNS #######################
  ######################################################################
  
  message("Removing all-grey runs")
  
  # count the grey
  list_vec_n_grey <- lapply(list_list_colors_matched, count_grey_in_list_of_vec)
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  list_logical_params_ok <- mapply(function(x,y) as.logical(mapply(function(a,b) a!=length(b), a=x, b=y, SIMPLIFY=F)), 
                                   x=list_vec_n_grey, 
                                   y=list_list_colors_matched, 
                                   SIMPLIFY=F)
  logical_subsets_ok <- sapply(list_logical_params_ok, any, simplify = T)
  
  # First, for each run, remove results of any parametrisation for which all the genes were assigned to grey
  # If for a subset all parameters gave only grey modules, take it out of the top level list.
  list_list_plot_label_ok <- mapply(function(x,y) x[y], x = list_list_plot_label, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_cutree_ok <- mapply(function(x,y) x[y], x = list_list_cutree, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_MEs_ok <- if (fuzzyModMembership=="kME")  mapply(function(x,y) x[y], x = list_list_MEs, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok] else NULL
  list_list_colors_matched_ok <- mapply(function(x,y) x[y], x=list_list_colors_matched, y=list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_pkMs_ok <- mapply(function(x,y) x[y], x=list_list_pkMs, y=list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_reassign_log <- if (!is.null(kM_reassign)) mapply(function(x,y) x[y], x=list_list_reassign_log, y=list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok] else NULL
  
  sNames_ok <- sNames[logical_subsets_ok]
  
  list_datExpr_ok <- list_datExpr_gg[logical_subsets_ok] # If for a subset all parameters gave only grey modules, take it out of the top level list.
  list_geneTree_ok <- list_geneTree[logical_subsets_ok]
  list_dissTOM_ok <- list_dissTOM[logical_subsets_ok]

  # Assign gene names to each color vector
  list_list_colors_matched_ok <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), 
                                                                                         dimension = NULL)), x = list_list_colors_matched_ok, y = list_datExpr_ok, SIMPLIFY=F)
  # Make a warning and report if any whole subsets produced no modules
  discarded <- setdiff(sNames, sNames_ok)
  if (length(discarded > 0)) {
    warning(paste0(discarded, " were dropped because the script didn't find any modules  "))
    fileConn <- file(description = sprintf("%s%s_ERROR_report_%s.txt", log_dir,  data_prefix, flag_date), open = 'a')
    writeLines(paste0(discarded, " were dropped because the script didn't find any modules"), fileConn)
    close(fileConn)
  }
  
  # Keep dissTOM for computing kIMs (just in case, also if fuzzyModMembership = kME)
  save(list_dissTOM_ok, file=sprintf("%s%s_list_dissTOM_ok.RData", scratch_dir, data_prefix))
  
  rm(list_datExpr, list_list_cutree, list_list_MEs, list_datExpr_gg, list_geneTree, list_list_plot_label, list_dissTOM)#, list_idx_mapped_gg)
  if (fuzzyModMembership != "kIM")  rm(list_dissTOM_ok) # remove from memory unless we need it soon to compute kIMs
  
  invisible(gc())
  
  ######################################################################
  ############################# CHECKPOINT #############################
  ######################################################################
  
  resume = "checkpoint_3"
  message("Reached checkpoint 3, saving session image")
  
  if (autosave==T) save.image( file=sprintf("%s%s_checkpoint_3_image.RData", scratch_dir, data_prefix))
  if (!is.null(quit_session)) if (quit_session=="checkpoint_3") quit(save="no")
  
} else if (resume == "checkpoint_3") {
  
  tryCatch({load(file=sprintf("%s%s_checkpoint_3_image.RData", scratch_dir, data_prefix))},
           error = function(x) {stop(paste0(resume, " session image file not found in ", scratch_dir))})
  
  # Load parameter values and utility functions anew (in case a bug was fixed)
  source(file = paste0(current.dir, "rwgcna_params.R"))
  source(file = paste0(current.dir, "rwgcna_functions.R"))
  options(stringsAsFactors = F)
}  
  
if (resume == "checkpoint_3") {
  
  ##########################################################################
  ############################ (UN)LOAD PACKAGES ############################
  ##########################################################################
  
  # Free up DLLs
  #try(detach("package:STRINGdb", unload=TRUE))
  invisible(R.utils::gcDLLs())
  # Load packages
  suppressPackageStartupMessages(library(STRINGdb))
  
  ######################################################################
  #################### CHECK MODULES FOR PPI ENRICHMENT ################
  ######################################################################
  # in:
  #   list_list_colors_matched_ok
  #   list_list_pkMs
  #
  # out:
  #   list_list_colors_matched_ok_PPI
  
  message("Checking modules for significant Protein-Protein Interactions through STRINGdb")
  
  if (data_organism == "hsapiens") STRINGdb_species <- 9606 else if (data_organism == "mmusculus") STRINGdb_species <- 10090
  
  if (fuzzyModMembership=="kME") {
    # Value already set in parameters; just make it into a vector of identical values, one entry per celltype 
    list_PPI_pkM_threshold <- as.list(rep(PPI_pkME_threshold, times = length(list_list_pkMs_ok)))
  } else if (fuzzyModMembership=="kIM") {
    # If using pkIMs to filter out genes with low module membership from the lists submitted to the PPI analysis, for each
    # cellType subset (upper level of the list_list_pkMs), compute the 0.05 quantile of gene primary module membership values (pkIMs)
    # and set that as a threshold for that datatype. Code below returns a numeric vector of quantile values.
    list_PPI_pkM_threshold <- lapply(list_list_pkMs_ok, function(x) quantile(unlist(x, use.names = F), probs = c(0.05), names=F))
  }
  
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_PPI_outer_for_vec.txt"))
  
  list_list_PPI <- mapply(function(a,b,c) clusterMap(cl, function(x,y) PPI_outer_for_vec(colors = x,
                                                                                         pkMs = y,
                                                                                         STRINGdb_species = STRINGdb_species,
                                                                                         PPI_pkM_threshold = c,
                                                                                         pvalThreshold = pvalThreshold,
                                                                                         project_dir = project_dir,
                                                                                         data_prefix = data_prefix,
                                                                                         flag_date = flag_date),
                                                         x = a,
                                                         y = b,
                                                         SIMPLIFY=F, .scheduling = c("dynamic")),
                              a = list_list_colors_matched_ok,
                              b = list_list_pkMs_ok,
                              c = list_PPI_pkM_threshold,
                              SIMPLIFY=F)
  
  list_list_colors_PPI <- if (checkPPI==T) lapply(list_list_PPI, function(x) lapply(x, function(y) y$colors_PPI)) else if (checkPPI==F) list_list_colors_matched_ok
  
  list_list_module_PPI <- lapply(list_list_PPI, function(x) lapply(x, function(y) data.frame(y$module_PPI, stringsAsFactors = F)))
  list_list_module_PPI_signif <- lapply(list_list_PPI, function(x) lapply(x, function(y) data.frame(y$module_PPI_signif, stringsAsFactors = F)))
  
  stopCluster(cl)
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  names(list_list_colors_PPI) <- names(list_list_module_PPI) <- names(list_list_module_PPI_signif) <- list_list_colors_matched_ok 
  
  # Name the PPI enrichment results as they will be output in this state; 
  # we will however soon remove a layer from the list_list_colors_PP double nested list of vector so no need to name them
  list_list_module_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension=NULL),
                                 x=list_list_module_PPI,
                                 y=list_list_colors_matched_ok, 
                                 SIMPLIFY = F)
  
  list_list_module_PPI_signif <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension=NULL),
                                        x=list_list_module_PPI_signif,
                                        y=list_list_colors_matched_ok, 
                                        SIMPLIFY = F)
  
  ######################################################################
  ############################# CHECKPOINT #############################
  ######################################################################
  
  resume = "checkpoint_4"
  message("Reached checkpoint 4, saving session image")
  
  if (autosave==T) save.image( file=sprintf("%s%s_checkpoint_4_image.RData", scratch_dir, data_prefix))
  if (!is.null(quit_session)) if (quit_session=="checkpoint_4") quit(save="no")
  
} else if (resume == "checkpoint_4") {
  
  tryCatch({load(file=sprintf("%s%s_checkpoint_4_image.RData", scratch_dir, data_prefix))},
           error = function(x) {stop(paste0(resume, " session image file not found in ", scratch_dir))})
  
  # Load parameter values and utility functions anew (in case a bug was fixed)
  source(file = paste0(current.dir, "rwgcna_params.R"))
  source(file = paste0(current.dir, "rwgcna_functions.R"))
  options(stringsAsFactors = F)
}

if (resume == "checkpoint_4") {
  
  ######################################################################
  ######################### LOAD PACKAGES ##############################
  ######################################################################
  
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(liger))
  suppressPackageStartupMessages(library(boot))
  
  ######################################################################
  ########### ORDER PARAMETER SETS BY PPI ENRICHMENT AND PLOT ##########
  ######################################################################
  
  message("Selecting parameters with the highest number of genes assigned to modules significantly enriched for Protein-Protein Interactions")
  
  # Count how many genes were assigned to grey under each parametrisation and order the parametrisations
  list_PPI_vec_n_grey <- lapply(list_list_colors_PPI, count_grey_in_list_of_vec)
  
  # Order all the outputs by how many genes were assigned to a (non-grey) module
  list_list_reassign_log_order <- if (!is.null(kM_reassign)) mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_reassign_log, y = list_PPI_vec_n_grey, SIMPLIFY=F) else NULL
  list_list_plot_label_ok_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_plot_label_ok, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_colors_matched_ok_order <- mapply(function(x,y) x[order(y, decreasing=F)], x  =  list_list_colors_matched_ok , y = list_PPI_vec_n_grey, SIMPLIFY = F )
  list_list_colors_PPI_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_colors_PPI, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_module_PPI_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_module_PPI, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_module_PPI_signif_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_module_PPI_signif, y = list_PPI_vec_n_grey, SIMPLIFY=F)
 
  ######################################################################
  ##### FOR EACH SUBSET SELECT PARAMETERS WITH BEST PPI ENRICHMENT #####
  ######################################################################
  
  # Eliminate a layer of nesting by selecting only the best parametrisation per celltype
  list_plot_label_final <- lapply(list_list_plot_label_ok_order, function(x) x[[1]])
  list_reassign_log <- if (!is.null(kM_reassign)) lapply(list_list_reassign_log_order, function(x) x[[1]]) else NULL
  list_colors_PPI <- lapply(list_list_colors_PPI_order, function(x) x[[1]])
  list_colors <- lapply(list_list_colors_matched_ok_order, function(x) x[[1]])
  list_module_PPI <- lapply(list_list_module_PPI_order, function(x) x[[1]])
  list_module_PPI_signif <- lapply(list_list_module_PPI_signif_order, function(x) x[[1]])
  
  # Name by cell clusters
  names(list_plot_label_final) <-  names(list_colors_PPI) <- names(list_colors) <- names(list_module_PPI) <- names(list_module_PPI_signif) <- sNames_ok
  if (!is.null(kM_reassign)) names(list_reassign_log) <-sNames_ok
    
  # Make list of list of final parameters
  param_names = c("minClusterSize", "deepSplit","pamStage", "moduleMergeCutHeight")
  list_list_cutree_params_final <- lapply(list_PPI_vec_n_grey, function(x) comb_list[order(x,decreasing = F)][[1]])
  list_list_cutree_params_final <- lapply(list_list_cutree_params_final, function(x) name_for_vec(to_be_named = x, given_names = param_names, dimension = NULL)) 
  
  ######################################################################
  ################### DELETE RUNS WITH ONLY GREY #######################
  ######################################################################
  
  PPI_vec_n_grey <- count_grey_in_list_of_vec(list_colors_PPI)
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  logical_subsets_PPI_ok <- as.logical(mapply(function(x,y) x != length(y), x = PPI_vec_n_grey, y = list_colors_PPI, SIMPLIFY = F))
  
  sNames_PPI <- sNames_ok[logical_subsets_PPI_ok]
  list_colors_PPI <- list_colors_PPI[logical_subsets_PPI_ok]
  list_module_PPI_ok <- list_module_PPI[logical_subsets_PPI_ok]
  list_module_PPI_signif_ok <- list_module_PPI_signif[logical_subsets_PPI_ok]
  list_datExpr_PPI <- list_datExpr_ok[logical_subsets_PPI_ok]
  
  list_colors <- list_colors[logical_subsets_PPI_ok]
  list_list_cutree_params_final_PPI <- list_list_cutree_params_final[logical_subsets_PPI_ok]
  list_plot_label_final_PPI  <- list_plot_label_final[logical_subsets_PPI_ok]
  list_geneTree_PPI <- list_geneTree_ok[logical_subsets_PPI_ok]
  
  # Get proportion of assigned genes
  diagnostic_stats$prop_genes_assign[sNames %in% sNames_ok] <- sapply(list_colors, function(x) round(sum(x!="grey")/length(x),2), simplify=T)
  diagnostic_stats$prop_genes_assign_PPI[sNames %in% sNames_PPI] <- sapply(list_colors_PPI, function(x) round(sum(x!="grey")/length(x),2), simplify=T)
  
  # Count module (minus 1 for grey)
  diagnostic_stats$n_modules[sNames %in% sNames_ok] <- sapply(list_colors, function(x) length(unique(as.character(x)))-1, simplify=T) 
  diagnostic_stats$n_modules_PPI[sNames %in% sNames_PPI] <- sapply(list_colors_PPI, function(x) length(unique(as.character(x)))-1, simplify=T) 
  
  ######################################################################
  ##################### MAKE MODULE COLORS UNIQUE ######################
  ######################################################################
  
  message("Making module colors unique across cell clusters")
  
  # Get nested list of modules
  list_mods <- lapply(list_colors_PPI, function(x) names(table(x))[-grep("^grey$", names(table(x)))])
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
    fakeCols <- paste0(all_cols_nogrey_uniq, "_", 1:(length(mods) - length(all_cols_nogrey)))
    mods_uniq <- mods
    mods_uniq[1:length(all_cols_nogrey)] <- all_cols_nogrey[sample(x=1:length(all_cols_nogrey), size=length(all_cols_nogrey), replace=F)]
    mods_uniq[(length(all_cols_nogrey)+1):length(mods_uniq)] <- fakeCols[sample(x=1:length(fakeCols), size=length(mods_uniq)-length(all_cols_nogrey), replace=F)]
  }
  
  # Nest the modules by celltype
  k = 0
  list_mods_uniq <- vector(mode="list", length=length(list_mods))
  for (i in 1:length(list_mods)) {
    list_mods_uniq[[i]] <- mods_uniq[(k+1):(k+length(list_mods[[i]]))]
    k <- k+length(list_mods[[i]])
  }
  
  names(list_mods_uniq) <- names(list_colors_PPI)
  
  # Update color vectors
  # TODO: review
  
  list_colors_PPI_uniq <- mapply(function(x,y,z) z[match(x, y)],
                                 x = list_colors_PPI, 
                                 y = list_mods, 
                                 z = list_mods_uniq, SIMPLIFY = F)
  
  # Update also the old module assignment colors, but only change those that were validated by PPI (and therefore among those replaced with unique colors)
  list_colors_uniq <-  mapply(function(x,y,z) ifelse(x==y, yes=z, no = x),
                              x = list_colors, 
                              y = list_colors_PPI, 
                              z = list_colors_PPI_uniq,
                              SIMPLIFY = F)
  
  # Replace NAs (i.e. no match in non-grey modules) with "grey"
  list_colors_uniq <- lapply(list_colors_uniq, function(x) ifelse(is.na(x), yes="grey", no = x))
  list_colors_PPI_uniq <- lapply(list_colors_PPI_uniq, function(x) ifelse(is.na(x), yes="grey", no = x))
  
  # Give gene names to vectors
  list_colors_uniq <- mapply(function(x,y) name_for_vec(to_be_named= x, given_names = names(y), dimension=NULL), x = list_colors_uniq, y = list_colors, SIMPLIFY=F)
  list_colors_PPI_uniq <- mapply(function(x,y) name_for_vec(to_be_named= x, given_names = names(y), dimension=NULL), x = list_colors_PPI_uniq, y = list_colors_PPI, SIMPLIFY=F)
  
  ######################################################################
  ########## RECOMPUTE MEs, kMs, pkMs AFTER PPI FILTER #########
  ######################################################################
  # in:
  #   list_colors_PPI_ok
  #   list_datExpr_PPI_ok
  #
  # out:
  #   list_Ms_PPI
  #   list_kMs_PPI
  
  # TODO: Surely we can just remove the grey ones?
  
  if (checkPPI==T) message("Computing kMs after filtering modules for significant Protein-Protein Interactions")
  
  if (fuzzyModMembership=="kIM" | scale_MEs_by_kIMs == T) {
  list_dissTOM_ok_path <- dir(path = scratch_dir, pattern = paste0(data_prefix, "_list_dissTOM_ok"), full.names = T)
  list_dissTOM_ok <- load_obj(list_dissTOM_ok_path)
  names(list_dissTOM_ok) <- sNames_ok
  list_dissTOM_PPI <- list_dissTOM_ok[sNames_ok %in% sNames_PPI]
  }
  
  if (fuzzyModMembership == "kME") {
    
    invisible(gc())
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_kMEs_PPI.txt"))
    
    list_MEs_PPI <- clusterMap(cl, function(x,y,z) moduleEigengenes_kIM_scale(expr = as.data.frame(x, col.names=col.names(x)),
                                                                  colors=y,
                                                                  excludeGrey=T,
                                                                  scale_MEs_by_kIMs=scale_MEs_by_kIMs,
                                                                  dissTOM=z), 
                               x = list_datExpr_PPI, 
                               y = list_colors_PPI_uniq,
                               z = list_dissTOM_PPI,
                               SIMPLIFY = F,
                               .scheduling = c("dynamic"))
    
    list_kMs_PPI <- clusterMap(cl, function(x,y) signedKME(as.matrix(x),
                                                           y$eigengenes,
                                                           outputColumnName = "",
                                                           corFnc = corFnc), 
                               x=list_datExpr_PPI, 
                               y=list_MEs_PPI,
                               SIMPLIFY=F,
                               .scheduling = c("dynamic"))
    stopCluster(cl)
    invisible(gc()); invisible(R.utils::gcDLLs())
    
  } else if (fuzzyModMembership == "kIM"){
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_kIMs_PPI.txt"))
    
    list_kMs_PPI <- clusterMap(cl, function(x,y) kIM_eachMod_norm(dissTOM = x, 
                                                                  colors = y,
                                                                  excludeGrey=T),
                               x = list_dissTOM_PPI, 
                               y = list_colors_PPI_uniq, 
                               SIMPLIFY=F, 
                               .scheduling = c("dynamic"))
    
    stopCluster(cl)
    
    names(list_kMs_PPI) <- sNames_PPI
    
    rm(list_dissTOM_PPI, list_dissTOM_ok)
    
    list_MEs_PPI <- NULL
    
    invisible(gc()); invisible(R.utils::gcDLLs())
  } 
  
  ##########################################################################
  ################## COMPUTE RARE VARIANTS ENRICHMENT ######################
  ##########################################################################
  
  message("Checking modules for enrichment with rare variants/mendelian genes")
  
  variants <- load_obj(f=variants_path)
  
  # convert kMs from dataframes to lists
  list_list_kMs_PPI <- lapply(list_kMs_PPI, function(x) as.list(x)) 
  
  # give gene names to each vector of kMs
  list_list_kMs_PPI <- mapply(function(a,b) lapply(a, function(x) name_for_vec(to_be_named = x, given_names = rownames(b), dimension = NULL)),
                              a = list_list_kMs_PPI,
                              b = list_kMs_PPI,
                              SIMPLIFY = F)
  
  names(list_list_kMs_PPI) = sNames_PPI 
  # prepend the celltype to the module name
  # list_list_kMs_PPI <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names=paste0(y, "_", names(x)), dimension=NULL), 
  #                             x=list_list_kMs_PPI, 
  #                             y=sNames_PPI,
  #                             SIMPLIFY = F)
  
  invisible(gc()); invisible(R.utils::gcDLLs())
  cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, data_prefix, "_log_variants_GSEA.txt"))
  
  list_list_variants_GSEA <- parLapplyLB(cl, list_list_kMs_PPI, 
                                         function(x) lapply(x, 
                                                            function(y) gsea(values = y, geneset = variants, plot = F, return.details=T)))
  # Returns a list (of modules) of dataframes with enrichment terms as rownames and columns: p.val, q.val, sscore and edge
  stopCluster(cl)
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  # dissolve the top level celltype list so we have a list of module enrichment vectors
  list_variants_GSEA <- unlist(list_list_variants_GSEA, recursive = F, use.names = T) 
  #variants_GSEA_mat <- matrix(nrow=length(list_variants_GSEA), ncol=6) 
  
  celltype_col <- rep(sNames_PPI, times=sapply(list_list_variants_GSEA, length))
  mods_col <- gsub("^.*\\.", "", names(list_variants_GSEA))
  GSEA_cols <- Reduce( f = rbind, list_variants_GSEA)
  
  variants_GSEA_df  <- data.frame(celltype=celltype_col, module=mods_col, GSEA_cols, row.names=NULL, stringsAsFactors = F)
  # adjust p.vals for multiple testing
  p.adjust(variants_GSEA_df[['p.val']], method = "fdr") %>% -log10(.) -> variants_GSEA_df$q.val 
  
  ##########################################################################
  ############# COMPUTE MAGMA COMMON VARIANTS GWAS ENRICHMENT ##############
  ##########################################################################
  
  # MAGMA docs: http://ctg.cncr.nl/software/MAGMA/doc/manual_v1.06.pdf
  
  if (!is.null(magma_gwas_dir)) {
    
    message("Scoring enriched modules for enrichment with genes linked by GWAS to phenotypes of interest")
    
    if (data_organism == "mmusculus") {
      invisible(gc()); invisible(R.utils::gcDLLs())
      # map mouse to human gene orthologs  
      mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
      
      cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "mapMMtoHs_par.txt"))
      list_kMs_hs <- parLapplyLB(cl, list_kMs_PPI, function(x) mapMMtoHs(modulekM = x, 
                                                                         log_dir = log_dir, 
                                                                         flag_date = flag_date, 
                                                                         data_prefix = data_prefix, 
                                                                         mapping_orthology = mapping_orthology))
      stopCluster(cl)
      invisible(gc()); invisible(R.utils::gcDLLs())
      
      names(list_kMs_hs) <- names(list_kMs_PPI)
      
    } else if (data_organism == "hsapiens") {
      list_kMs_hs <- list_kMs_PPI
      names(list_kMs_hs) <- names(list_kMs_PPI)
    }
    
    # Load gene mapping and annotation files
    
    # columns: ensembl_gene_id, entrezgene, hgnc_symbol. Use this for mapping entrezgene to ensembl
    mapping_hs_entrez2ensembl = read.csv(gzfile(mapping_hs_entrez2ensembl_filepath),sep="\t",header=T, stringsAsFactors = F)
    
    # Load MAGMA genes and remap to Ensembl gene IDs
    
    d = dir(path=magma_gwas_dir, pattern="[.]genes.out", recursive = T)
    gwas = vector(mode="list")
    for(i in 1:length(d)) {
      gwas[[i]] = read.table(paste(magma_gwas_dir, d[[i]],sep=""),head=T, check.names = FALSE, stringsAsFactors = F)
    }
    names(gwas) = gsub(".genes.out", "", d)
    
    # Remap from human Entrez to human Ensembl gene IDs
    for (i in 1:length(gwas)) {
      idx = match(gwas[[i]]$GENE, mapping_hs_entrez2ensembl$entrezgene)
      mapping = data.frame(entrez=gwas[[i]]$GENE, ensembl=mapping_hs_entrez2ensembl$ensembl_gene_id[idx])
      gwas[[i]]$gene_name = mapping$ensembl
    }
    invisible(gc()); invisible(R.utils::gcDLLs())
    cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "kM_magma_par.txt"))
    magma_results <- clusterMap(cl, function(x,y) kM_magma(cellType = x, 
                                                           modulekM = y,
                                                           gwas = gwas),
                                x = names(list_kMs_hs),
                                y = list_kMs_hs,
                                SIMPLIFY=F,
                                .scheduling = c("dynamic"))
    stopCluster(cl)
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    # Prepare tables for results
    list_magma.r <- list_magma.p <- list_magma.emp.p <- vector(mode="list", length = length(list_kMs_hs))
    
    # Extract the coefficient, analytical p-value and empirical p-value dataframes for each module, put them into lists
    
    for (i in 1:length(magma_results)) { 
      list_magma.r[[i]] <- magma_results[[i]][['corrCoef']]
      list_magma.p[[i]] <- magma_results[[i]][['p.val']]
      list_magma.emp.p[[i]] <- magma_results[[i]][['emp.p.val']]
    }
    
    # Merge each list into a single table
    for (t in 1:length(list_magma.p)) {
      if (t == 1) {
        magma.r.all <- list_magma.r[[t]]
        magma.p.all <- list_magma.p[[t]]
        magma.emp.p.all <- list_magma.emp.p[[t]]
      } else {
        magma.r.all <- rbind(magma.r.all, list_magma.r[[t]] )
        magma.p.all <- rbind(magma.p.all, list_magma.p[[t]] )
        magma.emp.p.all <- rbind(magma.emp.p.all, list_magma.emp.p[[t]] )
      }
    }
    # We will plot the analytical and empirical p-values side-by-side later as a diagnostic
    
    # adjust analytical p-values for multiple testing
    magma.p.fdr.all = p.adjust(magma.p.all, method="fdr")
    dim(magma.p.fdr.all) = dim(magma.p.all);  dimnames(magma.p.fdr.all) = dimnames(magma.p.all)
    
    # Transform analytical p-values to retain those associated with positive correlation coefficients only
    magma.p.fdr.signed = as.data.frame(magma.p.fdr.all*sign(magma.r.all))
    magma.p.fdr.signed[magma.p.fdr.signed<0]=1 #Only look for positive enrichment
    magma.p.fdr.log = -log10(magma.p.fdr.signed)
    
    # Convert to dataframes
    magma.p.fdr.log <- magma.p.fdr.log %>% as.data.frame
    magma.p.all <- magma.p.all %>% as.data.frame
    magma.emp.p.all <- magma.emp.p.all %>% as.data.frame
    magma.r.all <- magma.r.all %>% as.data.frame
    
    celltype <- gsub("__.+","",rownames(magma.p.fdr.log))    
    module <- gsub(".+__","",rownames(magma.p.fdr.log))
    
    magma.p.fdr.log <- cbind(celltype, module, magma.p.fdr.log)
    magma.p.all <- cbind(celltype, module, magma.p.all)
    magma.emp.p.all <- cbind(celltype, module, magma.emp.p.all)
    magma.r.all <- cbind(celltype, module, magma.r.all)
    
  } else if (is.null(magma_gwas_dir)) {
    
    magma.r.all <- NULL
    magma.p.all <- NULL
    magma.emp.p.all <- NULL 
    magma.p.fdr.log <- NULL
    
  }
  
  ##########################################################################
  ############## FILTER MODULES ON GENES WITH GWAS ENRICHMENT ##############
  ##########################################################################
  
  # TODO: check if MAGMA needs correcting for gene length
  # TODO: Need to add non-protein coding genes to MAGMA
  
  if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) {
    
    # If filtering out non-gwas enriched modules, also retain modules with rare variant / mendelian gene enrichment 
    idx_row_variants <- variants_GSEA_df$q.val[match(variants_GSEA_df$module, magma.p.fdr.log$module)] > -log10(pvalThreshold)
    
    # Check if the gwas_filter_traits strings provided by the user can be found in the magma_gwas_dir files
    
    if (sapply(gwas_filter_traits, function(x) any(grepl(x, colnames(magma.p.fdr.log), ignore.case=T)), simplify = T) %>% all) {
      #if (sapply(gwas_filter_traits, function(x) any(grepl(x, colnames(magma.p.fdr.log[,-grep("module", colnames(magma.p.fdr.log))]), ignore.case=T)), simplify = T) %>% all) {
      
      # Identify columns to keep
      idx_col_keep <- sapply(gwas_filter_traits, function(x) grep(paste0(x,"|module|celltype"), colnames(magma.p.fdr.log), ignore.case=T), simplify = T) %>% Reduce(union,.)# %>% as.logical 
      
      # fdr-corrected log-transformed analytical p-values
      magma.p.fdr.log.sig <- magma.p.fdr.log[,idx_col_keep]
      
      # correlation coefficients
      magma.r.sig <- magma.r.all[,idx_col_keep]
      
      # Identify rows to keep
      idx_row_gwas <- apply(magma.p.fdr.log.sig[,-grep("module|celltype", colnames(magma.p.fdr.log.sig))], MARGIN = 1, max) > -log10(pvalThreshold)
      
      # rows enriched for gwas common variants and/or for rare variants
      idx_row_keep <- idx_row_gwas | idx_row_variants
      
      if (sum(idx_row_keep)>0) {
        
        magma.p.fdr.log.sig <- magma.p.fdr.log.sig[idx_row_keep,]
        magma.r.sig <- magma.r.sig[idx_row_keep,]
        
      } else if (sum(idx_row_keep)==0) { # Warn if there are no significant enrichments
        
        gwas_filter_traits <- NULL
        warning("No modules enriched for gwas_filter_traits or rare variants - skipping gwas significance filtering step")
        
      }
    } else {
      
      gwas_filter_traits <- NULL 
      message("gwas_filter_traits not found in magma_gwas_dir files - skipping gwas significance filtering step")
      
    }
  }  
  
  if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
    magma.p.fdr.log.sig <- magma.p.fdr.log
    magma.r.sig <- magma.r.all
  }
  
  ######################################################################
  ####### FILTER MODULES ON GWAS SIGNIFICANCE ALSO IN OTHER FILES ######
  ######################################################################
  
  if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) { # NB: if no significant gwas correlations found in previous section, gwas_filter_traits <- NULL
    
    # get cell types with modules with gwas enrichment
    sNames_gwas <- names(table(magma.p.fdr.log.sig$celltype))
    
    list_module_gwas <- vector(mode="list", length=length(sNames_gwas))
    names(list_module_gwas) <- sNames_gwas
    
    # Get a list of vectors with modules per celltype
    for (name in sNames_gwas) {
      idx <- grep(pattern = name, x = magma.p.fdr.log.sig$celltype)
      list_module_gwas[[name]] <- sapply(idx, function(j) magma.p.fdr.log.sig$module[j], simplify = T)
    }
    
    # Filter out cell types with no significant modules
    list_colors_gwas <- list_colors_PPI_uniq[sNames_PPI  %in% sNames_gwas]
    
    # relabel non-GWAS modules 'grey'
    list_colors_gwas <- mapply(function(x,y) ifelse(x %in% y, yes = x, no = "grey"),
                               ### 180611
                               # x = list_colors_gwas,
                               x = list_colors_gwas, 
                               y = list_module_gwas,
                               SIMPLIFY = F)
    
    # give gene names to color assignment vectors
    list_colors_gwas <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(y), dimension = NULL), 
                               x = list_colors_gwas,
                               y = list_colors_PPI_uniq[sNames_PPI %in% sNames_gwas, drop=F],
                               SIMPLIFY = F)
    
    # remove any cell clusters without gwas enrichment: expression matrices
    list_datExpr_gwas <- list_datExpr_PPI[sNames_PPI %in% sNames_gwas, drop=F] 
    
  } else if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
    
    sNames_gwas <- sNames_PPI
    list_module_gwas <- lapply(list_kMs_PPI, function(x) colnames(x)) 
    names(list_module_gwas) <- sNames_gwas
    list_colors_gwas <- list_colors_PPI_uniq
    list_datExpr_gwas <- list_datExpr_PPI
  }
  
  ######################################################################
  #################### RECOMPUTE kMs AFTER MAGMA #######################
  ######################################################################
  
  if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) {
    
    message("Re-computing ", fuzzyModMembership, "s after filtering modules for MAGMA enrichment")
    
    if (fuzzyModMembership == "kIM" | scale_MEs_by_kIMs==T) {
      list_dissTOM_ok_path <- dir(path = scratch_dir, pattern = paste0(data_prefix, "_list_dissTOM_ok"), full.names = T)
      list_dissTOM_ok <- load_obj(list_dissTOM_ok_path)
      names(list_dissTOM_ok) <- sNames_ok
      list_dissTOM_gwas <- list_dissTOM_ok[names(list_dissTOM_ok) %in% sNames_gwas]
      rm(list_dissTOM_ok)  
    }
    
    if (fuzzyModMembership=="kME") {
      
      invisible(gc()); invisible(R.utils::gcDLLs())
      cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_moduleEigengenes_kMEs_magma.txt"))
      
      list_MEs_gwas <- clusterMap(cl, function(x,y,z) moduleEigengenes_kIM_scale(expr=as.data.frame(x,col.names=col.names(x)),
                                                                     colors=y,
                                                                     excludeGrey=T,
                                                                     scale_MEs_by_kIMs=scale_MEs_by_kIMs,
                                                                     dissTOM=z), 
                                  x = list_datExpr_gwas, 
                                  y = list_colors_gwas, 
                                  z = list_dissTOM_gwas,
                                  SIMPLIFY = F,
                                  .scheduling = c("dynamic"))
      
      list_kMs_gwas <- clusterMap(cl, function(x,y) signedKME(as.matrix(x),
                                                              y$eigengenes,
                                                              outputColumnName = "",
                                                              corFnc = corFnc), 
                                  x=list_datExpr_gwas, 
                                  y=list_MEs_gwas,
                                  SIMPLIFY=F,
                                  .scheduling = c("dynamic"))
      
      stopCluster(cl)
      invisible(gc()); invisible(R.utils::gcDLLs())
      
    } else if (fuzzyModMembership=="kIM") {
      

      invisible(gc()); invisible(R.utils::gcDLLs())
      
      cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_kIMs_magma.txt"))
      
      list_kMs_gwas <- clusterMap(cl, function(x,y) kIM_eachMod_norm(dissTOM = x, 
                                                                     colors = y),
                                  x = list_dissTOM_gwas, 
                                  y = list_colors_gwas, 
                                  SIMPLIFY=F,
                                  .scheduling = c("dynamic"))
      
      stopCluster(cl)
      invisible(gc()); invisible(R.utils::gcDLLs())
      names(list_kMs_gwas) <- sNames_gwas
      
      rm(list_dissTOM_gwas)
      
      list_MEs_gwas <- NULL
      
      invisible(gc())
      
    }
    
    # Remove grey modules (if any and filter out runs with only grey (probably not possible but still))
    list_MEs_gwas <- deleteGrey(list_MEs_gwas) %>% Filter(f=length)
    list_kMs_gwas <- deleteGrey(list_kMs_gwas) %>% Filter(f=length)
    
  } else if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
    
    list_MEs_gwas <- list_MEs_PPI
    list_kMs_gwas <- list_kMs_PPI
    
  }
  
  
  ######################################################################
  ####### COMPUTE MODULE - METADATA CORRELATION IN EACH CELL CLUSTER ###
  ######################################################################
  

  if (!is.null(metadata_corr_col)) {
    
    if (!is.null(metadata)) {
      
      message("Computing module-metadata correlation in each celltype")
      
      list_metadata <- lapply(list_datExpr_gwas, function(x) metadata[match(rownames(x), rownames(metadata)), , drop=F]) # get list of cell * metadata
      
      if (fuzzyModMembership=="kME") {
        # Compute correlation between metadata (columns) and eigengenes (columns). 
        list_mod_metadata_corr_rho <- mapply(function(x,y) cor(x=as.matrix(x), 
                                                               y=as.matrix(y$eigengenes), 
                                                               method = c("pearson"), 
                                                               use = 'pairwise.complete.obs'), 
                                             x=list_metadata, 
                                             y=list_MEs_gwas,  
                                             SIMPLIFY = F) 
        
        
        
        # Remove 'ME' from eigengene names
        for (j in 1:length(list_mod_metadata_corr_rho)) {
          colnames(list_mod_metadata_corr_rho[[j]]) <- gsub("^ME", "", colnames(list_mod_metadata_corr_rho[[j]]), ignore.case=F)
        }
        
      } else if (fuzzyModMembership == "kIM") {
        
        # compute the equivalent to eigengenes (i.e. embeddings) but using kIMs as eigenvectors on which to project each cell's gene-length vector
        #cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "list_cell_IM_embed_mat_for_meta_corr.txt"))
        list_embed_mat <- mapply(function(a,b,c,d) cellModEmbed(datExpr=a, 
                                                     colors=b, 
                                                     latentGeneType="IM",
                                                     cellType=c,
                                                     kMs=d),
                                     a = list_datExpr_gwas,
                                     b = list_colors_gwas,
                                     c = names(list_datExpr_gwas),
                                     d = list_kMs_gwas,
                                     SIMPLIFY=F)
        #stopCluster(cl)
        
        # Get correlations
        list_mod_metadata_corr_rho <- mapply(function(x,y) cor(x=as.matrix(x), 
                                                               y=as.matrix(y), 
                                                               method = c("pearson"), 
                                                               use = 'pairwise.complete.obs'), 
                                             x=list_metadata, 
                                             y=list_embed_mat,  
                                             SIMPLIFY = F) 
        
        # Clear space
        rm(list_embed_mat)
        
        list_mod_metadata_corr_rho <- lapply(list_mod_metadata_corr_rho, function(x) name_for_vec(to_be_named = x, given_names = colnames(metadata), dimension = 1)) 
        
      }
      
      # Name the rows of the correlation matrix with the metadata columns
      list_mod_metadata_corr_rho <- lapply(list_mod_metadata_corr_rho, function(x) name_for_vec(to_be_named = x, given_names = colnames(metadata), dimension = 1)) 
      
      # Compute p values
      list_mod_metadata_corr_pval <- mapply(function(x,y) WGCNA::corPvalueStudent(x, 
                                                                                  n = nrow(y)), 
                                            x=list_mod_metadata_corr_rho, 
                                            y=list_metadata, 
                                            SIMPLIFY=F)
      
      # Convert p values into one big matrix in order to adjust p-values for the number of modules tested across *all* celltypes
      for (j in 1:length(list_mod_metadata_corr_pval)) {
        if (j==1) {
          corr_pval <- as.data.frame(list_mod_metadata_corr_pval[[j]])
          colnames(corr_pval) <- paste0(names(list_mod_metadata_corr_pval)[j], "__", colnames(corr_pval)) 
        } else {
          tmp <- list_mod_metadata_corr_pval[[j]]
          colnames(tmp) <- paste0(names(list_mod_metadata_corr_pval)[j], "__", colnames(tmp)) 
          corr_pval <- cbind(corr_pval, tmp)
        }
      }
      
      # set NAs to 1
      corr_pval[is.na(corr_pval)] <- 1
      
      # Compute the false discovery rates
      corr_fdr <- t(apply(t(corr_pval), MARGIN=2, FUN = function(x) p.adjust(x, method = "fdr")))
      corr_fdr.log <- -log10(corr_fdr)     
      
      # Split single dataframe back into a list of celltypes
      list_mod_metadata_corr_fdr <- list_mod_metadata_corr_fdr.log <- vector(mode = "list", length = length(list_mod_metadata_corr_pval))
      names(list_mod_metadata_corr_fdr) <- names(list_mod_metadata_corr_fdr.log) <- names(list_mod_metadata_corr_rho)
      
      k = 0
      for (j in 1:length(list_mod_metadata_corr_fdr.log)) {
        list_mod_metadata_corr_fdr[[j]] <- as.data.frame(corr_fdr[,(k+1):(k+ncol(list_mod_metadata_corr_pval[[j]])), drop=F])
        list_mod_metadata_corr_fdr.log[[j]] <- as.data.frame(corr_fdr.log[,(k+1):(k+ncol(list_mod_metadata_corr_pval[[j]])), drop=F])
        k <- k + ncol(list_mod_metadata_corr_pval[[j]])
      }
      
      # Also make a single dataframe for saving 
      dim(corr_fdr.log) <- c(ncol(metadata), length(unlist(list_module_gwas)))
      rownames(corr_fdr.log) <- rownames(corr_pval)
      colnames(corr_fdr.log) <- colnames(corr_pval)
      corr_fdr.log <- as.data.frame(corr_fdr.log)
      
    } 
    
  } else {
    metadata = NULL
  }
  
  ##########################################################################
  #### FILTER COLORS VECS, GENE AND KME LISTS FOR METADATA CORRELATIONS ####
  ##########################################################################
  
  if (!is.null(metadata_corr_col) & !is.null(metadata_corr_filter_vals)) {
    
    if (!is.null(metadata)) {
      
      # Get a list of logical vectors indicating significant correlations
      list_idx_module_sig <- lapply(list_mod_metadata_corr_fdr.log, function(x) apply(x, 2, function(y) any(y>-log10(pvalThreshold)))) # logical
      
      # Only keep significantly correlated modules within each celltype
      list_module_meta <- mapply(function(x,y) colnames(x[,y, drop=F]),x=list_mod_metadata_corr_fdr.log, y=list_idx_module_sig, SIMPLIFY=F) 
      list_module_meta <- lapply(list_module_meta, function(x) gsub("^.*__", "", x)) 
      
      sNames_meta <- sNames_gwas[sapply(list_idx_module_sig, any, simplify = T)]
      
      # Filter out celltypes with no modules significantly correlated with metadata
      list_module_meta <- list_module_meta[sNames_gwas %in% sNames_meta]
      
      # reassign genes of filtered out modules to grey and remove any empty cell clusters
      list_colors_meta <- mapply(function(x,y) ifelse(x %in% y, yes = x, no = "grey"),
                                 x = list_colors_gwas[names(list_colors_gwas) %in% sNames_meta, drop=F],
                                 y = list_module_meta,
                                 SIMPLIFY = F)
      
      # give gene names to color assignment vectors
      list_colors_meta <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(y), dimension = NULL), 
                                 x = list_colors_meta,
                                 y = list_colors_gwas[names(list_colors_gwas) %in% sNames_meta, drop=F],
                                 SIMPLIFY = F)
      
      # Filter kM lists
      list_kMs_meta <- list_kMs_gwas[names(list_kMs_gwas) %in% sNames_meta, drop = F]
      list_kMs_meta <- mapply(function(x,y) x[,colnames(x)%in%y],x=list_kMs_meta, y = list_module_meta) %>% Filter(f=length)
      
      
    } else 
      
      sNames_meta <- sNames_gwas
    list_colors_meta <- list_colors_gwas
    list_module_meta <- list_module_gwas
    list_kMs_meta <- list_kMs_gwas
    
  } else if (is.null(metadata_corr_col) | is.null(metadata_corr_filter_vals)) {
    
    sNames_meta <- sNames_gwas
    list_colors_meta <- list_colors_gwas
    list_module_meta <- list_module_gwas
    list_kMs_meta <- list_kMs_gwas 
    # if (!is.null(magma_gwas_dir)) {
    #   magma.p.fdr.log.sig_meta <- magma.p.fdr.log.sig
    #   magma.r.sig_meta <- magma.r.sig
    # }
  }
  
  ######################################################################
  ##################### COMPUTE CELL x EIGENGENE MATRIX ################
  ######################################################################
  
  message("Computing all cell embeddings on all modules, across celltypes")
  
  scale_data_path <- dir(path = scratch_dir, pattern = paste0(data_prefix, "_scale_regr_data_ensembl"), full.names = T)
  load_obj(scale_data_path) %>% t -> datExpr 
  
  invisible(gc())
  
  if (fuzzyModMembership == "kIM" | scale_MEs_by_kIMs==T) {
    list_dissTOM_ok_path <- dir(path = scratch_dir, pattern = paste0(data_prefix, "_list_dissTOM_ok"), full.names = T)
    list_dissTOM_ok <- load_obj(list_dissTOM_ok_path)
    names(list_dissTOM_ok) <- sNames_ok
    list_dissTOM_meta <- list_dissTOM_ok[names(list_dissTOM_ok) %in% sNames_meta]
    rm(list_dissTOM_ok)  
  }
  
  latentGeneType <- if (fuzzyModMembership == "kME") "ME" else if (fuzzyModMembership=="kIM") "IM"
  
  cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "make_cell_embed_mat.txt"))
  
  list_cellModEmbed_mat <- clusterMap(cl, function(x,y,z,a) cellModEmbed(datExpr = datExpr, 
                                                                       colors = x,
                                                                       latentGeneType = latentGeneType,
                                                                       cellType = y,
                                                                       kMs = if (latentGeneType== "IM") z else NULL,
                                                                       scale_MEs_by_kIMs = scale_MEs_by_kIMs,
                                                                       dissTOM = if (scale_MEs_by_kIMs==T) a else NULL), 
                                      x = list_colors_meta, 
                                      y = names(list_colors_meta),
                                      z = list_kMs_meta,
                                      a = if (scale_MEs_by_kIMs==T) list_dissTOM_meta else numeric(length=length(sNames_meta)),
                                      SIMPLIFY=F,
                                      .scheduling = c("dynamic"))
  
  stopCluster(cl)
  
  invisible(gc()); invisible(R.utils::gcDLLs())
  
  list_cellModEmbed_mat %>% Reduce(function(mat1, mat2) cbind(mat1, mat2), .) -> cellModEmbed_mat
  
  rownames(cellModEmbed_mat) <- rownames(datExpr)
  
  rm(list_dissTOM_meta)
  ##########################################################################
  ######### PREPARE GENES LISTS AND DATAFRAME WITH MODULES, GENES ##########
  ##########################################################################
  
  # prepare nested lists of module genes
  list_list_module_meta_genes <- mapply(function(a,b) lapply(b, function(x) names(a)[a==x]), a=list_colors_meta, b=list_module_meta, SIMPLIFY=F)
  list_list_module_meta_genes <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = NULL), x=list_list_module_meta_genes, y=list_module_meta, SIMPLIFY=F)
  
  # order the gene lists by kME
  tmp <- list_list_module_meta_genes
  
  # iterate over celltypes
  for (i in 1:length(list_list_module_meta_genes)) {
    # iterate over modules
    for (j in 1:length(list_list_module_meta_genes[[i]])) {
      kMs <- list_kMs_meta[[i]] 
      genes <- list_list_module_meta_genes[[i]][[j]]
      tmp[[i]][[j]] <- genes[order(kMs[rownames(kMs) %in% genes, j], decreasing=T)]
    }
  }
  
  list_list_module_meta_genes <- tmp
  
  # Prepare module genes dataframe
  cell_cluster <- rep(sNames_meta, times=unlist(sapply(list_list_module_meta_genes, FUN=function(x) sum(sapply(x, function(y) length(y), simplify=T)), simplify=T)))
  module <- unlist(sapply(list_list_module_meta_genes, function(x) rep(names(x), sapply(x, function(y) length(y), simplify = T)), simplify=T), use.names = F)
  ensembl <- unlist(list_list_module_meta_genes, recursive = T, use.names = F)
  
  # if we mapped from hgnc to ensembl, retrieve the hgnc symbols
  if (file.exists(sprintf("%s%s_%s_hgnc_to_ensembl_mapping_df.RData", tables_dir, data_prefix, data_organism))) {
    mapping <- read.csv(file=sprintf("%s%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, data_organism), stringsAsFactors = F)
    list_list_module_meta_genes_hgnc <- lapply(list_list_module_meta_genes, function(x) lapply(x, function(y) mapping$symbol[match(y, mapping$ensembl)]))
    list_list_module_meta_genes_hgnc <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = NULL), x=list_list_module_meta_genes_hgnc, y=list_module_meta)
    hgnc <- unlist(list_list_module_meta_genes_hgnc, recursive = T, use.names=F)
    df_meta_module_genes <- data.frame(cell_cluster, module, ensembl, hgnc, row.names = NULL)
  } else {
    df_meta_module_genes <- data.frame(cell_cluster, module, ensembl, row.names = NULL)
  }
  
  ##########################################################################
  ############################# OUTPUT TABLES ##############################
  ##########################################################################
  
  ##################### SAVE LOGS OF REASSIGNED GENES #######################
  ###########################################################################
  
  
  if (!is.null(kM_reassign)) invisible(mapply(function(x,y) write.csv(x, file =  paste0(RObjects_dir, "_", data_prefix, "_", y, "_list_genes_reassigned.csv"), quote = F, row.names = F), 
                   x = list_reassign_log, 
                   y = names(list_reassign_log), 
                   SIMPLIFY=F))
  
  ##################### WRITE OUT TABLES OF MODULE GENES ####################
  ###########################################################################
  
  invisible(write.csv(df_meta_module_genes, file = sprintf("%s%s_cell_cluster_module_genes_%s.csv",tables_dir, data_prefix, flag_date), row.names = F, quote = F))
  
  # also output ordered gene lists
  for (i in names(list_list_module_meta_genes)) {
    for (j in names(list_list_module_meta_genes[[i]])) {
      invisible(write.csv(list_list_module_meta_genes[[i]][[j]], file = sprintf("%s%s_%s_%s_module_genes.csv",tables_dir, data_prefix, i, j), row.names = F, quote = F))
    }
  }
  
  save(list_list_module_meta_genes, file=sprintf("%s%s_list_list_module_genes.RData", RObjects_dir, data_prefix))
  
  ####################### SAVE STRINGDB PPI OUTPUT #########################
  ##########################################################################
  
  # convert the module PPI dataframe columns from list to numeric 
  
  list_module_PPI_ok <- lapply(list_module_PPI_ok, function(x) {
    out = sapply(X = x, FUN = as.numeric, simplify = T)
    rownames(out) = rownames(x)
    return(out)
    })
  
  list_module_PPI_signif_ok <- lapply(list_module_PPI_signif_ok, function(x) {
    out = sapply(X = x, FUN = as.numeric, simplify = T)
    rownames(out) = rownames(x)
    return(out)
  })
  
  # output as PPI
  
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_STRINGdb_output_all_%s.csv", tables_dir, data_prefix, y,  flag_date), row.names=T, quote = F), 
                   x=list_module_PPI_ok, 
                   y= sNames_PPI, 
                   SIMPLIFY = F))
  
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_STRINGdb_output_signif_%s.csv", tables_dir, data_prefix, y,  flag_date), row.names=T, quote = F), 
                   x=list_module_PPI_signif_ok, 
                   y= sNames_PPI, 
                   SIMPLIFY = F))
  
  ################## SAVE KMs FOR BMI-ENRICHED MODULES ####################
  ##########################################################################
  
  # Prepare dfs with a gene column followed by kMEs / kIMs 
  list_kMs_meta_out <- lapply(list_kMs_meta, function(x) cbind(genes=rownames(x), x))
  rownames(list_kMs_meta) <- NULL
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_%ss_meta_%s.csv", tables_dir, data_prefix, y, fuzzyModMembership, flag_date), row.names=F, quote = F), list_kMs_meta_out, sNames_meta, SIMPLIFY = F))
  
  ############### OUTPUT MENDELIAN/RARE VARIANT RESULTS #####################
  ###########################################################################
  
  variants_GSEA_df <- arrange(variants_GSEA_df, desc(scaled.score))  
  
  # Full set of results
  invisible(write.csv(variants_GSEA_df, file=sprintf("%s%s_mendelian_rareVariants_GSEA_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
  
  ######################## OUTPUT MAGMA RESULTS #############################
  ###########################################################################
  
  if (!is.null(magma_gwas_dir)) {
    invisible(write.csv(magma.p.fdr.log.sig, file=sprintf("%s%s_magma.fdr.log.sig_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
    invisible(write.csv(magma.r.sig, file=sprintf("%s%s_magma.r.all_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
    invisible(write.csv(magma.p.all, file=sprintf("%s%s_magma.p.all_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
    invisible(write.csv(magma.emp.p.all, file=sprintf("%s%s_magma.emp.p.all_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
  }
  
  ####################### OUTPUT METADATA CORRELATION RESULTS ###############
  ###########################################################################
  
  if (!is.null(metadata_corr_col) & !is.null(metadata)) {
    invisible(write.csv(corr_fdr.log, file=sprintf("%s%s_all_metadata_corr_logfdr_%s.csv", tables_dir, data_prefix, flag_date), row.names=T, quote = F))
    invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_metadata_corr_rho_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=T, quote = F), list_mod_metadata_corr_rho, sNames_meta, SIMPLIFY = F))
    invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_metadata_corr_logfdr_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=T, quote = F), list_mod_metadata_corr_fdr.log, sNames_meta, SIMPLIFY = F))
  }
  ################## OUTPUT CELL MODULE EMBEDDINGS MATRIX ###################
  ###########################################################################
  
  try(invisible(write.csv(cellModEmbed_mat, file=sprintf("%s%s_%s_cellModEmbed_%s.csv", tables_dir, data_prefix, fuzzyModMembership, flag_date), row.names=T, quote = F)))
  
  ######################## SAVE DIAGNOSTIC STATS ############################
  ###########################################################################
  
  write.csv(diagnostic_stats, file=sprintf("%s%s_INFO_diagnostic_stats_%s.csv", log_dir, data_prefix, flag_date), quote = F, row.names=F)
  
}
######################################################################
########################### CLEAR UP #################################
######################################################################

# Delete consensus TOMs from disk
for (subsetName in sNames) {
  if (file.exists(sprintf("%s%s_consensusTOM-block.1.RData", scratch_dir, subsetName))) {
    file.remove(sprintf("%s%s_consensusTOM-block.1.RData", scratch_dir, subsetName))
  }
}

message("Saving final session image")

save.image(file=sprintf("%s%s_final_session_image.RData", RObjects_dir, data_prefix, flag_date))

rm(list=ls())
invisible(gc())

##########################################################################
################################ FINISH ##################################
##########################################################################

message("Script DONE!")
  