## Script to filter modules from WGCNA runs 

# min/max Ngenes in module
# PPI enrichment
# 'Feature selection' of kME / clustering: hierarchical clustering vs caret feature selection

######################################################################
########################## FUNCTIONS #################################
######################################################################

source(file="/projects/jonatan/tools/functions-src/utility_functions.R")

######################################################################
########################### OPTPARSE #################################
######################################################################

ipak(c("optparse"))

option_list <- list(
  make_option("--path_df_NWA", type="character",
              help = "path to dataframe from a gene network analysis run, in long format (i.e. one row per gene per celltype), containing 'cell_cluster', 'module', one or two gene name columns, and a column of numeric scores. I.e. one row per gene, e.g. rwgcna cell_cluster_module_genes.csv files"),  
  make_option("--colMod", type="character", default="module",
              help ="string or regex for grepping dfNWA e.g. 'module' or 'module_merged', [default %default]"),
  make_option("--colCellCluster", type="character", default="cell_cluster",
              help ="string or regex for grepping nwa_df e.g. 'cell_cluster' or 'cell_cluster_merged', [default %default]"),
  make_option("--colGeneWeights", type="character", default ="pk.*M",
              help = "dfNWA column with gene weights, [default %default]"),  
  # make_option("--colGeneNames", type="character", default="hgnc",
  #             help ="string or regex for grepping dfNWA e.g. 'hgnc|symbol|gene_name' or 'ensembl', [default %default]"),
  make_option("--path_df_PPI", type="character", 
              help = "Path to STRINGdb results dataframe"),  
  make_option("--originCellClusters", type="character", default=NULL,
              help="Optionally, a quoted list of origin cell clusters for filtering modueles, defaults to using all [default %default]."),
  make_option("--minGeneClusterSize", type="integer", default=5L,
              help="What is the minimum number of genes in a module to continue? Integer, [default %default]."),
  make_option("--maxGeneClusterSize", type="integer", default=1000L,
              help="What is the max number of genes in a module to continue? Integer, [default %default]."),
  make_option("--dirOut", type="character",
              help = "Outputs go to /tables and /RObjects subdirectories"),  
  make_option("--prefixOut", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"),
  make_option("--RAMGbMax", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("dplyr", "Matrix", "parallel", "readr"))

######################################################################
############################ PARAMS ##################################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

path_df_NWA <- opt$path_df_NWA 
colMod <- opt$colMod
colCellCluster <- opt$colCellCluster
colGeneWeights <- opt$colGeneWeights
#colGeneNames <- opt$colGeneNames
path_df_PPI <- opt$path_df_PPI
originCellClusters <- opt$originCellClusters
if (!is.null(originCellClusters)) originCellClusters <- eval(parse(text=originCellClusters))
minGeneClusterSize <- opt$minGeneClusterSize
maxGeneClusterSize <- opt$maxGeneClusterSize
dirOut <- opt$dirOut
prefixOut <- opt$prefixOut
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
########################### LOAD DATA ################################
######################################################################

message("Loading gene module data")

df_geneModule <- load_obj(path_df_NWA)

#colGeneNames <- grep(colGeneNames, colnames(df_geneModule), value=T)
colGeneWeights <- grep(colGeneWeights, colnames(df_geneModule), value=T)#colnames(df_geneModule)[which(sapply(X=colnames(df_geneModule), FUN =function(colname) class(df_geneModule[[colname]]))=="numeric")]

######################################################################
######################## 0. initialise  ##############################
######################################################################

modules <- df_geneModule[[colMod]]  %>% unique %>% sort
df_modFilter <- data.frame("module"=modules, row.names=NULL)

######################################################################
######################## 1. FILTER ON CELL CLUSTER ###################
######################################################################

if (!is.null(originCellClusters)) { 
  df_modFilter[["cellClusterOK"]] <- modules %in% df_geneModule[[colMod]][df_geneModule[[colCellCluster]] %in% originCellClusters]
} else { 
  df_modFilter[["cellClusterOK"]] <-  rep(T, length(modules)) 
}  

######################################################################
######################## 2. FILTER ON SIZE ###########################
######################################################################

fun = function(module){
  sum(df_geneModule[[colMod]]==module)>=minGeneClusterSize & sum(df_geneModule[[colMod]]==module)<=maxGeneClusterSize
}

args = list("X"=modules)

logical_nGeneOK <- safeParallel(fun=fun, args=args, simplify=T)

df_modFilter[["nGeneOK"]] <- logical_nGeneOK
  
######################################################################
################# 3. FILTER ON STRINGdb PPI enrichment ###############
######################################################################

STRINGdbOuts <- load_obj(path_df_PPI)

fun <- function(module) {
  STRINGdbOuts %>% dplyr::filter(colors==module) %>% '[['("q.value") < 0.05 
}

args = list("X"=modules)

logical_modulesSTRINGenriched <- safeParallel(fun=fun, args=args, simplify=T)
  
df_modFilter[["STRINGdbOK"]] <- logical_modulesSTRINGenriched

######################################################################
############################## OUTPUTS ###############################
######################################################################

df_modFilter[["filterOK"]] <- with(data = df_modFilter, 
                                      expr = cellClusterOK & nGeneOK & STRINGdbOK)

write_csv(x=df_modFilter, path = paste0(dirTables, prefixOut, "_modFilter.csv"))

message("Script done!")