## script to merge gene modules from different WGCNA runs based on the overlap in genes

######################################################################
########################## FUNCTIONS #################################
######################################################################

source(file="/projects/jonatan/tools/functions-src/utility_functions.R")

######################################################################
########################### OPTPARSE #################################
######################################################################

ipak(c("optparse"))

option_list <- list(
  make_option("--df_pathNWA", type="character",
              help = "path to dataframe from a gene network analysis run, in long format (i.e. one row per gene per celltype), containing 'cell_cluster', 'module', one or two gene name columns, and a column of numeric scores. I.e. one row per gene, e.g. rwgcna cell_cluster_module_genes.csv files"),  
  make_option("--geneWeightsCol", type="character", default ="pk.*M",
              help = "nwa_df column with gene weights"),  
  make_option("--geneCol", type="character", default="hgnc",
              help ="string or regex for grepping nwa_df e.g. 'hgnc|symbol|gene_name' or 'ensembl', [default %default]"),
  make_option("--propIntersectThreshold", type="character", default = "pearson",
              help = "pearson, spearman or kendall"),
  make_option("--corMethod", type="character", default = "pearson",
              help = "pearson, spearman or kendall"),
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

ipak(c("dplyr", "ggplot2", "Biobase", "Matrix", "parallel", "WGCNA", "corrplot"))

######################################################################
############################ PARAMS ##################################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

df_pathNWA <- opt$df_pathNWA 
geneWeightsCol <- opt$geneWeightsCol
geneCol <- opt$geneCol
dirOut <- opt$dirOut
prefixOut <- opt$prefixOut
corMethod <- opt$corMethod
propIntersectThreshold <- opt$propIntersectThreshold
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
######################## LOAD MODULE DATA ############################
######################################################################

message("Loading gene module data")

df_geneModule <- load_obj(df_pathNWA)

geneCol <- grep(geneCol, colnames(df_geneModule), value=T)
geneWeightsCol <- grep(geneWeightsCol, colnames(df_geneModule), value=T)#colnames(df_geneModule)[which(sapply(X=colnames(df_geneModule), FUN =function(colname) class(df_geneModule[[colname]]))=="numeric")]
######################################################################
###################### INITIALISE MERGE LOOP #########################
######################################################################

df_geneModule[["module_merged"]] <- df_geneModule[["module"]]
df_geneModule[[paste0(geneWeightsCol, "_merged")]] <- df_geneModule[[geneWeightsCol]] 
k=1

while (T) {
  
  ######################################################################
  ######################## GET GENE WEIGHT VECTORS #####################
  ######################################################################
  # get a vector of gene weights, with gene as name, for each module
  
  modules <- sort(unique(df_geneModule[["module"]]))
  fun <- function(module) {
    df_geneModule[[geneWeightsCol]] %>% '['(df_geneModule[["module_merged"]]==module) %>% na.omit %>% as.numeric ->geneWeights 
    df_geneModule[[geneCol]] %>% '['(df_geneModule[["module_merged"]]==module) %>% na.omit %>% as.character -> names(geneWeights) 
    return(geneWeights)
  }
  
  args = list("X"=modules)
  list_pkIMs <- safeParallel(fun=fun,args=args)
  
  ######################################################################
  ####################### Compute the overlap matrix ###################
  ######################################################################
  
  # module * module matrix
  # for each column j, each row i gives the proportion of module j that overlaps with module i 
  
  fun = function(pkIMs) {
    sapply(list_pkIMs, function(pkIMsOther) {
      base::intersect(x=names(pkIMs), y=names(pkIMsOther)) %>% length %>% '/'(length(pkIMs))
    }, simplify=T)
  }
  
  args = list("X"=list_pkIMs)
  mat_moduleGeneIntersect <- safeParallel(fun=fun,args=args, simplify = T)
  
  ######################################################################
  ####################### Plot the overlap matrix ######################
  ######################################################################
  
  pdf(sprintf("%s%s_iter%s_mod_intersectpdf", dirPlots, prefixOut, k), 
      width=max(20,ncol(mat_moduleGeneIntersect) %/% 3),
      height=max(20,ncol(module_gene_intersect) %/% 3))
  corrplot(corr = mat_moduleGeneIntersect, 
           method = "color",
           col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")), bias = 1)(200),
           diag = F,
           is.corr=F,
           title=paste0(prefixOut, ": prop. of mod j intersecting mod i, iteration k"),
           order = "hclust",
           hclust.method = hclustMethod,
           addCoef.col = "black",
           tl.srt = 45,
           number.digits = 2L,
           number.cex = 0.5)
  #tl.cex = 2,
  #cl.cex = 1)
  
  invisible(dev.off())
  
  ######################################################################
  ########### Compute the overlap weighted correlation matrix ##########
  ######################################################################
  
  # the outer loop goes in the columns, the inner loop goes in the rows

  fun = function(pkIMs) {
    sapply(list_pkIMs, function(pkIMsother) {
      genesintersect <- base::intersect(names(pkIMs), names(pkIMsother))
      mod_mod_intersect_prop <- sum(names(pkIMs) %in% genesintersect)/ length(pkIMs) 
      if (mod_mod_intersect_prop>propIntersectThreshold) {
        mod_mod_intersect_corr <-   WGCNA::cor(x = pkIMs[match(genesintersect, names(pkIMs))], 
                                               weights.x = pkIMs[match(genesintersect, names(pkIMs))],
                                               y= pkIMsother[match(genesintersect, names(pkIMsother))],
                                               weights.y = pkIMsother[match(genesintersect, names(pkIMsother))],
                                               method = corMethod, 
                                               quick=0.25, 
                                               verbose=T, 
                                               nThreads=30)
        return(mod_mod_intersect_corr)
        
      } else {
        return(0)
      }
    }, simplify=T)
  }
  
  args = list("X"=list_pkIMs)
  
  mat_modIntersectCorrWeighted <- safeParallel(fun=fun,args=args, simplify = T)

  
  ######################################################################
  #########################  ##########################
  ######################################################################
  

  ######################################################################
  #########################  ##########################
  ######################################################################
  
  
  ######################################################################
  #########################  ##########################
  ######################################################################
  
  
  ######################################################################
  #########################  ##########################
  ######################################################################
  
  
  
  
  
  
}






######################################################################
######################### SCALE MERGED kIMS ##########################
######################################################################

#TODO

######################################################################
############################ WRAP UP #################################
######################################################################

message("Script done!")
