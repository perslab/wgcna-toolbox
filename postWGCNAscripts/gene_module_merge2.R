## script to merge gene modules from different WGCNA runs based on the overlap in genes

# Usage e.g.
# time Rscript /projects/jonatan/tools/wgcna-serc/wgcna-toolbox/gene_module_merge2.R path_dfNWA /projects/jonatan/tmp-epilepsy/tables/ep_ex_4_10000_genes_cell_cluster_module_genes.csv --geneWeightsCol pkIM   geneCol ="hgnc" --minPropIntersect  0.75 --minWeightedCor  0.75 --minPropIntersect_iter_increase  1.01 --maxIter25 --corMethod  "pearson" --dirOut "/projects/jonatan/tmp-epilepsy/" --prefixOut "modMerge2_test" --RAMGbMax 250

# Algorithm:
# Do until k == maxIter or until there are no modules to merge:
# 1. for each module j, find the module i with which it has the highest gene overlap as 
#     a proportion of the total genes in j (number in [0:1])
# 2. if minPropIntersect and minWeightedCor conditions are met, 
#     before merging module j into i, check if module i itself might be 
#     merged into a third module m, and the intersecting genes between 
#     i and m are more correlated, weighted by the gene weights and 
#     the overall intersect size, than that of j and i

#TODO
# * renormalize gene weights?
# * compute correlations?

######################################################################
########################## FUNCTIONS #################################
######################################################################

source(file="/projects/jonatan/tools/functions-src/utility_functions.R")

######################################################################
########################### OPTPARSE #################################
######################################################################

ipak(c("optparse"))

option_list <- list(
  make_option("--path_dfNWA", type="character",
              help = "path to dataframe from a gene network analysis run, in long format (i.e. one row per gene per celltype), containing 'cell_cluster', 'module', one or two gene name columns, and a column of numeric scores. I.e. one row per gene, e.g. rwgcna cell_cluster_module_genes.csv files"),  
  make_option("--geneWeightsCol", type="character", default ="pk.*M",
              help = "nwa_df column with gene weights, , [default %default]"),  
  make_option("--geneCol", type="character", default="hgnc",
              help ="string or regex for grepping nwa_df e.g. 'hgnc|symbol|gene_name' or 'ensembl', [default %default]"),
    make_option("--minPropIntersect", type="numeric", default = 0.65,
              help = "minimum proportion of module j intersecting with module i to merge j into i, [default %default]"),
  make_option("--minPropIntersect_iter_increase", type="numeric", default = 1.02,
              help = "At the end of each iteration, multiply the minPropIntersect threshold by some numeric constant >= 1 to help ensure convergence given that merging modules increases the probability that others will have a large overlap with them, [default %default]"),
  make_option("--minWeightedCor", type="numeric", default = 0.8,
              help = "if module j has at least minPropIntersect genes also in module i, set minimum correlation weighted by gene weight, to merge j into i, [default %default]"),
  make_option("--corMethod", type="character", default = "pearson",
              help = "pearson, spearman or kendall, [default %default]"),
  make_option("--maxIter", type="integer", default = 20,
              help = "maximum number of merge iterations , [default %default]"),
  make_option("--dirOut", type="character",
              help = "Outputs go to /tables and /RObjects subdirectories"),  
  make_option("--prefixOut", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"),
  make_option("--doPlot", type="logical", default = F,
              help = "Whether or not to plot intersect and correlation matrices at each iteration, [default %default]"),
  make_option("--RAMGbMax", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("dplyr", "Matrix", "parallel", "WGCNA", "corrplot", "readr", "RColorBrewer"))

######################################################################
############################ PARAMS ##################################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

path_dfNWA <- opt$path_dfNWA 
geneWeightsCol <- opt$geneWeightsCol
geneCol <- opt$geneCol
maxIter <- opt$maxIter
dirOut <- opt$dirOut
prefixOut <- opt$prefixOut
corMethod <- opt$corMethod
minPropIntersect <- opt$minPropIntersect
minPropIntersect_iter_increase <- opt$minPropIntersect_iter_increase
minWeightedCor <- opt$minWeightedCor
doPlot <- opt$doPlot
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

df_geneModule <- load_obj(path_dfNWA)

geneCol <- grep(geneCol, colnames(df_geneModule), value=T)
geneWeightsCol <- grep(geneWeightsCol, colnames(df_geneModule), value=T)#colnames(df_geneModule)[which(sapply(X=colnames(df_geneModule), FUN =function(colname) class(df_geneModule[[colname]]))=="numeric")]

######################################################################
###################### INITIALISE MERGE LOOP #########################
######################################################################

df_geneModule[["module_merged"]] <- df_geneModule[["module"]]
df_geneModule[["cell_cluster_merged"]] <- df_geneModule[["cell_cluster"]]

#df_geneModule[[paste0(geneWeightsCol, "_merged")]] <- df_geneModule[[geneWeightsCol]] 
pathLog <- paste0(dirLog, prefixOut, "_mod_merge_log.txt")
iteration=1
modules_to_exclude <- c()
while (T) {
  
  message(paste0(prefixOut, " merge iteration ", iteration))  
  ######################################################################
  ########################### GET MODULES ##############################
  ######################################################################
  
  unique(df_geneModule[["module_merged"]]) %>% sort %>% na.omit %>% as.character -> modules  
  
  modules <- modules[!modules %in% modules_to_exclude]
  
  ######################################################################
  ######################## GET GENE WEIGHT VECTORS #####################
  ######################################################################
  # get a vector of gene weights, with gene as name, for each module
  
  message("Getting gene weight vectors")
  
  fun <- function(module) {
    df_geneModule[[geneWeightsCol]] %>% '['(df_geneModule[["module_merged"]]==module) %>% na.omit %>% as.numeric ->geneWeights 
    df_geneModule[[geneCol]] %>% '['(df_geneModule[["module_merged"]]==module) %>% na.omit %>% as.character -> names(geneWeights) 
    return(geneWeights)
  }
  
  args = list("X"=modules)
  list_geneWeights <- safeParallel(fun=fun,args=args)
  
  ######################################################################
  ############## Compute the mod j - mod i prop. overlap matrix ########
  ######################################################################
  
  # module * module matrix
  # for each column j, each row i gives the proportion of module j that overlaps with module i 
  
  message("Computing module-module gene intersect matrix")
  
  fun = function(geneWeights) {
    sapply(list_geneWeights, function(geneWeightsOther) {
      base::intersect(x=names(geneWeights), y=names(geneWeightsOther)) %>% length %>% '/'(length(geneWeights))
    }, simplify=T)
  }
  
  args = list("X"=list_geneWeights)
  mat_moduleGeneIntersect <- safeParallel(fun=fun,args=args, simplify = T)
  # set diagonal to zero
  mat_moduleGeneIntersect[col(mat_moduleGeneIntersect)==row(mat_moduleGeneIntersect)] <- 0

  ### Compute modules to exclude in the coming iterations ###
  # They have no big intersections with other modules as a prop of their genes
  logical_colNoBigIntersect <- apply(mat_moduleGeneIntersect, MARGIN=2, FUN=function(j) max(j) < 0.75*minPropIntersect)  
  # No other modules have a big intersection with them
  logical_rowNoBigIntersect <- apply(mat_moduleGeneIntersect, MARGIN=1, FUN=function(j) max(j) < 0.75*minPropIntersect)  
  logical_noBigIntersect <- logical_colNoBigIntersect | logical_rowNoBigIntersect 
    
  modules_to_exclude <- c(modules_to_exclude, modules[!logical_noBigIntersect])
  
  ######################################################################
  ####################### Plot the overlap matrix ######################
  ######################################################################
  if (doPlot){
    message("Plotting gene intersect matrix")
    
    pdf(sprintf("%s%s_iter_%s_mod_intersect.pdf", dirPlots, prefixOut, k), 
        width=max(20,ncol(mat_moduleGeneIntersect) %/% 3),
        height=max(20,ncol(mat_moduleGeneIntersect) %/% 3))
    corrplot(corr = mat_moduleGeneIntersect, 
             method = "color",
             col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")), bias = 1)(200),
             diag = F,
             is.corr=F,
             title=paste0(prefixOut, ": prop. of mod j intersecting mod i, iteration k"),
             order = "original",#hclust",
             hclust.method = "average",
             addCoef.col = "black",
             tl.srt = 45,
             number.digits = 2L,
             number.cex = 0.5)
    
    invisible(dev.off())
  }
  ######################################################################
  ########### Compute the overlap weighted correlation matrix ##########
  ######################################################################
  
  message("Computing weighted intersect gene weight correlations for highly overlapping modules")
  
  mat_modIntersectCorrWeighted <- matrix(data=0, nrow=nrow(mat_moduleGeneIntersect), ncol=ncol(mat_moduleGeneIntersect))
  dimnames(mat_modIntersectCorrWeighted) <- dimnames(mat_moduleGeneIntersect)
  
  # Find columns of the overlap matrix with a big intersect
  idx_colBigIntersect <- apply(mat_moduleGeneIntersect, MARGIN=2, FUN=function(j) max(j) >= minPropIntersect)  %>% which
  idx_rowBigIntersect <- apply(mat_moduleGeneIntersect, MARGIN=2, which.max) %>% '['(idx_colBigIntersect)
   
  if (length(idx_colBigIntersect)==0) {
    message("no more modules to merge")
    break
  }
  
  # For columns and rows with big intersect, calculate the correlations
  for (k in 1:length(idx_colBigIntersect)) {
    j=idx_colBigIntersect[k]
    i=idx_rowBigIntersect[k]
    genesIntersect <- base::intersect(names(list_geneWeights[[j]]), names(list_geneWeights[[i]]))
    mat_modIntersectCorrWeighted[i,j] <- WGCNA::cor(x = list_geneWeights[[j]][match(genesIntersect, names(list_geneWeights[[j]]))], 
                                                     weights.x = list_geneWeights[[j]][match(genesIntersect, names(list_geneWeights[[j]]))],
                                                     y= list_geneWeights[[i]][match(genesIntersect, names(list_geneWeights[[i]]))],
                                                     weights.y = list_geneWeights[[i]][match(genesIntersect, names(list_geneWeights[[i]]))],
                                                     method = corMethod, 
                                                     quick=0.25, 
                                                     verbose=F, 
                                                     nThreads=30)
      
    
  }
  
  # for (k in 1:length(idx_colBigIntersect)) {
  #   j=idx_colBigIntersect[k]
  #   i=idx_rowBigIntersect[k]
  #   print(mat_modIntersectCorrWeighted[i,j])
  # }
  
  #args = list("X"=list_geneWeights)
  
  # the outer loop goes in the columns, the inner loop goes in the rows
  

  # fun = function(geneWeights) {
  #   sapply(list_geneWeights, function(geneWeightsother) {
  #     genesintersect <- base::intersect(names(geneWeights), names(geneWeightsother))
  #     propIntersect <- length(genesintersect)/length(geneWeights) 
  #     if (propIntersect>=minPropIntersect) {
  #       intersectCorrWeighted <-   
  #       return(intersectCorrWeighted)
  #       
  #     } else {
  #       return(0)
  #     }
  #   }, simplify=T)
  # }
  
  #mat_modIntersectCorrWeighted <- safeParallel(fun=fun,args=args, simplify = T)
  
  # set diagonal to zero
  #mat_modIntersectCorrWeighted[col(mat_modIntersectCorrWeighted)==row(mat_modIntersectCorrWeighted)] <- 0
  
  ######################################################################
  ############### Plot the overlap weighted correlation matrix #########
  ######################################################################
  # if(doPlot) {
  #   message("Plotting gene intersect weighted correlation matrix")
  #   
  #   pdf(sprintf("%s%s_iter%s_mat_modIntersectCorrWeighted.pdf", dirPlots, prefixOut, k), 
  #       width=max(20,ncol(mat_modIntersectCorrWeighted) %/% 3),
  #       height=max(20,ncol(mat_modIntersectCorrWeighted) %/% 3))
  #   corrplot(corr = mat_modIntersectCorrWeighted, 
  #            method = "color",
  #            col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")), bias = 1)(200),
  #            diag = F,
  #            is.corr=F,
  #            title=paste0(prefixOut, ": weighted correlation of mod j's intersection with mod i, iteration k"),
  #            order = "original",#"hclust",
  #            hclust.method = "average",
  #            addCoef.col = "black",
  #            tl.srt = 45,
  #            number.digits = 2L,
  #            number.cex = 0.5)
  #   
  #   invisible(dev.off())
  # }
  ######################################################################
  ################### IDENTIFY MODULES TO MERGE ########################
  ######################################################################
  
  message("Identifying modules to merge")
  # 1. For each column in the intersection matrix, identify the maximum row value. If no values are over threshold, STOP
  
  # For each column j, find the index of the max row i of the intersect matrix, i!=j
  #idx_rowmax <- apply(mat_moduleGeneIntersect, MARGIN=2, which.max)
  
  # Find the corresponding values of the intersect
  #toprow_vals <- apply(mat_moduleGeneIntersect, MARGIN=2, max)
  
  # 2. For each maximum row value, check whether the column correspond to that row has any higher values.
  # Then check if the intersect meets the threshold
  # logical_colmerge <- toprow_vals >= sapply(idx_rowmax, function(j) {
  #   max(mat_moduleGeneIntersect[,j])
  # }) & toprow_vals >= minPropIntersect
  
  #   3. For each marked j:i pair, verify that the corresponding weighted correlation of the 
  # intersecting genes is > minWeightedCor. If not, do not merge
  # logical_highcorr <- sapply(which(logical_colmerge), function(j){
  #   mat_modIntersectCorrWeighted[idx_rowmax[j],j]>=minWeightedCor
  # })
  # 
  # logical_colmerge[which(logical_colmerge)][!logical_highcorr]

  
  logical2_colBigIntersectCorr <- sapply(1:length(idx_colBigIntersect), function(k){
    j = idx_colBigIntersect[k]
    i = idx_rowBigIntersect[k]
    # For each module pair j,i, before merging module j into i,
    # check if module i itself might be merged into a third module m, and 
    # the intersecting genes between i and m are more correlated, weighted
    # by the gene weights and the overall intersect size, than that of j and i
    
    idx_matchedMax <- which.max(mat_modIntersectCorrWeighted[,i])
    bool <- mat_modIntersectCorrWeighted[i,j]*mat_moduleGeneIntersect[i,j] >= max(mat_modIntersectCorrWeighted[idx_matchedMax,i])*mat_moduleGeneIntersect[idx_matchedMax,i]
    # is the intersect correlation sufficient?
    bool2 <- bool & mat_modIntersectCorrWeighted[i,j] >= minWeightedCor
  })

  if (!any(logical2_colBigIntersectCorr)) {
    message("no more modules to merge")
    break
  }
  
  ### Merge module j into i in df_geneModule dataframe:
  
  # ii. for duplicated genes, add NA in df_geneModule[["module_merged"]] for the j copy
  
  mods_to_merge <- modules[idx_colBigIntersect[logical2_colBigIntersectCorr]]
  mods_to_merge_into <- modules[idx_rowBigIntersect[logical2_colBigIntersectCorr]]

  ######################################################################
  ###################### UPDATE DF_GENEMODULE ##########################
  ######################################################################
  
  message("Updating module dataframe")
  
  for (i in 1:length(mods_to_merge)) {
    cell_cluster_to_merge_into <- df_geneModule[["cell_cluster"]][df_geneModule[["module"]]==mods_to_merge_into[i]][1]
    # module_merged column
    logical_to_merge <- df_geneModule[["module_merged"]] %in% mods_to_merge[i] 
    logical_to_merge_into <- df_geneModule[["module_merged"]] %in% mods_to_merge_into[i]

    df_geneModule[["module_merged"]][logical_to_merge] <- mods_to_merge_into[i]
    df_geneModule[["cell_cluster_merged"]][logical_to_merge] <- cell_cluster_to_merge_into
    
    # in module_merged column set duplicate genes in merged module to NA_character 
    logical_dup <- duplicated(df_geneModule[["hgnc"]]) 
    logical_dup_merged <- logical_dup & logical_to_merge
    
    #logical2_isna <- is.na(df_geneModule[["module_merged"]][logical_to_merge | logical_to_merge_into][logical2_dup])
    df_geneModule[["module_merged"]][logical_dup_merged] <- NA_character_
    df_geneModule[["cell_cluster_merged"]][logical_dup_merged] <- NA_character_ 
    
  }
  
  # write out log 
  log_entry <- paste0(prefixOut, " merge iteration ", iteration, " with minPropIntersect = ", round(minPropIntersect,3), " and minWeightedCor = ", minWeightedCor)
  message(log_entry)
  cat(log_entry, file = pathLog, append=T, sep = "\n")

  log_entry <- paste0("    merging ", mods_to_merge, " into ", mods_to_merge_into, collapse="\n")
  message(log_entry)
  cat(log_entry, file = pathLog, append=T, sep = "\n")
  
  iteration = iteration+1
 
  if (k>maxIter) {
    message(paste0(maxIter, " iterations reached, stopping"))
  }
  
  #increment the minPropIntersect
  if (!is.null(minPropIntersect_iter_increase)) minPropIntersect <- min(1, minPropIntersect*minPropIntersect_iter_increase)
  
}

######################################################################
############################ WRAP UP #################################
######################################################################

readr::write_csv(x=df_geneModule, path = gzfile(paste0(dirTables, prefixOut, "_cell_cluster_module_genes.csv.gz")))

message("Script done!")
