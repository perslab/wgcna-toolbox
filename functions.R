# Title: sc functions
# Author: Jonatan Thompson, Pers Lab


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

FilterGenes <- function(seurat_obj_sub, min.cells) {
  #unchanged <- seurat_obj_sub
  # EDIT 180501_4 remove tryCatch in parallelised functions
  #tryCatch({
  if (min.cells > 0) {
    num.cells <- rowSums(seurat_obj_sub@data > 0)
    genes.use <- names(x = num.cells[which(x = num.cells >= min.cells)])
    seurat_obj_sub@raw.data <- seurat_obj_sub@raw.data[genes.use, ]
    seurat_obj_sub@data <- seurat_obj_sub@data[genes.use, ]
  }
      # EDIT 180501_4 remove tryCatch in parallelised functions
      #   }}, error = function(c) { warning("parFilterGenes produced an error - returning unchanged seurat subset")
  #     return(seurat_obj_sub_backup)
  # })
  return(seurat_obj_sub)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

wrapJackStraw = function(seurat_obj_sub, n_cores, jackstraw.num.replicate) {
  # gene.criterion: 'p.val' means selecting genes (used for PCA) with significant empirical p-val
  #                 'PC.loadings' means projecting all genes onto the PCs to get loadings and selecting 
  #                 genes that have a high absolute loading on a significant PC
  # Schematic:
  # 1. pcs.compute = ncol

  ### 180507
  prop.freq <- max(0.012, round(4/length(seurat_obj_sub@var.genes),3)) # to ensure we have at least 3 samples so the algorithm works well
  # see https://github.com/satijalab/seurat/issues/5
  ###
  
  pcs.compute = ncol(seurat_obj_sub@dr$pca@gene.loadings)
  
  if (jackstraw.num.replicate > 0) {
    
    seurat_obj_sub <- JackStraw(object = seurat_obj_sub,
                                num.pc = pcs.compute,
                                num.replicate = jackstraw.num.replicate, 
                                display.progress = T,
                                do.par = T,
                                num.cores = n_cores,
                                prop.freq = prop.freq) # https://github.com/satijalab/seurat/issues/5
  
    score.thresh = (5e-2)/pcs.compute # TODO: Is this a suitable multiple testing adjustment
  
    ### EDIT_180426_2
    #pAll <- GetDimReduction(seurat_obj_sub, reduction.type = "pca", slot = "jackstraw")@emperical.p.value.full
    pAll <- GetDimReduction(seurat_obj_sub, reduction.type = "pca", slot = "jackstraw")@emperical.p.value
    ###
    pAll <- pAll[, 1:pcs.compute, drop = FALSE]
    pAll <- as.data.frame(pAll)
    pAll$Contig <- rownames(x = pAll)
    pAll.l <- melt(data = pAll, id.vars = "Contig")
    colnames(x = pAll.l) <- c("Contig", "PC", "Value")
    
    score.df <- NULL
    
    for (i in (1:pcs.compute)) {
      pc.score <- suppressWarnings(prop.test( 
        x = c(length(x = which(x = pAll[, i] <= score.thresh)), floor(x = nrow(x = pAll) * score.thresh)),
        n = c(nrow(pAll), nrow(pAll)))$p.val)
      if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
        pc.score <- 1
      }
      if (is.null(x = score.df)) {
        score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
      } else {
        score.df <- rbind(score.df, data.frame(PC = paste0("PC",i), Score = pc.score))
      }
    }
    
    PC_select_idx <- which(score.df$Score < score.thresh)
    
    if (nrow(pAll) == length(seurat_obj_sub@var.genes)) {
      
      seurat_obj_sub <- ProjectPCA(seurat_obj_sub, 
                                   do.print = F, 
                                   pcs.print = NULL, 
                                   pcs.store = pcs.compute, 
                                   genes.print = NULL, 
                                   replace.pc=F, 
                                   do.center=T)
        
      loadings <- abs(seurat_obj_sub@dr$pca@gene.loadings.full[,PC_select_idx])
      max_loadings <- apply(loadings, 1, function(x) max(x))
      names_genes_use <- names(max_loadings[order(max_loadings, decreasing = T)])[1:5000]
      
    } else if (nrow(pAll) > length(seurat_obj_sub@var.genes)) {
      
      pAll[,sapply(pAll, function(x) class(x)!="numeric")] <- NULL # remove the column of gene names
      row_min <- apply(pAll[,PC_select_idx], MARGIN = 1, FUN = function(x) min(x))
      names_genes_use <- rownames(pAll)[row_min < score.thresh]
      
      if (length(names_genes_use) < 1000) names_genes_use <- rownames(pAll)[row_min < score.thresh*2]
    }
    
  } else if (jackstraw.num.replicate == 0) {
    
    seurat_obj_sub <- ProjectPCA(seurat_obj_sub, 
                                 do.print = F, 
                                 pcs.print = NULL, 
                                 pcs.store = pcs.compute, 
                                 genes.print = NULL, 
                                 replace.pc=F, 
                                 do.center=T)
    
    loadings <- abs(seurat_obj_sub@dr$pca@gene.loadings.full[,1:(min(13,pcs.compute))])
    max_loadings <- apply(loadings, 1, function(x) max(x))
    names_genes_use <- names(max_loadings[order(max_loadings, decreasing = T)])[1:5000]
    
  }
    
    ### EDIT_180426_2
    ### EDIT_180429_3
    # seurat_obj_sub@dr$pca@gene.loadings.full <- seurat_obj_sub@dr$pca@gene.loadings.full[,PC_select_idx]
    # max_loads <- apply(abs(seurat_obj_sub@dr$pca@gene.loadings.full), MARGIN=1, max, T)
    
    #max_loads <- apply(abs(seurat_obj_sub@dr$pca@gene.loadings.full[,PC_select_idx]), MARGIN=1, max, T)
    
    
    ###
    ### EDIT 180430_3
    # we were getting poor results on Campbell neurons..
    # threshold <- quantile(max_loads, 0.5)
    # names_genes_use = names(max_loads[max_loads > threshold])
    #threshold <- quantile(max_loads, 0.5)
    #names_genes_use = names(max_loads[order(max_loads, decreasing=T)][1:5000])
    
    ### 180504_v1.7_5
    # datExpr <- seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,] %>% t()
    # colnames(datExpr) <- rownames(seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,])
    # ###
    datExpr <- seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,] %>% t() #%>% as.matrix() 
    colnames(datExpr) <- rownames(seurat_obj_sub@scale.data[rownames(seurat_obj_sub@scale.data) %in% names_genes_use,])
    
    
  return(datExpr)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

seurat_to_datExpr = function(seurat_obj_sub, idx_genes_use) {
  # Only if we don't use PCA loading genes
  datExpr <- seurat_obj_sub@scale.data[idx_genes_use,] %>% t() 
  colnames(datExpr) <- rownames(seurat_obj_sub@scale.data[idx_genes_use,])  
  return(datExpr)  
}



############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

sft_for_par <- function(datExpr, subsetName) { 
  softPower <- 8 # Set a default value as fall back
  # EDIT 180501_4 remove tryCatch in parallelised functions
  #tryCatch({
    
    disableWGCNAThreads()
    
    powers = c(1:30)
    
    sft = pickSoftThreshold(data=datExpr,
                            powerVector = powers,
                            blockSize = min(2500, ncol(datExpr)), #try to prevent crashing
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
    fitIndices_filter <- fitIndices %>% dplyr::filter(median.k. <= quantile(median.k.,0.95, na.rm=T))
    if (sum(fitIndices_filter$SFT.R.sq >= 0.9) > 1) {
      softPower = min(fitIndices_filter$Power[fitIndices_filter$SFT.R.sq>=0.9])
    } else {
      softPower <- (fitIndices_filter %>% dplyr::arrange(desc(SFT.R.sq)) %>% dplyr::select(Power))[1,]
    }
  # EDIT 180501_4 remove tryCatch in parallelised functions
  # }, error = function(c) {
  #   write.csv("ERROR", file=sprintf("%s%s_%s_pickSoftThreshold_ERROR_set_to_8_as_default_%s.csv", project_dir, data_prefix, subsetName, flag_date), row.names = F)
  # })
  return(softPower)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

bootstrap <- function(datExpr,
                      nPermutations,
                      replace,
                      fraction,
                      randomSeed)
  # @Usage: Resample a dataset.
  # @args:
  #       datExpr: Dataset with samples in the rows, is coerced to matrix
  #       nPermutations: runs
  #       replace: sample with replacement?
  #       fraction: sample what fraction of the total each iteration?
  #       randomSeed: initial random seed
  # @return:
  #       result: list of resampled datasets in matrix format. If replace = T, first item is the unpermuted dataset.
  # @author: Jonatan Thompson jjt3f2188@gmail.com
  # @date: 180222

{
  startRunIndex = 1
  endRunIndex = startRunIndex + nPermutations 
  result = vector("list", length= if (replace==T) nPermutations+1 else nPermutations);
  nSamples = nrow(datExpr);
  nGenes = ncol(datExpr);
  
  for (run in startRunIndex:endRunIndex)
    
  {
    set.seed(randomSeed + 2*run + 1);
    
    if (run == startRunIndex & replace == T) {
      useSamples = c(1:nSamples) # The first run just returns the full dataset
    } else if (run>startRunIndex | replace == F) 
    {  useSamples = sample(nSamples, as.integer(nSamples * fraction), replace = replace)
    } 
    
    
    samExpr = as.matrix(datExpr[useSamples, ]);
    
    result[[run]]$data <- samExpr
  }
  
  return(result)
  
}
  

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################


TOM_for_par = function(datExpr, subsetName, softPower) {
  
  ### 180503 v1.7
  disableWGCNAThreads()
  ###
  
  adjacency = adjacency(datExpr, 
                        type=type, 
                        power = softPower, 
                        corFnc = corFnc, 
                        corOptions = corOptions)
  
  consTomDS = TOMsimilarity(adjMat=adjacency,
                            TOMType=TOMType,
                            TOMDenom=TOMDenom,
                            verbose=verbose,
                            indent = indent)
  
  save(consTomDS, file=sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, subsetName)) # Save TOM the way consensusTOM would have done
  goodGenesTOM_idx <- rep("TRUE", ncol(datExpr))
  return(goodGenesTOM_idx)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

dissTOM_for_par = function(subsetName) {
  # We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes.
  # We use Dynamic Tree Cut from the package dynamicTreeCut.    
  load(sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, subsetName)) # Load the consensus TOM generated by blockwiseConsensusModules
  
  dissTOM <- 1-as.dist(consTomDS) # Convert proximity to distance
  rm(consTomDS)

  return(dissTOM)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

cutreeHybrid_for_vec <- function(comb, geneTree, dissTOM) {
  # Utility function for more easily parallellising the cutreeHybrid function
  tree = cutreeHybrid(dendro = geneTree, 
                      cutHeight = NULL,
                      minClusterSize = comb[[1]], 
                      distM=as.matrix(dissTOM),
                      deepSplit=comb[[2]], 
                      pamStage=comb[[3]],
                      pamRespectsDendro=comb[[3]]) #TODO: would it make sense at all to have pamStage = T w/o pamRD = T as well ?
  # Gandal et al 2018:  cutHeight = 0.999
  return(tree)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

mergeCloseModules_for_vec <- function(cutree,comb, datExpr, excludeGrey) {
  # Utility function for more easily parallelising mergeCloseModules
  if (any(as.logical(cutree$labels))) { # Are there any modules at all? Avoids an error
    merged = mergeCloseModules(exprData=as.matrix(datExpr), 
                               colors = cutree$labels, 
                               cutHeight=comb[[4]])
    colors = labels2colors(merged$colors)
    #MEs = merged$newMEs 
    MEs = moduleEigengenes(expr = as.matrix(datExpr),
                           colors,
                           excludeGrey = excludeGrey)
  }
  else {# only grey modules
    colors = labels2colors(cutree$labels)
    MEs = NULL
  }
  return(list("cols" = colors, "MEGs"= MEs))
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# parvecMergeCloseModules = function(list_cutree, datExpr) {
#   # Inherits comb_list from the parent environment, same for every subsetName, no need to pass as an argument
#   list_merged <- mapply(function(x,y) vecMergeCloseModules(cutree=x,comb=y, datExpr_filter=datExpr_filter), list_cutree, comb_list, SIMPLIFY=FALSE) # Merge modules whose eigengenes are highly correlated
#   return(list_merged)
# }

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

parGetMEs = function(list_merged) {
  list_MEs = lapply(list_merged, function(x) x$MEGs) # list of merged color eigengenes
  return(list_MEs)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

extract_and_name_colors <- function(merged, datExpr) {
  colors <- merged[[1]]
  names(colors) = colnames(datExpr)
  return(colors)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

parGetColors = function(list_merged, datExpr) {
  list_colors = lapply(list_merged, function(x) extract_and_name_colors(merged=x, datExpr=datExpr)) # list of merged colors
  return(list_colors)
}  

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

parMatchColors <- function(list_colors) {
  # Inherits n_combs from the parent environment, same for every subsetName, no need to pass as an argument
  if (length(list_colors)>1) {
    # Match colors between iterations by reference to the set of colors with the highest number of unique colors
    max_colors_idx = which.max(lapply(list_colors, function(x) length(table(x))))
    
    diffParams_colors = matrix(0, nrow=length(list_colors[[1]]), ncol=n_combs)
    diffParams_colors[, max_colors_idx] = list_colors[[max_colors_idx]] # use the set of colors with the highest number of modules as 'reference'
    
    for (j in (1:n_combs)[-max_colors_idx]) {
      
      if (length(table(list_colors[[j]])) > 1) {
        diffParams_colors[,j] <- matchLabels(source = list_colors[[j]], # traverse gene colour assignments
                                             reference = diffParams_colors[,max_colors_idx],
                                             pThreshold = 5e-2,
                                             na.rm = TRUE,
                                             extraLabels = standardColors())#paste0("Extra_", 0:500)) # labels for modules in source that cannot be matched to any modules in reference
      
        } else if (length(table(list_colors[[j]])) == 1) {
        diffParams_colors[,j] <- "grey"
      }
    }
    # Split the matrix of matched colors back into a list 
    list_colors_matched <- split(diffParams_colors, rep(1:ncol(diffParams_colors), each = nrow(diffParams_colors))) 
    
  } else if (length(list_colors)==1)  {
    diffParams_colors <- as.matrix(list_colors[[1]], nrow=length(list_colors[[1]]), ncol=1)
    list_colors_matched <- list_colors
  } 
  
  return(list_colors_matched)
}



############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

parkMEs = function(list_MEs, datExpr) {
  
  # Compute kMEs and pkMEs for each set of colors corresponding to distinct parameters
  list_kMEs <- NULL
  list_kMEs <- vector(mode="list", length=length(list_MEs))
 
  for (j in 1:length(list_MEs)) {
    list_kMEs[[j]] = signedKME(datExpr=as.matrix(datExpr),
                               datME=list_MEs[[j]]$eigengenes,
                               #list_MEs[[j]],#$eigengenes,
                               #outputColumnName="",
                               corFnc = corFnc)
  }
  return(list_kMEs)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

replaceNA = function(replace_in, replace_from) {
  replace_in[is.na(replace_in)] <- replace_from[is.na(replace_in)]
  return(replace_in)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

parPkMEs = function(list_kMEs, list_colors) {
  
  list_pkMEs <- list()
  
  for (j in 1:length(list_colors)) {
    # Get the 'principal kMEs', i.e. kME of each gene to the module to which it is allocated
    # EDIT_180421_8
    #pkMEs <- vector(mode="numeric",length=nrow(diffParams_colors))
    
    pkMEs <- vector(mode="numeric",length=length(list_colors[[1]]))
    ###
    # Loop over every gene and get its pkME
    for (i in 1:length(pkMEs)) {
      #pkMEs[i] <- list_kMEs[[j]][diffParams_colors[i,j]][i,]
      pkMEs[i] <- list_kMEs[[j]][paste0("kME", list_colors[[j]][[i]])] [i,]    
    }
    
    list_pkMEs[[j]] <- pkMEs 
    
  }
  return(list_pkMEs)
}



parPkMEs_2 = function(list_kMEs, list_colors) {
  
  list_pkMEs <- list()
  
  for (j in 1:length(list_colors)) {
    # Get the 'principal kMEs', i.e. kME of each gene to the module to which it is allocated
    # EDIT_180421_8
    #pkMEs <- vector(mode="numeric",length=nrow(diffParams_colors))
    
    pkMEs <- vector(mode="numeric",length=length(list_colors[[1]]))
    ###
    # Loop over every gene and get its pkME
    for (i in 1:length(pkMEs)) {
      #pkMEs[i] <- list_kMEs[[j]][diffParams_colors[i,j]][i,]
      pkMEs[i] <- list_kMEs[[j]][list_colors[[j]][[i]]] [i,]    
    }
    
    list_pkMEs[[j]] <- pkMEs 
    
  }
  return(list_pkMEs)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# Delete grey kMEs and kIMs
# deleteGrey <- function(list_kMs) {
#   for (k in 1:length(list_k)) {
#     if (any(grepl("grey", colnames(list_kMs[[k]])))) list_kMs[[k]][['grey']] <- NULL
#     if (any(grepl("kMEgrey", colnames(list_kMs[[k]])))) list_kMs[[k]][['kMEgrey']] <- NULL
#   }
#   return(list_kMs)
# }

deleteGrey <- function(list_kMs) {
  for (k in 1:length(list_kMs)) {
    if (any(colnames(list_kMs[[k]]) == "grey")) list_kMs[[k]][['grey']] <- NULL
    if (any(colnames(list_kMs[[k]]) == "kMEgrey")) list_kMs[[k]][['kMEgrey']] <- NULL
  }
  return(list_kMs)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

kME_reassign_fnc = function(MEs, kMEs, pkMEs, kME_reassign_threshold, colors, filter, corFnc, excludeGrey) {
  # Status: hiatus - currently unused
  # Returns a list with new colors,l MEs, kMEs, pkMEs, and dataframe with information on the reassigned genes
  
  results <- list()
  old_colors = colors 
  reassigned = NULL
  
  col_most_negkME = max.col(-1*kMEs) # indices of most negative kMEs
  
  rows_to_reassign <- logical(length=length(colors))
  
  for (i in nrow(kMEs)) {
    rows_to_reassign[i] <- kMEs[i, col_most_negkME[i]] < -kME_reassign_threshold*(pkMEs[i])
  }
  
  if (sum(rows_to_reassign) > 0) {
    
    reassigned = colors[rows_to_reassign]
    colors[rows_to_reassign] <- colors[col_most_negkME][rows_to_reassign]
    
    # Recompute Module Eigengenes, kME, pkME
    MEs = moduleEigengenes(expr = as.matrix(datExpr_filter),
                           colors,
                           excludeGrey = excludeGrey)
    
    kMEs = signedKME(as.matrix(datExpr_filter),
                     MEs$eigengenes,
                     outputColumnName = "",
                     corFnc = corFnc)
    
    pkMEs <- vector(length=length(colors))
    
    
    for (i in 1:length(colors)) {
      pkMEs[i] <- as.numeric(unlist(kMEs %>% dplyr::select(matches(colors[[i]]))))[[i]]
    }
    
    #reassigned = data.frame(gene=names(colors[rows_to_reassign]), old_module = old_colors, new_module = colors[rows_to_reassign], old_pkME = old_pkMEs, new_pkME = pkMEs[rows_to_reassign], row.names = NULL)
  }
  
  results$colors <- colors
  results$old_colors <- old_colors
  results$MEs <- MEs
  results$kMEs <- kMEs
  results$pkMEs <- pkMEs
  results$reassigned <- names(colors[rows_to_reassign]) 
  
  return(results)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

kIM_eachMod_norm = function(dissTOM, colors, genes) {
  # compute intramodularConnectivity for every gene with regard to every cluster
  unique_colors = sort(names(table(colors)))#[-which(sort(names(table(colors))) == "grey")]
  
  mTOM <- as.matrix(1-dissTOM)
  
  results = matrix(nrow=length(genes), ncol=length(unique_colors))
  normTerms <- numeric(length(unique_colors))
  
  for (j in 1:length(unique_colors)) {
    for (i in 1:length(genes)) {
      results[i,j] <- sum(mTOM[i,colors %in% unique_colors[[j]]])
    }
    
    #normTerms[j] <- max(colSums(mTOM[,colors %in% unique_colors[[j]]]))
    normTerms[j] <- sum(colSums(mTOM[,colors %in% unique_colors[[j]]]))
  }
  
  results_norm <- as.data.frame(t(t(results)/normTerms), col.names = unique_colors, row.names=genes)

  #colnames(results_norm) <- unique_colors
  #rownames(results_norm) <- genes

  return(results_norm)
  
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

PrePPIfilter_for_vec <- function(list_pkMEs, list_colors) {
  # Filter the genes so we only submit to STRINGdb if abs(pkME) > 0.1
  list_idx_pkME_ok <- lapply(list_pkMEs, function(x) abs(x) > 0.1)
  list_colors_pkME_ok <- mapply(function(x,y) x[y], x=list_colors, y=list_idx_pkME_ok, SIMPLIFY = F)
  return(list_colors_pkME_ok)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# parVecPPI_outer <- function(list_colors_pkME_ok) {
#   list_colors_PPI <- lapply(list_colors_pkME_ok, 
#                             function(x) vecPPI_outer(colors = x,
#                                                      STRINGdb_species = STRINGdb_species,
#                                                      PPI_pval_threshold = PPI_pval_threshold,
#                                                      project_dir = project_dir, 
#                                                      data_prefix = data_prefix, 
#                                                      flag_date = flag_date))
#   
#   return(list_colors_PPI)
# }

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

count_grey_in_list_of_vec = function(list_colors) {
  vec_n_grey <- sapply(list_colors, function(x) sum(x=="grey"), simplify=T)
  return(vec_n_grey)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# TODO: why does this differ from the version above? Can they be merged?

getPkMEs <- function(colors, kMEs) {
  pkMEs <- vector(length=length(colors))
  
  for (i in 1:length(colors)) {
    pkMEs[i] <- kMEs[colors[i]][i,]
  }
  
  return(pkMEs)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# replace_unplottable_colors <- function(colors) {
#   # takes a vector of colors, i.e. module assignments
#   # if all are recognised R colors then it returns the vector as is
#   # otherwise it replaces them and returns the vector
#   
#   idx_not_real_colors = !(colors %in% colors())
#   n_replace <- sum(idx_not_real_colors)
#   
#   if (n_replace>0) {
#     unused = setdiff(colors(), colors) 
#     colors[idx_not_real_colors] <- unused[1:n_replace]
#   }
#   return(colors)
# }

### 180507_v1.8_dev2
replace_unplottable_colors <- function(colors) {
  # takes a vector of colors, i.e. module assignments
  # if all are recognised R colors then it returns the vector as is
  # otherwise, first try to remove full stops and numbers. Then try to replace them and returns the vector
  
  idx_not_real_colors = !(colors %in% colors())
  n_replace <- sum(idx_not_real_colors)
  
  if (n_replace>0) {
    colors[idx_not_real_colors] <- gsub("[\\.123456789]", "", colors[idx_not_real_colors])
  }

  idx_not_real_colors = !(colors %in% colors())
  n_replace <- sum(idx_not_real_colors)

  if (n_replace>0) {
    unused = setdiff(colors(), colors) 
    colors[idx_not_real_colors] <- unused[1:n_replace]
  }
  
  return(colors)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

checkGrey <- function(kMEs) {
  if (any(grepl("grey", colnames(kMEs)))) kMEs[['grey']] <- NULL 
  return(kMEs)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

parPlotDiffParams <- function(list_colors, subsetName) {
  
  list_colors <- lapply(list_colors, replace_unplottable_colors)
  
  diffParams_colors_order = matrix(unlist(list_colors), nrow=length(list_colors[[1]]), ncol=length(list_colors))
  
  pdf(sprintf("%s%s_%s_diffParams_colors_order_%s.pdf", plots_dir, data_prefix, subsetName, flag_date),width=9,height=6+length(list_colors) %/% 2)
  plotDendroAndColors(list_geneTree_ok[[subsetName]], 
                      #labels2colors(list_colors_PPI_order),
                      #matrix(unlist(list_colors_PPI_order), ncol = length(list_colors_PPI_order), byrow = F), 
                      diffParams_colors_order,
                      groupLabels = list_list_plot_label_ok_order[[subsetName]], 
                      addGuide= TRUE, 
                      dendroLabels=FALSE, 
                      main=sprintf("%s modules under different params, pre-PPI test, ranked by assigned genes in PPI test.", subsetName), 
                      cex.colorLabels=0.4)

  dev.off()
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

parPlotDiffParams_PPI <- function(list_colors, subsetName) {
  
  list_colors <- lapply(list_colors, replace_unplottable_colors)
  
  diffParams_colors_PPI_order = matrix(unlist(list_colors), nrow=length(list_colors[[1]]), ncol=length(list_colors))
  
  pdf(sprintf("%s%s_%s_diffParams_colors_PPI_order_%s.pdf", plots_dir, data_prefix, subsetName, flag_date),width=9,height=6+length(list_colors) %/% 2)
  plotDendroAndColors(list_geneTree_ok[[subsetName]], 
                      #labels2colors(list_colors_PPI_order),
                      #matrix(unlist(list_colors_PPI_order), ncol = length(list_colors_PPI_order), byrow = F), 
                      diffParams_colors_PPI_order,
                      groupLabels = list_list_plot_label_ok_order[[subsetName]], 
                      addGuide= TRUE, 
                      dendroLabels=FALSE, 
                      main=sprintf("%s modules under diff. params., PPI, ranked by assigned genes.", subsetName), 
                      cex.colorLabels=0.4)
  dev.off()
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

plotFinalColors <- function(subsetName) {
  
  ### EDIT180430
  list_colors_both_to_plot <- lapply(list_colors_both, function(x) as.data.frame(apply(x, MARGIN=2, FUN=function(y) replace_unplottable_colors(y))))
  ###
  
  pdf(sprintf("%s%s_%s_final_colors_%s.pdf", plots_dir, data_prefix, subsetName, flag_date),width=9,height=5)
  plotDendroAndColors(list_geneTree_ok[[subsetName]],
                      list_colors_both_to_plot[[subsetName]],
                      groupLabels = if (!is.null(dim(list_colors_both_to_plot[[subsetName]]))) c(list_plot_label[[subsetName]], "PPI enriched") else list_plot_label[[subsetName]],
                      addGuide=T,
                      dendroLabels=F,
                      main=sprintf("%s final colors, pre & post PPI test", subsetName),
                      cex.colorLabels=0.4)
  dev.off()
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

plotCorr_for_par = function(corrMatrix, vectorNames, subsetName, diag, is.corr, order, hclust.method) {
  pdf(sprintf("%s%s_%s_%s_corr_%s.pdf", plots_dir, data_prefix, subsetName, vectorNames, flag_date))#,width=nrow(corrMatrix) %/% 4, height=(nrow(corrMatrix) %/% 4)+1)
  corrplot(corr = corrMatrix,
           method = "color",
           add=F,
           #col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           title=sprintf("%s %s correlations", subsetName, vectorNames),
           is.corr = is.corr,
           diag=diag,
           #mar = c(0, 0, 0, 0),
           #addCoef.col = T, # this adds an annoying number label to each square
           order = order,
           hclust.method = hclust.method,
           number.digits = NULL)
  #clustering_distance_rows = dist(resR.G_G@listData$pvalue[select], method="euclidian"))
  dev.off()
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

### 180516 Not a good idea - rather use the eigen_mat

# get_kMEs_all_genes = function(kMEs, all_genes) {
#   # works equally well for kIMs
#   # kMEs should be a named dataframe or matrix
#   # all_genes should be a vector
#   
#   out <- matrix(0, nrow=length(all_genes), ncol = ncol(kMEs))
#   rownames(out) <- all_genes
#   colnames(out) <- colnames(kMEs)
#   row_idx <- all_genes %in% rownames(kMEs)
#   
#   ### EDIT_180429_5
#   # out[row_idx,] <- as.matrix(kMEs)[row_idx,]
#   out[row_idx,] <- as.matrix(kMEs) 
#   ###
#   return(out)
#   
# }

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

wrapModulePreservation <- function(listDatExpr,
                                   listColors,
                                   labels = if (is.null(names(listColors))) 1:length(listColors) else names(listColors),
                                   dataIsExpr,
                                   networkType, 
                                   corFnc,
                                   corOptions,
                                   # We actually just want to get quality stats on the reference network,
                                   # but the function requires testNetworks..
                                   nPermutations, 
                                   includekMEallInSummary,
                                   restrictSummaryForGeneralNetworks,
                                   calculateQvalue,
                                   randomSeed, 
                                   maxGoldModuleSize, 
                                   maxModuleSize, 
                                   quickCor, 
                                   ccTupletSize, 
                                   #calculateCor.kIMall,
                                   calculateClusterCoeff,
                                   useInterpolation, 
                                   checkData, 
                                   #greyName, 
                                   savePermutedStatistics , 
                                   loadPermutedStatistics, 
                                   permutedStatisticsFile, 
                                   plotInterpolation, 
                                   interpolationPlotFile, 
                                   discardInvalidOutput,
                                   parallelCalculation,
                                   verbose, 
                                   indent) {
  
  # @Usage: Wrapper for the WGCNA modulePreservation function to make it easier to use for single and multiple datasets.
  #        modulePreservation performs two distinct tasks:
  #           1. evaluates the quality of one or more sets of module assignments in a single dataset
  #           2. evaluates the preservation of a single set of colors in multiple datasets 
  #         However, in each case, modulePreservation function  requires 'reference' and 'test datasets regardless.
  #         The function also requires the data to be presented in a multiExpr format (see WGCNA function checkSets).
  #         In the first case, this wrapper function makes 'test' datasets just by copying the original to satisfy the required modulePreservation args.
  #         In the second case, this wrapper assumes that the first dataset is the reference and the rest are test sets, 
  #         and calls modulePreservation do do pairwise ref-test comparisons.
  #         The wrapper puts datasets into the correct multiExpr format.
  #
  # @args: 
  #         listDatExpr: a vector(list) of expression data (genes in columns, cells in rows). If the list contains a single dataset the function outputs the 
  #                     quality statistics only. Assumes that the first entry holds the reference dataset and the rest hold test sets, if relevant.
  #                     If the list entries are names the function uses them as labels by default unless a list of labels is provided.
  #         listColors: a vector(list) of color assignment vectors? TODO If given a single set of colors for multiple datasets the function 
  #                     cannot compute hypergeometric stats ? TODO
  #         labels:     a vector of character entries, one per set of colors. If not given, defaults to names(listColors); if these are NULL, defaults to integers
  #         ...         (see modulePreservation)
  #
  # @return: 
  #         result:     a list of list of dataframes? TODO, each nested list corresponds to a pairwise evaluation     
  # @author: Jonatan Thompson jjt3f2188@gmail.com
  # @date: 180316
  
  stopifnot(is.null(dim(listDatExpr)) & typeof(listDatExpr) == "list" & length(listDatExpr) >= 1 & is.null(dim(listColors)) & typeof(listColors)=="list" & length(listColors) >= 1) # basic input format checks
  stopifnot(length(listDatExpr) == 1 | length(listColors) == 1) # Cannot have many datasets and many colourings so at least one list must have length 1
  
  referenceNetworks = c(1:length(listColors)) # 
  testNetworks = as.list(rep(length(listColors)+1, length(listColors)))
  
  if (verbose >= 1) {
    if (length(listDatExpr) == 1) {
      if (length(listColors) == 1) {
        print("Measuring the quality of a module assignment on a dataset..")
      }
      else if (length(listColors) > 1) {
        print("Measuring the quality of several module assignments on a dataset..")
      }
    }
    
    # Set up the multi-set format (see WGCNA checkSets() ). The modulePreservation function requires test datasets even when evaluating module assignments in a single dataset,
    # so we create a multi-set format list with length(colors)+1 copies of the dataset, where the extra acts as test set.
    # the reason for this memory-inefficient copying of datasets is that modulePreservation requires a dataset to match each set of colors.
    
    multiExpr <- vector("list", length=length(listColors)+1)  # make a list of copies of the expression data, number of colourings to evaluate + test set
    for (i in 1:length(multiExpr)) {
      multiExpr[[i]]$data  <- as.matrix(listDatExpr[[1]]) # NB: this whole section only runs if there is just one dataset
    }
    
    names(multiExpr) <- c(labels, "test")
    
    multiColor <- vector("list", length = length(listColors)+1) # make a matching list with the colors to evaluate. 
    
    multiColor[[length(multiColor)]] <- vector("character", length = length(listColors[[1]])) # set last set of colors to empty strings 
    
    for (i in 1:(length(multiColor)-1)) {
      multiColor[[i]] <- as.array(listColors[[i]])
    }
    #  Individual expression sets and their module labels are matched using names of the corresponding components in multiExpr and multiColor
    names(multiColor) <- names(multiExpr) # i.e. =  c(labels, "test")
    
  } else if (length(listDatExpr) > 1) {
    
    if (verbose >= 1) {
      print("Measuring the preservation of a module assignment between a reference and one or more test datasets..")
    }
    
    # Set up the multi-set format for the modulePreservation function
    multiExpr <- vector(mode="list", length=length(listDatExpr))  # convert the list of expression data sets to the required format
    
    for (i in 1:length(multiExpr)) {
      multiExpr[[i]]$data  <- as.matrix(listDatExpr[[i]])
    }
    
    names(multiExpr) <- if (!is.null(names(listDatExpr))) names(listDatExpr) else 1:length(multiExpr) 
    
    multiColor <- listColors # rename for consistency but otherwise leave unchanged
    
    names(multiColor) <- names(multiExpr)[1] # listColors and hence multiColor is only allowed to have length 1
    
    referenceNetworks = rep(1, length(multiExpr)-1) # The reference dataset is always the first
    testNetworks = as.list(2:length(multiExpr)) # the test sets are the second, third, ..., last 
  }
  
  result <- modulePreservation(multiData = multiExpr,
                               multiColor = multiColor,
                               dataIsExpr = dataIsExpr,
                               networkType = networkType, 
                               corFnc = corFnc,
                               corOptions = corOptions,
                               referenceNetworks = referenceNetworks, 
                               testNetworks = testNetworks, 
                               # We actuall just want to get quality stats on the reference network,
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
  return(result)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# load_obj <- function(f) {
#   # Utility function for loading an object inside a new environment and returning it so it can
#   # stored in a variable
#   # EDIT_180421_12
#   #env <- new.env()
#   env <- new.env(parent = emptyenv())
#   ###
#   nm <- load(f, env)[1]
#   # EDIT_180421_12
#   #env[[nm]]
#   out <- env[[nm]]
#   #rm(list=ls(env), envir=env)
#   return(out)
#   ###
# }

load_obj <- function(f) {
  # Utility function for loading an object inside a new environment and returning it so it can
  # stored in a variable
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# EDIT_20180420_6
#parLabel <- function(comb) {
plotLabel_for_vec <- function(comb) {
###
  # Utility function for more easily parallelising making labels
  label = paste0("MMS=", comb[[1]], ",DS=", comb[[2]],",PAM=",comb[[3]], ",CUT=",comb[[4]],sep="") 
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

### 180524 v1.8_dev4
# Now it only converts character or factor columns to dummy

SeuratFactorToIndicator <- function(obj, meta.colnames) {
  # @Usage:   Converts factor/character metadata to new dummy variable columns
  # @args:    obj: Seurat object
  #           meta.colnames: list of metadata column names - character
  # @return:  dataframe with numeric metadata columns, column bound to 'dummy' variable columns, one for each level of each factor or character meta data column in obj
  # @depends: Seurat
  # @author:  Jonatan Thompson jjt3f2188@gmail.com
  # @date:    180524
  # @TODO: 
  idx_factor = sapply(obj@meta.data, function(x) class(x) %in% c("factor", "character"))
  
  meta.data.factor <- obj@meta.data[idx_factor]
  meta.data.factor <- meta.data.factor[colnames(meta.data.factor) %in% meta.colnames]
  meta.colnames.factor <- meta.colnames[meta.colnames %in% colnames(meta.data.factor)]
  
  meta.data.numeric <- obj@meta.data[!idx_factor]
  meta.data.numeric <- meta.data.numeric[colnames(meta.data.numeric) %in% meta.colnames]

  
  # For factor metadata columns, get logical vectors for each level
  list.list.idx = list()
  for (j in seq(length(meta.colnames.factor))){
    col <- grep(meta.colnames.factor[[j]], names(meta.data.factor)) # Find the meta data column index
    list.idx <- lapply(names(table(meta.data.factor[col])), function(x) meta.data.factor[col]==x) # find, for each unique value, a logical vector of occurences
    names(list.idx) = names(table(meta.data.factor[col]))
    list.list.idx[[j]] <- list.idx
  }
  
  flatlist.idx <- unlist(list.list.idx, recursive=F, use.names = T) # Flatten to list of logical vectors  
  #names(flatlist.idx) <- unlist(values, recursive=F) # Assign a feature name to each vector
  new_cols <- cbind(meta.data.numeric, data.frame(lapply(flatlist.idx, function(x) as.numeric(x)), row.names=row.names(obj@meta.data))) # Make numeric dataframe
  #return(AddMetaData(object=obj, new_cols)) # Add dataframe to Seurat object and return
  return(new_cols) # return dataframe, not seurat object
}
###
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

wrapSubset <- function(obj, ident, do.scale=T, do.center=T) {
  # @Status:  OK
  # @Usage:   Pseudo-wrapper for Seurat's SubsetData, ScaleData, FindVariableGenes and runPCA functions
  #           Main purpose: Allows for parallelised subsetting by putting all the steps into one function
  # @args: 
  #           obj: Seurat object from which to take a subset
  #           ident: Vector of characters within ident column on which to subset 
  #           do.scale: handled using ScaleData() rather than SubsetData()
  #           do.center: handled using ScaleData() rather than SubsetData()
  # @return: 
  #           Seurat subset object 
  # @depends: Seurat
  # @author:  Jonatan Thompson jjt3f2188@gmail.com
  # @date:    180326
  
  # Subset the data based on new identity; scale and center the data, restore the old identity
  obj_sub = SubsetData(obj, 
                       ident.use=ident, 
                       do.scale=F, 
                       do.center=F, 
                       subset.raw = T)
  
  obj_sub <- ScaleData(object = obj_sub, 
                       genes.use = NULL, # default: genes.use = all genes in @data
                       vars.to.regress = c("nUMI", "percent.mito", "percent.ribo"), 
                       model.use="linear", 
                       do.scale=do.scale, 
                       do.center=do.center) 
  
  obj_sub <- FindVariableGenes(object = obj_sub, 
                               mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, 
                               x.high.cutoff = 3, 
                               y.cutoff = 0.5, 
                               do.plot=F)
  
  obj_sub <- RunPCA(object = obj_sub, 
                    pc.genes = obj_sub@var.genes, 
                    pcs.compute = nPC_seurat, 
                    do.print = F)#, pcs.print = 1:5, genes.print = 5)
  
  
  return(obj_sub)
  
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

name_for_vec = function(to_be_named, given_names, dimension=NULL){
  # dimension = 1 names the rows, 2 names columns 
  
  if (!is.null(dim(to_be_named))) { # it's a matrix or data frame
    if (dimension == 1) {
      rownames(to_be_named) <- given_names
    } else if (dimension ==2) {
      colnames(to_be_named) <- given_names
    }
  } else if (is.null(dim(to_be_named))) {
    names(to_be_named) <- given_names
  } 
  return(to_be_named)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

PPI_outer_for_vec = function(colors, pkMEs, STRINGdb_species, PPI_pval_threshold, project_dir, data_prefix, flag_date) {
  # Rather than for parallelising over modules within a set of modules, parallelise over several sets of colors (produced by comparing parameters)
  # It calls PPI_innver_for_vec as a subroutine
  # Args:
  #   STRINGdb_species
  #   pkMEs
  #   colors       
  # Returns: 
  #   colors_PPI: colors where modules that did not satisfy the PPI threshold are set to grey
  
  suppressPackageStartupMessages(library(STRINGdb)) 
  
  module_PPI <- NULL
  unique_colors <- NULL
  unique_colors_PPI <- NULL
  colors_PPI <- NULL
  MEs_PPI <- NULL
  kMEs_PPI <- NULL
  pkMEs_PPI <- NULL 
  
  unique_colors <- unique(colors[abs(pkMEs)>PPI_pval_threshold])
  
  string_db <- STRINGdb$new(version="10", species = STRINGdb_species, score_threshold=0, input_directory="") # Default: 10090 (Mus musculus)
  module_PPI <- sapply(unique_colors, function(x) PPI_inner_for_vec(color=x, unique_colors = unique_colors, colors = colors, string_db = string_db), simplify=T) %>% t()
  module_PPI <- as.data.frame(module_PPI, row.names = unique_colors)
  
  # FILTER MODULES ON PPI ENRICHMENT  

  if (!is.null(module_PPI)) {
    unique_colors_PPI = unique_colors[module_PPI$'p-value' < PPI_pval_threshold]
    genes_PPI_idx <- colors %in% unique_colors_PPI
    colors_PPI <- colors
    colors_PPI[!genes_PPI_idx] <- "grey"
  } else {
    colors_PPI <- rep("grey", length(colors))
  }  
  return(colors_PPI)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################


PPI_inner_for_vec <- function(color, unique_colors, colors, string_db) {
  # Parallelise over modules within a set of modules
  # utility function to parallelise the PPI STRINGdb call
  # args: color should be one of the unique_colors
  
  ### 180504_v1.8_dev1
  # Not used
  #i = which(color == unique_colors)
  ###
  ppi <- data.frame(gene = names(colors[colors==color])) # extract the genes with the corresponding color to dataframe
  module_PPI <- list('p-value'=1, 'expected interactions'=NA)
  
  
  # EDIT 180504_v1.8_dev1 
  # We need trycatch here
  tryCatch({
   
  example1_mapped <- string_db$map(ppi, 'gene', removeUnmappedRows = TRUE ) # Check the dataframe genes' PPI. May produce error: we couldn't map to STRING 100% of your identifiers
  # (ctd.) .. Error in if (hitList[i] %in% V(ppi_network)$name) { : argument is of length zero
  hits <- example1_mapped$STRING_id
  module_PPI['p-value'] = p.adjust(string_db$get_ppi_enrichment(hits)$enrichment,method = 'fdr',n = length(names(colors)))
  module_PPI['expected interactions'] = string_db$get_ppi_enrichment(hits)$lambda
  
  
  }, error = function(c) {
     module_PPI <- list('p-value'=1, 'expected interactions' = NA) # probably unnecessary since defined above but just in case it has been changed in place
  }) 
  
  return(module_PPI)
  
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

mapMMtoHs = function(modulekME,
                     log_dir,
                     flag_date,
                     data_prefix,
                     mapping_orthology) {
  
  if (!is.null(modulekME$genes)) modulekME$genes <- NULL
  
  log_not_mapped_filepath = paste0(log_dir,flag_date,"_genes_orthology_not_mapped_",data_prefix,"_", ".tab")
  
  mapping = data.frame(ensembl.mouse=row.names(modulekME))
  # orthology mapping
  mapping$ensembl.human = mapping_orthology$ensembl_gene_id[ match(mapping$ensembl.mouse, mapping_orthology$mmusculus_homolog_ensembl_gene) ]
  #mapping$ensembl.human[mapping$ensembl.human == ""] = NA
  df_not_mapped = mapping[is.na(mapping$ensembl.human),]
  append = !file.exists(log_not_mapped_filepath)
  col.names = !append 
  write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t", col.names = col.names, row.names=F, append=append)
  #modulekME$symbol = mapping$symbol
  modulekME$ensembl = mapping$ensembl.human
  modulekME = na.omit(modulekME)
  
  ### 180508_v.18_dev2
  #tmp = within(modulekME, rm("symbol","ensembl"))
  tmp = within(modulekME, rm("ensembl"))
  ###
  
  # Average duplicated gene IDs
  modulekME_ens <-aggregate(tmp, by=list(modulekME$ensembl),FUN=mean, na.rm=TRUE)
  rownames(modulekME_ens) = modulekME_ens$Group.1
  modulekME_ens = within(modulekME_ens, rm("Group.1"))
  
  modulekME_ens$genes <- NULL
  
  return(modulekME_ens)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

kME_magma <- function(modulekME_name,
                      modulekME,
                      magma_gwas_dir,
                      sub_dir,
                      mapping_hs_entrez2ensembl) {
  
  ### 180522_v1.8_dev3
  # moved to separate function so as to map kME file
  # Initialise log file
  #log_not_mapped_filepath = paste0(log_dir,flag_date,"_genes_orthology_not_mapped_",data_prefix,"_", subsetNames[k], sub_dir,".tab")
  ###
  
  # Load MAGMA genes and remap to Ensembl gene IDs
  d = dir(path=paste0(magma_gwas_dir, sub_dir, "/"), pattern="[.]genes.out", recursive = T)
  gwas = vector(mode="list")
  for(i in 1:length(d)) {
    gwas[[i]] = read.table(paste(magma_gwas_dir, sub_dir, "/", d[[i]],sep=""),head=T, check.names = FALSE)
  }
  names(gwas) = gsub(".genes.out", "", d)
  
  # Match and replace with ENSG 
  #genes = union(gwas[[1]]$GENE, gwas[[2]]$GENE)
  #for(i in 3:length(gwas)) genes = union(genes, gwas[[i]]$GENE)
  
  # Remapping from human Entrez to human Ensembl gene IDs
  for(i in 1:length(gwas)) {
    idx = match(gwas[[i]]$GENE, mapping_hs_entrez2ensembl$entrezgene)
    mapping = data.frame(entrez=gwas[[i]]$GENE, ensembl=mapping_hs_entrez2ensembl$ensembl_gene_id[idx])
    gwas[[i]]$gene_name = mapping$ensembl
  }
  
  
  # Calculate spearman's correlation between gene module membership and GWAS gene significance
  ### 180522
  #colors = colnames(modulekME_ens)
  colors = colnames(modulekME)
  ###
  table.kme.cor.p = table.kme.cor.r <- matrix(NA,nrow=length(unique(colors)),ncol=length(gwas)) 
  rownames(table.kme.cor.r) = rownames(table.kme.cor.p) = unique(colors)
  colnames(table.kme.cor.r) = colnames(table.kme.cor.p) = names(gwas) 
  
  for(m in unique(colors)) {
    for(i in 1:length(gwas)) {
      #col = paste("kME", m, sep="")
      col = m
      genes = intersect(rownames(modulekME),gwas[[i]]$gene_name)
      x = -log10(gwas[[i]]$P[match(genes, gwas[[i]]$gene_name)])
      y = modulekME[match(genes,rownames(modulekME)), col]
      cor = cor.test(x,y,method="spearman", exact=F)
      table.kme.cor.r[m,i] <- cor$estimate
      table.kme.cor.p[m,i] <- cor$p.value
    }
  }
  
  rownames(table.kme.cor.r) <- rownames(table.kme.cor.p) <- paste0(modulekME_name, "__", rownames(table.kme.cor.p))
  #rownames(table.kme.cor.r) <- rownames(table.kme.cor.p) <- paste0("subset_", k, "_", rownames(table.kme.cor.p))

  return(list('p.val'= table.kme.cor.p, 'corrCoef' = table.kme.cor.r))  

}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

makeColorsUniqueForPlot = function(cols) {
  while (all(!duplicated(cols)) == F) {
    for (i in which(duplicated(cols))) {
      newColor <- "newcolor"
      k=1
      while (!newColor %in% colors()) {
        if (k<100) {
          newColor <- paste0(cols[i], sample.int(4, 1, replace=T))
        } else if (k == 100) { # stop trying and just plot with grey
          newColor <- paste0("grey", sample.int(50:80, 1, replace=T))
        }
        k = k+1
      }
      cols[i] <- newColor
      #if (!htca_BMI_fdr_colors[i] %in% colors()) htca_BMI_fdr_colors[i] <- colors()[agrep(pattern = htca_BMI_fdr_colors[i], x = colors(), max.distance = c("deletions"=1))[1]]
    }
  }
  return(cols)
}
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

magma_par = function(subsetNames, 
                     list_kMEs,
                     project_dir, 
                     plots_dir, 
                     tables_dir, 
                     log_dir, 
                     n_cores,
                     magma_gwas_dir, 
                     data_prefix, 
                     file_suffix, 
                     flag_date, 
                     organism,
                     gwas_filter_traits = NULL) {
  
  # Usage: Loop over immediate subdirectories in magma_gwas_dir, scoring modules on each
  # Args: 
  # Returns:
  #   For each magma_gwas_dir immediate subdirectory, a table with module scores across all cell types 

  # MAP MOUSE TO HUMAN GENES
  if (organism == "mmusculus") {
    
    mapping_hs_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz" # to map mouse to human ensembl
    mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T)
    
    cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "mapMMMtoHs_par.txt"))
    list_kMEs_hs <- parLapply(cl,list_kMEs, function(x) mapMMtoHs(modulekME = x, 
                                                            log_dir = log_dir, 
                                                            flag_date = flag_date, 
                                                            data_prefix = data_prefix, 
                                                            mapping_orthology = mapping_orthology))
    stopCluster(cl)
    invisible(gc())
    names(list_kMEs_hs) <- names(list_kMEs)
    #rm(list_kMEs_PPI_ok)
    
  } else if (organism == "hsapiens") {
    list_kMEs_hs <- list_kMEs
    names(list_kMEs_hs) <- names(list_kMEs)
    #rm(list_kMEs_PPI_ok)
  }
  
  options(stringsAsFactors = F)
  file_sep = ','
  
  # Set paths to mapping files to variables (these are independent of arguments)
  mapping_hs_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz" # columns: ensembl_gene_id, entrezgene, hgcn_symbol. Use this for mapping entrezgene to ensembl
  mapping_hs_entrez2ensembl = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T)
  
  ### 180522_v1.8_dev3
  # moved to separate function so as to map kME file
  # mapping_hs_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz" # to map mouse to human ensembl
  # mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T)
  ###

  gwas_sub_dirs = list.dirs(path = magma_gwas_dir, full.names = FALSE, recursive = FALSE)
  
  for (sub_dir in gwas_sub_dirs) {
    
    cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "kME_magma_par.txt"))
    magma_results <- clusterMap(cl, function(x,y) kME_magma(modulekME_name = x, 
                                                            modulekME = y,
                                                            magma_gwas_dir = magma_gwas_dir,
                                                            sub_dir = sub_dir,
                                                            mapping_hs_entrez2ensembl = mapping_hs_entrez2ensembl),
                                x = names(list_kMEs_hs),
                                y = list_kMEs_hs,
                                SIMPLIFY=F)
    stopCluster(cl)
    invisible(gc())
    # for (k in 1:length(list_kMEs_hs)) { 
    #   modulekME <- list_kMEs_hs[[k]]
    # 
    #   ### 180522_v1.8_dev3
    #   # moved to separate function so as to map kME file
    #   # Initialise log file
    #   #log_not_mapped_filepath = paste0(log_dir,flag_date,"_genes_orthology_not_mapped_",data_prefix,"_", subsetNames[k], sub_dir,".tab")
    #   ###
    #   
    #   # Load MAGMA genes and remap to Ensembl gene IDs
    #   d = dir(path=paste0(magma_gwas_dir, sub_dir, "/"), pattern="[.]genes.out", recursive = T)
    #   gwas = vector(mode="list")
    #   for(i in 1:length(d)) {
    #     gwas[[i]] = read.table(paste(magma_gwas_dir, sub_dir, "/", d[[i]],sep=""),head=T, check.names = FALSE)
    #   }
    #   names(gwas) = gsub(".genes.out", "", d)
    #   
    #   # Match and replace with ENSG 
    #   #genes = union(gwas[[1]]$GENE, gwas[[2]]$GENE)
    #   #for(i in 3:length(gwas)) genes = union(genes, gwas[[i]]$GENE)
    #   
    #   # Remapping from human Entrez to human Ensembl gene IDs
    #   for(i in 1:length(gwas)) {
    #     idx = match(gwas[[i]]$GENE, mapping_hs_entrez2ensembl$entrezgene)
    #     mapping = data.frame(entrez=gwas[[i]]$GENE, ensembl=mapping_hs_entrez2ensembl$ensembl_gene_id[idx])
    #     gwas[[i]]$gene_name = mapping$ensembl
    #   }
    #   
    #   
    #   # Calculate spearman's correlation between gene module membership and GWAS gene significance
    #   ### 180522
    #   #colors = colnames(modulekME_ens)
    #   colors = colnames(modulekME)
    #   ###
    #   table.kme.cor.p = table.kme.cor.r <- matrix(NA,nrow=length(unique(colors)),ncol=length(gwas)) 
    #   rownames(table.kme.cor.r) = rownames(table.kme.cor.p) = unique(colors)
    #   colnames(table.kme.cor.r) = colnames(table.kme.cor.p) = names(gwas) 
    #   
    #   for(m in unique(colors)) {
    #     for(i in 1:length(gwas)) {
    #       #col = paste("kME", m, sep="")
    #       col = m
    #       genes = intersect(rownames(modulekME),gwas[[i]]$gene_name)
    #       x = -log10(gwas[[i]]$P[match(genes, gwas[[i]]$gene_name)])
    #       y = modulekME[match(genes,rownames(modulekME)), col]
    #       cor = cor.test(x,y,method="spearman", exact=F)
    #       table.kme.cor.r[m,i] <- cor$estimate
    #       table.kme.cor.p[m,i] <- cor$p.value
    #     }
    #   }
    # 
    #   rownames(table.kme.cor.r) <- rownames(table.kme.cor.p) <- paste0(names(list_kMEs_hs)[k], "_", rownames(table.kme.cor.p))
    #   #rownames(table.kme.cor.r) <- rownames(table.kme.cor.p) <- paste0("subset_", k, "_", rownames(table.kme.cor.p))
    #   list_table.kme.cor.r[[k]] <- table.kme.cor.r
    #   list_table.kme.cor.p[[k]] <- table.kme.cor.p
    #   
    # } # end of list_kMEs_hs loop
    
    # merge each list of tables into one big table spanning all cell clusters
    # [deleted]
    
    list_table.kme.cor.r <- list_table.kme.cor.p <- vector(mode="list", length = length(list_kMEs_hs))
    
    # Extract the coefficient and p-value dataframes for each module, put them into lists
    for (i in 1:length(magma_results)) { 
      list_table.kme.cor.p[[i]] <- magma_results[[i]][['p.val']]
      list_table.kme.cor.r[[i]] <- magma_results[[i]][['corrCoef']]
    }
    
    # Merge each list
    for (t in 1:length(list_table.kme.cor.r)) {
      if (t == 1) {
        table.kme.cor.r_all <- list_table.kme.cor.r[[t]]
      } else {
        table.kme.cor.r_all <- rbind(table.kme.cor.r_all, list_table.kme.cor.r[[t]] )
      }
    }
    
    for (t in 1:length(list_table.kme.cor.p)) {
      if (t == 1) {
        table.kme.cor.p_all <- list_table.kme.cor.p[[t]]
      } else {
        table.kme.cor.p_all <- rbind(table.kme.cor.p_all, list_table.kme.cor.p[[t]] )
      }
    }
    
    table.kme.cor.p.fdr_all = p.adjust(table.kme.cor.p_all, method="fdr")
    dim(table.kme.cor.p.fdr_all) = dim(table.kme.cor.p_all);  dimnames(table.kme.cor.p.fdr_all) = dimnames(table.kme.cor.p_all)
    
    #d = -log10(table.kme.cor.p.fdr) * sign(table.kme.cor.r)  # we don't use this!
    #pdf("SampleGraph.pdf",width=7,height=5)
    #sizeGrWindow(9,7)
    #par(mfrow = c(2,2))
    #par(mar = c(4, 5, 4, 6));
    #labeledHeatmap(d,textMatrix = signif(table.kme.cor.r,1), xLabels = colnames(d), yLabels = rownames(d),invertColors = T, colors = blueWhiteRed(1000), main="GWAS - kME correlation", cex.text = 0.6)
    #dev.off()
    
    dat = as.data.frame(table.kme.cor.p.fdr_all*sign(table.kme.cor.r_all))
    dat[dat<0]=1 #Only look for positive enrichment
    dat = -log10(dat)
    dat$module = gsub("kME|kME\\.","",rownames(dat))

    ### 180523 v1.8_dev4
    #This was missing
    #table.kme.cor.r_all <- cbind(module=dat$module, table.kme.cor.r_all)
    ###
    # Write full fdr and coefficient tables to csv
    write.csv(dat, file=paste0(tables_dir, data_prefix, "_", file_suffix, "_magma_GWAS_fdr_all_modules_all_traits", sub_dir, "_", flag_date, ".csv"), quote=F, row.names = F)
    write.csv(table.kme.cor.r_all, file=paste0(tables_dir, data_prefix, "_", file_suffix, "_magma_GWAS_corrCoef_all_modules_all_traits", sub_dir, "_", flag_date,".csv"), quote=F, row.names = F)

    #table.kme.cor.r_all <- as.data.frame(table.kme.cor.r_all,row.names = rownames(table.kme.cor.r_all))
    #table.kme.cor.p_all <- as.data.frame(table.kme.cor.p_all, row.names = rownames(table.kme.cor.p_all))
  
    ### v1.8_dev4
    # filter on user specified gwas studies
    # dat_BMI <-  dat[c("bmi_height_yengo_2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz", "module")] [dat[["bmi_height_yengo_2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz"]] > -log10(0.05),]
    # rownames(dat_BMI) <- NULL
    # r_BMI <- cbind(table.kme.cor.r_all["bmi_height_yengo_2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz"], module=dat$module)
    # r_BMI <-r_BMI[dat[["bmi_height_yengo_2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz"]] > -log10(0.05),]
    # rownames(r_BMI) <- NULL
    
    if (!is.null(gwas_filter_traits)) {
      if (sapply(gwas_filter_traits, function(x) any(grepl(x, colnames(dat[,-grep("module", colnames(dat))]), ignore.case=T)), simplify = T) %>% all) {
      # filter columns
      idx_col_include <- sapply(gwas_filter_traits, function(x) grepl(x, colnames(table.kme.cor.r_all), ignore.case=T), simplify = T) %>% rowSums %>% as.logical 
      dat_filter_traits <- dat[,-grep("module", colnames(dat))][,idx_col_include]
      dat_filter_traits <- cbind(module = dat$module, dat_filter_traits)
      # Filter rows
      idx_row_include <- apply(dat_filter_traits[,-grep("module", colnames(dat_filter_traits))], MARGIN = 1, max) > -log10(0.05)
      dat_filter_traits <- dat_filter_traits[idx_row_include,]
      # Do the same for the correlation coefficients table
      r_filter_traits <- table.kme.cor.r_all[idx_row_include,idx_col_include]
      rownames(r_filter_traits) <- NULL
      } else {
        dat_filter_traits <- dat
        rownames(dat_filter_traits) <- NULL
        r_filter_traits <- cbind(table.kme.cor.r_all, module=dat$module)
        rownames(r_filter_traits) <- NULL
      }
    } else {
      dat_filter_traits <- dat
      rownames(dat_filter_traits) <- NULL
      r_filter_traits <- cbind(table.kme.cor.r_all, module=dat$module)
      rownames(r_filter_traits) <- NULL
    }
    ###
    
    #dat2 = melt(dat)
    #dat2$variable = as.character(dat2$variable)
    # #p=ggplot(melt(dat),aes(x=variable,y=value,fill=colors)) + 

    # prepare colors for plot: replace duplicate colors 
    ### 180528_v.18_dev4
    # No point in doing this for the gwas plot for all modules from all clusters since the user won't be able to see the colors anyway.
    # There is a high risk of unresolvable duplicates.
    
    # barColorsAll <- gsub(".*__", "", dat$module)
    # 
    # while (all(!duplicated(barColorsAll)) == F) {
    #   for (i in which(duplicated(barColorsAll))) {
    #     newColor <- "newcolor"
    #     k=1
    #     while (!newColor %in% colors()) {
    #       if (k<100) {
    #         newColor <- paste0(barColorsAll[i], sample.int(4, 1, replace=T))
    #       } else if (k == 100) { # stop trying and just plot with grey
    #         newColor <- paste0("grey", sample.int(50:80, 1, replace=T))
    #       }
    #       k = k+1
    #     }
    #     barColorsAll[i] <- newColor
    #     #if (!htca_BMI_fdr_colors[i] %in% colors()) htca_BMI_fdr_colors[i] <- colors()[agrep(pattern = htca_BMI_fdr_colors[i], x = colors(), max.distance = c("deletions"=1))[1]]
    #   }
    # }

    ###
    
    barColorsFilter <- gsub(".*_", "", dat_filter_traits$module)
    
    barColorsFilter <- makeColorsUniqueForPlot(barColorsFilter)
    ### 180528 moved into a function
    # while (all(!duplicated(barColorsFilter)) == F) {
    #   for (i in which(duplicated(barColorsFilter))) {
    #     newColor <- "newcolor"
    #     k=1
    #     while (!newColor %in% colors()) {
    #       if (k<100) {
    #         newColor <- paste0(barColorsFilter[i], sample.int(4, 1, replace=T))
    #       } else if (k == 100) { # stop trying and just plot with grey
    #         newColor <- paste0("grey", sample.int(50:80, 1, replace=T))
    #       }
    #       k = k+1
    #     }
    #     barColorsFilter[i] <- newColor
    #     #if (!htca_BMI_fdr_colors[i] %in% colors()) htca_BMI_fdr_colors[i] <- colors()[agrep(pattern = htca_BMI_fdr_colors[i], x = colors(), max.distance = c("deletions"=1))[1]]
    #   }
    # }
    ### 180528 Don't really need this
    
    # Plot all modules, all GWAS
    # invisible(p=ggplot(melt(dat),aes(x=variable,y=value, fill = module)) +
    #   geom_bar(stat="identity",position=position_dodge(),color= "black") +
    #   scale_fill_manual(values=unique(barColorsAll)) + theme_classic() +
    #   geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
    #   theme(axis.text.x=element_text(angle=50, size=10, hjust=1)))
    # p=ggplot(melt(dat),aes(x=variable,y=value, fill = module)) +
    #             geom_bar(stat="identity",position=position_dodge(),color= "black") +
    #             #scale_fill_manual(values=unique(barColorsAll)) + theme_classic() +
    #             scale_fill_manual(values=colors()) + theme_classic() +
    #             geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
    #             theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
    # ggsave(p, filename = paste0(plots_dir, data_prefix, "_", file_suffix, "_magma_all_modules_all_GWAS_", sub_dir, "_", flag_date, ".pdf") ,width=45,height=12)
     
    if (dim(dat_filter_traits)[1] > 0) {
      # Plot significant modules, only on selected gwas
      p=ggplot(melt(dat_filter_traits, value.name="value"),aes(x=variable,y=value, fill = module)) + 
        geom_bar(stat="identity",position=position_dodge(),color= "black") +
        scale_fill_manual(values=unique(barColorsFilter)) + theme_classic() +
        geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
        theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
      
      ggsave(p, filename = paste0(plots_dir, data_prefix, "_", file_suffix, "_magma_modules_gwas_filter_", sub_dir, "_", flag_date, ".pdf") ,width=45,height=12)
      
      return(list("fdr"=dat_filter_traits, "corrCoef"=r_filter_traits)) 
      
      } else {
      
      return(NULL)
    }  
    
  }
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# TODO: Should we test the gwas signal modules or all PPI modules?

mendelianGenes <- function() {
  # Load coding variant and Mendelian genes
  coding_variants <- read.table("/data/pub-others/turcot-biorxiv-2018/s21/s21.hgnc")
  mendelian <- read.table("/data/pub-others/turcot-biorxiv-2018/table1/table1.hgnc")
  
  # Load table of kME values
  modulekME = read.csv("~/ygg-projects/mludwig/WGCNA/Campbell/Results/agrp_kmemodule.csv", row.names=1, check.names = FALSE, sep = ",")
  modulekME <- kme 
  # Load table of module assignments
  load("~/ygg-projects/mludwig/WGCNA/Genetics/Data/Campbell_AgRP_WGCNA.RData")
  
  # Gene lists
  gene.lists = vector(mode="list")
  gene.lists[[1]] = coding_variants[,1]
  gene.lists[[2]] = mendelian[,1]
  names(gene.lists)= c("coding_variants", "mendelian")
  
  
  # Change gene names to upper case
  rownames(c_modules) <- toupper(rownames(c_modules))
  
  colors <- c_modules[,1]
  
  table.p = matrix(NA, nrow=length(unique(colors)), ncol=length(gene.lists))
  rownames(table.p) = unique(colors)
  colnames(table.p) = names(gene.lists) 
  table.or = logit.p = logit.or = table.p
  
  for(i in 1:ncol(table.or)) {
    for(j in 1:nrow(table.or)) {
      col = rownames(table.or)[j]
      
      binaryMat = as.data.frame(cbind(as.numeric(colors==col), 0))
      colnames(binaryMat) = c("Module", "GeneSet")
      idx = match(gene.lists[[i]], rownames(c_modules))
      binaryMat$GeneSet[idx] = 1
      
      glm.out <- glm(binaryMat$GeneSet~binaryMat$Module, family=binomial())
      summary(glm.out)
      
      logit.or[j,i] = exp(coefficients(glm.out)[2])  # Calculate odds ratio from logistic regression
      logit.p[j,i] = summary(glm.out)$coefficients[2,4]
    }
  }
  
  table.p.fdr = p.adjust(logit.p,method="fdr")
  dim(table.p.fdr) = dim(logit.p)
  dimnames(table.p.fdr) = dimnames(logit.p)
  table.p.fdr
  logit.p
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# NOT IN USE 
# 180522: /!\ can't install DOSE

# GO_for_par = function(list_kMEs_gwas) {
#   
#   suppressPackageStartupMessages(library(DOSE))
#   suppressPackageStartupMessages(library(GO.db))
#   suppressPackageStartupMessages(library(org.Hs.eg.db))
#   suppressPackageStartupMessages(library(GSEABase))
#   suppressPackageStartupMessages(library(clusterProfiler))
#   
#   ### Only for PPI
#   
#   list_GSEA_gwas <- NULL
#   
#   # # select only kME columns of gwas enriched module
#   # list_kMEs_gwas <- mapply(function(x,y) x[,colnames(x) %in% y], x = list_kMEs_PPI_ok_gwas, y = list_modules_PPI_gwas)
#    
#   # Order genes by kME to their own module. This allows us to submit them as ranked queries to gprofiler for GSEA style p-values (https://cran.r-project.org/web/packages/gProfileR/gProfileR.pdf)
#   list_list_module_gwas_genes_order <- lapply(list_kMEs_gwas, function(x) lapply(colnames(x), function(y) rownames(x)[order(x[[y]], decreasing=T)]))
#   
#   ### 190509_v1.8.dev2
#   ### # Got to here!
#   invisible(gc())
#   cl <- makeCluster(n_cores, type = "FORK", 
#                     outfile = paste0(log_dir, "log_gprofiler.txt"))
#   
#   list_list_ggo <- parLapply(list_list_module_gwas_genes_order, function(x) lapply(x, function(y) groupGO(gene = y,
#                                                                                                           OrgDb    = org.Mm.eg.db,
#                                                                                                           ont      = "CC",
#                                                                                                           #level    = 3,
#                                                                                                           readable = T)))
#   
#   
#   ego3 <- parLapply(list_list_module_gwas_genes_order, function(x) lapply(x, function(y) gseGO(geneList = y,
#                                                                                                OrgDb        = org.Mm.eg.db,
#                                                                                                ont          = "CC",
#                                                                                                nPerm        = 1000,
#                                                                                                minGSSize    = 100,
#                                                                                                maxGSSize    = 500,
#                                                                                                pvalueCutoff = 0.05,
#                                                                                                verbose      = FALSE)))
#   stopCluster(cl)
#   invisible(gc())
# }

# list_list_gprofiles <- lapply( list_list_module_gwas_genes_order, function(x) gprofiler(query=x, 
#                                                                                        organism = organism, 
#                                                                                        sort_by_structure = T,
#                                                                                        ordered_query = T, 
#                                                                                        significant = T, 
#                                                                                        exclude_iea = T, # TODO check with Dylan
#                                                                                        underrep = F,
#                                                                                        evcodes = F, 
#                                                                                        region_query = F, 
#                                                                                        max_p_value = 1, 
#                                                                                        min_set_size = 0,
#                                                                                        max_set_size = 0, 
#                                                                                        min_isect_size = 0, 
#                                                                                        correction_method = "analytical",
#                                                                                        hier_filtering = "none", 
#                                                                                        domain_size = "annotated", 
#                                                                                        custom_bg = "",
#                                                                                        numeric_ns = "", 
#                                                                                        png_fn = NULL, 
#                                                                                        include_graph = F, 
#                                                                                        src_filter = NULL))


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
# NOT IN USE
# GSEA_for_par = function() {
# 
#   list_GSEA_PPI_gwas <- NULL
#   
#   # select only kME columns of gwas enriched moduled 
#   list_kMEs_PPI_ok_gwas <- mapply(function(x,y) x[,colnames(x) %in% y], x = list_kMEs_ok[names(list_kMEs_PPI_ok) %in% sNames_PPI_gwas], y = list_modules_PPI_gwas)
#   
#   # Order genes by kME to their own module. This allows us to submit them as ranked queries to gprofiler for GSEA style p-values (https://cran.r-project.org/web/packages/gProfileR/gProfileR.pdf)
#   list_list_module_PPI_gwas_genes_order <- lapply(list_kMEs_PPI_ok_gwas, function(x) lapply(colnames(x), function(y) rownames(x)[order(x[[y]], decreasing=T)]))
#   
#   # name them
# 
#   invisible(gc())
#   cl <- makeCluster(n_cores, type = "FORK", 
#                     outfile = paste0(log_dir, "log_gprofiler_PPI.txt"))
#   
#   ### 180607_v1.8_dev2
#   list_list_gprofiles_PPI <- parLapply(cl, list_list_module_PPI_genes_order, function(x) gprofiler(query=x, 
#                                                                                                    organism = organism, 
#                                                                                                    sort_by_structure = T,
#                                                                                                    ordered_query = T, 
#                                                                                                    significant = T, 
#                                                                                                    exclude_iea = F, # TODO check with Dylan
#                                                                                                    underrep = F,
#                                                                                                    evcodes = F, 
#                                                                                                    region_query = F, 
#                                                                                                    max_p_value = 1, 
#                                                                                                    min_set_size = 0,
#                                                                                                    max_set_size = 0, 
#                                                                                                    min_isect_size = 0, 
#                                                                                                    correction_method = "analytical",
#                                                                                                    hier_filtering = "none", 
#                                                                                                    domain_size = "annotated", 
#                                                                                                    custom_bg = "",
#                                                                                                    numeric_ns = "", 
#                                                                                                    png_fn = NULL, 
#                                                                                                    include_graph = F, 
#                                                                                                    src_filter = NULL))
#   
#   stopCluster(cl)
#   invisible(gc())
#   
#   
#   ### edit_180504_v1.7_dev1
#   # Assign names to each list with gprofiler dataframes
#   #list_list_gprofiles_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names= colnames(y), dimension = NULL),  x=list_list_gprofiles_PPI,  y=list_kIMs_PPI, SIMPLIFY=F)
#   list_list_gprofiles_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = NULL),  x=list_list_gprofiles_PPI,  y=list_list_module_PPI_genes_order, SIMPLIFY=F)
#   ###
#   #names(list_list_gprofiles_PPI) <- sNames_ok_PPIfilter
# }
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

make_eigen_mat <- function(RObjects_dir,
                        list_list_module_genes,
                        n_cores,
                        log_dir) {

  sigmods<-list_list_module_genes
  
  sigmod_list<-lapply(sigmods,function(x) melt.list(x))
  sigmod_list %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="value"), .)->gene_modules.df
  gene_modules.df<-as.data.frame(gene_modules.df)
  #gene_modules.df["L1"] <- NULL
  gene_modules.df<-sapply(gene_modules.df, as.character)
  gene_modules.df<-as.data.frame(gene_modules.df)
  colnames(gene_modules.df)<-c("genes", names(sigmod_list))
  gene_modules.df[is.na(gene_modules.df)]<-'random'
  row.names(gene_modules.df) <-  gene_modules.df$genes
  gene_modules.df<-gene_modules.df[,-1, drop=F]
  
  # FILTER SEURATOBJ 
  
  # load seurat object
  seurat_obj <- load_obj(f=sprintf("%sseurat_obj_ensembl.RData",RObjects_dir))
  # Get ident
  ident <- seurat_obj@ident
  # filter
  s.coexp <- t(seurat_obj@data[rownames(seurat_obj@data) %in% row.names(gene_modules.df),]) 
  s.coexp <- as.matrix(s.coexp)
  # Free up memory
  rm(seurat_obj)
  # order genes and modules same
  s.coexp <- s.coexp[,order(colnames(s.coexp))]
  gene_modules.df<-gene_modules.df[order(row.names(gene_modules.df)),]
  identical(colnames(s.coexp),row.names(gene_modules.df))
  eiglist<-as.list(as.data.frame(gene_modules.df))
  
  # Score cells on eigengenes 
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_eigengene_allscore.txt"))
  eigenge_allscore.list<-parLapply(cl, eiglist, function(x) moduleEigengenes((s.coexp), x, impute = T)$eigengenes) #Score cell matrix on each eigengene found in each celltype
  stopCluster(cl)
  invisible(gc())
  
  #Rename Columns in each eigengene matrix
  for (i in 1:length(eigenge_allscore.list)){ 
    colnames(eigenge_allscore.list[[i]])<-paste0(names(eigenge_allscore.list)[[i]],'_',colnames(eigenge_allscore.list[[i]])) 
  }
  
  eigen_mat <- bind_cols(eigenge_allscore.list) # Merge Eigengene Matrices
  eigen_mat <- eigen_mat[,!grepl('.*merandom$',colnames(eigen_mat), ignore.case = T)] #Remove eigengene scores on genes that were not found in modules in a cell type
  row.names(eigen_mat)<-row.names(s.coexp) 

  # Reorder rows in eigenmat by cell type
  for (i in 1:length(unique(ident))) {
    if (i==1) {
      eigen_mat_order <- eigen_mat[ident== names(table(ident))[i],]
    } else {
      eigen_mat_order <- rbind(eigen_mat_order, eigen_mat[ident== names(table(ident))[i],])
    }
  }
  ### 180528 prefixing the cell cluster is redundant
  #prefix <- rep(x=names(list_list_module_genes), times=sapply(list_list_module_genes, length))
  colnames(eigen_mat) <- gsub("gene_modules.df|ME", "", colnames(eigen_mat))
  #colnames(eigen_mat) <- paste0(prefix,"_", colnames(eigen_mat))
  return(eigen_mat)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# WORK IN PROGRESS

plot_eigen_mat = function(eigen_mat,
                          RObjects_dir,
                          list_modules,
                          #dat_fdr,
                          #dat_corrCoef,
                          plots_dir, 
                          data_prefix, 
                          flag_date) {
  
  # load libraries within the function environment to reduce the risk of 'maximum DLLs reached' error later
  suppressPackageStartupMessages(library(ComplexHeatmap))

  # load seurat object
  seurat_obj <- load_obj(f=sprintf("%sseurat_obj_ensembl.RData",RObjects_dir))
  # Get ident
  ident <- seurat_obj@ident
  
  # make row annotation
  htra_colors <- sample(colors()[c(1:151,362:length(colors()))], size=length(table(ident)), replace=F) 
  names(htra_colors) <- names(table(ident))
  
  htra<-rowAnnotation(cellcluster = ident,
                      annotation_legend_param = list(cellcluster = list(nrow = 7, title = "Cell cluster", title_position = "topcenter")),
                      col=list(cellcluster = htra_colors),
                      width = unit(5, "mm"))

  # Make column annotations
  # make non-duplicate colors corresponding to modules
  # htca_colors <- unlist(list_modules, use.names =T)
  # names(htca_colors) <- colnames(eigen_mat)
  # 
  # # replace duplicate colors with related ones
  # while (all(!duplicated(htca_colors)) == F) {
  #   for (i in which(duplicated(htca_colors))) {
  #     newColor <- "newcolor"
  #     k=1
  #     while (!newColor %in% colors()) {
  #       if (k<100) {
  #         newColor <- paste0(htca_colors[i], sample.int(4, 1, replace=T))
  #       } else if (k==100) {
  #         newColor <- paste0("grey", sample.int(99, 1, replace=T))        
  #       }
  #       k = k+1
  #     }
  #     htca_colors[i] <- newColor
  #  }
  # }
  

  #colAnnodf <- data.frame(dat_fdr[-grep("module", colnames(dat_fdr))])
                          #dat_corrCoef[-grep("module", colnames(dat_corrCoef))])

  #module = htca_colors,
  
  #rownames(colAnnodf) <- dat_fdr$module
  
  ### 180524
  
  #$colnames(colAnnodf[-grep("module", colnames(dat_fdr))])
  # paste0(colnames(colAnnodf[-grep("module", colnames(dat_fdr))])[j], "_fdr")
  
  # idx_col <- 1:(ncol(colAnnodf)-1)
  # suffix <- c(rep(",", ncol(colAnnodf)-2),"")
  # htca <- HeatmapAnnotation(eval(parse(text = paste0("fdr_", idx_col , "= anno_barplot(colAnnodf[[", idx_col, "]], axis=T, axis_side='right', baseline = 0, gp = gpar(fill = htca_colors)), height = unit(50, 'mm'), show_annotation_name = TRUE,  annotation_name_offset = unit(10, 'mm'), annotation_name_side = 'right', annotation_name_rot = c(0, 0, 90))",suffix))))

  # plot heatmap without column annotations
  
  tryCatch({   
  # Without clustering rows
  ht1 <- Heatmap(eigen_mat[,], 
               cluster_rows = F,
               #row_order = row.names(eigen_mat_order),
               cluster_columns = T, 
               #show_row_dend = F, 
               show_column_dend = F, 
               show_heatmap_legend = FALSE, 
               show_row_names = F, 
               show_column_names = T)
               #top_annotation = eval(parse(text = paste0("list_htca_fdr", "[[", 1:length(list_htca_fdr),"]]"))))
  pdf(sprintf("%s%s_AllEigenHeatmap_rows_celltype_order_%s.pdf", plots_dir, data_prefix, flag_date ),h=20,w=20)
  draw(ht1+htra)
  dev.off()
  
  # With clustering rows
  ht1 <- Heatmap(eigen_mat[,], 
                 cluster_rows = T,
                 #row_order = row.names(eigen_mat_order),
                 cluster_columns = T, 
                 #show_row_dend = F, 
                 show_column_dend = F, 
                 show_heatmap_legend = FALSE, 
                 show_row_names = F, 
                 show_column_names = T)
                 #top_annotation = eval(parse(text = paste0("list_htca_fdr", "[[", 1:length(list_htca_fdr),"]]"))))
  pdf(sprintf("%s%s_AllEigenHeatmap_rows_cluster_%s.pdf", plots_dir, data_prefix, flag_date ),h=20,w=20)
  draw(ht1+htra)
  dev.off()
  
  }, error = function(x) warning("Heatmap plot failed. Maybe only one module?"))
 return(htra_colors)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# TODO: Reintroduce average cell type gene expression heatmap?
if (FALSE) {
  mean.celltype.eig.expr <- t(sapply(names(table(meta$neuron_groups)), function(x) rowMeans(t(eigen_mat[meta$neuron_groups==x,])), simplify=T))


  ht.mean.celltype.eig.expr <- Heatmap(mean.celltype.eig.expr, 
               cluster_rows = T, 
               cluster_columns = T, 
               show_row_dend = F, 
               show_column_dend = F, 
               show_heatmap_legend = FALSE, 
               show_row_names = T,
               row_names_side = "left",
               show_column_names = T)
               #top_annotation = htca)
  #bottom_annotation = htca)
  pdf(sprintf("%s%s_mean.celltype.eig.expr_%s.pdf", plots_dir, data_prefix, flag_date ),h=12,w=12)
  draw(ht.mean.celltype.eig.expr)
  dev.off()
}