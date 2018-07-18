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

wrapJackStraw = function(seurat_obj_sub, n_cores, jackstrawnReplicate, pvalThreshold) {
  # gene.criterion: 'p.val' means selecting genes (used for PCA) with significant empirical p-val
  #                 'PC.loadings' means projecting all genes onto the PCs to get loadings and selecting 
  #                 genes that have a high absolute loading on a significant PC
  # Schematic:
  # 1. pcs.compute = ncol

  ### 180507
  prop.freq <- max(0.016, round(4/length(seurat_obj_sub@var.genes),3)) # to ensure we have at least 3 samples so the algorithm works well
  # see https://github.com/satijalab/seurat/issues/5
  ###
  
  pcs.compute = ncol(seurat_obj_sub@dr$pca@gene.loadings)
  
  if (jackstrawnReplicate > 0) {
    
    seurat_obj_sub <- JackStraw(object = seurat_obj_sub,
                                num.pc = pcs.compute,
                                num.replicate = jackstrawnReplicate, 
                                display.progress = T,
                                do.par = T,
                                num.cores = n_cores,
                                prop.freq = prop.freq) # https://github.com/satijalab/seurat/issues/5
  
  
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
        x = c(length(x = which(x = pAll[, i] <= pvalThreshold)), floor(x = nrow(x = pAll) * pvalThreshold)),
        n = c(nrow(pAll), nrow(pAll)))$p.val)
      if (length(x = which(x = pAll[, i] <= pvalThreshold)) == 0) {
        pc.score <- 1
      }
      if (is.null(x = score.df)) {
        score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
      } else {
        score.df <- rbind(score.df, data.frame(PC = paste0("PC",i), Score = pc.score))
      }
    }
    
    PC_select_idx <- which(score.df$Score < pvalThreshold)
    
    if (nrow(pAll) == length(seurat_obj_sub@var.genes)) {
      
      seurat_obj_sub <- ProjectPCA(seurat_obj_sub, 
                                   do.print = F, 
                                   pcs.print = NULL, 
                                   pcs.store = pcs.compute, 
                                   genes.print = NULL, 
                                   replace.pc=F, 
                                   do.center=T)
        
      loadings <- abs(seurat_obj_sub@dr$pca@gene.loadings.full[,PC_select_idx, drop=F])
      max_loadings <- apply(loadings, 1, function(x) max(x))
      names_genes_use <- names(max_loadings[order(max_loadings, decreasing = T)])[1:5000]
      
    } else if (nrow(pAll) > length(seurat_obj_sub@var.genes)) {
      
      pAll[,sapply(pAll, function(x) class(x)!="numeric")] <- NULL # remove the column of gene names
      row_min <- apply(pAll[,PC_select_idx, drop=F], MARGIN = 1, FUN = function(x) min(x))
      names_genes_use <- rownames(pAll)[row_min < pvalThreshold]
      
      if (length(names_genes_use) < 1000) names_genes_use <- rownames(pAll)[row_min < pvalThreshold*2]
    }
    
  } else if (jackstrawnReplicate == 0) {
    
    seurat_obj_sub <- ProjectPCA(seurat_obj_sub, 
                                 do.print = F, 
                                 pcs.print = NULL, 
                                 pcs.store = pcs.compute, 
                                 genes.print = NULL, 
                                 replace.pc=F, 
                                 do.center=T)
    
    loadings <- abs(seurat_obj_sub@dr$pca@gene.loadings.full[,1:(min(13,pcs.compute))])
    max_loadings <- apply(loadings, MARGIN=1, FUN=function(x) max(x))
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
  endRunIndex = if (replace==T) nPermutations+1 else nPermutations 
  result = vector("list", length= if (replace==T) nPermutations+1 else nPermutations);
  nSamples = nrow(datExpr);
  nGenes = ncol(datExpr);
  
  try(
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
  })
  
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

cutreeHybrid_for_vec <- function(comb, geneTree, dissTOM, maxPamDist, useMedoids) {
  # Utility function for more easily parallellising the cutreeHybrid function
  tree = cutreeHybrid(dendro = geneTree, 
                      cutHeight = NULL,
                      minClusterSize = comb[[1]], 
                      distM=as.matrix(dissTOM),
                      deepSplit=comb[[2]], 
                      pamStage=comb[[3]],
                      pamRespectsDendro=comb[[3]],
                      maxPamDist = maxPamDist,
                      useMedoids = useMedoids) 
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

parPkMs = function(list_kMs, 
                   list_colors) {
  
  list_pkMs <- list()
  list_kMs <- lapply(list_kMs, function(x) name_for_vec(to_be_named = x, given_names = gsub("kME", "", colnames(x), ignore.case =F), dimension=2))
  
  for (j in 1:length(list_colors)) {
    # Get the 'principal kMEs', i.e. kME of each gene to the module to which it is allocated
    pkMs <- vector(mode="numeric",length=length(list_colors[[1]]))
    ###
    # Loop over every gene and get its pkME
    for (i in 1:length(pkMs)) {
      #pkMEs[i] <- list_kMEs[[j]][diffParams_colors[i,j]][i,]
      pkMs[i] <- list_kMs[[j]] [list_colors[[j]][i]] [i,]    
    }
    list_pkMs[[j]] <- pkMs 
    names(list_pkMs[[j]]) <- names(list_colors[[j]])
  }
  return(list_pkMs)
}


parPkMs_2 = function(list_kMs, list_colors) {
  
  list_pkMs <- list()
  
  for (j in 1:length(list_colors)) {
    # Get the 'principal kMEs', i.e. kME of each gene to the module to which it is allocated
    # EDIT_180421_8
    #pkMEs <- vector(mode="numeric",length=nrow(diffParams_colors))

    pkMs <- vector(mode="numeric",length=length(list_colors[[1]]))
    ###
    # Loop over every gene and get its pkME
    for (i in 1:length(pkMs)) {
      #pkMEs[i] <- list_kMEs[[j]][diffParams_colors[i,j]][i,]
      pkMs[i] <- list_kMs[[j]][list_colors[[j]][[i]]] [i,]    
    }
    
    list_pkMs[[j]] <- pkMs 
    
  }
  return(list_pkMs)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################


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

kIM_eachMod_norm = function(dissTOM, colors) {
  # compute intramodularConnectivity for every gene with regard to every cluster
  unique_colors = sort(names(table(colors)))#[-which(sort(names(table(colors))) == "grey")]
  # Convert the distance matrix to a proximity matrix. Nb: It is squary with 0 diagonal
  mTOM <- as.matrix(1-dissTOM)
  # make a gene * unique_colors matrix for results: each gene's degree w.r.t each color (module)
  results = matrix(nrow=length(colors), ncol=length(unique_colors))

  for (j in 1:length(unique_colors)) { # modules in columns
    message(paste0("Computing kIM for ", unique_colors[j]))
    for (i in 1:length(colors)) { # genes in rows
      results[i,j] <- sum(mTOM[i,colors %in% unique_colors[[j]]]) / sum(colors == unique_colors[[j]])  # for that gene, sum connections to other genes assigned to module j
    }
    #normTerms[j] <- sum(colSums(mTOM[,colors %in% unique_colors[[j]]])) 
    # 180706: Rather than normalize by the sum of degrees of a module's genes, instead normalize
    # just by the number of genes in the module, making the kIM = average degree between a gene
    # and genes in the module = "average" clustering parame
  }
  # Output a dataframe with genes in rows and modules (colors) as columns
  results <- as.data.frame(results, row.names= names(colors))
  colnames(results) = unique_colors
  return(results)
}

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


getPkMs <- function(colors, kMs) {
  pkMs <- vector(length=length(colors))
  
  for (i in 1:length(colors)) {
    pkMs[i] <- kMs[colors[i]][i,]
  }
  
  return(pkMs)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

plottable_colors <- function(colors) {
  # takes a vector of colors, i.e. module assignments
  # if all are recognised R colors then it returns the vector as is
  # otherwise, first try to remove full stops and numbers. Then try to replace them and returns the vector
  
  idx_not_real_colors = !(colors %in% colors())

  if (sum(idx_not_real_colors)>0) {
    colors[idx_not_real_colors] <- gsub("\\.\\d+|\\_\\d+", "", colors[idx_not_real_colors])
  }
  
  # If there are still any 'not real' colors, make them black
  idx_not_real_colors = !(colors %in% colors())
  n_replace <- sum(idx_not_real_colors)

  if (n_replace>0) {
    #unused = setdiff(colors(), colors) 
    colors[idx_not_real_colors] <- "black" #unused[1:n_replace]
  }
  
  return(colors)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

checkGrey <- function(kMs) {
  if (any(grepl("grey", colnames(kMs)))) kMs[['grey']] <- NULL 
  return(kMs)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

plotDendro_for_vec <- function(list_colors, 
                                geneTree,
                                list_labels,
                                subsetName,
                                title,
                                flag_file) {

  list_colors <- lapply(list_colors, plottable_colors)
  
  colors_mat = matrix(unlist(list_colors), nrow=length(list_colors[[1]]), ncol=length(list_colors))
  
  pdf(sprintf("%s%s_%s_%s_diffParams_colors_order_%s.pdf", plots_dir, data_prefix, subsetName, flag_file, flag_date), width = 9, height = 5+length(list_colors) %/% 2)
  plotDendroAndColors(geneTree, 
                      #labels2colors(list_colors_PPI_order),
                      #matrix(unlist(list_colors_PPI_order), ncol = length(list_colors_PPI_order), byrow = F), 
                      colors_mat,
                      groupLabels = list_labels, 
                      addGuide= TRUE, 
                      dendroLabels=FALSE, 
                      main=paste0(title, " - ", subsetName), 
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

plotLabel_for_vec <- function(comb) {
  # Utility function for more easily parallelising making labels
  label = paste0("MMS=", comb[[1]], ",DS=", comb[[2]],",PAM=",comb[[3]], ",CUT=",comb[[4]],sep="") 
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

factorToIndicator <- function(vec) {
  if (!class(vec) %in% c("character", "factor")) return(vec)
  vec <- as.character(x = vec)
  mat <- matrix(data=0, nrow=length(vec), ncol=length(unique(vec)))
  for (l in 1:length(unique(vec))) {
    mat[,l] <- ifelse(vec==unique(vec)[l], yes=1, no=0)
  }
  #class(mat) <- "numeric"
  dim(mat) <- c(length(vec), length(unique(vec)))
  colnames(mat) <- unique(vec)
  return(mat)
}
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# SeuratFactorToIndicator <- function(obj, meta.colnames) {
#   # @Usage:   Converts factor/character metadata to new dummy variable columns
#   # @args:    obj: Seurat object
#   #           meta.colnames: list of metadata column names - character
#   # @return:  dataframe with numeric metadata columns, column bound to 'dummy' variable columns, one for each level of each factor or character meta data column in obj
#   # @depends: Seurat
#   # @author:  Jonatan Thompson jjt3f2188@gmail.com
#   # @date:    180524
#   # @TODO: 
#   
#   idx_factor = sapply(obj@meta.data, function(x) class(x) %in% c("factor", "character"))
#   
#   metadata.factor <- obj@meta.data[idx_factor]
#   metadata.factor <- metadata.factor[colnames(metadata.factor) %in% meta.colnames]
#   meta.colnames.factor <- meta.colnames[meta.colnames %in% colnames(metadata.factor)]
#   
#   metadata.numeric <- obj@meta.data[!idx_factor]
#   metadata.numeric <- metadata.numeric[colnames(metadata.numeric) %in% meta.colnames]
# 
#   
#   # For factor metadata columns, get logical vectors for each level
#   list.list.idx = list()
#   for (j in 1:length(meta.colnames.factor)) {
#     col <- grep(meta.colnames.factor[j], names(metadata.factor)) # Find the meta data column index
#     list.idx <- lapply(names(table(metadata.factor[col])), function(x) metadata.factor[col]==x) # find, for each unique value, a logical vector of occurences
#     names(list.idx) = names(table(metadata.factor[col]))
#     list.list.idx[[j]] <- list.idx
#   }
#   
#   flatlist.idx <- unlist(list.list.idx, recursive=F, use.names = T) # Flatten to list of logical vectors  
#   #names(flatlist.idx) <- unlist(values, recursive=F) # Assign a feature name to each vector
#   new_cols <- cbind(metadata.numeric, data.frame(lapply(flatlist.idx, function(x) as.numeric(x)), row.names=row.names(obj@meta.data))) # Make numeric dataframe
#   #return(AddMetaData(object=obj, new_cols)) # Add dataframe to Seurat object and return
#   return(new_cols) # return dataframe, not seurat object
# }
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

PPI_outer_for_vec = function(colors, 
                             pkMs, 
                             STRINGdb_species, 
                             PPI_pkM_threshold, 
                             pvalThreshold, 
                             project_dir, 
                             data_prefix, 
                             flag_date) {
  
  # Rather than for parallelising over modules within a set of modules, parallelise over several sets of colors (produced by comparing parameters)
  # It calls PPI_innver_for_vec as a subroutine
  # Args:
  #   colors
  #   pkMs
  #   STRINGdb_species
  #   PPI_pkM_threshold: numeric scalar
  #   pvalThreshold: numeric scalar
  #   project_dir
  #   data_prefix
  #   flag_date
  #   
  # Returns: 
  #   colors_PPI: colors where modules that did not satisfy the PPI threshold are set to grey
  
  #suppressPackageStartupMessages(library(STRINGdb)) 
  
  module_PPI <- NULL
  unique_colors <- NULL
  unique_colors_PPI <- NULL
  colors_PPI <- NULL
  
  unique_colors <- unique(colors[pkMs>PPI_pkM_threshold])
  
  string_db <- STRINGdb$new(version="10", 
                            species = STRINGdb_species, 
                            score_threshold=0, 
                            input_directory="") # Default: 10090 (Mus musculus)
  sapply(unique_colors, function(x) PPI_inner_for_vec(color=x, 
                                                      unique_colors = unique_colors, 
                                                      colors = colors, 
                                                      string_db = string_db), simplify=T) %>% t() -> module_PPI 
  
  module_PPI <- as.data.frame(module_PPI, row.names = unique_colors)
  
  # FILTER MODULES ON PPI ENRICHMENT  

  if (!is.null(module_PPI)) {
    unique_colors_PPI = unique_colors[module_PPI$'p-value' < pvalThreshold]
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
  
  ppi <- data.frame(gene = names(colors[colors==color])) # extract the genes with the corresponding color to dataframe
  module_PPI <- list('p-value'=1, 'expected interactions'=NA)

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

mapMMtoHs = function(modulekM,
                     log_dir,
                     flag_date,
                     data_prefix,
                     mapping_orthology) {
  
  if (!is.null(modulekM$genes)) modulekM$genes <- NULL
  
  log_not_mapped_filepath = paste0(log_dir,flag_date,"_genes_orthology_not_mapped_",data_prefix,"_", ".tab")
  
  mapping = data.frame(ensembl.mouse=row.names(modulekM))
  # orthology mapping
  mapping$ensembl.human = mapping_orthology$ensembl_gene_id[ match(mapping$ensembl.mouse, mapping_orthology$mmusculus_homolog_ensembl_gene) ]
  #mapping$ensembl.human[mapping$ensembl.human == ""] = NA
  df_not_mapped = mapping[is.na(mapping$ensembl.human),]
  append = !file.exists(log_not_mapped_filepath)
  col.names = !append 
  write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t", col.names = col.names, row.names=F, append=append)
  #modulekM$symbol = mapping$symbol
  modulekM$ensembl = mapping$ensembl.human
  modulekM = na.omit(modulekM)
  
  ### 180508_v.18_dev2
  #tmp = within(modulekM, rm("symbol","ensembl"))
  tmp = within(modulekM, rm("ensembl"))
  ###
  
  # Average duplicated gene IDs
  modulekM_ens <-aggregate(tmp, by=list(modulekM$ensembl),FUN=mean, na.rm=TRUE)
  rownames(modulekM_ens) = modulekM_ens$Group.1
  modulekM_ens = within(modulekM_ens, rm("Group.1"))
  
  modulekM_ens$genes <- NULL
  
  return(modulekM_ens)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

cor_magma_pval <- function(gwas_pvals, 
                           data, 
                           indices) {
  x <- gwas_pvals
  y <- data[indices] # allows boot to select sample 
  cor = cor.test(x,y,method="spearman", exact=F)
  return(cor$estimate)
} 

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################


kM_magma <- function(cellType,
                      modulekM,
                      gwas) {
  # Usage: Subroutine to calculate spearman's correlation between gene module membership and GWAS gene significance
  # Args: 
  # Returns:
  
  colors = colnames(modulekM)
  table.kM.cor.p = table.kM.cor.r = table.kM.cor.emp.p <- matrix(NA,nrow=length(unique(colors)),ncol=length(gwas)) 
  rownames(table.kM.cor.r) = rownames(table.kM.cor.p) = rownames(table.kM.cor.emp.p) = unique(colors)
  colnames(table.kM.cor.r) = colnames(table.kM.cor.p) = colnames(table.kM.cor.emp.p) = names(gwas) 
  
  for (col in unique(colors)) {
    for (j in 1:length(gwas)) {
      #col = paste("kM", m, sep="")
      genes = intersect(rownames(modulekM),gwas[[j]]$gene_name)
      x = -log10(gwas[[j]]$P[match(genes, gwas[[j]]$gene_name)])
      y = modulekM[match(genes,rownames(modulekM)), col]
      
      cor = cor.test(x,y,method="spearman", exact=F)
      table.kM.cor.r[col,j] <- cor$estimate
      table.kM.cor.p[col,j] <- cor$p.value
      
      # Generate 1000 (10000) permutation distribution samples and their associated correlation coefficients
      # For a 0.05 significance threshold this gives an error of the p-value of sqrt(p(1-p)/k) ~ 0.00689202437 (0.0021794497)
      # see https://stats.stackexchange.com/questions/80025/required-number-of-permutations-for-a-permutation-based-p-value
      boot_out <- boot(data = y,
                       statistic=cor_magma_pval,
                       R=R,
                       sim = "permutation",
                       gwas_pvals=x,
                       parallel = "no")
      
      # compute the empirical probability of the p-value - should correspond to the p-value, ideally..
      table.kM.cor.emp.p[col,j] <- ecdf(boot_out$t)(boot_out$t0)
    }
  }
  
  rownames(table.kM.cor.r) <- rownames(table.kM.cor.p) <- rownames(table.kM.cor.emp.p) <- paste0(cellType, "__", rownames(table.kM.cor.p))
  
  return(list('p.val'= table.kM.cor.p, 'corrCoef' = table.kM.cor.r, 'emp.p.val' = table.kM.cor.emp.p))
  
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

makeColorsUniqueForPlot = function(cols) {
  # DEPRECATED
  
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

mendelianGenes <- function(cellType,
                           colors,
                           coding_variants,
                           mendelian_genes) {
  # Usage: compute rare variants enrichment for a set of modules
  # Arguments:
  #   kMs: table of kMs (WGCNA)
  #   colors: character vector of color (module) assignments (WGCNA), with names = genes (human or mouse?)
  #   coding_variants: vector of genes which...
  #   mendelian_genes: genes for which SNPs on a single gene is associated with phenotype (BMI)

  # Gene lists
  gene.lists = vector(mode="list")
  gene.lists[[1]] = coding_variants
  gene.lists[[2]] = mendelian_genes
  names(gene.lists)= c("coding_variants", "mendelian")
  
  # Change gene names to upper case
  names(colors) <- toupper(names(colors)) # change genes to uppercase 
  
  logit.or = matrix(NA, nrow=length(unique(colors)), ncol=length(gene.lists)) # make a n_modules * n_genelist matrix
  rownames(logit.or) = unique(colors) 
  logit.or <- logit.or[rownames(logit.or) != "grey", , drop=F] # remove the grey module
  colnames(logit.or) = names(gene.lists) 
  # Set up tables for odds ratios (table.or), p-value
  logit.p <- logit.or 
  
  for(j in 1:ncol(logit.or)) {
    for(i in 1:nrow(logit.or)) {
      
      #select the first color (i.e. first rowname)
      col = rownames(logit.or)[i]
      
      # make a n_genes * 2 dataframe; first column is 1 for genes in that color, 0 otherwise 
      binaryMat = as.data.frame(cbind(as.numeric(colors==col), 0))
      colnames(binaryMat) = c("Module", "GeneSet")
      # find idx in color vector of genes in sets and set those to '1' in the vector
      idx = match(gene.lists[[j]], names(colors))
      binaryMat$GeneSet[idx] = 1
      
      tryCatch({glm.out <- glm(binaryMat$GeneSet~binaryMat$Module, family=binomial())
      logit.or[i,j] = exp(coefficients(glm.out)[2])  # Calculate odds ratio from logistic regression
      logit.p[i,j] = summary(glm.out)$coefficients[2,4]  # get P value for coefficient
      
      }, 
      error = function(c) {
        logit.or[i,j] = 0.0
        logit.p[i,j] = 1.0  
      })
      
    }
  }
  
  rownames(logit.or) <- rownames(logit.p) <- paste0(cellType, "__", rownames(logit.or))
  
  
  
  
  return(list("or" = logit.or,  "p.val"=logit.p))
  
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# NOT IN USE 
# 180522: /!\ can't install DOSE

# GO_for_par = function(list_kMs_gwas) {
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
#   # # select only kM columns of gwas enriched module
#   # list_kMs_gwas <- mapply(function(x,y) x[,colnames(x) %in% y], x = list_kMs_PPI_ok_gwas, y = list_modules_PPI_gwas)
#    
#   # Order genes by kM to their own module. This allows us to submit them as ranked queries to gprofiler for GSEA style p-values (https://cran.r-project.org/web/packages/gProfileR/gProfileR.pdf)
#   list_list_module_gwas_genes_order <- lapply(list_kMs_gwas, function(x) lapply(colnames(x), function(y) rownames(x)[order(x[[y]], decreasing=T)]))
#   
#   ### 190509_v1.8.dev2
#   ### # Got to here!
#   invisible(gc())
#   cl <- makeCluster(n_cores, type = "FORK", 
#                     outfile = paste0(log_dir, "log_gprofiler.txt"))
#   
#   list_list_ggo <- parLapplyLB(list_list_module_gwas_genes_order, function(x) lapply(x, function(y) groupGO(gene = y,
#                                                                                                           OrgDb    = org.Mm.eg.db,
#                                                                                                           ont      = "CC",
#                                                                                                           #level    = 3,
#                                                                                                           readable = T)))
#   
#   
#   ego3 <- parLapplyLB(list_list_module_gwas_genes_order, function(x) lapply(x, function(y) gseGO(geneList = y,
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
#   # select only kM columns of gwas enriched moduled 
#   list_kMs_PPI_ok_gwas <- mapply(function(x,y) x[,colnames(x) %in% y], x = list_kMs_ok[names(list_kMs_PPI_ok) %in% sNames_PPI_gwas], y = list_modules_PPI_gwas)
#   
#   # Order genes by kM to their own module. This allows us to submit them as ranked queries to gprofiler for GSEA style p-values (https://cran.r-project.org/web/packages/gProfileR/gProfileR.pdf)
#   list_list_module_PPI_gwas_genes_order <- lapply(list_kMEs_PPI_ok_gwas, function(x) lapply(colnames(x), function(y) rownames(x)[order(x[[y]], decreasing=T)]))
#   
#   # name them
# 
#   invisible(gc())
#   cl <- makeCluster(n_cores, type = "FORK", 
#                     outfile = paste0(log_dir, "log_gprofiler_PPI.txt"))
#   
#   ### 180607_v1.8_dev2
#   list_list_gprofiles_PPI <- parLapplyLB(cl, list_list_module_PPI_genes_order, function(x) gprofiler(query=x, 
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


  # melt each nested list to a single dataframe with a gene column and a module (list name) column
  sigmod_list <- lapply(list_list_module_genes, function(x) melt.list(x))
  
  # Join list of dataframes into a single dataframe by the gene columns
  sigmod_list %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2, by = 'value'), .) -> gene_modules.df
  
  # sigmod_list %>%
  #   Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="value"), .)->gene_modules.df
  
  #gene_modules.df["L1"] <- NULL
  gene_modules.df <- sapply(gene_modules.df, as.character) %>% as.data.frame
  
  colnames(gene_modules.df) <- c("genes", names(sigmod_list))
  
  gene_modules.df[is.na(gene_modules.df)]<-'random'
  
  row.names(gene_modules.df) <-  gene_modules.df$genes
  
  gene_modules.df<-gene_modules.df[,-1, drop=F]
  
  # FILTER SEURATOBJ 
  
  # load seurat object
  seurat_obj <- load_obj(f=sprintf("%sseurat_obj_ensembl.RData",RObjects_dir))
  # Get ident
  ident <- seurat_obj@ident
  # filter
  s.coexp <- t(seurat_obj@data[rownames(seurat_obj@data) %in% row.names(gene_modules.df),, drop=F]) 
  s.coexp <- as.matrix(s.coexp)
  # Free up memory
  rm(seurat_obj)
  # order genes and modules same
  s.coexp <- s.coexp[,order(colnames(s.coexp)), drop=F]
  gene_modules.df<-gene_modules.df[order(row.names(gene_modules.df)),, drop=F]
  #identical(colnames(s.coexp),row.names(gene_modules.df))
  eiglist<-as.list(as.data.frame(gene_modules.df))
  
  # Score cells on eigengenes 
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_eigengene_allscore.txt"))
  eigenge_allscore.list<-parLapplyLB(cl, eiglist, function(x) moduleEigengenes((s.coexp), x, impute = T)$eigengenes) #Score cell matrix on each eigengene found in each celltype
  stopCluster(cl)
  invisible(gc())
  
  #Rename Columns in each eigengene matrix
  for (i in 1:length(eigenge_allscore.list)){ 
    colnames(eigenge_allscore.list[[i]])<-paste0(names(eigenge_allscore.list)[[i]],'_',colnames(eigenge_allscore.list[[i]])) 
  }
  
  eigen_mat <- bind_cols(eigenge_allscore.list) # Merge Eigengene Matrices
  eigen_mat <- eigen_mat[,!grepl('.*merandom$',colnames(eigen_mat), ignore.case = T), drop=F] #Remove eigengene scores on genes that were not found in modules in a cell type
  row.names(eigen_mat)<-row.names(s.coexp) 

  # Reorder rows in eigenmat by cell type
  for (i in 1:length(unique(ident))) {
    if (i==1) {
      eigen_mat_order <- eigen_mat[ident== names(table(ident))[i],,drop=F]
    } else {
      eigen_mat_order <- rbind(eigen_mat_order, eigen_mat[ident== names(table(ident))[i],,drop=F])
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

cellModEmbed <- function(datExpr, 
                        colours=NULL, 
                        latentGeneType,
                        cellType=NULL,
                        kMs=NULL) {
  # datExpr should be cell x gene matrix or data.frame with ALL genes and ALL cells in the analysis
  # colors is a character vector with gene names
  # the genes should be ordered identically
  
  # prepare celltype character for naming columns
  if (is.null(cellType)) cellType_prefix <- "" else cellType_prefix <- paste0(cellType, "__")
  if (latentGeneType == "ME") { 
    list_datExpr <- lapply(unique(colours)[unique(colours)!="grey"], function(x) datExpr[,match(names(colours)[colours==x], colnames(datExpr))])
    embed_mat <- as.matrix(sapply(list_datExpr, function(x) prcomp_irlba(t(x), n=1, retx=T)$rotation, simplify = T))
    colnames(embed_mat) <- paste0(cellType_prefix, unique(colours)[unique(colours)!="grey"]) 
  } else if (latentGeneType == "IM") {
    list_datExpr <- lapply(colnames(kMs), function(x) datExpr[,match(names(colours)[colours==x], colnames(datExpr))])
    list_kMs <- mapply(function(x,y) kMs[match(names(colours)[colours==y], rownames(kMs)),y],
                       x= kMs,
                       y=colnames(kMs),
                       SIMPLIFY=F)
    
    embed_mat <- as.matrix(mapply(function(x,y) x %*% as.matrix(y),
                                  x = list_datExpr,
                                  y = list_kMs,
                                  SIMPLIFY=T))
    # datExpr_1 = datExpr[, match(rownames(kMs), colnames(datExpr), nomatch=0) [match(rownames(kMs), colnames(datExpr),  nomatch=0)>0]]
    # embed_mat <- as.matrix(datExpr_1) %*% as.matrix(kMs)
    colnames(embed_mat) <- paste0(cellType_prefix, colnames(kMs)) 
  } 
  rownames(embed_mat) <- rownames(datExpr)
  return(embed_mat)
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
  mean.celltype.eig.expr <- t(sapply(names(table(meta$neuron_groups)), function(x) rowMeans(t(eigen_mat[meta$neuron_groups==x,,drop=F])), simplify=T))


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

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

ipak <- function(pkg, bioconductor = F){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    if (bioconductor==F) {
      try(install.packages(new.pkg, dependencies = TRUE))
    } else {
      source("https://bioconductor.org/biocLite.R")
      biocLite(pkg)
    }
  }
  sapply(pkg, require, character.only = TRUE)
}

# usage
# packages <- c("ggplot2", "plyr", "reshape2", "RColorBrewer", "scales", "grid")
# ipak(packages)

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
# Remove packages installed under 3.5.0

#pak3.5 <- installed.packages()[,"Package"][grep("3.5.0", installed.packages()[,"Built"])]

rpak <- function(pkg){
  installed.pkg <- pkg[(pkg %in% installed.packages()[, "Package"])]
  if (length(installed.pkg)) {
    try(remove.packages(installed.pkg))
  }
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

wrap_consensusTOM <- function(...) {
  # a function wrapper the sole purpose of which is to create a function environment which contains only the 
  # necessary variables (...) required by consensusTOM. These should be named.
  # This creates a suitable environment for using makeCluster to set up fork clusters while sending only the necessary variables
  # from the global environment to the workers, thereby saving memory.
  
  cl <- makeCluster(n_cores_consTOM, type = "FORK", outfile = paste0(log_dir, "log_consensusTOM.txt"))
  
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
                                                              cacheDir = RObjects_dir,
                                                              cacheBase = ".blockConsModsCache",
                                                              verbose = verbose,
                                                              indent = indent), 
                               x=list_multiExpr, 
                               y=sNames, 
                               SIMPLIFY=F)
  stopCluster(cl)
  invisible(gc())
  
  return(list_consensus)
}