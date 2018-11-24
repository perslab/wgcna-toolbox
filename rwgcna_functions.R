# Title: sc functions
# Author: Jonatan Thompson, Pers Lab

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

FilterGenes <- function(seurat_obj_sub, min.cells) {
  
  if (min.cells > 0 & !is.null(dim(seurat_obj_sub@raw.data))) { 
    num.cells <- rowSums(seurat_obj_sub@raw.data > 0)
    genes.use <- names(x = num.cells[which(x = num.cells >= min.cells)])
    if (length(genes.use)>0) {
      seurat_obj_sub@raw.data <- seurat_obj_sub@raw.data[genes.use, ]
    } else {
      seurat_obj_sub <- NULL
    }
  }
  
  return(seurat_obj_sub)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

wrapJackStraw = function(seurat_obj_sub, n_cores, jackstrawnReplicate, pvalThreshold, n_genes_use=5000) {
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
    
    if (nrow(pAll) == length(seurat_obj_sub@var.genes) | length(PC_select_idx)<5) {
      
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
      
    } else if (nrow(pAll) > length(seurat_obj_sub@var.genes) & length(PC_select_idx)>=5) {
      
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
    names_genes_use <- names(max_loadings[order(max_loadings, decreasing = T)])[1:min(n_genes_use, length(max_loadings))]
    
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
# refactored - no need for function here
# TOM_for_par = function(datExpr, subsetName, softPower, data_prefix) {
#   
#   ### 180503 v1.7
#   disableWGCNAThreads()
#   ###
#   
#   adjacency = adjacency(datExpr, 
#                         type=type, 
#                         power = softPower, 
#                         corFnc = corFnc, 
#                         corOptions = corOptions)
#   
#   consTomDS = TOMsimilarity(adjMat=adjacency,
#                             TOMType=TOMType,
#                             TOMDenom=TOMDenom,
#                             verbose=verbose,
#                             indent = indent)
#   
#   save(consTomDS, file=sprintf("%s%s_%s_consensusTOM-block.1.RData", scratch_dir, data_prefix, subsetName)) # Save TOM the way consensusTOM would have done
#   goodGenesTOM_idx <- rep("TRUE", ncol(datExpr))
#   return(goodGenesTOM_idx)
# }

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

dissTOM_for_par = function(subsetName, data_prefix) {
  # We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes.
  # We use Dynamic Tree Cut from the package dynamicTreeCut.    
  load(sprintf("%s%s_%s_consensusTOM-block.1.RData", scratch_dir, data_prefix, subsetName)) # Load the consensus TOM generated by blockwiseConsensusModules
  
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
  # deprecated
  colors <- NULL
  MEs <- NULL 
  tryCatch({
    merged = mergeCloseModules(exprData=as.matrix(datExpr), 
                               colors = cutree$labels, 
                               impute =T,
                               corFnc = corFnc,
                               corOptions = corOptions,
                               cutHeight=comb[[4]],
                               iterate = T,
                               getNewMEs = F,
                               getNewUnassdME = F)
    
    
    colors = labels2colors(merged$colors)
    #MEs = merged$newMEs 
    MEs = moduleEigengenes_uv(expr = as.matrix(datExpr),
                                     colors=colors,
                                     excludeGrey = excludeGrey)
    
    
  }, error = function(c) {
    warning(paste0("MergeCloseModules failed"))
  })
  return(list("cols" = colors, "MEs"= MEs))
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# extract_and_name_colors <- function(merged, datExpr) {
#   colors <- merged[[1]]
#   names(colors) = colnames(datExpr)
#   return(colors)
# }

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

mergeCloseModskIM = function(datExpr,
                             colors,
                             kIMs,
                             dissTOM = NULL,
                             moduleMergeCutHeight,
                             verbose=2,
                             cellType) {
  
  # Usage: compute correlations between gene module cell embeddings. 
  #         merge modules with a positive correlation > 1-moduleMergeCutHeight
  #
  # Args: datExpr: a cell * gene expression matrix
  #       colors: vector module assignments with gene names
  #       kIMs: generalised IntraModular connectivity, dataframe of gene row and module column. Should not include the grey module
  #       dissTOM: only needed if iterate=T, in order to recompute kIMs
  #       cellModEmbed_mat: embedding matrix of cell rows and module columbs. Passing the matrix saves time; otherwise will be computed
  #       moduleMergeCutHeight: max distance (1-correlation coefficient) for merging modules
  #       verbose: verbosity of WGCNA::cor function, 1-5
  #
  # Returns: list with entries colors=colors_merged, kIMs = kIMs_merged
  
  
  corr_clust <- character(length = length(unique(colors)))
  colors_original <- colors
  # Remove 'grey' (unassigned) module 
  #if (any(colnames(kIMs) == "grey")) kIMs[["grey"]] <- NULL
  
  if (length(unique(colors)) > 1) { # if there is more than one non-grey module
    while (TRUE) {
      
      message(paste0(cellType, ": computing cell-module embedding matrix"))
      cellModEmbed_mat <- cellModEmbed(datExpr=datExpr, 
                                       colors=colors, 
                                       latentGeneType = "IM",
                                       cellType=NULL,
                                       kMs=kIMs)
      
      #cellModEmbed_mat <- cellModEmbed_mat[,-grep("^grey$", colnames(cellModEmbed_mat))]
      # Cluster modules using the Pearson correlation between module cell embeddings
      mod_corr <- WGCNA::cor(x=cellModEmbed_mat, method=c("pearson"), verbose=verbose)
      mod_corr[mod_corr<0] <- 0 #only keep positive correlations
      corr_dist <- as.dist(m = (1-mod_corr), diag = F, upper = F) # convert to (1-corr) distance matrix
      corr_dendro <- hclust(d = corr_dist, method = "average")
      
      # use a simple cut to determine clusters
      corr_clust = cutreeStatic(dendro = corr_dendro, 
                                cutHeight = moduleMergeCutHeight, 
                                minSize=1)
      names(corr_clust) =  corr_dendro$labels
      
      # At this point if corr_clust has as many unique modules as the original partition, the while loop will exit
      if(length(unique(corr_clust)) == length(unique(colors))) {  
        message(paste0(cellType, " done"))
        break
      }
      
      # if not we merge modules
      n_merge <-  length(unique(colors))  - length(unique(corr_clust)) 
      if (verbose>0) message(paste0(cellType, ": ", n_merge, " modules to be merged with others"))
      
      merge_idx <- sapply(unique(corr_clust), function(x) sum(corr_clust==x)>1)
      
      # merge modules
      for (clust in unique(corr_clust)[merge_idx]) { # loop over numeric cluster assignments 
        new_color = names(corr_clust)[corr_clust==clust][1] # use the colors name of the first module in the cluster
        for (color in names(corr_clust[corr_clust==clust])) { # loop over all module colors in the cluster
          colors[colors==color] <- new_color # give them all the new colour
        }
      }
      
      # compute merged module kIMs 
      message(paste0("Computing merged module kIMs for ", cellType))
      # Compute new kIMs
      kIMs <- kIM_eachMod_norm(dissTOM = dissTOM, 
                               colors = colors,
                               verbose = 1,
                               excludeGrey = F)
      
    }
    names(colors) = names(colors_original)
  }
  return(list(colors=colors, kIMs = kIMs))
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

parkMEs = function(list_MEs, datExpr, cellCluster) {
  list_kMEs <- NULL
  tryCatch({
    # Compute kMEs and pkMEs for each set of colors corresponding to distinct parameters
    list_kMEs <- vector(mode="list", length=length(list_MEs))
    
    for (j in 1:length(list_MEs)) {
      list_kMEs[[j]] = signedKME(datExpr=as.matrix(datExpr),
                                 datME=list_MEs[[j]],
                                 #list_MEs[[j]],#$eigengenes,
                                 outputColumnName="",
                                 corFnc = corFnc )
    }
  }, error = function(cellCluster) warning(paste0("signedkME failed for ", cellCluster)))
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
  
  #if (!all(sapply(list_kMs, function(x) any("grey" %in% colnames(x)) ))) stop("'grey' module kMEs not found") 
  
  list_pkMs <- list()
  #if(!is.null(list_kMs)) list_kMs <- lapply(list_kMs, function(x) name_for_vec(to_be_named = x, given_names = gsub("kME", "", colnames(x), ignore.case =F), dimension=2))
  
  for (j in 1:length(list_colors)) {
    pkMs <- NULL
    # Get the 'principal kMEs', i.e. kME of each gene to the module to which it is allocated
    if (!is.null(list_kMs[[j]]) & length(unique(list_colors[[j]]))>1) {
      pkMs <- vector(mode="numeric",length=length(list_colors[[j]]))
      ###
      # Loop over every gene and get its pkME
      for (i in 1:length(pkMs)) {
        #pkMEs[i] <- list_kMEs[[j]][diffParams_colors[i,j]][i,]
        pkMs[i] <- if (list_colors[[j]][i]=="grey") 0 else list_kMs[[j]] [list_colors[[j]][i]] [i,]    
      }
      
      list_pkMs[[j]] <- pkMs 
      names(list_pkMs[[j]]) <- names(list_colors[[j]])
      
    } else {
      list_pkMs[[j]] = rep(0, length(list_colors[[j]]))
      names(list_pkMs[[j]]) <- names(list_colors[[j]])
    }
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

pkMs_fnc = function(kMs, colors) {
  pkMs <- NULL
  if (!is.null(kMs) & length(unique(colors))>1) {
    for (i in 1:length(colors)) {
      pkMs[i] <- if (colors[i]=="grey") 0 else kMs[colors[i]] [i,]
    }
    names(pkMs) <- names(colors)
  } else {
    pkMs = rep(0, length(colors))
    names(pkMs) <- names(colors)
  }
  return(pkMs)
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

kM_reassign_fnc = function(colors,
                           fuzzyModMembership,
                           dissTOM = NULL,
                           datExpr = NULL,
                           verbose=2,
                           corFnc=NULL,
                           max_iter=3,
                           cellType) {
  
  # Usage: iteratively reassign genes to modules based on kIM / kME until the number to reassign <= stop_condition
  # Args: 
  #   fuzzyModMembership: "kME" or "kIM"
  #   dissTOM: list of TOM distance matrices. Only required if fuzzyModMembership = "kIM"
  #   datExpr: cell x gene expression matrix. Only required if fuzzyModMembership = "kME"
  #   colors: vector of module assignments with gene names
  #   corFnc: WGCNA correlation function; only needed if fuzzyModMembership = "kME"
  #   verbose: 0-5, controls message output
  #   max_iter: integer; max iterations.
  # Value:
  #   list with entries "colors", "kMs" and "log". 
  
  print(paste0(cellType, ": Reassigning genes to modules with higher ", fuzzyModMembership))
  
  # initialise
  tryCatch({
    colors_original <- colors_new <- colors
    reassign_total  <- reassign_t_1 <- logical(length(colors))
    reassign_t <- ! logical(length(colors))
    min_change = 5 # the minimum change from one iteration to the next required to continue
    t = 1
    MEs <- NULL
    kMs <- NULL 
    log=NULL 
    
    if (!is.null(colors) & length(unique(colors))>1) {
      while(TRUE) {
        if (fuzzyModMembership == "kME") {
          message(paste0(cellType, ": Computing Module Eigengenes"))
          MEs <- moduleEigengenes_uv(expr=datExpr, # TODO do we need to make it into a dataframe with names?..
                                            colors = colors,
                                            excludeGrey=T)$eigengenes
          
          message(paste0(cellType, ": Computing ", fuzzyModMembership, "s"))
          kMs <- signedKME(as.matrix(datExpr),
                           MEs,
                           outputColumnName = "",
                           corFnc = corFnc)
          
          message(paste0(cellType, ": Computing primary "), fuzzyModMembership, "s")
          
        }  else if (fuzzyModMembership == "kIM") {
          
          kMs <- kIM_eachMod_norm(dissTOM=dissTOM, 
                                  colors=colors, 
                                  verbose=verbose,
                                  excludeGrey=T)
        }
        
        pkMs <- pkMs_fnc(kMs=kMs, colors=colors)
        
        maxkMs <- max.col(kMs, ties.method = "random") #  integer vector of length nrow(kMs)
        # Reassign if there is a kM value to another module which is more than 1.05 times the current pkM value
        colors_new <- mapply(function(i, j, pkM) if (kMs[i,j] > 1.05*pkM) colnames(kMs)[j] else colors[i], i = 1:length(maxkMs), j = maxkMs, pkM=pkMs, SIMPLIFY=T) # get new module assignment vector
        colors_new[colors=="grey"] <- "grey" # Do not reallocate genes previously allocated to grey, since they are likely to get reallocated as "grey" is not a real cohesive module
        names(colors_new) <- names(colors)
        
        reassign_t <- colors_new != colors
        
        if ((t>1 & sum(reassign_t_1) - sum(reassign_t) < min_change) | sum(reassign_t) < min_change | t >= max_iter) break #else message(paste0(cellType, ", iteration ", t, ": reassigning ",  sum(reassign_t), " genes to new modules based on their ", fuzzyModMembership))
        
        colors <- colors_new
        
        reassign_total <- reassign_total | reassign_t
        reassign_t_1 <- reassign_t
        
        t = t+1
      }
      
      log <- data.frame("gene" = names(colors)[reassign_total], 
                        "original_module" = colors_original[reassign_total], 
                        "new_module" = colors_new[reassign_total], stringsAsFactors=F, row.names=NULL)
      
      if (verbose > 0) print(paste0(cellType, ": a total of ", sum(reassign_total), " genes reassigned to new modules"))
      
    } else {
      warning(paste0(cellType, ": one or no modules, nothing to reassign"))
    }
    return(list("colors" = colors_new, "kMs" = kMs, "log" = log))
  }, 
  error = function(c) {
    warning(paste0("kM_reassign_fnc failed for ", cellType, " with the following error: ", c))
    return(list("colors" = colors_original, 
                "kMs" = NULL,
                "log" = NULL))}
  )
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
# DEPRECATED 
kM_reassign_neg_fnc = function(Ms, 
                               kMs, 
                               pkMs, 
                               kM_reassign_threshold, 
                               colors, 
                               datExpr, 
                               corFnc, 
                               excludeGrey)  {
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
    MEs = moduleEigengenes_uv(expr = as.matrix(datExpr),
                                     colors,
                                     excludeGrey = excludeGrey)$eigengenes
    
    kMEs = signedKME(as.matrix(datExpr),
                     MEs,
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

kIM_eachMod_norm = function(dissTOM, colors, verbose=2, excludeGrey=T, do.par=F, n_cores=3) {
  # compute intramodularConnectivity for every gene with regard to every module
  # Args:
  #   dissTOM: distance matrix from WGCNA in sparse matrix form
  #   colors: vector of module assignment 'color' labels with gene names
  #   verbose: controls messages
  #   excludeGrey: compute "grey" kME?
  # Value:
  #   kIMs: gene x module dataframe of normalised intramodular gene-module connectivity scores,
  #         which are just the mean distance between the gene and every gene in that module

  unique_colors = sort(names(table(colors)))#[-which(sort(names(table(colors))) == "grey")]
  if (excludeGrey ==T) unique_colors <- unique_colors[!unique_colors=="grey"]
  # Convert the distance matrix to a proximity matrix. Nb: It is square with 0 diagonal
  mTOM <- as.matrix(1-dissTOM)
  
  rm(dissTOM)
  
  colnames(mTOM) <- rownames(mTOM) <- names(colors)
  
  # make a gene * unique_colors matrix for results: each gene's degree w.r.t each color (module)
  results = matrix(nrow=nrow(mTOM), ncol=length(unique_colors))
  
  for (j in 1:length(unique_colors)) { # modules in columns
    
    idx_mTOM_col <- match(names(colors)[colors %in% unique_colors[[j]]], colnames(mTOM))
    idx_mTOM_col <- idx_mTOM_col[!is.na(idx_mTOM_col)]
    norm_term <- max(1,length(idx_mTOM_col)-1)
    #norm_term <- (sum(colors == unique_colors[[j]])-1)
    mTOM_mod_sub <- mTOM[, idx_mTOM_col, drop=F]
    cl <- NULL
    if (do.par) {
      invisible(gc()); invisible(R.utils::gcDLLs())
      cl <- try(makeCluster(spec=n_cores, type="FORK", timeout=30))
      if (!"try-error" %in% class(cl)) {
        if (verbose>2) message(paste0("Computing kIM for ", unique_colors[j], " using ", n_cores, " cores"))
        results[,j] <- parSapply(cl=cl, X = 1:nrow(mTOM_mod_sub), FUN=function(i) {
          sum(mTOM_mod_sub[i, ]) / norm_term
        })
        stopCluster(cl)
      }
    }
    if ("try-error" %in% class(cl) | !do.par) {
      do.par <- F
      if (verbose>2) message(paste0("Computing kIM for ", unique_colors[j], " serially"))
      results[,j] <- sapply( X = 1:nrow(mTOM_mod_sub), FUN=function(i) {
        sum(mTOM_mod_sub[i, ]) / norm_term
      })
    }
    rm(mTOM_mod_sub)
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    if (verbose>2) message(paste0("Done computing kIM for ", unique_colors[j]))
  }
  
  if (verbose>0) message("Done computing kIMs for set of colors")
  
  # Output a dataframe with genes in rows and modules (colors) as columns
  
  results <- as.data.frame(results, row.names= rownames(mTOM))
  colnames(results) = unique_colors
  return(results)
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

pkIM_norm_var <- function(dissTOM, 
                          colors,
                          cellType,
                          verbose=2, 
                          excludeGrey=T) {
  
  if (verbose>0) message(paste0(cellType, ": computing pkIM variances"))
  # compute variance of intramodularConnectivity for every gene with regard to its own
  # Args:
  #   dissTOM: distance matrix from WGCNA in sparse matrix form
  #   colors: vector of module assignment 'color' labels with gene names
  #   verbose: controls messages
  #   excludeGrey: compute "grey" kME?
  # Value:
  #   kIMs: gene x module dataframe of normalised intramodular gene-module connectivity scores,
  #         which are just the mean distance between the gene and every gene in that module
  
  # Convert the distance matrix to a proximity matrix. Nb: It is square with 0 diagonal
  mTOM <- as.matrix(1-dissTOM)
  # make a gene * unique_colors matrix for results: each gene's degree w.r.t each color (module)
  results <- mapply(function(gene,color) var(mTOM[names(colors)==gene, colors == color & names(colors) != gene] / (sum(colors == color)-1)), gene=names(colors), color=colors, SIMPLIFY=T)
  if (verbose>0) message(paste0(cellType, ": Done computing pkIM variances"))
  return(results)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

kEM_norm = function(dissTOM, colors, cellType, verbose=2, excludeGrey=T) {
  # compute normalised extramodularConnectivity for every gene with regard to all genes outside its module (a single number)
  # Args:
  #   dissTOM: distance matrix from WGCNA in sparse matrix form
  #   colors: vector of module assignment 'color' labels with gene names
  #   verbose: controls messages
  #   excludeGrey: compute "grey" kME?
  # Value:
  #   kEMs: gene-named vector of normalised extramodular connectivity scores, 
  #   i.e. the average link between each gene and all genes outside its module
  
  # Convert the distance matrix to a proximity matrix. Nb: It is square with 0 diagonal
  message(paste0(cellType, ": Computing kEMs"))
  mTOM <- as.matrix(1-dissTOM)
  results <- mapply(function(gene,color) sum(mTOM[names(colors)==gene, colors != color]) / sum(colors != color), gene=names(colors), color=colors, SIMPLIFY=T)
  if (verbose>0) message(paste0(cellType, ": Done computing kEMs"))
  names(results) <- names(colors)
  return(results)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

kEM_norm_var = function(dissTOM, colors, cellType) {
  # compute variance of extramodular Connectivity for each gene
  message(paste0(cellType, ": Computing kEM variances"))
  mTOM <- as.matrix(1-dissTOM)
  tryCatch({
    results = mapply(function(color, gene) stats::var(x=(mTOM[names(colors)==gene, !colors==color])/sum(colors!=color), na.rm=T), 
                     color=colors, 
                     gene=names(colors), 
                     SIMPLIFY=T) # for gene i, compute variance of connections to other modules
    names(results) = names(colors)
  }, error = function(c) {
    results = vector(mode="numeric", length=length(colors))
    warning(paste0(cellType, ": kEM variance calculation failed"))
  })
  if (verbose>0) message(paste0(cellType, ": Done computing kEM variances"))
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
                                   nThreads,
                                   verbose, 
                                   indent) {
  
  # STATUS: 181026: Fixed a bug in WGCNA::modulePreservation (https://www.biostars.org/p/339950/) (actually arises when calling signedKME)
  # @Usage: Wrapper for the WGCNA modulePreservation function to make it easier to use for single and multiple datasets.
  #         The wrapper puts datasets and module assignment (colors) into the correct multiData and multiColor formats.
  #         
  #         WGCNA::modulePreservation notes
  #         * This function calculates module preservation statistics pair-wise between given reference sets and all other sets in multiExpr.
  #         * Reference sets must have their corresponding module assignment specified in multiColor; module assignment is optional for test sets.
  #         * Individual expression sets and their module labels are matched using names of the corresponding components in multiExpr and multiColor.
  #         *
  #         *
  #
  #         WGCNA::modulePreservation performs two distinct tasks:
  #
  #           1. evaluates the quality of one or more sets of module assignments in a single dataset
  #           2. evaluates the preservation of a single set of colors in multiple datasets
  #           
  #           Multiple colors in multiple datasets is not supported.
  #
  #         In each case, modulePreservation function  requires 'reference' and 'test datasets regardless.
  #         The function also requires the data to be presented in a multiExpr format (see WGCNA function checkSets).
  #
  #         In case 1 (multiple module assignments, one dataset) this wrapper function makes 'test' datasets just by copying the original to satisfy the required modulePreservation args.
  #         As the colors may come from a different dataset, the function uses the gene intersect to match each set of colors to the dataset
  # 
  #         In the case 2 (one set of module assignments, multiple datasaets), this wrapper assumes that the first dataset is the reference and the rest are test sets, 
  #         and calls modulePreservation do do pairwise ref-test comparisons.
  #         
  #         TODO: How to match genes between 1 set of colors and multiple datasets?
  #
  # @args: 
  #         listDatExpr: a vector(list) of expression data (genes in columns, cells in rows). If the list contains a single dataset the function outputs the 
  #                     quality statistics only. Assumes that the first entry holds the reference dataset and the rest hold test sets, if relevant.
  #                     If the list entries are names the function uses them as labels by default unless a list of labels is provided.
  #         listColors: list of color assignment vectors with gene names
  #                     TODO If given a single set of colors for multiple datasets the function cannot compute hypergeometric stats ? 
  #         labels:     a vector of character, one per set of colors. If not given, defaults to names(listColors); if these are NULL, defaults to integers
  #         ...         (see modulePreservation)
  #
  # @return: 
  #         result:     a list of list of dataframes? 
  #                     TODO, each nested list corresponds to a pairwise evaluation     
  # @author: Jonatan Thompson jjt3f2188@gmail.com
  # @date: 180316
  
  if (parallelCalculation) require("doParallel")
  
  stopifnot(is.null(dim(listDatExpr)) & typeof(listDatExpr) == "list" & length(listDatExpr) >= 1 & is.null(dim(listColors)) & typeof(listColors)=="list" & length(listColors) >= 1) # basic input format checks
  stopifnot(length(listDatExpr) == 1 | length(listColors) == 1) # Cannot have many datasets and many colourings so at least one list must have length 1
  
  # Average duplicate genes
  # listDatExpr <- lapply(listDatExpr, function(datExpr) {
  #   if (any(duplicated(colnames(datExpr)))) datExpr <- aggregate(datExpr, by=list(colnames(datExpr)), FUN=mean, na.rm=TRUE)
  #   return(datExpr)
  #   })
  
  if (length(listDatExpr) == 1) {
    
    referenceNetworks = c(1:length(listColors)) # 
    testNetworks = as.list(rep(length(listColors)+1, length(listColors)))
    names(testNetworks) = referenceNetworks
    # Set up the multi-set format (see WGCNA checkSets() ). The modulePreservation function requires test datasets even when evaluating module assignments in a single dataset,
    # so we create a multi-set format list with length(colors)+1 copies of the dataset, where the extra acts as test set.
    # the reason for this memory-inefficient copying of datasets is that modulePreservation requires a dataset to match each set of colors.
    
    multiExpr <- vector("list", length=length(listColors)+1)  # make a list of copies of the expression data, length = number of colourings to evaluate + 1 for 'test set'
    
    for (i in 1:length(multiExpr)) {
      multiExpr[[i]]$data  <- as.matrix(listDatExpr[[1]]) # NB: this whole section only runs if there is just one dataset
    }
    
    names(multiExpr) <- c(labels, "test")
    
    multiColor <- vector("list", length = length(listColors)+1) # make a matching list with the colors to evaluate. 
    
    multiColor[[length(multiColor)]] <- vector("character", length = ncol(multiExpr[[1]]$data))# length(listColors[[1]])) # set last set of colors to empty strings 
    names(multiColor[[length(multiColor)]]) <- colnames(multiExpr[[1]]$data)
    
    # Cut module assignment genes down to the dataset genes
    for (i in 1:(length(multiColor)-1)) {
      colors_tmp <- vector(mode="character", length=ncol(multiExpr[[i]]$data))
      names(colors_tmp) <- colnames(multiExpr[[i]]$data)
      # TODO: can this handle NAs?
      
      colors_tmp[match(names(listColors[[i]]), names(colors_tmp))] <- listColors[[i]][names(listColors[[i]]) %in% names(colors_tmp)]
      colors_tmp[nchar(colors_tmp)==0] <- "grey"
      multiColor[[i]] <- colors_tmp
    }
    #  Individual expression sets and their module labels are matched using names of the corresponding components in multiExpr and multiColor
    names(multiColor) <- names(multiExpr) # i.e. =  c(labels, "test")
    
  } else if (length(listDatExpr) > 1) {
    
    # multicolor must have length 1
    
    if (verbose >= 1) {
      print("Measuring the preservation of a module assignment between a reference and one or more test datasets..")
    }
    
    # Set up the multi-set format for the modulePreservation function
    multiExpr <- vector(mode="list", length=length(listDatExpr))  # convert the list of expression data sets to the required format
    
    for (i in 1:length(multiExpr)) {
      multiExpr[[i]]$data  <- as.matrix(listDatExpr[[i]])
    }
    
    names(multiExpr) <- if (!is.null(names(listDatExpr))) names(listDatExpr) else 1:length(multiExpr) 
    
    #########
    
    multiColor <- listColors # rename for consistency but otherwise leave unchanged
    
    names(multiColor) <- names(multiExpr)[1] # listColors and hence multiColor is only allowed to have length 1
    
    ########
    
    referenceNetworks = rep(1, length(multiExpr)-1) # The reference dataset is always the first
    
    # TODO: check testNetworks	
    # a list with one component per each entry in referenceNetworks above, 
    # giving the test networks in which to evaluate module preservation 
    # for the corresponding reference network. If not given, 
    # preservation will be evaluated in all networks (except each reference network). 
    # If referenceNetworks is of length 1, testNetworks can also be a vector 
    # (instead of a list containing the single vector).
    testNetworks = as.list(2:length(multiExpr)) # the test sets are the second, third, ..., last 
  }
  
  if (parallelCalculation) enableWGCNAThreads(nThreads = nThreads)
  
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
  
  try(disableWGCNAThreads())
  
  invisible(gc())
  
  return(result)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

load_obj <- function(f) {
  # Usage: Loads (gzip compressed) file from .RData, .RDS, .loom, .csv, .txt, .tab, .delim
  #
  # Args: 
  #   f: path to file
  # returns: 
  #   RObject
  
  compressed = F
  
  if (grepl(pattern = "\\.gz|\\.gzip", x=f))  {
    compressed <- T
    f = paste0("gzfile('",f,"')")
  }
  
  if (grepl(pattern = "\\.RDS", x = f, ignore.case = T)) {
    out <- readRDS(file=if(compressed) eval(parse(text=f)) else f)
  } else if (grepl(pattern="\\.RData|\\.Rda", x=f, ignore.case = T)) { 
    env <- new.env()
    nm <- load(f, env)[1]
    out <- env[[nm]]
  } else if (grepl(pattern="\\.loom", x=f)) {
    out <- connect(filename=f, mode = "r+")
  } else if (grepl(pattern = "\\.csv", x=f)) {
    out <- read.csv(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
  } else if (grepl(pattern = "\\.tab|\\.tsv", x=f)) {
    out <- read.table(file=if(compressed) eval(parse(text=f)) else f, sep="\t", stringsAsFactors = F, quote="", header=T) 
  } else if (grepl(pattern = "\\.txt", x=f)) {
    out <- read.delim(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
  } 
  out
}
# load_obj <- function(f) {
#   # Utility function for loading an object stored in any of 
#   # .RData, .RDS, .loom, .csv, .txt, .tab, .delim
#   # File may be gzip compressed 
#   
#   compressed = F
#   
#   if (grepl(pattern = "\\.gz|\\.gzip", x=f))  {
#     compressed <- T
#     f = paste0("gzfile('",f,"')")
#   }
#   
#   if (grepl(pattern = "\\.RDS", x = f, ignore.case = T)) {
#     readRDS(file=if(compressed) eval(parse(text=f)) else f)
#   } else if (grepl(pattern="\\.RData|\\.rda", x=f, ignore.case = T)) { 
#     env <- new.env()
#     nm <- load(f, env)[1]
#     env[[nm]]
#   } else if (grepl(pattern="\\.loom", x=f)) {
#     connect(filename=f, mode = "r+")
#   } else if (grepl(pattern = "\\.csv", x=f)) {
#     read.csv(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
#   } else if (grepl(pattern = "\\.tab", x=f)) {
#     read.table(file=if(compressed) eval(parse(text=f)) else f, sep="\t", stringsAsFactors = F, quote="", header=T) 
#   } else if (grepl(pattern = "\\.txt", x=f)) {
#     read.delim(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
#   }
# }

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
                             pvalThreshold) {
  
  # Rather than for parallelising over modules within a set of modules, parallelise over several sets of colors (produced by comparing parameters)
  # It calls PPI_innver_for_vec as a subroutine
  # Args:
  #   colors
  #   pkMs
  #   STRINGdb_species
  #   PPI_pkM_threshold: numeric scalar
  #   pvalThreshold: numeric scalar
  #   
  # Returns: 
  #   colors_PPI: colors where modules that did not satisfy the PPI threshold are set to grey
  
  #suppressPackageStartupMessages(library(STRINGdb)) 
  
  module_PPI <- NULL
  module_PPI_signif <- NULL
  unique_colors <- NULL
  unique_colors_PPI <- NULL
  colors_PPI <- NULL
  
  unique_colors <- unique(colors[pkMs>PPI_pkM_threshold])
  
  if (length(unique_colors)>0){ # grey isn't included since the pkMs will not pass the threshold
    
    # generate a STRINGdb object instance
    string_db <- STRINGdb$new(version="10", 
                              species = STRINGdb_species, 
                              score_threshold=0, 
                              input_directory="") 
    
    # For each module, check STRING_db enrichment
    sapply(unique_colors, function(x) PPI_inner_for_vec(color=x, 
                                                        unique_colors = unique_colors, 
                                                        colors = colors, 
                                                        string_db = string_db), simplify=T) %>% t() -> module_PPI 
    
    module_PPI <- data.frame("colors"= unique_colors, module_PPI, stringsAsFactors=F, row.names = NULL)
    module_PPI$q.value <- as.numeric(module_PPI$q.value)
    module_PPI$expected.interactions <- as.numeric(module_PPI$expected.interactions)
    
  } else {
    module_PPI <- data.frame(colors= "grey", q.value=1, expected.interactions=0)
  }
  # FILTER MODULES ON PPI ENRICHMENT  
  
  if (!is.null(module_PPI)) {
    if (sum(module_PPI$q.value < pvalThreshold)>0) {
      module_PPI_signif <- module_PPI[module_PPI$q.value < pvalThreshold,]
      unique_colors_PPI = unique_colors[module_PPI$q.value < pvalThreshold]
      genes_PPI_idx <- colors %in% unique_colors_PPI
      colors_PPI <- colors
      colors_PPI[!genes_PPI_idx] <- "grey"
    } else {
      colors_PPI <- rep("grey", length(colors))
    }
  } else {
    module_PPI_signif <- NULL
    colors_PPI <- rep("grey", length(colors))
  }  
  
  return(list("colors_PPI" = colors_PPI, "module_PPI" = module_PPI, "module_PPI_signif" = module_PPI_signif))
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

PPI_inner_for_vec <- function(color, unique_colors, colors, string_db) {
  # Usave: 
  #   Parallelise over modules within a set of modules
  #   utility function to parallelise the PPI STRINGdb call
  # args: 
  #   color should be one of the unique_colors
  # value:
  #   module_PPI: a list with named entries 'p-value' and 'expected interactions'
  
  ppi <- data.frame(gene = names(colors[colors==color])) # extract the genes with the corresponding color to dataframe
  module_PPI <- list('q.value' = 1, 'expected.interactions' = NA)
  
  tryCatch({
    
    example1_mapped <- string_db$map(ppi, 'gene', removeUnmappedRows = TRUE ) # Check the dataframe genes' PPI. May produce error: we couldn't map to STRING 100% of your identifiers
    # (ctd.) .. Error in if (hitList[i] %in% V(ppi_network)$name) { : argument is of length zero
    hits <- example1_mapped$STRING_id
    module_PPI['q.value'] = p.adjust(string_db$get_ppi_enrichment(hits)$enrichment,method = 'fdr',n = length(names(colors)))
    module_PPI['expected.interactions'] = string_db$get_ppi_enrichment(hits)$lambda
    
  }, error = function(c) {
    module_PPI <- list('q.value'=1, 'expected.interactions' = NA) # probably unnecessary since defined above but just in case it has been changed in place
  }) 
  
  return(module_PPI)
  
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

mapMMtoHs = function(modulekM,
                     colors = NULL,
                     log_dir,
                     data_prefix,
                     run_prefix,
                     mapping_orthology) {
  
  if (!is.null(modulekM$genes)) modulekM$genes <- NULL
  
  log_not_mapped_filepath = paste0(log_dir,"genes_orthology_not_mapped_",data_prefix,"_", run_prefix, ".tab")
  
  mapping = data.frame(ensembl.mouse=row.names(modulekM))
  # orthology mapping
  mapping$ensembl.human = mapping_orthology$ensembl_gene_id[ match(mapping$ensembl.mouse, mapping_orthology$mmusculus_homolog_ensembl_gene) ]
  #mapping$ensembl.human[mapping$ensembl.human == ""] = NA
  df_not_mapped = mapping[is.na(mapping$ensembl.human),]
  append = !file.exists(log_not_mapped_filepath)
  col.names = !append 
  write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t", col.names = col.names, row.names=F, append=append)
  
  modulekM$ensembl = mapping$ensembl.human
  
  # Also name colors
  if (!is.null(colors)) {
    colors_ens <- colors[!is.na(modulekM$ensembl)]
    names(colors_ens) <- na.omit(mapping$ensembl.human)
  }
  
  modulekM <- na.omit(modulekM) # does this not produce a na.omit out object?
  
  ### 180508_v.18_dev2
  #tmp = within(modulekM, rm("symbol","ensembl"))
  tmp = within(modulekM, rm("ensembl"))
  ###
  
  # Average duplicated gene IDs
  modulekM_ens <-aggregate(tmp, by=list(modulekM$ensembl),FUN=mean, na.rm=TRUE)
  rownames(modulekM_ens) = modulekM_ens$Group.1
  modulekM_ens = within(modulekM_ens, rm("Group.1"))
  
  if (!is.null(colors)) colors_ens <- colors_ens[names(colors_ens) %in% rownames(modulekM_ens)] else NULL 
  
  modulekM_ens$genes <- NULL
  
  prop.mapped <- round(sum(!is.na(mapping$ensembl.human))/nrow(mapping),2)
  
  return(list(kM=modulekM_ens, colors=colors_ens, prop.mapped=prop.mapped))
  
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
                     gwas,
                     test_type = c("full_kME_spearman", "module_genes_spearman", "module_genes_t"),
                     genes_background_ensembl_hs,
                     colors_hs = NULL) {
  # Usage: Subroutine to calculate spearman's correlation between gene module membership and GWAS gene significance
  # Args: 
  # Returns:
  
  tryCatch({
    unique_colors = colnames(modulekM)
    table.kM.cor.p = table.kM.cor.r = table.kM.cor.emp.p <- matrix(NA,nrow=length(unique_colors),ncol=length(gwas)) 
    rownames(table.kM.cor.r) = rownames(table.kM.cor.p) = rownames(table.kM.cor.emp.p) = unique_colors
    colnames(table.kM.cor.r) = colnames(table.kM.cor.p) = colnames(table.kM.cor.emp.p) = names(gwas) 
    
    for (col in unique_colors) {
      for (j in 1:length(gwas)) {
        genes = if (test_type=="full_kME_spearman") intersect(rownames(modulekM), gwas[[j]]$gene_name) else if (test_type %in% c("module_genes_spearman", "module_genes_t")) intersect(names(colors_hs)[colors_hs==col], gwas[[j]]$gene_name)  
        
        if (test_type %in% c("full_kME_spearman", "module_genes_spearman")) {
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
          
        } else if (test_type == "module_genes_t") {
          
          x = -log10(gwas[[j]]$P[match(genes, gwas[[j]]$gene_name)])
          shared_genes_tmp <- intersect(genes_background_ensembl_hs, gwas[[j]]$gene_name)
          #shared_genes_tmp <- shared_genes_tmp[!shared_genes_tmp %in% genes]
          y = -log10(gwas[[j]]$P[match(shared_genes_tmp, gwas[[j]]$gene_name)])
          
          test = tryCatch({t.test(x=x,y=y, alternative = "g")}, error = function(c) {return(list(p.value=1))})
          
          table.kM.cor.p[col,j] <- test$p.value
          table.kM.cor.r[col,j] <- 1e-5 # to avoid 0, which will lead to problems when taking sign of the correlation to adjust p-values
          table.kM.cor.emp.p[col,j] <- 1        
        }
      }
    }
    
    rownames(table.kM.cor.r) <- rownames(table.kM.cor.p) <- rownames(table.kM.cor.emp.p) <- paste0(cellType, "__", rownames(table.kM.cor.p))
    
    return(list('p.val'= table.kM.cor.p, 'corrCoef' = table.kM.cor.r, 'emp.p.val' = table.kM.cor.emp.p))}, 
    error = function(c) {warning(paste0("kM_magma failed for ", cellType))})
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

# disused

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
      
      tryCatch({
        glm.out <- glm(binaryMat$GeneSet~binaryMat$Module, family=binomial())
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

cellModEmbed <- function(datExpr, 
                         colors=NULL, 
                         latentGeneType,
                         cellType=NULL,
                         kMs=NULL,
                         excludeGrey=T,
                         dissTOM=NULL) {
  # datExpr should be cell x gene matrix or data.frame with ALL genes and ALL cells in the analysis
  # colors is a character vector with gene names
  # the genes should be ordered identically
  
  # prepare celltype character for naming columns
  tryCatch({
    if (is.null(cellType)) cellType_prefix <- "" else cellType_prefix <- paste0(cellType, "__")
    if (latentGeneType == "ME") { 
      
      colors_full <- rep("grey", times = ncol(datExpr))
      names(colors_full) <- colnames(datExpr)
      
      for (col in unique(colors)) {
        colors_full[names(colors_full) %in% names(colors[colors==col])] <- col
      }
      
      embed_mat <- moduleEigengenes_uv(expr=datExpr, 
                                              colors = colors_full, 
                                              impute = TRUE, 
                                              nPC = 1, 
                                              align = "along average", 
                                              excludeGrey = excludeGrey, 
                                              subHubs = TRUE,
                                              trapErrors = FALSE, 
                                              returnValidOnly = trapErrors, 
                                              scale = TRUE,
                                              verbose = 0, 
                                              indent = 0)$eigengenes
      
      # list_datExpr <- lapply(unique(colours)[unique(colours)!="grey"], function(x) datExpr[,match(names(colours)[colours==x], colnames(datExpr))])
      # embed_mat <- as.matrix(sapply(list_datExpr, function(x) prcomp_irlba(t(x), n=1, retx=T)$rotation, simplify = T))
      
      colnames(embed_mat) <- paste0(cellType_prefix, gsub("ME", "", colnames(embed_mat)))
      
    } else if (latentGeneType == "IM") {
      
      list_datExpr <- lapply(colnames(kMs), function(x) datExpr[ , match(names(colors[colors==x]), colnames(datExpr))])
      
      list_kMs <- lapply(colnames(kMs), function(module) kMs[match(names(colors)[colors==module], rownames(kMs)),module])
      
      embed_mat <- as.matrix(mapply(function(x,y) x %*% as.matrix(y),
                                    x = list_datExpr,
                                    y = list_kMs,
                                    SIMPLIFY=T))
      # datExpr_1 = datExpr[, match(rownames(kMs), colnames(datExpr), nomatch=0) [match(rownames(kMs), colnames(datExpr),  nomatch=0)>0]]
      # embed_mat <- as.matrix(datExpr_1) %*% as.matrix(kMs)
      colnames(embed_mat) <- paste0(cellType_prefix, colnames(kMs)) 
      rownames(embed_mat) <- rownames(datExpr)
    } 
    return(embed_mat)
  }, error = function(c) {warning(paste0("cellModEmbed failed for ", cellType))})
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
  pdf(sprintf("%s%s_%s_mean.celltype.eig.expr_%s.pdf", plots_dir, data_prefix, run_prefix, flag_date ),h=12,w=12)
  draw(ht.mean.celltype.eig.expr)
  dev.off()
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

moduleEigengenes_uv <- function (expr, colors, 
                                 impute = TRUE, 
                                 nPC = 1, 
                                 align = "along average", 
                                 excludeGrey = FALSE, 
                                 grey = if (is.numeric(colors)) 0 else "grey", 
                                 subHubs = TRUE, 
                                 trapErrors = FALSE, 
                                 returnValidOnly = trapErrors, 
                                 softPower = 8, 
                                 scale = TRUE, 
                                 verbose = 0, 
                                 indent = 0) 
{
  # moduleEigengenes function modified to also return left (u) eigenvectors of length n_gene_in_module
  
  spaces = indentSpaces(indent)
  if (verbose == 1) 
    printFlush(paste(spaces, "moduleEigengenes: Calculating", 
                     nlevels(as.factor(colors)), "module eigengenes in given set."))
  if (is.null(expr)) {
    stop("moduleEigengenes: Error: expr is NULL. ")
  }
  if (is.null(colors)) {
    stop("moduleEigengenes: Error: colors is NULL. ")
  }
  if (is.null(dim(expr)) || length(dim(expr)) != 2) 
    stop("moduleEigengenes: Error: expr must be two-dimensional.")
  if (dim(expr)[2] != length(colors)) 
    stop("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene).")
  if (is.factor(colors)) {
    nl = nlevels(colors)
    nlDrop = nlevels(colors[, drop = TRUE])
    if (nl > nlDrop) 
      stop(paste("Argument 'colors' contains unused levels (empty modules). ", 
                 "Use colors[, drop=TRUE] to get rid of them."))
  }
  if (softPower < 0) 
    stop("softPower must be non-negative")
  alignRecognizedValues = c("", "along average")
  if (!is.element(align, alignRecognizedValues)) {
    printFlush(paste("ModulePrincipalComponents: Error:", 
                     "parameter align has an unrecognised value:", align, 
                     "; Recognized values are ", alignRecognizedValues))
    stop()
  }
  maxVarExplained = 10
  if (nPC > maxVarExplained) 
    warning(paste("Given nPC is too large. Will use value", 
                  maxVarExplained))
  nVarExplained = min(nPC, maxVarExplained)
  modlevels = levels(factor(colors))
  if (excludeGrey) 
    if (sum(as.character(modlevels) != as.character(grey)) > 0) {
      modlevels = modlevels[as.character(modlevels) != 
                              as.character(grey)]
    } else {
    stop(paste("Color levels are empty. Possible reason: the only color is grey", 
               "and grey module is excluded from the calculation."))
    }
  PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
                                ncol = length(modlevels)))
  
  PrinComps_l = vector(mode="list", length= length(modlevels))  
  
  # PrinComps_l = data.frame(matrix(NA, nrow = dim(expr)[[2]], 
  #                               ncol = length(modlevels)))
  
  averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
                               ncol = length(modlevels)))
  
  averExpr_l = vector(mode="list", length= length(modlevels)) 
  # averExpr_l = data.frame(matrix(NA, nrow = dim(expr)[[2]], 
  #                              ncol = length(modlevels)))
  
  
  varExpl = data.frame(matrix(NA, nrow = nVarExplained, ncol = length(modlevels)))
  validMEs = rep(TRUE, length(modlevels))
  validAEs = rep(FALSE, length(modlevels))
  isPC = rep(TRUE, length(modlevels))
  isHub = rep(FALSE, length(modlevels))
  validColors = colors
  names(PrinComps) = paste(moduleColor.getMEprefix(), modlevels, 
                           sep = "")
  names(PrinComps_l) = modlevels
  
  names(averExpr) = paste("AE", modlevels, sep = "")
  names(averExpr_l) = paste("AE", modlevels, sep = "")
  
  for (i in c(1:length(modlevels))) {
    if (verbose > 1) 
      printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", 
                       modlevels[i]))
    modulename = modlevels[i]
    restrict1 = as.character(colors) == as.character(modulename)
    if (verbose > 2) 
      printFlush(paste(spaces, " ...", sum(restrict1), 
                       "genes"))
    datModule = as.matrix(t(expr[, restrict1])) # datModule is a genes_mod * cells = n * p matrix
    n = dim(datModule)[1]
    p = dim(datModule)[2]
    
    svd_1 = try({
      if (nrow(datModule) > 1 && impute) {
        seedSaved = FALSE
        if (exists(".Random.seed")) {
          saved.seed = .Random.seed
          seedSaved = TRUE
        }
        if (any(is.na(datModule))) {
          if (verbose > 5) 
            printFlush(paste(spaces, " ...imputing missing data"))
          datModule = impute.knn(datModule, k = min(10, 
                                                    nrow(datModule) - 1))
          try({
            if (!is.null(datModule$data)) 
              datModule = datModule$data
          }, silent = TRUE)
        }
        if (seedSaved) 
          .Random.seed <<- saved.seed
      }
      if (verbose > 5) 
        printFlush(paste(spaces, " ...scaling"))
      if (scale) 
        datModule = t(scale(t(datModule)))

      if (verbose > 5) 
        printFlush(paste(spaces, " ...calculating SVD"))
      
      svd_1 = svd(datModule, nu = min(n, p, nPC), nv = min(n,p, nPC))
      
      if (verbose > 5) 
        printFlush(paste(spaces, " ...calculating PVE"))
      veMat = cor(svd_1$v[, c(1:min(n, p, nVarExplained))], 
                  t(datModule), use = "p")
      varExpl[c(1:min(n, p, nVarExplained)), i] = rowMeans(veMat^2, 
                                                           na.rm = TRUE)
      svd_1
    }, silent = TRUE)
    
    pc <- svd_1$v[, 1]
    pc_l <- svd_1$u[, 1]
    
    ###############################
    
    if (class(pc) == "try-error") {
      if ((!subHubs) && (!trapErrors)) 
        stop(pc)
      if (subHubs) {
        if (verbose > 0) {
          printFlush(paste(spaces, " ..principal component calculation for module", 
                           modulename, "failed with the following error:"))
          printFlush(paste(spaces, "     ", pc, spaces, 
                           " ..hub genes will be used instead of principal components."))
        }
        isPC[i] = FALSE
        pc = try({
          scaledExpr = scale(t(datModule))
          covEx = cov(scaledExpr, use = "p")
          covEx[!is.finite(covEx)] = 0
          modAdj = abs(covEx)^softPower
          kIM = (rowMeans(modAdj, na.rm = TRUE))^3
          if (max(kIM, na.rm = TRUE) > 1) 
            kIM = kIM - 1
          kIM[is.na(kIM)] = 0
          hub = which.max(kIM)
          alignSign = sign(covEx[, hub])
          alignSign[is.na(alignSign)] = 0
          isHub[i] = TRUE
          pcxMat = scaledExpr * matrix(kIM * alignSign, 
                                       nrow = nrow(scaledExpr), ncol = ncol(scaledExpr), 
                                       byrow = TRUE)/sum(kIM)
          pcx = rowMeans(pcxMat, na.rm = TRUE)
          varExpl[1, i] = mean(cor(pcx, t(datModule), 
                                   use = "p")^2, na.rm = TRUE)
          pcx
        }, silent = TRUE)
      }
    }
    if (class(pc) == "try-error") {
      if (!trapErrors) 
        stop(pc)
      if (verbose > 0) {
        printFlush(paste(spaces, " ..ME calculation of module", 
                         modulename, "failed with the following error:"))
        printFlush(paste(spaces, "     ", pc, spaces, 
                         " ..the offending module has been removed."))
      }
      warning(paste("Eigengene calculation of module", 
                    modulename, "failed with the following error \n     ", 
                    pc, "The offending module has been removed.\n"))
      validMEs[i] = FALSE
      isPC[i] = FALSE
      isHub[i] = FALSE
      validColors[restrict1] = grey
    } else {
      PrinComps[, i] = pc
      ae = try({
        if (isPC[i]) 
          scaledExpr = scale(t(datModule))
        averExpr[, i] = rowMeans(scaledExpr, na.rm = TRUE)
        if (align == "along average") {
          if (verbose > 4) 
            printFlush(paste(spaces, " .. aligning module eigengene with average expression."))
          corAve = cor(averExpr[, i], PrinComps[, i], 
                       use = "p")
          if (!is.finite(corAve)) 
            corAve = 0
          if (corAve < 0) 
            PrinComps[, i] = -PrinComps[, i]
        }
        0
      }, silent = TRUE)
      if (class(ae) == "try-error") {
        if (!trapErrors) 
          stop(ae)
        if (verbose > 0) {
          printFlush(paste(spaces, " ..Average expression calculation of module", 
                           modulename, "failed with the following error:"))
          printFlush(paste(spaces, "     ", ae, spaces, 
                           " ..the returned average expression vector will be invalid."))
        }
        warning(paste("Average expression calculation of module", 
                      modulename, "failed with the following error \n     ", 
                      ae, "The returned average expression vector will be invalid.\n"))
      }
      validAEs[i] = !(class(ae) == "try-error")
      
      ########## ADDED ##########
      # Also align left singular components with average gene expression
      names(pc_l) <- rownames(datModule)
      PrinComps_l[[i]] = pc_l
      
      # if (isPC[i]) {
      #   scaledExpr_l = scale(t(datModule)) # a cell * gene_mod_i matrix
      # }
      averExpr_l[[i]] = colMeans(scaledExpr, na.rm = TRUE) # for each gene averaging over cells
      if (align == "along average") {
        if (verbose > 4) 
          printFlush(paste(spaces, " .. aligning reverse module eigengene with average gene expression."))
        corAve_l = cor(averExpr_l[[i]], PrinComps_l[[i]], 
                       use = "p")
        if (!is.finite(corAve_l)) 
          corAve_l = 0
        if (corAve_l < 0) 
          PrinComps_l[[i]] = -PrinComps_l[[i]]
      }
      
      #################
    }
  }
  allOK = (sum(!validMEs) == 0)
  if (returnValidOnly && sum(!validMEs) > 0) {
    PrinComps = PrinComps[, validMEs]
    #### ADDED #####
    PrinComps_l = PrinComps[validMEs]
    ################
    averExpr = averExpr[, validMEs]
    #### ADDED #####
    averExpr_l = averExpr_l[validMEs]
    ################
    varExpl = varExpl[, validMEs]
    validMEs = rep(TRUE, times = ncol(PrinComps))
    isPC = isPC[validMEs]
    isHub = isHub[validMEs]
    validAEs = validAEs[validMEs]
  }
  allPC = (sum(!isPC) == 0)
  allAEOK = (sum(!validAEs) == 0)
  list(eigengenes = PrinComps, 
       ##### ADDED ####
       u = PrinComps_l,
       #####
       averageExpr = averExpr, 
       ##### ADDED ####
       averageExpr_l = averExpr_l, 
       ################
       
       varExplained = varExpl, 
       nPC = nPC, 
       validMEs = validMEs, 
       validColors = validColors, 
       allOK = allOK, 
       allPC = allPC, 
       isPC = isPC, 
       isHub = isHub, 
       validAEs = validAEs, 
       allAEOK = allAEOK)
}

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# install and require packages using a function (allows for automation)
ipak <- function(pkgs){
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new.pkgs)) {
    sapply(new.pkgs, function(pkg) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
      }, warning = function(war) {
        tryCatch({
          source("https://bioconductor.org/biocLite.R")
          biocLite(pkg, suppressUpdates = T)
        }, error = function(err1) {
          warning(paste0(pkg, " encountered the error: ", err1))
          dependency <- gsub("\\W|ERROR: dependency | is not available for package.*", "", err)
          ipak(dependency)
        })
      } ,
      error = function(err)
      {
        dependency <- gsub("\\W|ERROR: dependency | is not available for package.*", "", err)
        ipak(dependency)
      })
    })
  }
  suppressPackageStartupMessages(sapply(pkgs, require, character.only = TRUE))
  failed <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(failed)>0) warning(paste0(paste0(failed, collapse = " "), " failed to install"))
}

############################################################################################################################################################
############################################################## SPECIFITY ###################################################################################
############################################################################################################################################################
# 
# specificity.index(pSI.in = t(datExpr), 
#                   pSI.in.filter, 
#                   bts = 50, p_max = 0.1, e_min = 0.3, hist = FALSE, SI = FALSE)
#  returns a 
# 
# 
# mod_specificity <- function(modEmbedMat, module) {
#   # usage: Get specificity values for a set of modules across a set of celltypes / tissues etc
#   # args
#   #   modEmbedMat: cellType (or tissue etc) x module embedding matrix
#   #   module: a vector of gene weights with gene names
#   # value
#   #   - s matrix of dim == dim(modEmbedMat) with the specificity scores 
#   #   
#   # dependencies: 
#   #   pSI package https://cran.r-project.org/web/packages/pSI/pSI.pdf
# }

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

ipak <- function(pkgs){
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new.pkgs)) {
    sapply(new.pkgs, function(pkg) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
      }, warning = function(war) {
        tryCatch({
          source("https://bioconductor.org/biocLite.R")
          biocLite(pkg, suppressUpdates = T)
        }, error = function(err1) {
          warning(paste0(pkg, " encountered the error: ", err1))
          dependency <- gsub("\\W|ERROR: dependency | is not available for package.*", "", err)
          ipak(dependency)
        })
      } ,
      error = function(err)
      {
        dependency <- gsub("\\W|ERROR: dependency | is not available for package.*", "", err)
        ipak(dependency)
      })
    })
  }
  suppressPackageStartupMessages(sapply(pkgs, require, character.only = TRUE))
  failed <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(failed)>0) warning(paste0(paste0(failed, collapse = " "), " failed to install"))
}

######################################################################
########################### GENE REMAP ###############################
######################################################################

gene_map <- function(df,
                     idx_gene_column = NULL,
                     mapping,
                     from="hgnc", 
                     to="ensembl",
                     replace = F,
                     na.rm = T) {
  # args:
  #   df: data.frame with a column or rownames containing genes to remap
  #   idx_gene_column: df gene column. NULL implies rownames (implement later)
  #   mapping: data.frame or matrix with columns corresponding to 'from' and 'to' arguments (using grep matching)
  #   replace: boolean; T to replace original gene names, F to add a new column to the data.frame
  # value:
  #   df with new gene names in place of or in addition to the original gene names
  
  # usage: 
  # df_test <- gene_map(df=load_obj("/projects/jonatan/tmp-mousebrain/tables/mousebrain_Vascular_ClusterName_1_PER3_kIM.csv"),
  #                     idx_gene_column =1,
  #                     mapping=load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"),
  #                     from="ensembl", 
  #                     to="gene_name_optimal",
  #                     replace=T,
  #                     na.rm=T)
  
  orig_class <- class(df)
  
  df <- as.data.frame(df)
  
  if (is.null(idx_gene_column) & replace==T) warning("Duplicate genes will be averaged and merged to keep row.names unique")
  
  if (!from %in% colnames(mapping)) from <- grep(pattern=from, x = colnames(mapping), ignore.case=T, value = T)
  if (!to %in% colnames(mapping)) to <- grep(pattern=to, x = colnames(mapping), ignore.case=T, value = T)
  
  stopifnot(length(from)>0 & length(to)>0)
  
  genes_from <- if (is.null(idx_gene_column)) rownames(df) else df[[idx_gene_column]]
  
  genes_to <- mapping[[to]][match(genes_from, mapping[[from]])]
  
  
  if (replace) { # remove NAs
    if (is.null(idx_gene_column)) {
      # average identical gene names to ensure unique row names 
      df_aggr <- aggregate(df, by= list(genes_to), FUN=mean, na.rm=T)
      rownames(df_aggr) <- df_aggr[["Group.1"]]
      df <- within(df_aggr, rm("Group.1"))
    } else {
      if (na.rm) {
        df <- df[!is.na(genes_to),]
        genes_to <- genes_to[!is.na(genes_to)]
      }
      df[[idx_gene_column]] <- genes_to
      colnames(df)[idx_gene_column] <- to
    }
  }  else {
    if (na.rm) {
      df <- df[!is.na(genes_to),]
      df[[to]] <- genes_to[!is.na(genes_to)]
    }
    df[[to]] <- genes_to  
  }
  
  if (orig_class != "data.frame") {
    class(df) <- orig_class
  }
  
  return(df)
}

#############

signedKME = function (datExpr, datME, outputColumnName = "kME", corFnc = "cor", 
                      corOptions = "use = 'p'") 
{
  datExpr = data.frame(datExpr)
  datME = data.frame(datME)
  output = list()
  if (dim(as.matrix(datME))[[1]] != dim(as.matrix(datExpr))[[1]]) 
    stop("Number of samples (rows) in 'datExpr' and 'datME' must be the same.")
  varianceZeroIndicatordatExpr = as.vector(apply(as.matrix(datExpr), 
                                                 2, var, na.rm = TRUE)) == 0
  varianceZeroIndicatordatME = as.vector(apply(as.matrix(datME), 
                                               2, var, na.rm = TRUE)) == 0
  if (sum(varianceZeroIndicatordatExpr, na.rm = TRUE) > 0) 
    warning("Some genes are constant. Hint: consider removing constant columns from datExpr.")
  if (sum(varianceZeroIndicatordatME, na.rm = TRUE) > 0) 
    warning(paste("Some module eigengenes are constant, which is suspicious.\n", 
                  "    Hint: consider removing constant columns from datME."))
  no.presentdatExpr = as.vector(apply(!is.na(as.matrix(datExpr)), 
                                      2, sum))
  if (min(no.presentdatExpr) < 4) 
    warning(paste("Some gene expressions have fewer than 4 observations.\n", 
                  "    Hint: consider removing genes with too many missing values or collect more arrays."))
  corExpr = parse(text = paste("data.frame(", corFnc, "(datExpr, datME ", 
                               prepComma(corOptions), "))"))
  output = eval(corExpr)
  output[no.presentdatExpr < 4, ] = NA
  names(output) = paste(outputColumnName, substring(names(datME), 
                                                    first = 3, last = 100), sep = "")
  
  dimnames(output)[[1]] = names(datExpr)
  output
}

#############

detectCores_plus <- function(Gb_max=250, 
                             additional_Gb=1) {
  # args: 
  #  Gb_max: ceiling on session memory usage in Gb, assuming that each worker duplicates the session memory
  #  additional_Gb: max additional memory requirement for new (temporary) objects created within a parallel session 
  # value:
  #   n_cores (integer)
  obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x))))) / 1024^3)
  max(1, min(detectCores(), Gb_max %/% (obj_size_Gb + additional_Gb))-1) 
}

############################################################################################################################################################
########################################################## safeParallel ####################################################################################
############################################################################################################################################################

safeParallel = function(fun, args, simplify=F, MARGIN=NULL, n_cores=NULL, Gb_max=NULL, outfile=NULL,  ...) {
  # summary: calls the appropriate parallel computing function, with load balancing, 
  #          and falls back on vectorised equivalent if makeCluster hangs or the parallelised computation fails.
  #          simplify defaults to FALSE
  # args:
  #   fun: function to run in parallel. The (unevaluated) object
  #   args: a named list of arguments for fun to iterate over
  #   n_cores: n_cores. If NULL, a safe estimate is made based on size of objects in the global environment 
  #        and the length of the iterables in args
  #   outfile: passed to the parallelising function
  #   .. : additional arguments for fun to evaluate, but not iterate over
  # value:
  #   list (or if args contains "SIMPLIFY" = T, a vector or matrix); in case of failure, NULL
  
  if (length(args)==1) names(args) <- "X"
  
  list_out <- NULL
  
  if (is.null(n_cores)) {
    if (is.null(Gb_max)) Gb_max=200
    additional_Gb = max(as.numeric(sapply(args, FUN = function(x) object.size(x), simplify = T)))/1024^3
    obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
    n_cores <- max(1, min(detectCores()%/%3, Gb_max %/% (obj_size_Gb + additional_Gb))-1)
    n_cores <- min(n_cores, max(sapply(args, length)))
  }
  
  cl <-  if (!is.null(outfile)) try(makeCluster(spec=max(1,n_cores), type="FORK", timeout=30, outfile = outfile)) else try(makeCluster(spec=max(1,n_cores), type="FORK", timeout=30))
  #cl <- "hi"
  #class(cl) <- "try-error"
  if (!"try-error" %in% class(cl)) {
    if (length(args)>1) {
      fnc <- "clusterMap"
      args[[".scheduling"]] = c("dynamic")
      args[["SIMPLIFY"]] <- simplify
    } else {
      fnc <- "parLapplyLB"
      if (simplify) {
        fnc <- "parSapplyLB"
        args[["simplify"]] <- simplify
      }
      if (!is.null(MARGIN)) {
        fnc <- "parApplyLB"
        args[["MARGIN"]] <- MARGIN
        args[["simplify"]] <- NULL
        }
    }
    args[["fun"]] <- fun
    args[["cl"]] = cl
    
  } else if ("try-error" %in% class(cl)) {
    
    if (length(args)>1) {
      
      fnc = "mapply"
      args[["SIMPLIFY"]] <- simplify
      
    } else {
      
      fnc <- "lapply"
      if (simplify) {
        fnc <- "sapply"
        args[["simplify"]] = simplify
      }
      if (!is.null(MARGIN)) {
        fnc <- "apply"
        args[["MARGIN"]] <- MARGIN
        args[["simplify"]] <- NULL
        }
    }
    
    args[["FUN"]] <- fun
    
  }
  
  list_out <- tryCatch({ 
    do.call(what=fnc, args = args)
  }, error = function(err) {
    invisible(gc())
    
    if (fnc=="clusterMap") {
      fnc <- "mapply" 
    } else if (fnc=="parLapplyLB") {
      fnc <- "lapply"
    } else if (fnc=="parSapplyLB") {
      fnc <- "sapply"
      args[["SIMPLIFY"]] <- NULL
      args[["simplify"]] <- simplify
    } else if (fnc=="parApplyLB") {
      fnc <- "apply"
    }
    
    args[["cl"]] <- args[["fun"]] <- args[[".scheduling"]] <- NULL
    args[["FUN"]] <- fun
    do.call(what=fnc, args = args)
  })
  
  if (!"try-error" %in% class(cl)) try(stopCluster(cl))
  invisible(gc())
  
  list_out 
  
}
