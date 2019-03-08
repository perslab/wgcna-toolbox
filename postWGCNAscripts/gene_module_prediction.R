### gene module uni- and multivariate model script

# take modules as features
# make metadata as outcome to predict

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
# input data

make_option("--metadata_corr_cols", type="character", 
            help = "Specify seurat_obj@meta.data$... column(s) for which to compute correlations with gene modules. Takes a character with a vector in single (double) quotes of seurat_obj@meta.data column names in double (single) quotes, without whitspace, e.g. 'nUMI' or 'c('Sex','Age')'. For factor or character metadata, each levels is analysed as a dummy variable, so exercise caution.  [default %default]"),
make_option("--RAM_Gb_max", type="integer", default=250,
            help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)


metadata_corr_col <- opt$metadata_corr_col
if (!is.null(metadata_corr_col)) metadata_corr_col <- eval(parse(text=metadata_corr_col))

metadata_corr_filter_vals <- opt$metadata_corr_filter_vals
if (!is.null(metadata_corr_filter_vals)) metadata_corr_filter_vals <- eval(parse(text=metadata_corr_filter_vals))

RAM_Gb_max <- opt$RAM_Gb_max

######################################################################
############################ CONSTANTS ###############################
######################################################################


randomSeed <- 12345
set.seed(randomSeed)
######################################################################
######## EXTRACT METADATA AND CONVERT FACTORS TO MODEL MATRIX ########
######################################################################

# metadat_names <- colnames(seurat_obj@meta.data)
# 
# # Convert any character or factor meta.data to numeric dummy variables each level with its own numeric column
# if (!is.null(metadata_corr_col)) {
#   if (any(colnames(seurat_obj@meta.data) %in% metadata_corr_col)) {
#     metadata <- matrix(NA, nrow=nrow(seurat_obj@meta.data), ncol=1)
#     include <- seurat_obj@meta.data[,colnames(seurat_obj@meta.data) %in% metadata_corr_col, drop=F]
#     for (i in 1:ncol(include)) {
#       if (class(include[,i]) %in% c("factor", "character")) {
#         metadata <- cbind(metadata, factorToIndicator(include[,i, drop=T]))        
#       } else {
#         metadata <- cbind(metadata, include[,i, drop=F])
#       }
#     }
#     #rownames(metadata) <- rownames(seurat_obj@meta.data)
#     metadata <- metadata[,-1, drop=F]
#     if (!is.null(metadata_corr_filter_vals)) metadata <- metadata[, toupper(colnames(metadata)) %in% toupper(metadata_corr_filter_vals), drop = F]
#     metadata <- as.data.frame(metadata)
#     
#     # Filter out any metadata columns where all the values are identical
#     metadata <- metadata[apply(metadata, MARGIN=2, FUN = function(x) length(unique(x))>1)]
#     if (ncol(metadata) == 0) metadata <- NULL    
#     rownames(metadata) = rownames(seurat_obj@meta.data)
#   } else metadata <- NULL
# } else metadata <- NULL





#####################################################################
####### COMPUTE MODULE - METADATA CORRELATION IN EACH CELL CLUSTER ###
######################################################################

if (fuzzyModMembership=="kIM") {
  list_dissTOM_path <- dir(path = dirScratch, pattern = paste0(prefixData, "_", prefixRun, "_list_dissTOM"), full.names = T)
  list_dissTOM <- load_obj(list_dissTOM_path)
  list_dissTOM_gwas <- list_dissTOM[match(sNames_gwas, names(list_dissTOM))]
  rm(list_dissTOM)  
}

if (!is.null(metadata_corr_col)) {
  
  if (!is.null(metadata)) {
    
    message("Computing module-metadata correlation in each celltype")
    
    list_metadata <- lapply(list_datExpr_gwas, function(x) metadata[match(rownames(x), rownames(metadata)), , drop=F]) # get list of cell * metadata. The space between the commas is intentional!
    
    if (fuzzyModMembership=="kME") {
      # Compute correlation between metadata (columns) and eigengenes (columns). 
      list_mod_metadata_corr_rho <- mapply(function(x,y) cor(x=as.matrix(x), 
                                                             y=as.matrix(y), 
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
      list_embed_mat <- mapply(function(a,b,c,d,e) cellModEmbed(datExpr=a, 
                                                                colors=b, 
                                                                latentGeneType="IM",
                                                                cellType=c,
                                                                kMs=d,
                                                                dissTOM=e),
                               a = list_datExpr_gwas,
                               b = list_colors_gwas,
                               c = names(list_datExpr_gwas),
                               d = list_kMs_gwas,
                               e = list_dissTOM_gwas,
                               SIMPLIFY=F)
      
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
    
    # Compute p values
    list_mod_metadata_corr_pval <- mapply(function(x,y) WGCNA::corPvalueStudent(x, 
                                                                                n = nrow(y)), 
                                          x=list_mod_metadata_corr_rho, 
                                          y=list_metadata, 
                                          SIMPLIFY=F)
    
    # Convert p values into one big matrix in order to adjust p-values for the number of modules tested across *all* celltypes
    corr_pval <- Reduce(x=list_mod_metadata_corr_pval, f = cbind)
    
    # set NAs to 1
    corr_pval[is.na(corr_pval)] <- 1
    
    # Compute the false discovery rates
    corr_fdr <- apply(corr_pval, MARGIN=1, FUN = function(x) p.adjust(x, method = "fdr")) %>% as.matrix # apply outputs the vectors as columns
    if (ncol(corr_pval)>1) corr_fdr <- t(corr_fdr)
    corr_fdr.log <- -log10(corr_fdr)     
    colnames(corr_fdr.log) <- colnames(corr_pval)
    
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
    
    # Also single rho df
    corr_rho <- Reduce(x=list_mod_metadata_corr_rho, f=cbind)
    rownames(corr_rho) <- rownames(corr_pval)
    colnames(corr_rho) <- colnames(corr_pval)
    
  } 
  
} else {
  metadata = NULL
}

if (fuzzyModMembership=="kIM" ) rm(list_dissTOM_gwas)

##########################################################################
#### FILTER COLORS VECS, GENE AND KME LISTS FOR METADATA CORRELATIONS ####
##########################################################################

if (!is.null(metadata_corr_col) & !is.null(metadata_corr_filter_vals)) {
  
  if (!is.null(metadata)) {
    
    # Get a list of logical vectors indicating significant correlations
    list_idx_module_meta_sig <- lapply(list_mod_metadata_corr_fdr.log, function(x) apply(x, 2, function(y) any(y>-log10(pvalThreshold)))) # logical
    
    if (sum(sapply(list_idx_module_meta_sig, sum))==0) {
      
      warning("After filtering for metadata correlations no modules remain! Skipping metadata filtering step")
      sNames_meta <- sNames_gwas
      list_colors_meta <- list_colors_gwas
      list_module_meta <- list_module_gwas
      list_kMs_meta <- list_kMs_gwas
      list_MEs_meta <- list_MEs_gwas
      list_u_meta <- list_u_gwas 
      list_datExpr_meta <- list_datExpr_gwas
      
      metadata_corr_filter_vals = NULL
      metadata_corr_col = NULL        
      metadata=NULL
      
    } else {
      
      # Only keep significantly correlated modules within each celltype
      list_module_meta <- mapply(function(x,y) colnames(x[,y, drop=F]),x=list_mod_metadata_corr_fdr.log, y=list_idx_module_meta_sig, SIMPLIFY=F) 
      list_module_meta <- lapply(list_module_meta, function(x) gsub("^.*__", "", x)) %>% Filter(f=length)
      
      sNames_meta <- sNames_gwas[match(names(list_module_meta),sNames_gwas)] # ordered correctly
      
      # Keep the order
      list_module_meta <- list_module_meta[match(sNames_meta, names(list_module_meta))]
      
      list_datExpr_meta <- list_datExpr_gwas[match(sNames_meta, names(list_module_meta))]
      
      # reassign genes of filtered out modules to grey and remove any empty cell clusters
      list_colors_meta <- mapply(function(x,y) ifelse(x %in% y, yes = x, no = "grey"),
                                 x = list_colors_gwas[match(sNames_meta,names(list_colors_gwas))],
                                 y = list_module_meta,
                                 SIMPLIFY = F)
      
      # give gene names to color assignment vectors
      list_colors_meta <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(y), dimension = NULL), 
                                 x = list_colors_meta,
                                 y = list_colors_gwas[match(sNames_meta,names(list_colors_gwas))],
                                 SIMPLIFY = F)
      
      # Filter kM lists
      list_kMs_meta <- list_kMs_gwas[match(sNames_meta,names(list_kMs_gwas))]
      list_kMs_meta <- mapply(function(x,y) x[,match(y,colnames(x)), drop=F],x=list_kMs_meta, y = list_module_meta, SIMPLIFY=F) %>% Filter(f=length)
      
      # Filter left eigenvectors / kIMs 
      list_u_meta <- list_u_gwas[match(sNames_meta,names(list_u_gwas))]
      list_u_meta <- mapply(function(x,y) x[match(y,names(x))],
                            x = list_u_meta, 
                            y = list_module_meta,
                            SIMPLIFY=F) %>% Filter(f=length)
      
      # Filter ME lists
      if (fuzzyModMembership=="kME") {
        list_MEs_meta <- list_MEs_gwas[match(sNames_meta,names(list_MEs_gwas))]
        list_MEs_meta <- mapply(function(x,y) x[match(y,colnames(x))],
                                x = list_MEs_meta, 
                                y = list_module_meta,
                                SIMPLIFY=F) %>% Filter(f=length)
        
        list_u_meta <- list_u_gwas[match(sNames_meta,names(list_u_gwas))]
        list_u_meta <- mapply(function(x,y) x[match(y, names(x))],
                              x = list_u_meta, 
                              y = list_module_meta,
                              SIMPLIFY=F) %>% Filter(f=length)
        
        
      }  else {
        list_MEs_meta <- list_MEs_gwas 
        list_u_meta <- list_u_gwas
      }
    }
    
  } else {
    
    sNames_meta <- sNames_gwas
    list_colors_meta <- list_colors_gwas
    list_module_meta <- list_module_gwas
    list_kMs_meta <- list_kMs_gwas
    list_MEs_meta <- list_MEs_gwas
    list_u_meta <- list_u_gwas
    list_datExpr_meta <- list_datExpr_gwas
  }
  
} else if (is.null(metadata_corr_col) | is.null(metadata_corr_filter_vals)) {
  
  sNames_meta <- sNames_gwas
  list_colors_meta <- list_colors_gwas
  list_module_meta <- list_module_gwas
  list_kMs_meta <- list_kMs_gwas 
  list_MEs_meta <- list_MEs_gwas
  list_u_meta <- list_u_gwas
  list_datExpr_meta <- list_datExpr_gwas
}

# count number of enriched module per celltype for summary stats
n_modules_meta_enriched <- rep(NA, times=length(sNames_0))
names(n_modules_meta_enriched) <-sNames_0
if (!is.null(metadata_corr_col) & !is.null(metadata_corr_filter_vals)) n_modules_meta_enriched[names(n_modules_meta_enriched) %in% sNames_0] <- sapply(list_idx_module_meta_sig, sum) 
