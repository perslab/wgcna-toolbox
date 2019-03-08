


randomSeed <- 12345
set.seed(seed=randomSeed)

##########################################################################
############# COMPUTE MAGMA COMMON VARIANTS GWAS ENRICHMENT ##############
##########################################################################

# MAGMA docs: http://ctg.cncr.nl/software/MAGMA/doc/manual_v1.06.pdf
# TODO: Need to add non-protein coding genes to MAGMA

if (!is.null(magma_gwas_dir)) {
  
  message("Scoring modules for enrichment with genes linked by GWAS to phenotypes of interest")
  
  magma_test_type = "module_genes_t" # TODO: make this an argument?
  
  if (dataOrganism == "mmusculus") {
    
    # map mouse to human gene orthologs  
    mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
    
    outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "mapMMtoHs_par.txt")
    # cl <- makeCluster(min(n_cores, detectCores_plus(Gb_max = RAMGbMax, additional_Gb = 2)-1), type="FORK", outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "mapMMtoHs_par.txt"))
    fun = function(modulekM,colors_PPI_uniq) mapMMtoHs(modulekM = modulekM,
                                                       colors = colors_PPI_uniq,
                                                       dirLog = dirLog, 
                                                       prefixData = prefixData, 
                                                       prefixRun = prefixRun,
                                                       mapping_orthology = mapping_orthology)
    args = list(modulekM = list_kMs_PPI,
                colors_PPI_uniq =list_colors_PPI_uniq)
    list_MMtoHsmapping <- safeParallel(fun=fun, args=args, outfile=outfile, dirLog=dirLog, prefixData=prefixData, prefixRun=prefixRun, mapping_orthology=mapping_orthology)
    
    
    list_kMs_hs <- lapply(list_MMtoHsmapping, function(x) x$kM) 
    list_colors_hs <- lapply(list_MMtoHsmapping, function(x) x$colors)  
    
    vec_MMtoHsmapping_prop.mapped <- sapply(list_MMtoHsmapping, function(x) x$prop.mapped, simplify = T)
    
    names(list_kMs_hs) <- names(list_colors_hs) <- names(list_kMs_PPI)
    
    # Also map all 10x genes to homo sapiens
    
    genes_background_ensembl_mm <- read.csv(sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", dirTables, prefixData, prefixRun, dataOrganism))[["ensembl"]]
    genes_background_ensembl_hs <- mapping_orthology$ensembl_gene_id[match(genes_background_ensembl_mm, mapping_orthology$mmusculus_homolog_ensembl_gene)]
    genes_background_ensembl_hs <- as.character(na.omit(genes_background_ensembl_hs)) 
    
  } else if (dataOrganism == "hsapiens") {
    
    list_kMs_hs <- list_kMs_PPI
    list_colors_hs <- list_colors_PPI_uniq
    names(list_kMs_hs) <- names(list_colors_hs) <- names(list_kMs_PPI)
    genes_background_ensembl_hs <- read.csv(sprintf("%s%s_%s_%s_hgnc_to_ensembl_mapping_df.csv", dirTables, prefixData, prefixRun, dataOrganism))[["ensembl"]]
    
  }
  
  # Load gene mapping and annotation files
  
  # columns: ensembl_gene_id, entrezgene, hgnc_symbol. Use this for mapping entrezgene to ensembl
  mapping_hs_entrez2ensembl = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T, stringsAsFactors = F)
  
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
  
  
  # invisible(gc()); invisible(R.utils::gcDLLs())
  # cl <- makeCluster(min(n_cores, detectCores_plus(Gb_max = RAMGbMax, additional_Gb = 2)-1), type="FORK", outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "kM_magma_par.txt"))
  outfile = paste0(dirLog, prefixData, "_", prefixRun, "_", "kM_magma_par.txt")
  fun = function(cellType,kMs_hs,colors_hs) kM_magma(cellType =cellType, 
                                                     modulekM = kMs_hs,
                                                     gwas = gwas,
                                                     test_type = magma_test_type,
                                                     colors_hs = colors_hs,
                                                     genes_background_ensembl_hs = genes_background_ensembl_hs)
  args = list(cellType = names(list_kMs_hs),
              kMs_hs = list_kMs_hs,
              colors_hs = list_colors_hs)
  magma_results <- safeParallel(fun=fun, args=args, outfile = outfile, gwas=gwas, magma_test_type=magma_test_type, genes_background_ensembl_hs=genes_background_ensembl_hs)
  
  
  # Prepare tables for results
  list_magma.r <- list_magma.p <- list_magma.emp.p <- vector(mode="list", length = length(list_kMs_hs))
  
  # Extract the coefficient, analytical p-value and empirical p-value dataframes for each module, put them into lists
  
  for (i in 1:length(magma_results)) { 
    list_magma.r[[i]] <- magma_results[[i]][['corrCoef']]
    list_magma.p[[i]] <- magma_results[[i]][['p.val']]
    list_magma.emp.p[[i]] <- magma_results[[i]][['emp.p.val']]
  }
  
  # rowbind each list into a single table
  magma.r.all <- Reduce(f=rbind, x = list_magma.r)
  magma.p.all <- Reduce(f=rbind, x = list_magma.p)
  magma.emp.p.all <- Reduce(f=rbind, x = list_magma.emp.p)
  
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
  vec_MMtoHsmapping_prop.mapped <- NULL
}

rm(gwas)

##########################################################################
############## FILTER MODULES ON GENES WITH GWAS ENRICHMENT ##############
##########################################################################

if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) {
  
  # Check if the gwas_filter_traits strings provided by the user can be found in the magma_gwas_dir files
  
  if (sapply(gwas_filter_traits, function(x) any(grepl(x, colnames(magma.p.fdr.log), ignore.case=T)), simplify = T) %>% any) {
    
    # Identify columns to keep
    idx_col_keep <- sapply(gwas_filter_traits, function(x) grep(paste0(x,"|module|celltype"), colnames(magma.p.fdr.log), ignore.case=T), simplify = T) %>% Reduce(union,.)# %>% as.logical 
    
    # fdr-corrected log-transformed analytical p-values
    magma.p.fdr.log.sig <- magma.p.fdr.log[,idx_col_keep]
    
    # correlation coefficients
    magma.r.sig <- magma.r.all[,idx_col_keep]
    
    # Identify rows to keep
    idx_row_gwas <- apply(magma.p.fdr.log.sig[,-grep("module|celltype", colnames(magma.p.fdr.log.sig))], MARGIN = 1, max) > -log10(pvalThreshold)
    
    # rows enriched for MAGMA GWAS or GSEA
    idx_row_keep <- idx_row_gwas | (if (!is.null(list_genesets_path)) idx_row_GSEA[match(GSEA_df$module, magma.p.fdr.log$module)] else rep(TRUE, length(idx_row_gwas)) )
    
    if (sum(idx_row_keep)>0) {
      
      magma.p.fdr.log.sig <- magma.p.fdr.log.sig[idx_row_keep,]
      magma.r.sig <- magma.r.sig[idx_row_keep,]
      
    } else if (sum(idx_row_keep)==0) { # Warn if there are no significant enrichments
      
      gwas_filter_traits <- NULL # this will activate the logic in the next if statement
      warning("After filtering for GWAS enrichment no modules remain! Skipping GWAS filtering step")
      
    }
  } else {
    
    gwas_filter_traits <- NULL # this will activate the logic in the next if statement 
    warning("None of gwas_filter_traits strings found in magma_gwas_dir files! Skipping GWAS filtering step")
    
  }
}  

if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
  magma.p.fdr.log.sig <- NULL # magma.p.fdr.log
  magma.r.sig <- NULL #  magma.r.all
}

# count number of enriched module per celltype for summary stats
n_modules_gwas_enriched <- rep(NA, times = length(sNames_0))  
names(n_modules_gwas_enriched) <- sNames_0
n_modules_gwas_enriched[names(n_modules_gwas_enriched) %in% sNames_PPI] <- if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) sapply(sNames_PPI, function(x) sum(magma.p.fdr.log$celltype[idx_row_gwas]==x)) else rep(NA, times=length(sNames_PPI))

######################################################################
####### FILTER MODULES ON GWAS SIGNIFICANCE ALSO IN OTHER FILES ######
######################################################################

if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) { # NB: if no significant gwas correlations found in previous section, gwas_filter_traits <- NULL
  
  # get cell types with modules with gwas enrichment
  sNames_gwas <- sNames_PPI[match(names(table(magma.p.fdr.log.sig$celltype)[table(magma.p.fdr.log.sig$celltype)>0]),sNames_PPI)]
  
  list_module_gwas <- vector(mode="list", length=length(sNames_gwas))
  
  names(list_module_gwas) <- sNames_gwas
  
  # Get a list of vectors with modules per celltype
  for (name in sNames_gwas) {
    idx <- grep(pattern = name, x = magma.p.fdr.log.sig$celltype)
    list_module_gwas[[name]] <- as.character(sapply(idx, function(j) magma.p.fdr.log.sig$module[j], simplify = T))
  }
  
  # Filter out cell types with no significant modules
  list_colors_gwas <- list_colors_PPI_uniq[match(sNames_gwas,sNames_PPI)]
  
  # relabel non-GWAS modules 'grey'
  list_colors_gwas <- mapply(function(x,y) ifelse(x %in% y, yes = x, no = "grey"),
                             x = list_colors_gwas, 
                             y = list_module_gwas,
                             SIMPLIFY = F)
  
  # give gene names to color assignment vectors
  list_colors_gwas <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(y), dimension = NULL), 
                             x = list_colors_gwas,
                             y = list_colors_PPI_uniq[match(sNames_gwas,sNames_PPI), drop=F],
                             SIMPLIFY = F)
  
  # remove any cell clusters without gwas enrichment: expression matrices
  list_datExpr_gwas <- list_datExpr_PPI[match(sNames_gwas,sNames_PPI), drop=F] 
  
  # filter kMs
  list_kMs_gwas <- list_kMs_PPI[match(sNames_gwas,sNames_PPI), drop = F]
  list_kMs_gwas <- mapply(function(x,y) x[colnames(x) %in% y],
                          x = list_kMs_gwas,
                          y = list_module_gwas,
                          SIMPLIFY=F)
  
  # Filter left eigenvectors / kIMs
  list_u_gwas <- list_u_PPI[match(sNames_gwas,sNames_PPI), drop = F]
  
  list_u_gwas <- mapply(function(x,y) x[names(x) %in% y],
                        x = list_u_gwas,
                        y = list_module_gwas,
                        SIMPLIFY=F)
  
  # filter MEs
  if (fuzzyModMembership=="kME") {
    
    list_MEs_gwas <- list_MEs_PPI[match(sNames_gwas,sNames_PPI), drop = F]
    list_MEs_gwas <- mapply(function(x,y) x[colnames(x) %in% y], 
                            x = list_MEs_gwas,
                            y = list_module_gwas,
                            SIMPLIFY=F)
    
  } else {
    list_MEs_gwas <- NULL
    list_u_gwas <- NULL
  }
  
} else if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
  
  sNames_gwas <- sNames_PPI
  list_module_gwas <- lapply(list_kMs_PPI, function(x) colnames(x)) 
  names(list_module_gwas) <- sNames_gwas
  list_colors_gwas <- list_colors_PPI_uniq
  list_datExpr_gwas <- list_datExpr_PPI
  list_kMs_gwas <- list_kMs_PPI 
  list_MEs_gwas <- list_MEs_PPI
  list_u_gwas <- list_u_PPI
}
