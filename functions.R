# Title: sc functions
# Author: Jonatan Thompson, Pers Lab
# Date: April 2018
# Version 1.2

######################################################################
#
# # Change log
#
# ## v1.2 unreleased
# ### Added
# * kME_reassign
# * extract_and_name_colors
# ### 
#
# ## Changed
# * PPI_for_par_inner
# * PPI_for_par_outer
# * parCutreeHybrid, parMergeCloseModules and parPlotLabel now take named parameters and do not rely on global variables
######################################################################



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


######################################################################


extract_and_name_colors <- function(merged, datExpr_filter) {
  colors <- merged[[1]]
  names(colors) = colnames(datExpr_filter)
  return(colors)
}


######################################################################



kME_reassign_fnc = function(MEs, kMEs, pkMEs, kME_reassign_threshold, colors, datExpr_filter, corFnc) {
  
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

    reassigned = color[rows_to_reassign]
    colors[rows_to_reassign] <- colors[col_most_negkME][rows_to_reassign]
    
    # Recompute Module Eigengenes, kME, pkME
    MEs = moduleEigengenes(expr = as.matrix(datExpr_filter),
                           colors,
                           excludeGrey = F)
    
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

######################################################################

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

######################################################################

load_obj <- function(f) {
  # Utility function for loading an object inside a new environment and returning it so it can
  # stored in a variable
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

######################################################################

parCutreeHybrid <- function(comb, geneTree, dissTOM) {
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

######################################################################

parMergeCloseModules <- function(cutree,comb, datExpr_filter) {
  # Utility function for more easily parallelising mergeCloseModules
  if (any(as.logical(cutree$labels))) { # Are there any modules at all? Avoids an error
    merged = mergeCloseModules(exprData=as.matrix(datExpr_filter), 
                               colors = cutree$labels, 
                               cutHeight=comb[[4]])
    colors = labels2colors(merged$colors)
    #MEs = merged$newMEs 
    MEs = moduleEigengenes(expr = as.matrix(datExpr_filter),
                           colors,
                           excludeGrey = F)
  }
  else {# only grey modules
    colors = labels2colors(cutree$labels)
    MEs = NULL
  }
  return(list("cols" = colors, "MEGs"= MEs))
}

######################################################################
# EDIT_20180420_6
#parLabel <- function(comb) {
parPlotLabel <- function(comb) {
###
  # Utility function for more easily parallelising making labels
  label = paste0("MMS=", comb[[1]], ",DS=", comb[[2]],",PAM=",comb[[3]], ",CUT=",comb[[4]],sep="") 
}

######################################################################

SeuratFactorToIndicator <- function(obj, colnames, features) {
  # @Usage:   Converts factor/character metadata to new dummy variable columns
  # @args:    obj: Seurat object
  #           colnames: list of metadata column names - character
  #           features: list of vectors of character names within each metadata column to extract as new dummy variable columns
  #                       length(colnames) := length(features)
  # @return:  Seurat subset object with additional new indicator variable metadata columns (old are preserved)
  # @depends: Seurat
  # @author:  Jonatan Thompson jjt3f2188@gmail.com
  # @date:    180223
  # @TODO: 
  list.list.idx = list()
  for (j in seq(length(colnames))){
    col <- grep(colnames[[j]], names(obj@meta.data)) # Find the meta data column index
    list.idx <- lapply(features[[j]], function(x) obj@meta.data[col]==x) # find, for each string, a logical vector of occurences
    list.list.idx[[j]] <- list.idx
  }
  flatlist.idx <- unlist(list.list.idx, recursive=F, use.names = T) # Flatten to list of logical vectors  
  names(flatlist.idx) <- unlist(features, recursive=F) # Assign a feature name to each vector
  new_cols <- data.frame(lapply(flatlist.idx, function(x) as.numeric(x)), row.names=row.names(obj@meta.data)) # Make dataframe
  return(AddMetaData(object=obj, new_cols)) # Add dataframe to Seurat object and return
}

######################################################################

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


# subSetF <- function(obj, meta, ident) {
# # @Status DEPRECATED. ONLY USEFUL IF MAKING SUBSETS ACROSS METADATA COLUMNS, OTHERWISE USE SubsetMeta
# # @Usage: Subset a seurat object using the 'ident.use' argument. 
# #         Allows for subsetting in parallel by putting all the steps into one function.
# # @args: 
# #       obj: Seurat object from which to take a subset
# #       meta: vector of metadata
# #       ident: character string(s) within metadata on which to subset
# # @return: 
# #       result: Seurat object with subset columns
# # @author: Jonatan Thompson jjt3f2188@gmail.com
# # @date: 180222
#   vec=character(length=length(obj@cell.names)) # make empty vector
#   idx <- ident %in% meta # Identify cells to subset based on metadata
#   vec[idx] <- as.character(meta[idx])
#   vec <- data.frame(vec, stringsAsFactors = T) # Set up the dataframe
#   row.names(vec) <- obj@cell.names
#   
#   # Store the current identity, add the dataframe as metadata, make it the cell identity
#   obj <- StashIdent(obj, save.name="oldIdent")
#   obj <- AddMetaData(obj, vec, col.name="vec")
#   obj <- SetAllIdent(obj, id="vec")
#   
#   # Subset the data based on new identity; scale and center the data, restore the old identity
#   obj_sub = SubsetData(obj, ident.use=ident, do.scale=T, do.center=T)
#   obj <- SetAllIdent(obj, id="oldIdent") # Revert cell identities in full dataset object
#   obj_sub <- SetAllIdent(obj_sub, id="oldIdent") # Revert cell identities in subset
#   return(obj_sub)
# }


######################################################################

PPI_for_par_outer = function(pkMEs, colors, STRINGdb_species, PPI_pval_threshold, project_dir, data_prefix, subsetName, flag_date) {
  # Rather than for parallelising over modules within a set of modules, parallelise over several sets of colors (produced by comparing parameters)
  # It calls PPI_for_par as a subroutine
  # Args:
  #   STRINGdb_species
  #   pkMEs
  #   colors       
  # Returns: 
  #   colors_PPI: colors where modules that did not satisfy the PPI threshold are set to grey
  
  module_PPI <- NULL
  unique_colors <- NULL
  unique_colors_PPI <- NULL
  colors_PPI <- NULL
  idx_pkME_ok <- NULL
  MEs_PPI <- NULL
  kMEs_PPI <- NULL
  pkMEs_PPI <- NULL 
  
  #idx_pkME_ok <- rep("TRUE", length(pkMEs))
  #idx_pkME_ok <- abs(pkMEs)>median(abs(pkMEs))#order(pkMEs, decreasing=T)[0:(length(pkMEs)/2)]
  idx_pkME_ok <- abs(pkMEs) > 0.1
  
  unique_colors <- unique(colors[idx_pkME_ok])
  
  string_db <- STRINGdb$new(version="10", species = STRINGdb_species, score_threshold=0, input_directory="") # Default: 10090 (Mus musculus)
  module_PPI <- sapply(unique_colors, function(x) PPI_for_par_inner(color=x, subsetName = subsetName, unique_colors = unique_colors, idx_pkME_ok = idx_pkME_ok, colors = colors, string_db = string_db), simplify=T) %>% t()
  module_PPI <- as.data.frame(module_PPI, row.names = unique_colors)
  
  # FILTER MODULES ON PPI ENRICHMENT  

  if (!is.null(module_PPI)) {
    unique_colors_PPI = unique_colors[module_PPI$'p-value' < PPI_pval_threshold]
    genes_PPI_idx <- colors %in% unique_colors_PPI
    colors_PPI <- colors
    colors_PPI[!genes_PPI_idx] <- "grey"
  } else {
    write.csv("WARNING" , file=sprintf("%s%s_%s_WARNING_sum_genes_PPI_0_%s.csv", project_dir, data_prefix, subsetName, flag_date))
    colors_PPI <- rep("grey", length(colors))
  }  
  return(colors_PPI)
}

######################################################################

PPI_for_par_inner <- function(color, subsetName, unique_colors, idx_pkME_ok, colors, string_db) {
  # Parallelise over modules within a set of modules
  # utility function to parallelise the PPI STRINGdb call
  # args: color should be one of the unique_colors
  
  i = which(color == unique_colors) 
  ppi <- data.frame(gene = names(colors[idx_pkME_ok][colors[idx_pkME_ok]==color])) # extract the genes with the corresponding color to dataframe
  module_PPI <- list('p-value'=1, 'expected interactions'=NA)
  
  tryCatch({
    
    example1_mapped <- string_db$map(ppi, 'gene', removeUnmappedRows = TRUE ) # Check the dataframe genes' PPI. May produce error: we couldn't map to STRING 100% of your identifiers
    # (ctd.) .. Error in if (hitList[i] %in% V(ppi_network)$name) { : argument is of length zero
    hits <- example1_mapped$STRING_id
    module_PPI['p-value'] = p.adjust(string_db$get_ppi_enrichment(hits)$enrichment,method = 'fdr',n = length(names(colors[idx_pkME_ok])))
    module_PPI['expected interactions'] = string_db$get_ppi_enrichment(hits)$lambda
    
  }, error = function(c) {
    write.csv("ERROR", file=sprintf("%s%s_%s_%s_STRINGdb_ERROR_%s.csv", project_dir, data_prefix, subsetName, color, flag_date), row.names = F)
    module_PPI <- list('p-value'=1, 'expected interactions' = NA) # probably unnecessary since defined above but just in case it has been changed in place
  }) 
  
  return(module_PPI)
  
}



############################################################################################################################################################


parMagma = function(subsetName, project_dir, plots_dir, tables_dir, log_dir, magma_gwas_dir, data_prefix, file_suffix, flag_date) {
  
  # Usage:
  # Args:
  #   subsetName: subsetname to loop over
  #   magma_gwas_dirs: take a vector of directories and loops over them, producing separate output for each
  # Returns: does not return an object but saves plots and .csv to the respective directories

  options(stringsAsFactors = F)
  file_sep = ','
  
  # Set paths to mapping files to variables (these are independent of arguments)
  mapping_hs_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz" # columns: ensembl_gene_id, entrezgene, hgcn_symbol
  mapping_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
  mapping_mm_synonyms_filepath = "/data/genetic-mapping/ncbi/Mus_musculus.gene_info_symbol2ensembl.gz"
  mapping_hs_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
  
  # Load WGCNA results
  modulekME = read.csv(file=sprintf("%s%s_%s_%s_%s.csv", tables_dir, data_prefix, subsetName, file_suffix, flag_date), row.names=1)#, check.names = FALSE, sep = file_sep)

  gwas_sub_dirs = list.dirs(path = magma_gwas_dir, full.names = FALSE, recursive = FALSE)
  
  for (sub_dir in gwas_sub_dirs) {
    
    # Initialise log file
    log_not_mapped_filepath = paste0(log_dir,flag_date,"_magma_wgcna_not_mapped_",data_prefix,"_",subsetName, sub_dir,".tab")
    
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
    mapping_hs_entrez2ensembl = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T)
    for(i in 1:length(gwas)) {
      idx = match(gwas[[i]]$GENE, mapping_hs_entrez2ensembl$entrezgene)
      mapping = data.frame(entrez=gwas[[i]]$GENE, ensembl=mapping_hs_entrez2ensembl$ensembl_gene_id[idx])
      gwas[[i]]$gene_name = mapping$ensembl
    }
    
    # Remapping the WGCNA data to human ensembl IDs (using synonyms)
    # Step 1: direct mapping
    mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T)
    mapping = data.frame(symbol=row.names(modulekME), ensembl.mouse=mapping_direct$ensembl_gene_id[ match(row.names(modulekME), mapping_direct$gene_name_optimal) ])
    
    # Step 2: map remaing using synonyms
    mapping_synonyms = read.csv(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T)
    mapping$ensembl.mouse[ which(is.na(mapping$ensembl.mouse)) ] = mapping_synonyms$ensembl[ match( mapping$symbol[which(is.na(mapping$ensembl.mouse)) ] ,mapping_synonyms$symbol) ]
    
    # Step 3: orthology mapping
    #mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    #mapping_mm_orthologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), mart=mart)
    mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T)
    mapping$ensembl.human = mapping_orthology$ensembl_gene_id[ match(mapping$ensembl.mouse,mapping_orthology$mmusculus_homolog_ensembl_gene) ]
    #mapping$ensembl.human[mapping$ensembl.human == ""] = NA
    df_not_mapped = mapping[is.na(mapping$ensembl.human),]
    write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
    modulekME$symbol = mapping$symbol
    modulekME$ensembl = mapping$ensembl.human
    modulekME = na.omit(modulekME)
    tmp = within(modulekME, rm("symbol","ensembl"))
    
    # Average duplicated gene IDs
    modulekME_ens <-aggregate(tmp, by=list(modulekME$ensembl),FUN=mean, na.rm=TRUE)
    rownames(modulekME_ens) = modulekME_ens$Group.1
    modulekME_ens = within(modulekME_ens, rm("Group.1"))
    
    # Calculate spearman's correlation between gene module membership and GWAS gene significance
    colors = colnames(modulekME_ens)
    table.kme.cor.p = table.kme.cor.r<- matrix(NA,nrow=length(unique(colors)),ncol=length(gwas))
    rownames(table.kme.cor.r) = rownames(table.kme.cor.p) = unique(colors)
    colnames(table.kme.cor.r) = colnames(table.kme.cor.p) = names(gwas) 
    
    for(m in unique(colors)) {
      for(i in 1:length(gwas)) {
        #col = paste("kME", m, sep="")
        col = m
        genes = intersect(rownames(modulekME_ens),gwas[[i]]$gene_name)
        x = -log10(gwas[[i]]$P[match(genes, gwas[[i]]$gene_name)])
        y = modulekME_ens[match(genes,rownames(modulekME_ens)), col]
        cor = cor.test(x,y,method="spearman", exact=F)
        table.kme.cor.r[m,i] = cor$estimate
        table.kme.cor.p[m,i] = cor$p.value
      }
    }
    
    table.kme.cor.p.fdr = p.adjust(table.kme.cor.p, method="fdr")
    dim(table.kme.cor.p.fdr) = dim(table.kme.cor.p);  dimnames(table.kme.cor.p.fdr) = dimnames(table.kme.cor.p)
    
    d = -log10(table.kme.cor.p.fdr) * sign(table.kme.cor.r) 
    #pdf("SampleGraph.pdf",width=7,height=5)
    #sizeGrWindow(9,7)
    #par(mfrow = c(2,2))
    #par(mar = c(4, 5, 4, 6));
    #labeledHeatmap(d,textMatrix = signif(table.kme.cor.r,1), xLabels = colnames(d), yLabels = rownames(d),invertColors = T, colors = blueWhiteRed(1000), main="GWAS - kME correlation", cex.text = 0.6)
    #dev.off()
    
    #dat = as.data.frame(table.kme.cor.p.fdr*sign(table.kme.cor.r))[c("tan","blue","yellow","purple","turquoise","green","greenyellow","salmon"),]
    dat = as.data.frame(table.kme.cor.p.fdr*sign(table.kme.cor.r))
    dat[dat<0]=1 #Only look for positive enrichment
    dat = -log10(dat)
    dat$module = gsub("kME","",rownames(dat))
    #dat$module = gsub("ME","",rownames(dat))
    dat2 = melt(dat)
    dat2$variable=as.character(dat2$variable)
    
    # #p=ggplot(melt(dat),aes(x=variable,y=value,fill=colors)) + 
    
    # p=ggplot(melt(dat),aes(x=variable,y=value,fill="blue")) + 
    #   geom_bar(stat="identity",position=position_dodge(),color="black") +
    #   scale_fill_manual(values=sort(unique(dat2$module)))+ theme_classic() +
    #   geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
    #   theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
    
    p=ggplot(melt(dat),aes(x=variable,y=value,fill=module)) + 
      geom_bar(stat="identity",position=position_dodge(),color="black") +
      scale_fill_manual(values=sort(unique(dat2$module)))+ theme_classic() +
      geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
      theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
    
    #p
    
    ggsave(p, filename = paste0(plots_dir, data_prefix,"_",subsetName, "_", file_suffix, "_magma_GWAS_", sub_dir, "_", flag_date, ".pdf") ,width=45,height=12)
    
  }
    
}




