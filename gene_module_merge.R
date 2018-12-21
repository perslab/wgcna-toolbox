# Cross-dataset gene network correlation and preservation analysis pipeline
# Usage: e.g.
# export PATH="/usr/local/R-3.5.1/bin/:$PATH"  # see https://stackoverflow.com/questions/26897335/how-can-i-load-a-specific-version-of-r-in-linux and https://gist.github.com/cheuerde/8fb9fd0dc8c0eca17c16#file-r_openblas-L64
# export R_MAX_NUM_DLLS=999
# time Rscript /projects/jonatan/tools/wgcna-src/wgcna-toolbox/gene_module_merge.R --dir_project_WGCNA /projects/jonatan/tmp-mousebrain/ --prefixes_WGCNA_run "c('mousebrain_Neurons_ClusterName_1b')" --dir_out /projects/jonatan/tmp-mousebrain/ --prefix_out mod_merge_1 --data_organism mmusculus  --corr_space cell --gene_names hgnc --fuzzyModMembership kIM --path_data /projects/jonatan/tmp-mousebrain/RObjects/mousebrain_hypoth.RDS --moduleMergeCutHeight 0.1 --n_cores 10

# args: see below
# value:
#   TODO

######################################################################
########################## FUNCTIONS #################################
######################################################################

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

source(file=paste0(current.dir,"rwgcna_functions.R"))

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("optparse", "Seurat", "dplyr", "Biobase", "Matrix", "parallel", "WGCNA"))

######################################################################
########################### OPTPARSE #################################
######################################################################
######################################################################
########################### OptParse #################################
######################################################################

option_list <- list(
  # inputs
  make_option("--dir_project_WGCNA", type="character",
              help = "Directory of WGCNA project, which must contain /tables and /RObjects subdirectories."), 
  make_option("--prefixes_WGCNA_run", type="character", default = NULL,
              help = "Search for matching files within the dir_project_WGCNA, if left as NULL will take all runs, [default %default]"),
  make_option("--dir_out", type="character", 
              help = "Directory for outputs, must contain /tables and /RObjects subdirectories"),
  make_option("--prefix_out", type="character", 
              help = "Prefix for output files"),
  make_option("--data_organism", type="character", default="mmusculus",
              help = "'hsapiens' or 'mmusculus', [default %default]"),
  make_option("--gene_names", type="character", default="hgnc|symbol|gene_name",
              help = "string or regex for grepping WGCNA outputs and gene mapping dataframe,s e.g. 'hgnc|symbol|gene_name' or 'ensembl', [default %default]"),
  make_option("--corr_space", type="character", default="cell",
              help = "Measure correlation in cell embedding ('cell') or gene loading ('gene') space, where the correlation takes the gene intersect? [default %default]"),
  make_option("--path_data", type="character", default = NULL,
              help = "If corr_space == cell, path to dataset with which to find module cell embeddings to measure module correlation in cell space."), 
  make_option("--fuzzyModMembership", type="character", default="kIM",
              help = "'kIM' or 'kME', [default %default]"), #TODO: implement for kME as well?
  make_option("--moduleMergeCutHeight", type="numeric", default=0.1,
              help = "Cut-off level for merging modules based on pairwise distance. [default %default]"),
  make_option("--n_cores", type="integer", default = 10L,
              help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]")
)

######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

if (FALSE) { 
  dir_project_WGCNA = c("/projects/jonatan/tmp-mousebrain/")
  prefixes_WGCNA_run = c("mousebrain_Neurons_ClusterName_1b")
  dir_out = "/projects/jonatan/tmp-mousebrain/"
  prefix_out = "mod_merge_1"
  data_organism = "mmusculus"
  corr_space = "cell"
  gene_names <- "hgnc"
  fuzzyModMembership="kIM"
  path_data = "/projects/jonatan/tmp-mousebrain/RObjects/mousebrain_hypoth.RDS"
  moduleMergeCutHeight = 0.1
  n_cores = 5
}

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

dir_project_WGCNA <- opt$dir_project_WGCNA 
prefixes_WGCNA_run <- opt$prefixes_WGCNA_run
if (!is.null(prefixes_WGCNA_run)) prefixes_WGCNA_run <- eval(parse(text=prefixes_WGCNA_run))
path_data <- opt$path_data
dir_out <- opt$dir_out
prefix_out <- opt$prefix_out
data_organism <- opt$data_organism
gene_names <- opt$gene_names
corr_space <- opt$corr_space
fuzzyModMembership <- opt$fuzzyModMembership
path_data <- opt$path_data
moduleMergeCutHeight <- opt$moduleMergeCutHeight
n_cores <- opt$n_cores

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F)
disableWGCNAThreads()

######################################################################
############################ VERIFY INPUT ############################
######################################################################

if (!file.exists(dir_project_WGCNA)) stop("dir_project_WGCNA not found")
if (!file.exists(path_data)) stop("path_data not found")
if (corr_space == "cell") if (is.null(path_data)) stop("if corr_space = cell,path_data is required")

######################################################################
############################ SET PATHS ###############################
######################################################################

# set up output directories 

dir_out_plots = paste0(dir_out,"plots/")
if (!file.exists(dir_out_plots)) dir.create(dir_out_plots) 

dir_out_tables = paste0(dir_out,"tables/")
if (!file.exists(dir_out_tables)) dir.create(dir_out_tables)

dir_out_RObjects = paste0(dir_out,"RObjects/")
if (!file.exists(dir_out_RObjects)) dir.create(dir_out_RObjects)

dir_out_log = paste0(dir_out, "log/")
if (!file.exists(dir_out_log)) dir.create(dir_out_log)

dir_scratch = "/scratch/tmp-wgcna/"

######################################################################
###################### CLUSTER AND MERGE MODULES #####################
######################################################################

if (corr_space=="cell") {
  
  # Set up directory
  dir_tables_WGCNA = paste0(dir_project_WGCNA,"tables/")
  if (!file.exists(dir_tables_WGCNA)) stop("dir_project_WGCNA must contain a /tables subdir")
  
  message("Loading gene module assignments")

  run_cellClusterModuleGenes <- list()
  for (prefix_run in prefixes_WGCNA_run) {
    cellClusterModuleGenes_path <- dir(path = dir_tables_WGCNA, pattern = paste0(prefix_run, "_cell_cluster_module_genes\\.csv"), full.names = T)
    run_cellClusterModuleGenes[[prefix_run]] <- load_obj(f = cellClusterModuleGenes_path)
    run_cellClusterModuleGenes[[prefix_run]][["module"]] <- paste0(prefix_run, "_",run_cellClusterModuleGenes[[prefix_run]][["module"]]) # ensure modules fromm different runs remain distinct
  }
  
  # merge into one long data frame
  runCellClusterModuleGenes <- Reduce(x = run_cellClusterModuleGenes,f = rbind)

  rm(run_cellClusterModuleGenes)
  
  # filter out very small modules. In future, this should have been done in WGCNA main
  mods_too_small <- names(table(runCellClusterModuleGenes[["module"]]))[table(runCellClusterModuleGenes[["module"]])<5]
  runCellClusterModuleGenes <- runCellClusterModuleGenes[!runCellClusterModuleGenes[["module"]] %in% mods_too_small,]
  
  pkM_col <- grep("^pk.?M$", colnames(runCellClusterModuleGenes), ignore.case=T, value=T)
  gene_col <- grep(gene_names, colnames(runCellClusterModuleGenes), value=T)
  ######################################################################
  ############################# LOAD DATA ##############################
  ######################################################################
  
  message("Loading expression matrix")
  
  data_obj <- NULL
  
  if (grepl(pattern = "\\.loom", path_data)) {
    suppressPackageStartupMessages(library("loomR"))
    data_obj <-  connect(filename = path_data, mode = "r+")
  } else {
    data_obj <- load_obj(path_data)
  } 
  
  if (!"loom" %in% class(data_obj) & class(data_obj)!="seurat") { # not a loom object, nor a seurat object
    data_obj <- CreateSeuratObject(raw.data=data_obj, project= prefix_out, min.cells = -Inf, min.genes = -Inf)
  } else if ("loom" %in% class(data_obj)) {
    data_obj <- Seurat::Convert(from=data_obj, to="seurat")
  } # if already a seurat object, do nothing

  ####### MAP EXPRESSION MATRIX TO ENSEMBL (AS USED IN WGNCA) #########
  
  mapping = if (data_organism=="hsapiens") {
    load_obj("/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz")
  } else if (data_organism=="mmusculus") {
    load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
  }  
  
  if (any(grepl("ENSG|ENSMUSG",rownames(data_obj@raw.data)))) {
    
    current_ensembl <- T
    
    tmp <- gene_map(df = data_obj@raw.data,
                    idx_gene_column = NULL,
                    mapping=mapping,
                    from="ensembl",
                    to="gene_name_optimal",
                    replace = T,
                    na.rm = T)
  } else {
    
    current_ensembl <- F
    tmp <- data_obj@raw.data
    
  }
  
  message("Counting percent ribo and percent mito")
  
  idx_mito.genes <- grepl(pattern = "^mt-", x = rownames(tmp), ignore.case=T)
  idx_ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(tmp), ignore.case=T)
  
  percent_mito <- Matrix::colSums(tmp[idx_mito.genes, ])/Matrix::colSums(tmp)
  percent_ribo <- Matrix::colSums(tmp[idx_ribo.genes, ])/Matrix::colSums(tmp)
  
  data_obj <- AddMetaData(object = data_obj, metadata = percent_mito, col.name = "percent_mito")
  data_obj <- AddMetaData(object = data_obj, metadata = percent_ribo, col.name = "percent_ribo")
  
  if (grepl("ensembl", gene_names) & !current_ensembl) {
    
    data_obj@raw.data <- gene_map(df = data_obj@raw.data,
                                  idx_gene_column = NULL,
                                  mapping=mapping,
                                  from="gene_name_optimal",
                                  to="ensembl",
                                  replace = T,
                                  na.rm = T)
    
  } else if (grepl("symbol|hgcn", gene_names) & current_ensembl) {
    
    data_obj@raw.data <- tmp
    
  }
  
  rm(tmp, mapping)
  
  ################# NORMALISE, SCALE AND REGRESS DATA ################
  
  data_obj <- NormalizeData(data_obj)
  vars.to.regress <- c("nUMI", "percent.mito", "percent.ribo")[c("nUMI", "percent.mito", "percent.ribo") %in% colnames(data_obj@meta.data)]
  data_obj <- ScaleData(object = data_obj, vars.to.regress = vars.to.regress, do.par = T, num.cores = n_cores)
  
  # take only scale.data
  # only take WGCNA genes
  # transpose to cell x gene

  datExpr <- data_obj@scale.data[rownames(data_obj@scale.data) %in% runCellClusterModuleGenes[[gene_col]],] %>% t

  rm(data_obj)

  ######################################################################
  ##################### CLUSTER AND MERGE ##############################
  ######################################################################
  
  ######### COMPUTE WGCNA DISTANCE TOPOLOGICAL OVERLAP MATRIX ##############
  
  powers = c(1:30)

  enableWGCNAThreads(nThreads=n_cores)
  sft = pickSoftThreshold(data=datExpr,
                          powerVector = powers,
                          blockSize = min(5000 , ncol(datExpr)), #try to prevent crashing
                          corFnc = "cor",
                          corOptions =  list(use = 'p'),
                          networkType = "signed hybrid",
                          verbose = 5)
  invisible(gc())
  disableWGCNAThreads()

  fitIndices <- as.data.frame(sft$fitIndices)
  fitIndices %>% dplyr::filter(median.k. <= quantile(median.k.,0.95, na.rm=T)) -> fitIndices_filter

  if (sum(fitIndices_filter$SFT.R.sq >= 0.93) > 1) {
    fitIndices_filter_2 <- fitIndices_filter[fitIndices_filter$SFT.R.sq>=0.93,]
    list_sft <- fitIndices_filter_2[which.min(fitIndices_filter_2$Power), c(1,2,6)]
    #softPower = min(fitIndices_filter$Power[fitIndices_filter$SFT.R.sq>=0.9])
  } else {
    #(fitIndices_filter %>% dplyr::arrange(desc(SFT.R.sq)) %>% dplyr::select(Power))[1,] -> list_sft
    (fitIndices_filter %>% dplyr::arrange(desc(SFT.R.sq)))[1,c(1,2,6)] -> list_sft
  }

  adjacency = adjacency(datExpr=datExpr,
                        type="signed hybrid",
                        power = list_sft$Power,
                        corFnc = "cor",
                        corOptions = list(use = 'p'))

  TOM = TOMsimilarity(adjMat=adjacency,
                            TOMType="signed",
                            TOMDenom="mean",
                            verbose=5,
                            indent = 0)

  rm(adjacency)
  
  dissTOM <- 1-as.dist(TOM) %>% as.matrix
  
  colnames(dissTOM) <- rownames(dissTOM) <- colnames(datExpr)
  
  rm(TOM) # Convert proximity to distance
  
  ############## MERGE CLOSE MODULES USING kIM #################
  
  corr_clust <- character(length = length(unique(runCellClusterModuleGenes[["module"]])))

  unique_modules <- unique(runCellClusterModuleGenes[["module"]])
  
  # check the intersect between module gene names and datExpr gene names
  # only try to merge modules where sufficient genes are found in the datExpr
  
  idx_datExpr_gene_intersect_ok <- sapply(unique_modules, function(module) {
    sum(runCellClusterModuleGenes[[gene_col]][runCellClusterModuleGenes[["module"]]==module] %in% colnames(datExpr))>10
  }, simplify=T)
  
  unique_modules_ok <- unique_modules[idx_datExpr_gene_intersect_ok]
  
  while (TRUE) {
      
    message("computing cell-module embedding matrix")

    unique_modules <- unique(runCellClusterModuleGenes[["module"]])
    unique_modules <- unique_modules[unique_modules %in% unique_modules_ok]
    
    cl <- try(makeCluster(spec=n_cores, type="FORK", outfile =paste0(dir_out_log,  prefix_out, "_parSapply_embed_mat.txt"), timeout=30))
    
    if (!"try-error" %in% class(cl)) {
      embed_mat <- as.matrix(parSapplyLB(cl, X=unique_modules, FUN=function(module) {
      datExpr_sub <- datExpr[,match(runCellClusterModuleGenes[[gene_col]][runCellClusterModuleGenes[["module"]]==module], colnames(datExpr))]
      pkIMs <- runCellClusterModuleGenes[[pkM_col]][runCellClusterModuleGenes[["module"]]==module]
      datExpr_sub %*% as.matrix(pkIMs)
      }))
    stopCluster(cl)
    
    } else {
      warning("makeCluster failed - probably due limited session memory. Computing serially instead")
      embed_mat <- as.matrix(sapply(X=unique_modules, FUN=function(module) {
        datExpr_sub <- datExpr[,match(runCellClusterModuleGenes[[gene_col]][runCellClusterModuleGenes[["module"]]==module], colnames(datExpr))]
        pkIMs <- runCellClusterModuleGenes[[pkM_col]][runCellClusterModuleGenes[["module"]]==module]
        datExpr_sub %*% as.matrix(pkIMs)
      }))
    }
    
    invisible(gc()); invisible(R.utils::gcDLLs())
    
    colnames(embed_mat) <- unique_modules
    rownames(embed_mat) <- rownames(datExpr)
    
    #cellModEmbed_mat <- cellModEmbed_mat[,-grep("^grey$", colnames(cellModEmbed_mat))]
    # Cluster modules using the Pearson correlation between module cell embeddings
    message("computing module-module correlation matrix")
    
    enableWGCNAThreads(nThreads=n_cores)
    mod_corr <- WGCNA::cor(x=embed_mat, method=c("pearson"), verbose=5)
    disableWGCNAThreads()
    rm(embed_mat)
    
    mod_corr[mod_corr<0] <- 0 #only keep positive correlations
    corr_dist <- as.dist(m = (1-mod_corr), diag = F, upper = F) # convert to (1-corr) distance matrix
    
    rm(mod_corr)
    message("clustering modules")
    corr_dendro <- hclust(d = corr_dist, method = "average")
      
    # use a simple cut to determine clusters
    clust_labels = cutreeStatic(dendro = corr_dendro,
                              cutHeight = moduleMergeCutHeight,
                              minSize=1)
    
    # corr_clust = cutreeHybrid(dendro = corr_dendro, 
    #                     #cutHeight = 1-moduleMergeCutHeight,
    #                     minClusterSize = 1, 
    #                     distM=as.matrix(corr_dist),
    #                     deepSplit=3, 
    #                     pamStage=T,
    #                     pamRespectsDendro=T,
    #                     verbose=5)#,
                        # maxPamDist = maxPamDist,
                        # useMedoids = useMedoids) 
    #clust_labels <- corr_clust$labels
    
    rm(corr_dist)
    
    names(clust_labels) =  corr_dendro$labels
    
    # At this point if corr_clust has as many unique modules as the original partition, the while loop will exit
    if(length(unique(corr_clust)) == length(unique_modules)) {  
      message(paste0("done merging modules"))
      break
    }
    # if not we merge modules
    merge_idx <- sapply(unique(clust_labels), function(x) sum(clust_labels==x)>1)
    
    n_merge <-  length(unique_modules)  - length(unique(clust_labels)) 
    message(paste0(n_merge, " modules to be merged with others"))
    
    # Th
    # merge modules
    new_mods <- c()
    for (clust in unique(clust_labels)[merge_idx]) { # loop over numeric cluster assignments 
      new_module = names(clust_labels)[clust_labels==clust][1] # use the colors name of the first module in the cluster
      new_mods <- c(new_mods, new_module)
      for (module in names(clust_labels)[clust_labels==clust][2:length(names(clust_labels)[clust_labels==clust])]) { # loop over all other gene modules in the cluster
        runCellClusterModuleGenes[["module"]][runCellClusterModuleGenes[["module"]]==module] <- new_module # give them all the new colour
        # Remove duplicates (going by gene). Doesn't matter which duplicates we remove as we'll recompute kIM anyway.
        #idx_dupl <- duplicated(runCellClusterModuleGenes[["ensembl"]][runCellClusterModuleGenes[["module"]]==new_module]) 
        #runCellClusterModuleGenes <- runCellClusterModuleGenes[!runCellClusterModuleGenes[["module"]]==new_module][!idx_dupl]
      }
    }
    ########## eliminate duplicate rows ###########
    
    idx_dupl <- duplicated(runCellClusterModuleGenes[[gene_col]]) & duplicated(runCellClusterModuleGenes[["module"]])
    # NB: either ensembl or symbol will do here.
    
    # NB: it appears that nearly all of the duplication of genes across different modules dissappears after merging module. only ~3500 genes are 
    # duplicated across different modules.
    # > sum(duplicated(runCellClusterModuleGenes[["ensembl"]]))
    # [1] 586344
    # > nrow(runCellClusterModuleGenes)
    # [1] 604875
    # > sum(duplicated(runCellClusterModuleGenes[["ensembl"]]) & duplicated(runCellClusterModuleGenes[["module"]]))
    # [1] 582776
    runCellClusterModuleGenes <- runCellClusterModuleGenes[!idx_dupl,]
    # > dim(runCellClusterModuleGenes)
    # [1] 22099     5
    # we are left with roughly the number of genes present in the cell ranger data. 
    
    ########### compute merged module kIMs ###########
    message(paste0("Computing merged module kIMs"))
    # Compute new kIMs
    colors <- runCellClusterModuleGenes[["module"]][runCellClusterModuleGenes[["module"]] %in% new_mods]
    names(colors) <- runCellClusterModuleGenes[[gene_col]][runCellClusterModuleGenes[["module"]] %in% new_mods]

    kIMs <- kIM_eachMod_norm(dissTOM = dissTOM, 
                             colors = colors,
                             verbose = 3,
                             excludeGrey = F,
                             do.par=F,
                             n_cores=n_cores)
    
    for (module in new_mods){
      idx_match <- match(runCellClusterModuleGenes[[gene_col]][runCellClusterModuleGenes[["module"]] == module], rownames(kIMs))
      if(length(idx_match[!is.na(idx_match)])>0) runCellClusterModuleGenes[[pkM_col]][runCellClusterModuleGenes[["module"]] == module][!is.na(idx_match)] <- kIMs[idx_match[!is.na(idx_match)], colnames(kIMs)==module] 
    }
    
  }

} #else if (corr_space=="gene") {
  
  # distance matrix:
  # for every pair, find intersect
  # find rank correlation between intersecting genes
  # do until no more modules to merge:
  #   hierarchical clustering and merging
  
  ######################################################################
  ######################### LOAD GENE LOADINGS VECTORS #################
  ######################################################################
  
  # message("Loading module gene loading vectors")
  # 
  # dir_RObjects_WGCNA = paste0(dir_project_WGCNA,"RObjects/")
  # if (!file.exists(dir_RObjects_WGCNA)) stop("dir_project_WGCNA must contain a /RObjects subdir")
  # 
  # run_cellType_module_u <- list()
# TODO
#   run_cellType_module_u_paths <- dir(path = dir_RObjects_WGCNA[i], pattern = paste0(prefixes_WGCNA_run, "_list_list_module_u.RDS", collapse = "|"), full.names = T)
#     for (cellType_module_u_path in run_cellType_module_u_paths) {
#       cellType_module_u_name <- gsub(".*/|_list_list_module_u.RDS", "", cellType_module_u_path)
#       WGCNAproj_cellType_module_u[[WGCNA_projects[i]]][[cellType_module_u_name]] <- load_obj(f = cellType_module_u_path)
#     }
#     names(run_cellType_module_u[[i]]) <- gsub(".*/|_list_list_module_u.RDS", "", run_cellType_module_u_paths)
#   }
  
#   names(WGCNAproj_run_cellType_module_u) = WGCNA_projects
#   
#   # Unlist 
#   cellType_module_u <- unlist(x = run_cellType_module_u, recursive = F, use.names = T)
#   module_u <- unlist(x = cellType_module_u, recursive = F, use.names = T)
#   names(module_u) <- paste0(prefix_out, ".", names(module_u))
#   
#   
# }
invisible(write.csv(runCellClusterModuleGenes, file = sprintf("%s%s_%s_cell_cluster_module_genes_merged.csv", tables_dir, prefix_data, prefix_out), row.names = F, quote = F))
