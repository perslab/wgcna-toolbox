# module specificity filter

# Usage: 
# export R_MAX_NUM_DLLS=999
# time RScript..

# TODO: delete all non-specificity measures
######################################################################
############################# SET OPTIONS ############################
######################################################################

options(stringsAsFactors = F, use="pairwise.complete.obs")

######################################################################
######################### UTILITY FUNCTIONS ##########################
######################################################################

source(file="/projects/jonatan/tools/functions-src/utility_functions.R")
source(file="/projects/jonatan/tools/functions-src/functions_sc.R")

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("optparse", "Matrix", "ggplot2", "dplyr", "parallel", "pSI", "scales", "liger", "ComplexHeatmap"))#, "reshape", "reshape2"))#, "loomR", "doubletFinder")

######################################################################
########################### OptParse #################################
######################################################################

option_list <- list(
  make_option("--path_df_NWA", type="character",
              help = "path to dataframe from a gene network analysis run, in long format (i.e. one row per gene per celltype), containing 'cell_cluster', 'module', one or two gene name columns, and a column of numeric scores. I.e. one row per gene, e.g. rwgcna cell_cluster_module_genes.csv files"),  
  make_option("--colCellCluster", type="character", default = "cell_cluster", 
              help = "The column for cell cluster annotation in the cell embedding matrix and cell_cluster_module_gene.csv files, [default %default]"), 
  make_option("--colMod", type="character", default = "module", 
              help = "The column for module cell_cluster_module_gene.csv files, [default %default]"), 
  make_option("--colGeneNames", type="character", default = "module", 
              help = "The column for gene names in cell_cluster_module_gene.csv files, [default %default]"), 
  make_option("--colGeneWeights", type="character", default = "module", 
              help = "The column for gene weights in cell_cluster_module_gene.csv files, [default %default]"), 
  make_option("--path_df_modFilter", type="character", default=NULL, 
              help = "Filter modules by values in this df [default %default]"), 
  make_option("--colFilter", type="character", default = "filterOK", 
              help = "Filter modules by this column of df_modFilter [default %default]"), 
  make_option("--path_df_kMs", type="character", default = NULL,
              help = "Path to gene module weights as dataframe. The first column is assumed to have the genes"),
  make_option("--path_mat_celltypeModEmbed", type="character", default = NULL,
              help = "Path to precomputed module celltype average expression or preservation matrix"),
  make_option("--scaleModExpr", type="logical", default = T, 
              help = "Convert module expression to [0:1]? Useful if cell module embedding is in Z-score format since specificity index requires non-negative values"), 
  make_option("--path_df_celltypeGenes", type="character", 
              help = "Path to celltype specific genes, e.g. SEM or differentially expressed genes, in a table with a 'gene' column and celltype columns with gene scores"),
  make_option("--SEMAbsCutoff", type="integer", default = 100L, 
              help = "Take top x SEM genes per celltype"), 
  make_option("--SEMPropCutoff", type="integer", default = 0.99, 
              help = "Take top x proportion of SEM genes"), 
  make_option("--dirOut", type="character",
              help = "Project directory to which to write files. Should include subdirectories /tables, /RObjects, /plots, /log"),
  make_option("--RAMGbMax", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)


######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

path_df_NWA <- opt$path_df_NWA
colCellCluster <- opt$colCellCluster
colMod <- opt$colMod
colGeneWeights <- opt$colGeneWeights
path_df_modFilter <- opt$path_df_modFilter
colFilter <- opt$colFilter
path_df_kMs <- opt$path_df_kMs
colGeneNames <- opt$colGeneNames
path_mat_celltypeModEmbed <- opt$path_mat_celltypeModEmbed
scaleModExpr <- opt$scaleModExpr
path_df_celltypeGenes <- opt$path_df_celltypeGenes
SEMAbsCutoff <- opt$SEMAbsCutoff
SEMPropCutoff <- opt$SEMPropCutoff
dirOut <- opt$dirOut
RAMGbMax <- opt$RAMGbMax

######################################################################
########################## DERIVED CONSTANTS #########################
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
############################ LOAD DATA ###############################
######################################################################

mat_celltypeModEmbed <- if (!is.null(path_mat_celltypeModEmbed))  load_obj(path_mat_celltypeModEmbed) else NULL

df_celltypeGenes <- load_obj(path_df_celltypeGenes)

df_NWA <- load_obj(path_df_NWA) #WGCNA outputs

df_kMs <- read.csv(path_df_kMs, quote="", header=T)

df_modFilter <- if (!is.null(path_df_modFilter)) load_obj(path_df_modFilter) else NULL
#cellModEmbedMat <- if (!is.null(pathCellModEmbedMat))  load_obj(pathCellModEmbedMat) else NULL 

#if (!xor(is.null(mat_CelltypeModEmbed), is.null(cellModEmbedMat))) stop("provide one and only one of pathmat_CelltypeModEmbed or pathCellModEmbedMat")

#kMs <- if (!is.null(pathkMs)) read.csv(pathkMs) else NULL

######################################################################
################# TODO: CONVERT TO HOMO SAPIENS ######################
######################################################################

df_NWA <- gene_map(df=df_NWA,
                   idx_gene_column = 5,
                   mapping= load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"),
                   from="mmusculus_homolog_ensembl_gene", 
                   to="ensembl_gene_id",
                   replace = F,
                   na.rm = F)

df_kMs <- gene_map(df=df_kMs,
                   idx_gene_column = 1,
                   mapping= load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"),
                   from="mmusculus_homolog_ensembl_gene", 
                   to="ensembl_gene_id",
                   replace = T,
                   na.rm = F)

# Get rid of NA genes in kMs
df_kMs <- df_kMs[!is.na(df_kMs[[1]]),]

######################################################################
################# GET MODULE WEIGHTS WITH GENE NAMES #################
######################################################################

modules <- df_NWA[[colMod]] %>% unique %>% sort 

fun <- function(module) {
  vec_geneWeights<-df_NWA[[colGeneWeights]] %>% '['(df_NWA[[colMod]]==module)
  geneNames <- df_NWA[[colGeneNames]] %>% '['(df_NWA[[colMod]]==module)
  vec_geneWeights <- vec_geneWeights[!is.na(geneNames)]  
  names(vec_geneWeights) <- geneNames[!is.na(geneNames)]  
  return(vec_geneWeights)
}
args = list("X"=modules)
list_vec_geneWeights <- safeParallel(fun=fun,args=args)
names(list_vec_geneWeights) <- modules

######################################################################
###################### DISTINGUISH ANNOTATION AND VALUES #############
######################################################################

# if (!is.null(cellModEmbedMat)){
#   logical_cellMod_numeric <- sapply(X=cellModEmbedMat, FUN = class) == "numeric"
#   if (all(logical_cellMod_numeric)) stop("no annotations in cellModEmbedMat")
#   
#   cellModEmbedMat_annot <- cellModEmbedMat[,!logical_cellMod_numeric] 
#   cellModEmbedMat_num <-  cellModEmbedMat[,logical_cellMod_numeric] 
#   rownames(cellModEmbedMat_annot) <- rownames(cellModEmbedMat_num) <- cellModEmbedMat_annot[[colCell]]
#   

######################################################################
######################### ALL SEM GENES  #############################
######################################################################

vec_genesAll <- df_celltypeGenes[,1]


if (is.null(SEMAbsCutoff)) {
  SEMAbsCutoff  <- as.integer(nrow(df_celltypeGenes)*(1-SEMPropCutoff))
}

######################################################################
##### Specificity Index test for mean celltype module expression #####
######################################################################

if (!is.null(mat_CelltypeModEmbed)) {
### scale the module expression to [0:1] score
  if (scaleModExpr ) {
    fun = function(module_vec) {
      rescale(x=module_vec, from = range(module_vec, na.rm = TRUE,
                                         finite = TRUE), to = c(0, 1))
      # module_vec %>% '-'(mean(.)) %>% '/'(stats::sd(.)) -> out
    }
    mat_CelltypeModEmbed <- apply("X"=mat_CelltypeModEmbed, MARGIN=2, FUN =fun)
  }
    
  
  #idx_neg <- apply(X=mat_CelltypeModEmbed_num, MARGIN=2, FUN = function(col) col<0)
  #mat_CelltypeModEmbed_num_pSI_tmp <- mat_CelltypeModEmbed_num
  #mat_CelltypeModEmbed_num_pSI_tmp[idx_neg] <- 0
  specIndex <- specificity.index(pSI.in <- as.data.frame(mat_CelltypeModEmbed_num),
                                  bts=50,
                                  p_max = 1,
                                  e_min = 0.0001) # what is difference between SI value and pSI value?
  
  specIndex %>% melt %>% "[["(2) %>% p.adjust(., method ="fdr") -> tmp # correct for n celltypes * n modules
  specIndex_adj <- matrix(tmp, nrow=nrow(specIndex)) # matrix builds the matrix column by column
  dimnames(specIndex_adj) <- dimnames(specIndex)
  
}

######################################################################
####################### Fisher's exact test ##########################
######################################################################

# Fisher’s exact test on SEM genes from PT. Define cut-off for top cell-type SEM genes

# iterate over celltypes
# iterate over modules

phyperStatFnc <- function(vec_logicalCelltypeGenes, 
                          vec_genesAll, 
                          vec_geneWeights) {
  vec_celltypeGenes <- vec_genesAll[vec_logicalCelltypeGenes]
  vec_genesModule = names(vec_geneWeights)
  m = length(vec_celltypeGenes) # number of white balls in the urn
  n = length(vec_genesAll)-m  # number of black balls in the urn
  k = length(vec_genesModule) # number of balls drawn from the urn
  q = length(intersect(vec_genesModule, vec_celltypeGenes)) # number of white balls drawn without replacement from the urn
  phyper(q=q, m=m, n=n, k=k, lower.tail=F)
}
  
R = 200

fun = function(vec_celltypeScores) {
  vec_celltypeScores %>% order(.,decreasing = T) %>% '['(1:SEMAbsCutoff) -> vec_idxCelltypeGenes 
  #vec_CelltypeGenes <- vec_genesAll[vec_idxCelltypeGenes]
  vec_logicalCelltypeGenes <- logical(length(vec_genesAll))
  vec_logicalCelltypeGenes[vec_idxCelltypeGenes] <- T
  sapply(list_vec_geneWeights, function(vec_geneWeights) {
    
    # vec_genesModule = names(vec_geneWeights)
    # m = length(vec_celltypeGenes) # number of white balls in the urn
    # n = length(vec_genesAll)-m  # number of black balls in the urn
    # k = length(vec_genesModule) # number of balls drawn from the urn
    # q = length(intersect(vec_genesModule, vec_celltypeGenes)) # number of white balls drawn without replacement from the urn
     
    boot_out <- boot(data = vec_genesAll,
                     statistic=phyperStatFnc,
                     R=R,
                     sim = "permutation",
                     vec_logicalCelltypeGenes=vec_logicalCelltypeGenes,
                     vec_geneWeights = vec_geneWeights,
                     parallel = "no")
    
    # compute the empirical probability of the p-value - should correspond to the p-value, ideally..
    ecdf(boot_out$t)(boot_out$t0)
    
    #phyper(q=q, m=m, n=n, k=k, lower.tail=F) # P[X > x]
  }, simplify=T)
}

args = list("X"= df_celltypeGenes[,-1])
mat_Fisher <- safeParallel(fun=fun, args=args, 
                          list_vec_geneWeights=list_vec_geneWeights,
                          vec_genesAll=vec_genesAll,
                          simplify=T)

mat_Fisher %>% as.data.frame %>% 
  melt %>% "[["(2) %>% p.adjust(., method ="fdr") -> tmp # correct for n celltypes * n modules
mat_FisherAdj <- matrix(tmp, nrow=nrow(mat_Fisher)) # matrix builds the matrix column by column
dimnames(mat_FisherAdj) <- dimnames(mat_Fisher)

######################################################################
############################ SEM-based test ##########################
######################################################################

# on SEM genes from PT. No need for cut-off.
# Groups = [gene module members] and [gene non-module members]
# variable = SEM values for a given cell-types
# Statistic: “GSEA-like” test

############################ t.test ##################################
# Test: do the genes in the module have higher cell-type SEM values than genes not in the module?

fun <- function(vec_celltypeScores) {
  sapply(list_vec_geneWeights, function(vec_geneWeights) {
    vec_genesModule = names(vec_geneWeights)
    vec_idxModGenes <- vec_genesAll%in%vec_genesModule
    t.test(x=vec_celltypeScores[vec_idxModGenes],
           y=vec_celltypeScores[!vec_idxModGenes],
           alternative= "greater",
           conf.level = 0.95)[["p.value"]]
  }, simplify=T)
}

args = list("X"=df_celltypeGenes[-1])

# mat_t.test <- safeParallel(fun=fun, args=args, 
#                            list_vec_geneWeights=list_vec_geneWeights,
#                            vec_genesAll = vec_genesAll,
#                            simplify=T)
mat_t.test <- sapply(X = df_celltypeGenes[-1], FUN = fun, simplify = T)

mat_t.test %>% as.data.frame %>% 
  melt %>% "[["(2) %>% p.adjust(., method ="fdr") -> tmp # correct for n celltypes * n modules
mat_t.testAdj <- matrix(tmp, nrow=nrow(mat_t.test)) # matrix builds the matrix column by column
dimnames(mat_t.testAdj) <- dimnames(mat_t.test)


########################### Wilcoxon test ##########################

# Test: do the genes in the module have higher cell-type SEM rank than genes not in the module?

fun <- function(vec_celltypeScores) {
  sapply(list_vec_geneWeights, function(vec_geneWeights) {
    vec_genesModule = names(vec_geneWeights)
    vec_idxModGenes <- vec_genesAll%in%vec_genesModule
    wilcox.test(x=vec_celltypeScores[vec_idxModGenes],
           y=vec_celltypeScores[!vec_idxModGenes],
           alternative= "greater",
           conf.level = 0.95)[["p.value"]]
  }, simplify=T)
}

args = list("X"=df_celltypeGenes[-1])

# mat_wilcox.test <- safeParallel(fun=fun, 
#                                 args=args, 
#                            list_vec_geneWeights=list_vec_geneWeights,
#                            vec_genesAll = vec_genesAll,
#                            simplify=T)

mat_wilcox.test <- sapply(X = df_celltypeGenes[-1], FUN = fun, simplify = T)

mat_wilcox.test %>% as.data.frame %>% 
  melt %>% "[["(2) %>% p.adjust(., method ="fdr") -> tmp # correct for n celltypes * n modules
mat_wilcox.testAdj <- matrix(tmp, nrow=nrow(mat_wilcox.test)) # matrix builds the matrix column by column
dimnames(mat_wilcox.testAdj) <- dimnames(mat_wilcox.test)

######################################################################
####################### full SEM-based test  #########################
######################################################################

############################### KS test ##############################

# Test: do the module genes rank higher in the SEM genes than by change?

# Get full celltype genes
fun = function(vec_celltypeScores) {
  vec_celltypeScores %>% order(.,decreasing = T) -> vec_idxCelltypeOrder 
  vec_celltypeScoresOrder <- vec_celltypeScores[vec_idxCelltypeOrder]
  names(vec_celltypeScoresOrder) <- vec_genesAll[vec_idxCelltypeOrder]
  vec_celltypeScoresOrder
}
args <- list("X"=df_celltypeGenes[,-1])
list_vec_celltypeGenesFull <- safeParallel(fun=fun, args=args, simplify=F)

# get module genes
list_vec_modGenes <- sapply(list_vec_geneWeights, names)

# Run GSEA
fun = function(vec_celltypeScoresOrder) { 
  try(iterative.bulk.gsea(
    values = vec_celltypeScoresOrder, 
    set.list = list_vec_modGenes))
}

args = list("X"=list_vec_celltypeGenesFull)
list_GSEAsem <- safeParallel(fun=fun,args=args, simplify=F)

# Get adjusted p-values
fun <- function(GSEAsem) {
  n = nrow(GSEAsem) * length(list_GSEAsem)  
  GSEAsem[["p.val"]] %>% p.adjust(., method="fdr",n=n) -> GSEAsemPadj
  data.frame("p.adj" = GSEAsemPadj,genes=rownames(GSEAsem), row.names=NULL)
}

args= list("X"=list_GSEAsem)
list_df_GSEAsemPadj <- safeParallel(fun=fun,args=args, simplify=F)

fun = function(df1, df2) dplyr::full_join(df1, df2)

#  make into a matrix
list_df_GSEAsemPadj %>% Reduce(f=fun, x=.) %>% as.matrix -> mat_GSEAsemPadj

######################################################################
####################### full kME-based test  #########################
######################################################################

############################### KS test ##############################

# Test: do the top cell-type SEM genes higher in the module kME than by chance?

# Get top celltype genes
fun = function(vec_celltypeScores) {
  vec_celltypeScores %>% order(.,decreasing = T) %>% '['(1:SEMAbsCutoff) -> vec_idxCelltypeGenes 
  vec_genesAll[vec_idxCelltypeGenes]
}
args <- list("X"=df_celltypeGenes[,-1])
list_vec_celltypeGenes <- safeParallel(fun=fun, args=args, simplify=F)

# Get full kMs as list of named vectors
fun = function(vec_kMs) {
  names(vec_kMs) <-  df_kMs[[1]]
  vec_kMs<-vec_kMs[!is.na(vec_kMs)]
}

args <- list("X"=df_kMs[,-1])
list_vec_kMs <- safeParallel(fun=fun, args=args, simplify=F)

# Run GSEA
fun = function(vec_kMs) { 
  try(iterative.bulk.gsea(
    values = vec_kMs, 
    set.list = list_vec_celltypeGenes))
}
 
args = list("X"=list_vec_kMs)
list_GSEAkME <- safeParallel(fun=fun,args=args, simplify=F)

# Get adjusted p-values
fun <- function(GSEAkME) {
  n = nrow(GSEAkME) * length(list_GSEAkME)  
  GSEAkME[["p.val"]] %>% p.adjust(., method="fdr",n=n) -> GSEAkMEpadj
  data.frame("p.adj" = GSEAkMEpadj,genes=rownames(GSEAkME), row.names=NULL)
}

args= list("X"=list_GSEAkME)
list_df_GSEAkMEpadj <- safeParallel(fun=fun,args=args, simplify=F)

fun = function(df1, df2) dplyr::full_join(df1, df2)

#  make into a matrix
list_df_GSEAkMEpadj %>% Reduce(f=fun, x=.) %>% as.matrix -> mat_GSEAkMEadj
  
######################################################################
######################## Plot and output results #####################
######################################################################

pvalThreshold <- 0.05

mat_consensus <- mat_t.testAdj < pvalThreshold & 
  mat_FisherAdj < pvalThreshold & 
  mat_wilcox.testAdj < pvalThreshold & 
  GSEAsemPadj < pvalThreshold & 
  mat_GSEAkMEpadj < pvalThreshold

ht1_consensus <- Heatmap(mat_consensus,
                        cluster_rows = T,
                        cluster_columns = F, 
                        show_row_dend = F, 
                        show_column_dend = F, 
                        show_heatmap_legend = T, 
                        show_row_names = T, 
                        show_column_names = T,
                        #bottom_annotation = htca_noLabel,
                        use_raster=T,
                        raster_device = c("png"),
                        raster_quality = 2,
                        heatmap_legend_param = list(title = "Module-celltype specificity"))#

