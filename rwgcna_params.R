# Title: Secondary script to set up WGCNA parameters
# Author: Jonatan Thompson, Pers Lab

############################## ADJACENCY #############################

#selectCols = NULL
# similarity = NULL # a (signed) similarity matrix: square, symmetric matrix with entries between -1 and 1.
type = networkType # a_{ij} = cor^Beta / "signed" # a_{ij} = ((cor+1)/2)^Beta / "signed hybrid" # a_{ij} = cor^Beta if cor>0; 0 otherwise / "distance" # Note that in some functions the parameter is called 'type', in others 'networktype'. We save the value under 'networkType'
# Contentious https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/TechnicalReports/signedTOM.pdf and http://www.peterlangfelder.com/signed-or-unsigned-which-network-type-is-preferable/.
# "signed" throws away negative correlations, while "unsigned" treats them equally (at least for the adjancy between two nodes) but this can be remediated at
# least for the part of the TOM value that takes into account neighbours by using signed TOM.
# The crux of the matter is that the TOM network has to have positive weights (?) meaning that at some point we'll have to jettison information
corOptions = if (corFnc == "cor") list(use = 'p') else NULL # "use = 'p', method = 'spearman'"  or list(use = 'p', method = 'spearman') to obtain Spearman correlation
#weights = NULL # 	 optional observation weights for datExpr to be used in correlation calculation. A matrix of the same dimensions as datExpr, containing non-negative weights. Only used with Pearson correlation.
#distFnc = "dist" #	character string specifying the function to be used to calculate co-expression similarity for distance networks. Defaults to the function dist. Any function returning non-negative values can be used. # TWEAKPARAM
#distOptions = "method = 'euclidean'" #	character string or a list specifying additional arguments to be passed to the function given by distFnc. For example, when the function dist is used, the argument method can be used to specify various ways of computing the distance.
#weightArgNames = c("weights.x", "weights.y") # character vector of length 2 giving the names of the arguments to corFnc that represent weights for variable x and y. Only used if weights are non-NULL.

############################## BICOR #################################

#x = NULL
#y = NULL
#robustX = TRUE
#robustY = TRUE
use = 'pairwise.complete.obs' #TWEAKPARAM
maxPOutliers = 0.1 # Suggested by package author in https://support.bioconductor.org/p/65124/
quick = 0
pearsonFallback = "individual" # in case of columns with zero MAD
#cosine = FALSE
#cosineX = cosine
#cosineY = cosine

############################## BLOCKWISEMODULES ######################

#blocks = NULL # optional specification of blocks in which hierarchical clustering and module detection should be performed. If given, must be a numeric vector with one entry per gene of multiExpr giving the number of the block to which the corresponding gene belongs.
maxBlockSize = 10e7 # integer giving maximum block size for module detection. Should be higher than the number of samples to avoid block-wise module detection. Ignored if blocks above is non-NULL. Otherwise, if the number of genes in datExpr exceeds maxBlockSize, genes will be pre-clustered into blocks whose size should not exceed maxBlockSize.
blockSizePenaltyPower = 10e-7 # number specifying how strongly blocks should be penalized for exceeding the maximum size. Set to a lrge number or Inf if not exceeding maximum block size is very important.
#nPreclusteringCenters = as.integer(min(ncol(datExpr0)/20, 100*ncol(datExpr0)/maxBlockSize))

loadTOM = FALSE # load TOM from previously saved file?

# Network construction arguments: correlation options

corType = if (corFnc == "cor") "pearson" else "bicor" # */!\* This same parameter is called 'corFnc' in adjancency().
# /!\: "bicor" gives: Error in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use = "all.obs",  : unused arguments (weights.x = NULL, weights.y = NULL) # TWEAKPARAM
quickCor = 0 # as defined for bicor. Real number between 0 and 1 that controls the handling of missing data in the calculation of correlations. 1 speeds up computation at the risk of large errors.
cosineCorrelation = FALSE # Logical: should the cosine version of the correlation calculation be used? The cosine calculation differs from the standard one in that it does not subtract the mean.

# Adjacency function options
#networkType = # Set value to same as in adjacency()
replaceMissingAdjacencies = TRUE # Should missing values in the calculation of adjacency be set to 0?
suppressTOMForZeroAdjacencies = FALSE

# Saving or returning TOM
getTOMs = NULL
saveTOMs = FALSE # should the networks (topological overlaps) be saved for each run?
#saveTOMFileBase = "TOM" # RObject name for a saved TOM

# Basic tree cut options
#detectCutHeight = 0.995 # Default value # TWEAKPARAM
#minModuleSize = minClusterSize # TWEAKPARAM

#useBranchEigennodeDissim = FALSE # should branch eigennode (eigengene) dissimilarity be considered when merging branches in Dynamic Tree Cut?
#minBranchEigennodeDissim = mergeCutHeight # Minimum consensus branch eigennode (eigengene) dissimilarity for branches to be considerd separate. The branch eigennode dissimilarity in individual sets is simly 1-correlation of the eigennodes; the consensus is defined as quantile with probability consensusQuantile.

#stabilityLabels = NULL # 	Optional matrix of cluster labels that are to be used for calculating branch dissimilarity based on split stability. The number of rows must equal the number of genes in multiExpr; the number of columns (clusterings) is arbitrary. See branchSplitFromStabilityLabels for details.
#stabilityCriterion = c("Individual fraction", "Common fraction") # One of c("Individual fraction", "Common fraction"), indicating which method for assessing stability similarity of two branches should be used. We recommend "Individual fraction" which appears to perform better; the "Common fraction" method is provided for backward compatibility since it was the (only) method available prior to WGCNA version 1.60.
#minStabilityDissim = NULL # Minimum stability dissimilarity criterion for two branches to be considered separate. Should be a number between 0 (essentially no dissimilarity required) and 1 (perfect dissimilarity or distinguishability based on stabilityLabels). See branchSplitFromStabilityLabels for details.

# Gene reassignment, module trimming, and module "significance" criteria

reassignThreshold = 1e-6
minCoreKME = 0.5 # a number between 0 and 1. If a detected module does not have at least minModuleKMESize genes with eigengene connectivity at least minCoreKME, the module is disbanded (its genes are unlabeled and returned to the pool of genes waiting for module detection). 0.5 is default.
#minCoreKMESize = minModuleSize/3 # see minCoreKME above. # TWEAKPARAM
minKMEtoStay = 0.2 # 0.3 # 	genes whose eigengene connectivity to their module eigengene is lower than minKMEtoStay are removed from the module. # TWEAKPARAM
kME_reassign_threshold = 1.25
# Module merging options

#mergeCutHeight = 0.2 # Dendrogram cut height for module merging. Note that the parameter moduleMergeCutHeight, used in mergeCloseModules, = mergeCutHeight. # TWEAKPARAM
# Value of 0.1 is used in Gandal,.., Geschwind et al 2018. Default 0.15.

############################## BOOTSTRAP #############################

# home-made - defined further down
#nPermutations = if (test_run==FALSE) 100 else 3 # TODO set to 100 as in Gandal,...,Geschwind et al 2018
startRunIndex = 1
#replace = T # Sample with replacement? Galdal et al use FALSE, default FALSE. Using replacement increases variance of the sampling mean. TWEAKPARAM
fraction = if (replace) 1.0 else 0.66 # TWEAKPARAM.

####################### CONSENSUSCALCULATION #########################

# Not (yet) used

############################ CONSENSUSKME ############################

signed = if (networkType == "signed" | networkType == "signed hybrid") TRUE else if (networkType == "unsigned") FALSE # In signed networks (TRUE), negative kME values are not considered significant and the corresponding p-values will be one-sided.
# In unsigned networks (FALSE), negative kME values are considered significant and the corresponding p-values will be two-sided.
useModules = NULL # Optional specification of module labels to which the analysis should be restricted. This could be useful if there are many modules, most of which are not interesting.
metaAnalysisWeights = NULL

# Correlation options

#corAndPvalueFnc = if (corFnc == "cor") corAndPvalue else bicorAndPvalue # TODO check that this works
corOptions = list()
corComponent = if (corFnc == "cor") "cor" else "bicor"
getQvalues = TRUE # should q-values (estimates of FDR) be calculated?
useRankPvalue = TRUE # Logical: should the rankPvalue function be used to obtain alternative meta-analysis statistics?
rankPvalueOptions = list(calculateQvalue = getQvalues, pValueMethod = "scale")
setNames = NULL # names for the input sets. If not given, will be taken from names(multiExpr). If those are NULL as well, the names will be "Set_1", "Set_2", ....
excludeGrey = TRUE

########### CONSENSUSTOM / BLOCKWISECONSENSUSMODULES #################

# Finds ConsensusModules by computing adjacenies and TOMs, aligning and 'merging' the TOM matrices.
# See https://rdrr.io/cran/WGCNA/man/consensusTOM.html

# TOM precalculation arguments, if available

individualTOMInfo = NULL
useIndivTOMSubset = NULL # If individualTOMInfo is given, this argument allows to only select a subset of the individual set networks contained in individualTOMInfo. It should be a numeric vector giving the indices of the individual sets to be used. Note that this argument is NOT applied to multiExpr.

# Save individual TOMs?

saveIndividualTOMs = if (plot_permuted == TRUE) TRUE else FALSE # unless we need to plot them
individualTOMFileNames = "individualTOM-Set%s-Block%b.RData"

# Consensus calculation options: network calibration

networkCalibration = "full quantile" # "single quantile",  "none" # network calibration method. One of "single quantile", "full quantile", "none" (or a unique abbreviation of one of them).
# Full quantile normalization, implemented in normalize.quantiles, adjusts the TOM matrices such that all quantiles equal each other (and equal to the quantiles of the component-wise average of the individual TOM matrices).

# Simple quantile calibration options

# calibrationQuantile = 0.95 # if networkCalibration is "single quantile", topological overlaps (or adjacencies if TOMs are not computed) will be scaled such that their calibrationQuantile quantiles will agree.

sampleForCalibration = TRUE # if TRUE, calibration quantiles will be determined from a sample of network similarities. Note that using all data can double the memory footprint of the function and the function may fail.
sampleForCalibrationFactor = 2000 # determines the number of samples for calibration: the number is 1/calibrationQuantile * sampleForCalibrationFactor. Should be set well above 1 to ensure accuracy of the sampled quantile.
getNetworkCalibrationSamples = FALSE # logical: should samples used for TOM calibration be saved for future analysis? This option is only available when sampleForCalibration is TRUE.

# Consensus definition

consensusQuantile = 0.5 # The desired quantile to use in the consensus similarity calculation. Lower is conservative. Gandal,..,Geschwind use 0.2, but this seems really low!  We use the median value.
# Set to quartile recommended for different datasets in http://www.genetics.ucla.edu/courses/statgene/networks/files/Langfelder-Thursday-ConsensusModules.pdf.
useMean = FALSE # should the consensus be determined from a (possibly weighted) mean across the data sets rather than a quantile?
#setWeights = NULL # Optional vector (one component per input set) of weights to be used for weighted mean consensus. Only used when useMean above is TRUE.

# Saving the consensus TOM

saveConsensusTOMs = TRUE # should the consensus topological overlap matrices for each block be saved and returned?
#consensusTOMFilePattern = "consensusTOM-block.%b.RData" # 	character string containing the file namefiles containing the consensus topological overlaps. The tag %b will be replaced by the block number.

# Internal handling of TOMs

useDiskCache = FALSE # should calculated network similarities in individual sets be temporarily saved to disk? Saving to disk is somewhat slower than keeping all data in memory, but for large blocks and/or many sets the memory footprint may be too big.
chunkSize = NULL # network similarities are saved in smaller chunks of size chunkSize.
cacheBase = ".blockConsModsCache" # character string containing the desired name for the cache files. The actual file names will consists of cacheBase and a suffix to make the file names unique.
cacheDir = "." # character string containing the desired path for the cache files.

# Alternative consensus TOM input from a previous calculation

consensusTOMInfo = NULL # optional list summarizing consensus TOM, output of consensusTOM. It contains information about pre-calculated consensus TOM. Supplying this argument replaces TOM calculation, so none of the individual or consensus TOM calculation arguments are taken into account.

# Basic tree cut options
checkMinModuleSize = TRUE # should sanity checks be performed on minModuleSize?

# Gene reassignment and trimming from a module, and module "significance" criteria

reassignThresholdPS = 1e-4 # per-set p-value ratio threshold for reassigning genes between modules.
trimmingConsensusQuantile = consensusQuantile # a number between 0 and 1 specifying the consensus quantile used for kME calculation that determines module trimming according to the arguments below.

equalizeQuantilesForModuleMerging = FALSE # specific to blockwiseConsensusModules
quantileSummaryForModuleMerging = "mean"  # specific to blockwiseConsensusModules

#Module merging options

#equalizeQuantilesForModuleMerging = FALSE # If TRUE, the quantiles of the eigengene correlation matrices (interpreted as a single vectors of non-redundant components) will be equalized across the input data sets. Note that although this seems like a reasonable option, it should be considered experimental and not necessarily recommended.
#quantileSummaryForModuleMerging = "median" # "mean" # One of "mean" or "median". If quantile equalization of the module eigengene networks is performed, the resulting "normal" quantiles will be given by this function of the corresponding quantiles across the input data sets.

mergeConsensusQuantile = consensusQuantile # consensus quantile for module merging. See mergeCloseModules for details.

# Output options

numericLabels = FALSE # should the returned modules be labeled by colors (FALSE), or by numbers (TRUE)?

# detectCutHeight	= cutHeight # dendrogram cut height for module detection
minHeight = 0.1 # TODO if we want to use blockwiseConsensusModules

########################### CUTREEHYBRID ############################
# See https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/BranchCutting/
# TODO: read up on this clustering algorithm

#minClusterSize = 10 # min(20, ncol(datExpr0)/2 ) # 40 is the value in Gandal et al, 2018 # Set value in notebook after evaluation
#deepSplit =  2 # integer value between 0 and 4. Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive. 
pamStage = T # Only used for method "hybrid". If TRUE, the second (PAM-like) stage will be performed. Default = TRUE. Increases specificity at the cost of sensitivity. Gandal et al 2018 set this to FALSE.
pamRespectsDendro = T # Only used if pamStage = TRUE. Logical, only used for method "hybrid". If TRUE, the PAM stage will respect the dendrogram in the sense that objects and small clusters will only be assigned to clusters that belong to the same branch that the objects or small clusters being assigned belong to. #TWEAKPARAM

# cutHeight = NULL # Maximum joining heights that will be considered. For method=="tree" it defaults to 0.99. For method=="hybrid" it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram. #TODO research this algorithm
# /!\ the parameter name is aliased wth the cutHeight parameter of the mergeCloseModules function which instead denotes maximum dissimilarity (i.e., 1-correlation) that qualifies modules for merging. If you need to set a non-default value, define it here as 'cutreeHybridCutHeight'

# Basic tree cut options - DEFINE THESE IN THE SCRIPT
method = "hybrid"
#distM = "dissTOM" # Only used for method "hybrid". The distance matrix used as input to hclust. If not given and method == "hybrid", the function will issue a warning and default to method = "tree".

# Advanced options /!\ TODO: check if we need to set these to non-default values
#maxCoreScatter = NULL #	 Only used for method "hybrid". Maximum scatter of the core for a branch to be a cluster, given as the fraction of cutHeight relative to the 5th percentile of joining heights.
#minGap = NULL # Only used for method "hybrid". Minimum cluster gap given as the fraction of the difference between cutHeight and the 5th percentile of joining heights.
#maxAbsCoreScatter = NULL # Only used for method "hybrid". Maximum scatter of the core for a branch to be a cluster given as absolute heights. If given, overrides maxCoreScatter.
#minAbsGap = NULL # 	Only used for method "hybrid". Minimum cluster gap given as absolute height difference. If given, overrides minGap.

#minSplitHeight = NULL # Minimum split height given as the fraction of the difference between cutHeight and the 5th percentile of joining heights. Branches merging below this height will automatically be merged. Defaults to zero but is used only if minAbsSplitHeight below is NULL.
#minAbsSplitHeight = NULL # Minimum split height given as an absolute height. Branches merging below this height will automatically be merged. If not given (default), will be determined from minSplitHeight above.

# External (user-supplied) measure of branch split
#externalBranchSplitFnc = NULL # Optional function to evaluate split (dissimilarity) between two branches. [...] This argument is only used for method "hybrid".
#minExternalSplit = NULL #	 Thresholds to decide whether two branches should be merged. It should be a numeric vector of the same length as the number of functions in externalBranchSplitFnc above. Only used for method "hybrid".
#externalSplitOptions = list() # Further arguments to function externalBranchSplitFnc. [If only one external function is specified in externalBranchSplitFnc above, externalSplitOptions can be a named list of arguments or a list with one component that is in turn the named list of further arguments for externalBranchSplitFnc[[1]]. [...] Only used for method "hybrid".
#externalSplitFncNeedsDistance = NULL # 	Optional specification of whether the external branch split functions need the distance matrix as one of their arguments. Either NULL or a logical vector with one element per branch split function that specifies whether the corresponding branch split function expects the distance matrix as one of its arguments. The default NULL is interpreted as a vector of TRUE. When dealing with a large number of objects, setting this argument to FALSE whenever possible can prevent unnecessary memory utilization.
#assumeSimpleExternalSpecification = TRUE # 	Logical: when minExternalSplit above is a scalar (has length 1), should the function assume a simple specification of externalBranchSplitFnc and externalSplitOptions? If TRUE, externalBranchSplitFnc is taken as the function specification and externalSplitOptions the named list of options. This is suitable for simple direct calls of this function. If FALSE, externalBranchSplitFnc is assumed to be a list with a single component which specifies the function, and externalSplitOptions is a list with one component that is in turn the named list of further arguments for externalBranchSplitFnc[[1]].

# Partitioning Around Medoids (PAM) stage options - these are set in the notebook script after evaluating the effects of different values
#useMedoids = FALSE # Only used for method "hybrid" and only if labelUnlabeled==TRUE. If TRUE, the second stage will be use object to medoid distance; if FALSE, it will use average object to cluster distance. Default FALSE recommended by Langfelder et al.
#maxPamDist = cutHeight # 	Only used for method "hybrid" and only if labelUnlabeled==TRUE. Maximum object distance to closest cluster that will result in the object assigned to that cluster. Defaults to cutHeight.
#respectSmallClusters = TRUE #	 Only used for method "hybrid" and only if labelUnlabeled==TRUE. If TRUE, branches that failed to be clusters in stage 1 only because of insufficient size will be assigned together in stage 2. If FALSE, all objects will be assigned individually.

############################ GENERAL WGCNA ###########################

# The following setting is important for WGCNA
options(stringsAsFactors = FALSE);
# Multiple threads possible in RStudio server but not in desktop
#enableWGCNAThreads(nThreads = n_cores)#
randomSeed = 12345

options(pairwise.complete.obs = T) # How to handle missing values in correlations
checkMissingData = TRUE # should data be checked for excessive missing values in genes and samples, and for genes with zero variance?
# Problem: in blockwiseConsensusModules, our consensus TOM seems to lose genes, which makes it impossible to find eigengenes

trapErrors = FALSE
verbose = 5 
indent = 0

############################## HCLUST ################################

# Takes the place of flashClust as of 2014

### 180611
#hclustMethod = "complete" # Gandal et al 2018 use "average".
###
# We call it hclustMethod because method is aliased with argument for cutreedynamic

################ HIERARCHICALCONSENSUSCALCULATION ####################

# Not implemented

################# HIERARCHICALCONSENSUSMODULES #######################

# Not implemented

###################### HIERARCHICALCONSENSUSTOM ######################

# Not implemented

######################### MERGECLOSEMODULES ##########################

# See also pquantile()

# Optional starting eigengenes
# MEs = NULL

# mergeCutHeight = 0.15 # defined in blockwise modules. value of 1-cor of eigengene correlation at which to merge modules. This is the default value

# Optional restriction to a subset of all sets
# useSets = NULL

# Input handling options
checkDataFormat = TRUE # If TRUE, the function will check exprData and MEs for correct multi-set structure. If single set data is given, it will be converted into a format usable for the function. If FALSE, incorrect structure of input data will trigger an error.
unassdColor = if (is.numeric(colors)) 0 else "grey"

# Options for eigengene network construction
useAbs = FALSE # 	Specifies whether absolute value of correlation or plain correlation (of module eigengenes) should be used in calculating module dissimilarity. Default: FALSE

# Options for constructing the consensus
equalizeQuantiles = TRUE # Should quantiles of the eigengene dissimilarity matrix be equalized ("quantile normalized")? The default is FALSE for reproducibility of old code; when there are many eigengenes (e.g., at least 50), better results may be achieved if quantile equalization is used. # TWEAKPARAM
quantileSummary = "mean" # One of "mean" or "median". Controls how a reference dissimilarity is computed from the input ones (using mean or median, respectively). # TWEAKPARAM

#TODO How does this actually work?

# Merging options
#moduleMergeCutHeight = mergeCutHeight # Maximum dissimilarity (i.e., 1-correlation) that qualifies modules for merging. /!\ This parameter is actually called cutHeight and hence is aliased with 'cutHeight' in cutreeHybrid.  Renamed to avoid confusion. # TWEAKPARAM
iterate = TRUE # Controls whether the merging procedure should be repeated until there is no change. If FALSE, only one iteration will be executed

# Output options
relabel = FALSE # Controls whether, after merging, color labels should be ordered by module size.
colorSeq = NULL # Color labels to be used for relabeling. Defaults to the standard color order used in this package if colors are not numeric, and to integers starting from 1 if colors is numeric.
getNewMEs = TRUE # 	Controls whether module eigengenes of merged modules should be calculated and returned.
getNewUnassdME = TRUE # When doing module eigengene manipulations, the function does not normally calculate the eigengene of the 'module' of unassigned ('grey') genes. Setting this option to TRUE will force the calculation of the unassigned eigengene in the returned newMEs, but not in the returned oldMEs.

############################ MODULEEIGENGENES ########################

# NB: If you want to call blockwiseModules where method="complete" rather than "average", you need to re-define the function!

impute = TRUE # should WGCNA call impute.knn to impute NA values for module eigengene calculation? Note that '0's will NOT be imputed and that this is only for calculating the eigengenes (which is sensitive to NAs)
nPC_wgcna = 10 # Number of principal components and variance explained entries to be calculated. Note that only the first principal component is returned; the rest are used only for the calculation of proportion of variance explained. The number of returned variance explained entries is currently min(nPC, 10)
align = "along average" # Controls whether eigengenes, whose orientation is undetermined, should be aligned with average expression (align = "along average", the default) or left as they are (align = ""). Any other value will trigger an error.
excludeGrey = FALSE # Should the improper module consisting of 'grey' genes be excluded from the eigengenes?
grey = if (is.numeric(colors)) 0 else "grey" # Value of colors designating the improper module. Note that if colors is a factor of numbers, the default value will be incorrect.
subHubs = TRUE #Controls whether hub genes should be substituted for missing eigengenes. If TRUE, each missing eigengene (i.e., eigengene whose calculation failed and the error was trapped) will be replaced by a weighted average of the most connected hub genes in the corresponding module. If this calculation fails, or if subHubs==FALSE, the value of trapErrors will determine whether the offending module will be removed or whether the function will issue an error and stop.
returnValidOnly = trapErrors  # logical; controls whether the returned data frame of module eigengenes contains columns corresponding only to modules whose eigengenes or hub genes could be calculated correctly (TRUE), or whether the data frame should have columns for each of the input color labels (FALSE).
#softPower = 6 # This is defined based on plotting the fit with scale-free topology
scale = TRUE # logical; can be used to turn off scaling of the expression data before calculating the singular value decomposition. The scaling should only be turned off if the data has been scaled previously, in which case the function can run a bit faster. Note however that the function first imputes, then scales the expression data in each module. If the expression contain missing data, scaling outside of the function and letting the function impute missing data may lead to slightly different results than if the data is scaled within the function. # TWEAKPARAM

######################### MODULEPRESERVATION #########################

# Compares test sets with a module assignments and test sets optionally with modules assigned.
# Computes the *preservation* of modules through statistics that compare adjacencies and overall module statistics.
# Also computes the *quality* of modules **in the reference set**. The quality statistics are calculated with respect
# to genes in common wiht the test set, hence the function returns *a set* of quality statistics for each reference-test pair.
# This may be counter-intuitive but it is convenient for returning the quality and preservation statistics together, even though
# the former only concern the reference network.

# multiData	  # expression data or adjacency data in the multi-set format (see checkSets). A vector of lists, one per set.
# Each set must contain a component data that contains the expression or adjacency data. If expression data are used,
# rows correspond to samples and columns to genes or probes. In case of adjacencies, each data matrix should be a symmetric
# matrix with entries between 0 and 1 and unit diagonal. Each component of the outermost list should be named.
# multiColor	# a list in which every component is a vector giving the module labels of genes in multiExpr.
# The components must be named using the same names that are used in multiExpr; these names are used top match labels to
# expression data sets.
dataIsExpr = TRUE # logical: if TRUE, multiData will be interpreted as expression data; if FALSE, multiData will be interpreted as adjacencies.
#referenceNetworks = c(2:length(multiExpr))
# A vector giving the indices of expression data to be used as reference networks.
# Reference networks must have their module labels given in multiColor.
#testNetworks = list(rep(1,length(referenceNetworks))) #	a list with one component per each entry in referenceNetworks above, giving the test networks in which
# to evaluate module preservation for the corresponding reference network. If not given, preservation will
# be evaluated in all networks (except each reference network). If referenceNetworks is of length 1,
# testNetworks can also be a vector (instead of a list containing the single vector).
# NB: We
#nPermutations = nPermutations # TODO nPermutations #	specifies the number of permutations that will be calculated in the permutation test. TODO: How costly is it?
includekMEallInSummary = FALSE  # logical: should cor.kMEall be included in the calculated summary statistics?
# Because kMEall takes into account all genes in the network, this statistic measures preservation
# of the full network with respect to the eigengene of the module.
# This may be undesirable, hence the default is FALSE.
restrictSummaryForGeneralNetworks = TRUE  # logical: should the summary statistics for general (not correlation) networks be restricted
# (density to meanAdj, connectivity to cor.kIM and cor.Adj)? The default TRUE corresponds to published work.
calculateQvalue = F # logical: should q-values (local FDR estimates) be calculated? Package qvalue must be installed for this calculation. Fails.
#  Note that q-values may not be meaningful when the number of modules is small and/or most modules are preserved.
maxGoldModuleSize = 1000 # maximum size of the "gold" module, i.e., the random sample of all network genes. TODO what is this?
maxModuleSize = 1000 # alternatively, max(levels(multiColor)) ? maximum module size used for calculations. Larger modules will be reduced
# by randomly sampling maxModuleSize genes.
ccTupletSize = 2 # 	tuplet size for co-clustering calculations. TODO ?
calculateCor.kIMall = FALSE # logical: should cor.kMEall be calculated? This option is only valid for adjacency input.
# If FALSE, cor.kIMall will not be calculated, potentially saving significant amount of time
# if the input adjacencies are large and contain many modules.
calculateClusterCoeff = FALSE # logical: should statistics based on the clustering coefficient be calculated? TOTO
# While these statistics may be interesting, the calculations are also computationally expensive.
useInterpolation = FALSE  # logical: should permutation statistics be calculated by interpolating an artificial set of evenly spaced modules?
# This option may potentially speed up the calculations, but it restricts calculations to density measures.
checkData = FALSE #TRUE  # logical: should data be checked for excessive number of missing entries? See goodSamplesGenesMS for details.
#greyName = NULL # label used for unassigned genes. Traditionally such genes are labeled by grey color or numeric label 0.
# These values are the default when multiColor contains character or numeric vectors, respectively.
savePermutedStatistics = TRUE
loadPermutedStatistics = FALSE # If interpolation is used (see useInterpolation above), the function can optionally generate
# diagnostic plots that can be used to assess whether the interpolation makes sense.
permutedStatisticsFile = if (useInterpolation) "permutedStats-intrModules.RData" else "permutedStats-actualModules.RData"
plotInterpolation = TRUE # file name to save the interpolation plots into.
interpolationPlotFile = "modulePreservationInterpolationPlots.pdf"
discardInvalidOutput = TRUE # should output columns containing no valid data be discarded? This option may be useful when input dataIsExpr is FALSE and some of the output statistics cannot be calculated. This option causes such statistics to be dropped from output.
parallelCalculation = FALSE  # Note that parallel calculations are turned off by default and will lead to somewhat DIFFERENT results
# than serial calculations because the random seed is set differently. For the calculation to actually run in parallel mode,
# a call to enableWGCNAThreads must be made before this function is called.
# Better control this locally

# out:
# The function returns a nested list of preservation statistics.
# At the top level, the list components are:
#
# quality
# observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally) q-values of quality statistics. All logarithms are in base 10.
# **relevant for evalating clustering parameters**
#
# preservation
# observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally) q-values of density and connectivity preservation statistics. All logarithms are in base 10.
#
# accuracy
# observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally) q-values of cross-tabulation statistics. All logarithms are in base 10.
#
# referenceSeparability
# observed values, Z scores, log p-values, Bonferoni-corrected log p-values, and (optionally) q-values of module separability in the reference network. All logarithms are in base 10.
# **relevant for evalating clustering parameters**
#
# testSeparability
# observed values, Z scores, p-values, Bonferoni-corrected p-values, and (optionally) q-values of module separability in the test network. All logarithms are in base 10.
# permutationDetails	results of individual permutations, useful for diagnostics

############# SAMPLEDHIERARCHICALCONSENSUSMODULES ####################

# Not implemented

############################ SEURAT ##################################

#min.cells = 5 # Filter out genes with few cells 
#do.center = if (corFnc == "bicor") FALSE else TRUE # for ScaleData() # update: just used scale_data==F
nPC_seurat = 120 # for RunPCA() and ProjectPCA
maxit = 1000 # for RunPCA() IRLBA package - default is 1000 https://cran.r-project.org/web/packages/irlba/irlba.pdf
fastpath = T # for RunPCA()
adj.p.val.threshold = 5e-2

############################### STRINGdb ##################################

PPI_pkME_threshold = 10e-2
#adj.p.val.threshold = 5e-2

############################### TOM ##################################

TOMType = "signed" # Takes into account the sign of the adjacency between neighbours
TOMDenom = "mean" # "min" gives the standard TOM described in Zhang and Horvath (2005), and "mean" in which the min
# function in the denominator is replaced by mean.
# The "mean" may produce better results but at this time should be considered experimental.