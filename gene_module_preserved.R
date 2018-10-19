## Gene module cross-dataset preservation analysis

# TODO: Copy structure from e.g. heatmap script
# TODO: Read preservation analysis description
# TODO: Pull out wrap preservation code from tmp_leftovers
# TODO: Change WGCNA_main script to output colors

# STATUS: down due to WGCNA::modulePreservation throwing a duplicate rownames error

# Usage: e.g.
# TODO 

# QUESTIONS
# Allow for multiple sets of modules? Yes, if no too complicated, output pairwise combinations. 
# Allow for multiple expression datasets? Could be useful and the function can handle it.


# outputs:
#   module preservation scores : n_wgcna_runs * n_datasets * n_modules vectors of preservation statistics

######################################################################
########################## DEFINE FUNCTIONS ##########################
######################################################################

source(file="/projects/jonatan/tools/functions-src/utility_functions.R")
# Get wrapmodulepreservation
source(file="/projects/jonatan/tools/wgcna-src/wgcna-toolbox/rwgcna_functions.R")
source(file="/projects/jonatan/tools/wgcna-src/wgcna-toolbox/tmp_modulePreservationFnc.R")
######################################################################
########################## SOURCE WGCNA PARAMS #######################
######################################################################

networkType="signed" # set by hand to what it was for run - needed for rwgcna_params.R
corFnc = "cor" # set by hand to what it was for run - needed for rwgcna_params.R
source(file="/projects/jonatan/tools/wgcna-src/wgcna-toolbox/rwgcna_params.R")

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("WGCNA", "dplyr", "Biobase", "Matrix", "parallel", "loomR", "ggplot2", "ComplexHeatmap"))

######################################################################
########################### OptParse #################################
######################################################################

option_list <- list(
  make_option("--path_data", type="character",
              help = "Path(s) to dataset(s) in which to evaluate preservation"), 
  make_option("--dir_project_WGCNA", type="character",
              help = "Directory of WGCNA projects which must contain /tables and /RObjects subdirectories. "), 
  make_option("--prefix_WGCNA_run", type="character", default = NULL,
              help = "Search for matching files within the dir_project_WGCNA. If left as NULL will take the first run [default %default]"),
  make_option("--dir_out", type="character", 
              help = "Directory for outputs, must contain /tables and /RObjects subdirectories"),
  make_option("--prefix_out", type="character", 
              help = "Prefix for output files"),
  make_option("--n_cores", type="integer", default = 10L,
              help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]")
)


######################################################################
############################ GET OPTIONS #############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

path_data <- opt$path_data
dir_project_WGCNA <- opt$dir_project_WGCNA
prefix_WGCNA_run <- eval(parse(text=opt$prefix_WGCNA_run))
dir_out <- opt$dir_out
prefix_out <- opt$prefix_out
n_cores <- opt$n_cores

######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

if (FALSE) { 
  path_data = c("/projects/jonatan/tmp-epilepsy/RObjects/EX_so2_UPF_QC.RData")
  dir_project_WGCNA = "/projects/jonatan/tmp-epilepsy/"
  prefix_WGCNA_run = c("UPF_5")
  dir_out = "/projects/jonatan/tmp-epilepsy/"
  prefix_out = "mod_preservation_1"
  n_cores = 10
}

######################################################################
############################# SET PARAMS #############################
######################################################################

options(stringsAsFactors = F)

######################################################################
###################### SET AND VERIFY PATHS ##########################
######################################################################

if (!file.exists(dir_project_WGCNA)) stop("dir_project_WGCNA not found")
if (!all(sapply(path_data,file.exists))) stop("one or more path_data not found")

# set up and verify WGCNA directory paths

dir_tables_WGCNA = paste0(dir_project_WGCNA,"tables/")
if (!file.exists(dir_tables_WGCNA)) stop("dir_project WGCNA must contain a /tables subdir")

dir_RObjects_WGCNA = paste0(dir_project_WGCNA,"RObjects/")
if (!file.exists(dir_RObjects_WGCNA)) stop("dir_project_WGCNA must contain a /RObjects subdir")

# set up expression data / output directories 
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
############################## LOAD DATA #############################
######################################################################

datExpr <- load_obj(path_data) 

if (class(datExpr)=="seurat") {
  datExpr@scale.data %>% t -> datExpr #NB: cell x gene
}

######################################################################
################# LOAD WGCNA MODULE ASSIGNMENT VECTORS ###############
######################################################################

WGCNA_session_image_path <- dir(path = dir_RObjects_WGCNA, pattern = paste0(prefix_WGCNA_run, "_final_session_image\\.RData"), full.names = T)

workspace_prior <- ls()

load(WGCNA_session_image_path)

#rm(list=ls()[!ls() %in% c(workspace_prior, "list_colors_meta")])

######################################################################
######################### DO MODULE PRESERVATION #####################
######################################################################

# testing a single dataset and multiple sets of colors

listDatExpr <- list(datExpr); rm(datExpr)
names(listDatExpr) <- paste0("datExpr_", prefix_out)

out <- wrapModulePreservation(listDatExpr=listDatExpr,
                              listColors=list_colors_meta,
                              labels = NULL,
                              dataIsExpr=T,
                              networkType=networkType, 
                              corFnc=corFnc,
                              corOptions=corOptions,
                              # We actually just want to get quality stats on the reference network,
                              # but the function requires testNetworks..
                              nPermutations=nPermutations, 
                              includekMEallInSummary=includekMEallInSummary,
                              restrictSummaryForGeneralNetworks=restrictSummaryForGeneralNetworks,
                              calculateQvalue=calculateQvalue,
                              randomSeed=randomSeed, 
                              maxGoldModuleSize=maxGoldModuleSize, 
                              maxModuleSize=maxModuleSize, 
                              quickCor=quickCor, 
                              ccTupletSize=ccTupletSize, 
                              #calculateCor.kIMall,
                              calculateClusterCoeff=calculateClusterCoeff,
                              useInterpolation=useInterpolation, 
                              checkData=checkData, 
                              #greyName, 
                              savePermutedStatistics=savePermutedStatistics , 
                              loadPermutedStatistics=loadPermutedStatistics, 
                              permutedStatisticsFile=permutedStatisticsFile, 
                              plotInterpolation=plotInterpolation, 
                              interpolationPlotFile=interpolationPlotFile, 
                              discardInvalidOutput=discardInvalidOutput,
                              parallelCalculation=parallelCalculation,
                              verbose=verbose, 
                              indent=indent)

