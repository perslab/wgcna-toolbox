# complex heatmap plotting

# Usage: todo

######################################################################
############################# SET OPTIONS ############################
######################################################################

options(stringsAsFactors = F, use="pairwise.complete.obs")

######################################################################
######################### UTILITY FUNCTIONS ##########################
######################################################################

source(file="/projects/jonatan/tools/functions-src/utility_functions.R")

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("optparse", "dplyr", "Matrix", "parallel", "ComplexHeatmap", "readr"))

######################################################################
########################### OptParse #################################
######################################################################

option_list <- list(
  make_option("--pathMatModEmbed", type="character", default = NULL,
              help = "Path to precomputed module embedding matrix, assumed to be cell x gene module"), 
  make_option("--pathDfEmbedMetadata", type="character", default = NULL,
              help = "Path to metadata for the embedding matrix stored in a dataframe with nrow equal to nrow of the matrix"), 
  make_option("--vec_colAnnotEmbed", type="character", default = NULL,
              help = "Colnames for indexing metadata dataframe for sample annotations"), 
  make_option("--vec_pathsDfNWA", type="character", 
              help = "Named vector of paths to dataframes containing WGCNA module genes and scores in long format, i.e. one row per gene e.g. rwgcna cell_cluster_module_genes.csv files, e.g.''c('mousebrain='/projects/jonatan/tmp-mousebrain/tables/example1.csv', 'tab_muris' = '/projects/jonatan/tmp-maca/example2.csv)''"),  
  make_option("--colModule", type="character", default = "module", 
              help = "The module name column in the cell_cluster_module_gene.csv files, [default %default]"), 
  make_option("--vec_colAnnotNWA", type="character", default=NULL,
              help = "Colnames for other module annotations such as cell cluster of origin, tissue, dataset, etc [default %default]"),
  make_option("--prefixOut", type="character", default = paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_", sample(x = 999, size = 1)),
              help = "Unique prefix for output files, [default %default]"), 
  make_option("--dirOut", type="character",
              help = "Project directory to which to write files. Should include subdirectories /tables, /RObjects, /plots, /log"),
  make_option("--RAMGbMax", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

pathMatModEmbed <- opt$pathMatModEmbed
pathDfEmbedMetadata <- opt$pathDfEmbedMetadata
vec_colAnnotEmbed <- opt$vec_colAnnotEmbed
if (!is.null(vec_colAnnotEmbed)) vec_colAnnotEmbed <- eval(parse(text=vec_colAnnotEmbed))
vec_pathsDfNWA <- eval(parse(text=opt$vec_pathsDfNWA))
colModule <- opt$colModule
vec_colAnnotNWA <- opt$vec_colAnnotNWA
if (!is.null(vec_colAnnotNWA)) vec_colAnnotNWA <- eval(parse(text=vec_colAnnotNWA))
prefixOut <- opt$prefixOut
dirOut <- opt$dirOut
RAMGbMax <- opt$RAMGbMax

# TODO: Got to here


list_subsetModEmbed <- NULL
if (FALSE){
  modDist <- dist(x=t(cellModEmbed_loom[["matrix"]][,]), method="euclidean", diag=F, upper=F) 
  modDendro <- hclust(d=modDist, method="average")
  modClust_order <- modDendro$order
  
  # cellDist <- dist(x=cellModEmbed_loom[["matrix"]][,], method="euclidean", diag=F, upper=F) 
  # cellDendro <- hclust(d=cellDist, method="average")
  # cellClust_order <- cellDendro$order
  
  cl <- makeCluster(spec=min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = as.numeric(max(sapply(list_subsetModEmbed, object.size())))/1024^3)-1), type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", ident_col, "_mean_mod_expr_clust.txt"))
  
  list_subsetModEmbed_modClustOrder <- parLapplyLB(cl, list_subsetModEmbed, function(subsetModEmbed) {
    cellType_modDist <- dist(x=t(subsetModEmbed), method="euclidean", diag=F, upper=F) 
    cellType_modDendro <- hclust(d=cellType_modDist, method="average")
    cellType_modDendro$order
  })
  
  list_subsetModEmbed_cellClustOrder <- parLapplyLB(cl, list_subsetModEmbed, function(subsetModEmbed) {
    cellType_cellDist <- dist(x=subsetModEmbed, method="euclidean", diag=F, upper=F) 
    cellType_cellDendro <- hclust(d=cellType_modDist, method="average")
    cellType_cellDendro$order
  })
  
  stopCluster(cl)
}



######################################################################
#################### PLOT CELL*MODULE MATRIX #########################
######################################################################

if (do_plot) {
  # generate colors for annotation
  colors_uniq <- gsub("\\d", "", colors()) %>% unique
  
  ######################## ROW ANNOTATION ##############################
  list_htra <- list_htra_link <- NULL 
  
  if (!is.null(datExpr_ident_cols)) {
    list_htra <- lapply(datExpr_ident_cols, function(ident_col) {
      ident <- metadata[[ident_col]]
      htra_colors <-   sample(if(length(unique(ident))<length(colors_uniq)) colors_uniq else colors(), size=length(unique(ident)), replace=F) 
      names(htra_colors) <- sort(unique(ident), decreasing=F)  
      
      htra <- rowAnnotation(ident_col = sort(ident, decreasing=F),
                            #annotation_legend_param = list(cellCluster = list(nrow = length(unique(ident)), title = "Cell cluster", title_position = "topcenter")),
                            col=list(ident_col = htra_colors),
                            width = unit(5, "mm"),
                            show_legend=F,
                            #annotation_legend_param = list(show_legend = F), 
                            show_annotation_name = F)#,
      htra
    })
    
    # make row annotation labels as a separate row annotation
    list_htra_link <- lapply(datExpr_ident_cols, function(ident_col)  {
      
      ident <- metadata[[ident_col]]
      htra_labels <- sort(unique(ident), decreasing=F)
      
      fontsize <- min(20, max(8, as.integer(500/length(unique(ident)))))
      
      ident <- metadata[[ident_col]]
      htra_label_pos <- cumsum(table(ident)) - as.integer(round(table(ident)/2,0))
      
      rowAnnotation(link = anno_link(which= "row",
                                     at = htra_label_pos, 
                                     labels = htra_labels,
                                     link_width = unit(3, "mm"),
                                     labels_gp=gpar(fontsize=fontsize), padding=0.5), 
                    width = unit(5, "mm") + max_text_width(labels, gp=gpar(fontsize=fontsize)))
    })
  }
  ######################## COLUMN ANNOTATION ###########################
  
  #  WGCNA run column annotation
  WGCNA_run <- unlist(mapply(function(names_run, cellType_module_u) rep(names_run, times= length(unlist(cellType_module_u, recursive = F))), names_run=prefixes_WGCNA_run, cellType_module_u=run_cellType_module_u, SIMPLIFY=T))
  htca_WGCNA_run_colors <- sample(if(length(prefixes_WGCNA_run)<length(colors_uniq)) colors_uniq else colors(), size=length(prefixes_WGCNA_run), replace=F) 
  names(htca_WGCNA_run_colors) <- sort(unique(prefixes_WGCNA_run, decreasing=F))
  
  # WGCNA celltype column annotation
  WGCNA_cellCluster <- rep(WGCNA_cellTypes,times= sapply(cellType_module_u, length))
  
  htca_WGCNA_cellCluster_colors <- sample(if (length(WGCNA_cellTypes)<length(colors_uniq)) colors_uniq else colors(), size=length(WGCNA_cellTypes), replace=F) 
  names(htca_WGCNA_cellCluster_colors) <- sort(WGCNA_cellTypes, decreasing=F) # no need to sort?
  
  # Prepare WGCNA celltype column annotation labels as separate column annotation
  htca_labels <- sort(WGCNA_cellTypes)
  htca_label_pos <- cumsum(table(WGCNA_cellCluster)) - as.integer(round(table(WGCNA_cellCluster)/2,0))
  
  # Add PPI enrichment column annotation
  if (plot_PPI_enrichment) {
    run_cellType_STRINGdb_output_all <- list()
    for (run in prefixes_WGCNA_run) {
      cellType_STRINGdb_output_all_paths <- dir(path = dir_tables_WGCNA, pattern = paste0(run, ".*STRINGsb_output_all\\.csv$"), full.names = T)
      run_cellType_STRINGdb_output_all[[run]] <- lapply(cellType_STRINGdb_output_all_paths, load_obj) # list of list of data.frame
    }
    names(run_cellType_STRINGdb_output_all) <- prefixes_WGCNA_run
  }
  
  cl <- makeCluster(spec=min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = 0.2)-1), type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "make_PPI_column_annotation.txt"))
  
  run_cellType_modPPIenrichment <- clusterMap(cl, function(cellType_mod_u, cellType_STRINGdb_output_all) {
    mapply(function(modname, STRINGdb_output_all) {
      STRINGdb_output_all[["q.value"]][grep(modname, STRINGdb_output_all[["colors"]])] %>% -log10
    },
    modname = names(cellType_mod_u), STRINGdb_output_all=cellType_STRINGdb_output_all, SIMPLIFY=F)
  }, cellType_mod_u = run_cellType_module_u, cellType_STRINGdb_output_all= run_cellType_STRINGdb_output_all, SIMPLIFY=F, .scheduling=c("dynamic"))
  
  stopCluster(cl)
  
  cellType_modPPIenrichment <- unlist(run_cellType_modPPIenrichment, recursive = F, use.names = T)
  modPPIenrichment <- unlist(cellType_modPPIenrichment, recursive = F, use.names = T)
  htca_modPPIenrichment_colors <- colorRamp2(c(0, 20), c("white", "red"))
  
  # Add MAGMA enrichment column annotation
  if (!is.null(plot_magma_enrichment)) {
    run_magma <- list()
    for (run in prefixes_WGCNA_run) {
      path_magma <- dir(pattern = paste0(run, "_magma\\.fdr\\.log\\.csv$"), path = dir_tables_WGCNA, full.names = T)
      run_magma[[run]] <- read.csv(file = path_magma, header = T, quote = "")
    }
    names(run_magma) <- prefixes_WGCNA_run
    
    cl <- makeCluster(spec=min(n_cores, detectCores_plus(Gb_max = RAM_Gb_max, additional_Gb = 0.2)-1), type="FORK", outfile = paste0(dir_log_data, prefix_out,"_", "make_magma_column_annotation.txt"))
    
    run_cellType_modMagmaEnrichment <- clusterMap(cl, function(cellType_mod_u, magma) {
      lapply(names(cellType_mod_u), function(modname){
        gwas_columns <- unique(sapply(plot_magma_enrichment, function(gwas_pattern) {
          grep(pattern=gwas_pattern, x = colnames(magma))
        }
        , SIMPLIFY=T))
        out <- magma[[grep(modname, magma[["module"]]), c(grep("module", colnames(magma)),gwas_columns), drop=T]] 
        out
      })
    }, cellType_mod_u = run_cellType_module_u, magma = run_magma, SIMPLIFY=F, .scheduling=c("dynamic"))
    
    stopCluster(cl)
  }
  
  cellType_modMagmaEnrichment <- unlist(run_cellType_modMagmaEnrichment, recursive = F, use.names = T)
  
  # deal with multiple GWAS magma columns 
  # also check module order 
  list_modMagmaEnrichment <- list()
  for (cname in colnames(cellType_modMagmaEnrichment[[1]])[-grep("module", colnames(cellType_modMagmaEnrichment[[1]]))]) { #assumes the GWAS are always the same
    list_modMagmaEnrichment[[cname]] <- sapply(cellType_modMagmaEnrichment, function(modMagmaEnrichment) modMagmaEnrichment[[cname]])
    names(list_modMagmaEnrichment[[cname]])  <- modMagmaEnrichment[["module"]] # list of named annotation vectors
  }
  
  
  # TODO: check that the magma vectors (of length module) are ordered the same way as the rest of the annotaiton
  # TODO: because we omit the neext unlist step, does the name scheme work out?
  #modMagmaEnrichment <- unlist(list_modMagmaEnrichment, recursive = F, use.names = T)
  # TODO: generate a vector of 'up' colors, one per gwas phenotype. Then generate a list of colours scheme vectors
  #htca_modMagmaEnrichment_colors <- colorRamp2(c(0, 20), c("white", "blue"))
  
  # Add metadata correlation column annotation
  #TODO 
  
  # Add WGCNA ident group column annotation
  htca_WGCNA_ident_group_df  <- NULL 
  if (!is.null(WGCNA_ident_cats)) {
    
    sapply(X=WGCNA_ident_cats, FUN = function(ident_cat) {
      level_vec <- rep(NA_character_, length=length(module_u))
      for (lvl in unique(ident_cat)){
        idx <- grepl(lvl, names(module_u))
        level_vec[idx] <- lvl
      }
    }, simplify = T) %>% data.frame -> htca_WGCNA_ident_df_tmp
    
    # Prepare WGCNA ident group colors
    list_htcs_WGCNA_ident_group_colors <- apply(X = htca_WGCNA_ident_df_tmp, MARGIN = 2, FUN = function(level_vec){
      out <- sample(if (length(unique(level_vec))<length(colors_uniq)) colors_uniq else colors(), size=length(unique(level_vec)), replace=F) 
      names(out) <- sort(unique(level_vec), decreasing=F)
      out
    })
    
  }
  
  # Finally gathercolumn annotations in a dataframe
  # TODO htca_df <- data.frame()
  
  # If there are common celltypes between the data and the WGCNA module celltypes, they get the same heatmap annotation colors
  if (!TODO) { 
    match_cellType_idx <- match(names(htca_WGCNA_cellCluster_colors), names(htra_cellCluster_colors), nomatch = NA)
    new_colors_tmp <- htra_cellCluster_colors[na.omit(match_cellType_idx)]
    if (length(new_colors_tmp) > 0) htca_WGCNA_cellCluster_colors[names(htca_WGCNA_cellCluster_colors) %in% names(htra_cellCluster_colors)] <- new_colors_tmp
  }
  
  # Make column annotation objects
  htca <- columnAnnotation(#WGCNA_run=WGCNA_run,
    WGCNA_cellCluster=WGCNA_cellCluster,
    # htca <- columnAnnotation(WGCNA_run=sort(WGCNA_run,  1:length(WGCNA_run), decreasing=F),
    #                          WGCNA_cellCluster=sort(WGCNA_cellCluster, decreasing=F),
    show_legend=F,
    link = anno_link(which = c("column"),
                     at = htca_label_pos, 
                     labels = htca_labels,
                     link_width = unit(0, "mm"),
                     labels_gp=gpar(fontsize=14), 
                     padding=0.3),  
    #col = list(WGCNA_run = htca_WGCNA_run_colors, WGCNA_cellCluster = htca_WGCNA_cellCluster_colors),
    col = list(WGCNA_cellCluster = htca_WGCNA_cellCluster_colors),
    annotation_height = unit.c(unit(1, "cm"), unit(1,"cm"), unit(1, "cm")))
  
  # When clustering module columns omit the rowannotation labels
  
  htca_noLabel <- columnAnnotation(WGCNA_run=sort(WGCNA_run,  1:length(WGCNA_run), decreasing=F),
                                   WGCNA_cellCluster=sort(WGCNA_cellCluster, decreasing=F),
                                   show_legend=F,
                                   col = list(WGCNA_run = htca_WGCNA_run_colors, WGCNA_cellCluster = htca_WGCNA_cellCluster_colors),
                                   annotation_height = unit.c(unit(1, "cm"), unit(1,"cm")))
  
  # Make plots
  ht1 <- Heatmap(cellModEmbed_loom[["matrix"]][order(ident, 1:nrow(cellModEmbed_loom[["matrix"]][,]), decreasing=F),order(WGCNA_run,  1:length(WGCNA_run), decreasing=F)], 
                 cluster_rows = F,
                 cluster_columns = F, 
                 show_heatmap_legend = T,
                 show_row_names = F, 
                 show_column_names = F,
                 bottom_annotation = htca,
                 use_raster=T,
                 raster_device = c("png"),
                 raster_quality = 1,
                 heatmap_legend_param = list(title = "Expression"))
  pdf(sprintf("%s%s_cellModEmbed.pdf", dir_plots_data, prefix_out), h=max(25, min(40, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2000)), w=max(30, min(40, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 50)))
  draw(ht1+list_htra[[1]]+list_htra_link[[1]])
  dev.off()
  
  ht1_cellClust <- Heatmap(cellModEmbed_loom[["matrix"]][order(ident, 1:nrow(cellModEmbed_loom[["matrix"]][,]), decreasing=F),order(WGCNA_run,  1:length(WGCNA_run), decreasing=F)], 
                           cluster_rows = T,
                           cluster_columns = F, 
                           show_row_dend = F, 
                           show_heatmap_legend = T,
                           show_row_names = F, 
                           show_column_names = F,
                           bottom_annotation = htca,
                           use_raster=T,
                           raster_device = c("png"),
                           raster_quality = 1,
                           heatmap_legend_param = list(title = "Expression"))
  pdf(sprintf("%s%s_cellModEmbed_cell_clust.pdf", dir_plots_data, prefix_out), h=max(25, min(40, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2000)), w=max(25, min(40, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 50)))
  draw(ht1_cellClust+htra)#+htra_link)
  dev.off()
  
  ht1_modClust <- Heatmap(cellModEmbed_loom[["matrix"]][order(ident, 1:nrow(cellModEmbed_loom[["matrix"]][,]), decreasing=F),order(WGCNA_run,  1:length(WGCNA_run), decreasing=F)], 
                          cluster_rows = F,
                          cluster_columns = T, 
                          show_column_dend = F, 
                          show_heatmap_legend = T,
                          show_row_names = F, 
                          show_column_names = F,
                          bottom_annotation = htca_noLabel,
                          use_raster=T,
                          raster_device = c("png"),
                          raster_quality = 1,
                          heatmap_legend_param = list(title = "Expression"))#,
  pdf(sprintf("%s%s_cellModEmbed_mod_clust.pdf", dir_plots_data, prefix_out), h=max(25, min(40, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2000)), w=max(25, min(40, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 50)))
  draw(ht1_modClust+htra+htra_link)
  dev.off()
  
  ht1_cellModClust <- Heatmap(cellModEmbed_loom[["matrix"]][order(ident, 1:nrow(cellModEmbed_loom[["matrix"]][,]), decreasing=F),order(WGCNA_run,  1:length(WGCNA_run), decreasing=F)], 
                              cluster_rows = T,
                              cluster_columns = T, 
                              show_row_dend = F, 
                              show_column_dend = F, 
                              show_heatmap_legend = T,
                              show_row_names = F, 
                              show_column_names = F,
                              bottom_annotation = htca_noLabel,
                              use_raster=T,
                              raster_device = c("png"),
                              raster_quality = 1,
                              heatmap_legend_param = list(title = "Expression"))
  pdf(sprintf("%s%s_cellModEmbed_cell_mod_clust.pdf", dir_plots_data, prefix_out), h=max(25, min(40, nrow(cellModEmbed_loom[["matrix"]][,]) %/% 2000)), w=max(25, min(40, ncol(cellModEmbed_loom[["matrix"]][,]) %/% 50)))
  draw(ht1_cellModClust+htra)#+htra_link)
  dev.off()
  
  ht2 <- Heatmap(cellTypeModEmbed[,], 
                 cluster_rows = F,
                 cluster_columns = F, 
                 show_heatmap_legend = T, 
                 show_row_names = T, 
                 show_column_names = F,
                 bottom_annotation = htca,
                 use_raster=T,
                 raster_device = c("png"),
                 raster_quality = 2,
                 heatmap_legend_param = list(title = "Expression"))
  pdf(sprintf("%s%s_cellTypeModEmbed.pdf", dir_plots_data, prefix_out), h=max(10, min(40, nrow(cellTypeModEmbed) %/% 4)), w=max(10, min(40, ncol(cellTypeModEmbed) %/% 50)))
  draw(ht2)
  dev.off()
  
  ht2_cellTypeClust <- Heatmap(cellTypeModEmbed[,], 
                               cluster_rows = T,
                               cluster_columns = F, 
                               show_row_dend = F, 
                               show_heatmap_legend = T, 
                               show_row_names = T, 
                               show_column_names = F,
                               bottom_annotation = htca,
                               use_raster=T,
                               raster_device = c("png"),
                               raster_quality = 2,
                               heatmap_legend_param = list(title = "Expression"))#
  pdf(sprintf("%s%s_cellTypeModEmbed_cellType_clust.pdf", dir_plots_data, prefix_out), h=max(10, min(40, nrow(cellTypeModEmbed) %/% 4)), w=max(10, min(40, ncol(cellTypeModEmbed) %/% 50)))
  draw(ht2_cellTypeClust)
  dev.off()
  
  ht2_modClust <- Heatmap(cellTypeModEmbed[,], 
                          cluster_rows = F,
                          cluster_columns = T, 
                          show_column_dend = F, 
                          show_heatmap_legend = T, 
                          show_row_names = T, 
                          show_column_names = F,
                          bottom_annotation = htca_noLabel,
                          use_raster=T,
                          raster_device = c("png"),
                          raster_quality = 2,
                          heatmap_legend_param = list(title = "Expression"))
  
  pdf(sprintf("%s%s_cellTypeModEmbed_mod_clust.pdf", dir_plots_data, prefix_out), h=max(10, min(40, nrow(cellTypeModEmbed) %/% 4)), w=max(10, min(40, ncol(cellTypeModEmbed) %/% 50)))
  draw(ht2_modClust)
  dev.off()
  
  ht2_cellModClust <- Heatmap(cellTypeModEmbed[,], 
                              cluster_rows = T,
                              cluster_columns = T, 
                              show_row_dend = F, 
                              show_column_dend = F, 
                              show_heatmap_legend = T, 
                              show_row_names = T, 
                              show_column_names = F,
                              bottom_annotation = htca_noLabel,
                              use_raster=T,
                              raster_device = c("png"),
                              raster_quality = 2,
                              heatmap_legend_param = list(title = "Expression"))#
  pdf(sprintf("%s%s_cellTypeModEmbed_cellType_mod_clust.pdf", dir_plots_data, prefix_out), h=max(10, min(40, nrow(cellTypeModEmbed) %/% 4)), w=max(10, min(40, ncol(cellTypeModEmbed) %/% 50)))
  draw(ht2_cellModClust)
  dev.off()
}
######################################################################
############################# WRAP UP ################################
######################################################################

message("Script done!")