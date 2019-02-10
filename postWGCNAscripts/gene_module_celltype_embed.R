# Make celltype mod embed
# Use gene loadings from one data (sub)set and possibly expression matrix from another

# AVERAGE EMBEDDINGS FOR datExpr_ident_cols

# if (FALSE) {
#   if (!is.null(datExpr_ident_cols)) {
#     
#     message(paste0("Computing expression averaged over ",paste0(datExpr_ident_cols, collapse=", ")))
#     
#     # remove NAs from ident columns
#     
#     for (ident_col in datExpr_ident_cols) {
#       metadata[[ident_col]][is.na(metadata[[ident_col]])] <- "Not_attributed"
#     }
#     
#     list_subsetModEmbed <- lapply(datExpr_ident_cols, function(ident_col) {
#       
#       cl <- makeCluster(spec=min(n_cores, detectCores_plus(Gb_max = RAMGbMax, additional_Gb = .5)-1), type="FORK", outfile = paste0(dirLog, prefixOut,"_", ident_col, "_mean_mod_expr.txt"))
#       
#       for (i in 1:length(sort(unique(ident)))) {
#         # cellModEmbed_loom$apply(name = paste0("row_attrs/",unique(ident)[i]), FUN = colMeans, MARGIN = 1, index.use = which(ident==unique(ident)[i]),
#         #             dataset.use = "matrix", display.progress = FALSE, overwrite = TRUE)
#         ident <- metadata[[ident_col]]
#         subEmbed <- if (scaleToZScore) cellModEmbed_loom[["layers/scale.data"]][ident==sort(unique(ident))[i],] else cellModEmbed_loom[["matrix"]][ident==sort(unique(ident))[i],] 
#         meanExpr <- matrix(parApply(cl, X=subEmbed, FUN=mean, MARGIN=2), nrow=1)
#         rownames(meanExpr) <- sort(unique(ident))[i]
#         colnames(meanExpr) <- cellModEmbed_loom[["row_attrs/module"]][]
#         # it might be ok to do this w/o parallelising
#         if (i==1) subsetModEmbed <- meanExpr else subsetModEmbed <- rbind(subsetModEmbed, meanExpr)
#       }
#       stopCluster(cl)
#       
#       subsetModEmbed
#       
#     })
#     
#     names(list_subsetModEmbed) <- datExpr_ident_cols
#   }
# }
# cellTypeModEmbed <- t(parSapply(cl, X=paste0("row_attrs/",unique(ident)), FUN = function(index) cellModEmbed_loom[[index]][], simplify = T))
# rownames(cellTypeModEmbed) <- gsub("col_attrs/", "", rownames(cellTypeModEmbed))
# colnames(cellTypeModEmbed) <- 

#rownames(cellTypeModEmbed) = sort(unique(ident)); colnames(cellTypeModEmbed) <- names(module_u)




# if (FALSE) {
#   if (!is.null(datExpr_ident_cols)) {
#     for (datExpr_ident in datExpr_ident_cols) {
#       write_csv(list_subsetModEmbed[[datExpr_ident]], file =  paste0(dirTables, prefixOut, "_", datExpr_ident, "_cellModEmbed.csv") )
#     }
#   }
# }