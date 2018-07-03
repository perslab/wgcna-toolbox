# Title: Script to find robust WGCNA modules

######################################################################
############################## USAGE #################################
#######################o###############################################

# e.g.
# export R_MAX_NUM_DLLS=999

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-bmi-brain/data/hypo_topBMI.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/hypo_topBMI_3/  --data_prefix hypo_topBMI_3 --use.imputed FALSE --min.cells 6 --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod average --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(F)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --n_cores 9 --autosave TRUE --checkPPI TRUE --resume checkpoint_4

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-L5.top10_cts.seurat_obj.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/mousebrain-top10-3/  --data_prefix mousebrain-top10-3 --metadata_subset_col ClusterName --metadata_corr_col 'c("Age", "AnalysisPool", "AnalysisProject",  "CellConc", "Species", "Strain", "Subclass", "CellConc", "Class", "Description", "Developmental_compartment", "DonorID", "Neurotransmitter", "OriginalClusters", "PCRCycles", "Probable_location", "Project", "Region", "Sex", "Species", "Strain", "Subclass", "TaxonomyRank1", "TaxonomyRank2", "TaxonomyRank3", "TaxonomyRank4", "TaxonomyGroup", "Tissue")' --use.imputed FALSE --min.cells 6 --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod average --minClusterSize "c(15)" --deepSplit "c(4)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --n_cores 10  --autosave T

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-bmi-brain/data/hypo_topBMI.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/hypo_topBMI_2/  --data_prefix hypo_topBMI_2 --use.imputed FALSE --min.cells 6 --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod average --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(F)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --n_cores 18 --autosave TRUE --checkPPI TRUE --resume checkpoint_1


# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter_ensembl.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-9/  --data_prefix EP-9 --metadata_corr_col 'c("dataset")' --metadata_corr_filter_vals 'c("upf")' --use.imputed FALSE --min.cells 8 --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod complete --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.1)" --jackstrawnReplicate 500 --TOMnReplicate 0  --organism hsapiens --magma_gwas_dir /projects/jonatan/tmp-epilepsy/data/magma/ilae-lancet-2014/ --n_cores 20 --autosave TRUE --resume checkpoint_3 


# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-L5.top10_cts.seurat_obj.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/mousebrain-top10-2/  --data_prefix mousebrain-top10-2 --metadata_subset_col ClusterName --metadata_corr_col 'c("Age", "AnalysisPool", "AnalysisProject",  "CellConc", "Species", "Strain", "Subclass", "CellConc", "Class", "Description", "Developmental_compartment", "DonorID", "Neurotransmitter", "OriginalClusters", "PCRCycles", "Probable_location", "Project", "Region", "Sex", "Species", "Strain", "Subclass", "TaxonomyRank1", "TaxonomyRank2", "TaxonomyRank3", "TaxonomyRank4", "TaxonomyGroup", "Tissue")' --use.imputed FALSE --min.cells 10 --genes_use PCA --pca_genes all --corFnc bicor --networkType signed --hclustMethod average --minClusterSize "c(20)" --deepSplit "c(2)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits "c('BMI', 't1d', 't2d')" --n_cores 20  --autosave T

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter_ensembl_proteinCoding.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-10/  --data_prefix EP-10-protein --metadata_corr_col 'c("dataset")' --metadata_corr_filter_vals 'c("UPF")' --use.imputed FALSE --min.cells 8 --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod complete --minClusterSize "c(20)" --deepSplit "c(4)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0  --organism hsapiens --magma_gwas_dir /projects/jonatan/tmp-epilepsy/data/magma/ilae-lancet-2014/ --n_cores 10 --autosave TRUE --checkPPI T --resume checkpoint_1

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_pancreas.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-pancreas-12/  --data_prefix maca-pancr-12  --use.imputed FALSE --min.cells 5 --genes_use PCA --pca_genes all --corFnc bicor --networkType signed --hclustMethod single --minClusterSize "c(20)" --deepSplit "c(3)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits "c('BMI', 't1d', 't2d')" --n_cores 14 --autosave T --resume checkpoint_1

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_pancreas.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-pancreas-12/  --data_prefix maca-pancr-12  --use.imputed FALSE --min.cells 5 --genes_use PCA --pca_genes all --corFnc bicor --networkType signed --hclustMethod single --minClusterSize "c(20)" --deepSplit "c(3)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits "c('BMI', 't1d', 't2d')" --n_cores 14 --autosave T --resume checkpoint_1

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter_ensembl.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-9/  --data_prefix EP-9-noPPI --metadata_corr_col 'c("dataset")' --metadata_corr_filter_vals 'c("upf")' --use.imputed FALSE --min.cells 8 --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod complete --minClusterSize "c(20)" --deepSplit "c(1)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.1)" --jackstrawnReplicate 500 --TOMnReplicate 0  --organism hsapiens --magma_gwas_dir /projects/jonatan/tmp-epilepsy/data/magma/ilae-lancet-2014/ --n_cores 20 --autosave TRUE --resume checkpoint_4 --checkPPI F

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter_ensembl.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-9/  --data_prefix EP-9 --metadata_corr_col 'c("dataset")' --metadata_corr_filter_vals 'c("upf")' --use.imputed FALSE --min.cells 8 --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod complete --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.1)" --jackstrawnReplicate 500 --TOMnReplicate 0  --organism hsapiens --magma_gwas_dir /projects/jonatan/tmp-epilepsy/data/magma/ilae-lancet-2014/ --n_cores 20 --autosave TRUE --resume checkpoint_2 

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-bmi-brain/data/hypo_topBMI.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/hypo_topBMI_2/  --data_prefix hypo_topBMI_2 --use.imputed FALSE --min.cells 6 --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod average --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(F)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --n_cores 18 --autosave TRUE --checkPPI TRUE --resume checkpoint_1


# TMUX 67 time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-bmi-brain/data/hypo_topBMI.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/hypo_topBMI_1/  --data_prefix hypo_topBMI_1 --metadata_corr_col "c('major_cell_type')" --use.imputed FALSE --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod average --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.1)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --n_cores 18 --autosave TRUE 
# TMUX 68 time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-bmi-brain/data/hypo_topBMI.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/hypo_topBMI_2/  --data_prefix hypo_topBMI_2 --metadata_corr_col "c('major_cell_type')" --use.imputed FALSE --min.cells 5 --do.center FALSE --genes_use PCA --pca_genes all --corFnc bicor --networkType signed --hclustMethod average --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.1)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --n_cores 18 --autosave TRUE 
# TMUX 69 time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-bmi-brain/data/hypo_topBMI.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/hypo_topBMI_3/  --data_prefix hypo_topBMI_3 --metadata_corr_col "c('major_cell_type')" --use.imputed FALSE --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --hclustMethod single --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(FALSE)" --moduleMergeCutHeight "c(0.1)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --n_cores 18 --autosave TRUE 




# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-L5.top10_cts.seurat_obj.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/mousebrain-top10-/  --data_prefix mousebrain-top10-1 --metadata_subset_col ClusterName --metadata_corr_col 'c("Age", "AnalysisPool", "AnalysisProject",  "CellConc", "Species", "Strain", "Subclass", "CellConc", "Class", "Description", "Developmental_compartment", "DonorID", "Neurotransmitter", "OriginalClusters", "PCRCycles", "Probable_location", "Project", "Region", "Sex", "Species", "Strain", "Subclass", "TaxonomyRank1", "TaxonomyRank2", "TaxonomyRank3", "TaxonomyRank4", "TaxonomyGroup", "Tissue")' --use.imputed FALSE --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --hclustMethod complete --minClusterSize "c(20)" --deepSplit "c(1,2,3,4)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.1)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits "c('BMI', 't1d', 't2d')" --n_cores 10 --resume checkpoint_3 

# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-L5.top10_cts.seurat_obj.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/mousebrain-top10-4/  --data_prefix mousebrain-top10-4 --metadata_subset_col ClusterName  --use.imputed FALSE --min.cells 5 --do.center FALSE --genes_use PCA --pca_genes all --corFnc bicor --networkType signed --hclustMethod single --minClusterSize "c(20)" --deepSplit "c(1)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.1)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits "c('BMI', 't1d', 't2d')" --n_cores 10 --resume checkpoint_3 



# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-22/ --metadata_corr_col 'c("nGene", "nUMI", "percent.mito", "percent.ribo", "X2.group", "X3.batches", "X4.sex", "X5.Diet", "X6.FvF")' --metadata_corr_filter_vals "c('X2.group', 'X3.batches', 'X4.sex', 'X5.Diet', 'X6.FvF')" --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits 'c("BMI", "T1D", "T2D")' --data_prefix campbell-AgRP-22 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 0 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/EX_so2_mito0.1.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-6/ --metadata_subset_col "res.0.6" --magma_gwas_dir /projects/jonatan/tmp-epilepsy/data/magma/BMI-brain/ --metadata_corr_col 'c("dataset", "percent.mito", "nUMI")' --metadata_corr_filter_vals 'c("lake", "upf")' --data_prefix epilepsy-6 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 100 --replace T --organism hsapiens --n_cores 30 --resume checkpoint_3
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_cell_ontology_class_filterRibo0.05_min20cell.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-5/ --metadata_corr_col 'c("nGene", "nUMI")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --gwas_filter_traits 'c("BMI", "T1D", "T2D")'  --data_prefix maca-5 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 20 --resume checkpoint_3
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_n01_to_n04.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell_n01_to_n04-11/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-n01_to_n04-11 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.20)" --jackstrawnReplicate 0 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_clust_all_sub_n.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-clust-all-sub-n-1/ --metadata_corr_col 'c("nGene", "nUMI", "percent.mito", "percent.ribo", "X2.group", "X3.batches", "X4.sex", "X5.Diet", "X6.FvF")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-clust-all-sub-n-1 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 400 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 12 --resume checkpoint_2
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_cell_ontology_class_filterRibo0.05_min20cell.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-5/ --metadata_corr_col 'c("nGene", "nUMI")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix maca-5 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 20 --resume checkpoint_3
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_pancreas.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-pancreas-10/ --metadata_corr_col 'c("percent.mito", "percent.ribo")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix maca-pancreas-10 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 14 --resume checkpoint_2
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_clust_all.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-clust-all-2/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-clust-all-1 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 400 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 21 --resume checkpoint_4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_clust_all_sub_n.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-clust-all-sub-n-2/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-clust-all-sub-n-2 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 300 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 12
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_pancreas.Rdata --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-pancreas-9/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix maca-pancreas-9 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 400 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 14 --resume checkpoint_4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-4/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix epilepsy-4 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(3,4)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism hsapiens --n_cores 10
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_n01_to_n04.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell_n01_to_n04-10/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-n01_to_n04-10 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.20)" --jackstrawnReplicate 0 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_n05_to_n10.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/campbell_n05_to_n10-9/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-n05_to_n10-9 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.20)" --jackstrawnReplicate 0 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 6
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_AgRP_neurons.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-AgRP-6/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-AgRP-6 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(1,2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 0 --TOMnReplicate 0 --replace T --organism mmusculus--plot_permuted F --n_cores 4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_neurons_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-neurons-17/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-neurons-17 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus--plot_permuted F --n_cores 14 --resume checkpoint_4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-holst-hsl/RObjects/campbell_s_sub.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-campbell-glia-17/ --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix campbell-glia-17 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(20)" --deepSplit "c(2,3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 500 --TOMnReplicate 0 --replace T --organism mmusculus--plot_permuted F --n_cores 11 --resume checkpoint_4
# time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/maca/RObjects/maca_seurat_cell_ontology_class_filterRibo0.05_min20cell.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-maca-6/ --metadata_corr_col 'c("nGene", "nUMI")' --magma_gwas_dir /projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix maca-6 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(1,2)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 400 --TOMnReplicate 0 --replace T --organism mmusculus --n_cores 10 
#  time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter_ensembl_lake.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-lake-3/ --metadata_subset_col "res.0.6" --magma_gwas_dir projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/ --data_prefix epilepsy-lake-3 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes var.genes --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 400 --TOMnReplicate 0 --replace T --organism hsapiens --n_cores 13
#  time Rscript /projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_main.R --seurat_path /projects/jonatan/tmp-epilepsy/RObjects/seurat_obj_filter_ensembl_UPF.RData --project_dir /projects/jonatan/tmp-rwgcna-tests/tmp-epilepsy-UPF-3/ --metadata_subset_col "res.0.6" --magma_gwas_dir NULL --data_prefix epilepsy-UPF-3 --min.cells 5 --do.center TRUE --genes_use PCA --pca_genes all --corFnc cor --networkType signed --minClusterSize "c(15)" --deepSplit "c(3)" --pamStage "c(TRUE)" --moduleMergeCutHeight "c(0.2)" --jackstrawnReplicate 400 --TOMnReplicate 100 --replace T --organism hsapiens --n_cores 13

######################################################################
################# TEST PARAMS FOR MANUAL RUNS ########################
######################################################################

if (FALSE) { 
  seurat_path = "/projects/jonatan/tmp-holst-hsl/RObjects/campbell_n05_to_n10.RData"
  project_dir = "/projects/jonatan/tmp-rwgcna-tests/tmp-campbell-n05_to_n10-10/"
  data_prefix = "campbell-n10_to_n10-10"
  resume = NULL 
  autosave = T
  metadata_subset_col = NULL
  metadata_corr_col = NULL #c("nUMI", "nGene", "percent.mito", "percent.ribo")
  metadata_corr_filter_vals = "UPF"
  use.imputed = F
  min.cells = 5
  genes_use = "PCA"
  pca_genes = 'all'
  corFnc = "cor"
  networkType = "signed"
  anti_cor_action = NULL
  minClusterSize = c(20)
  deepSplit = c(2)
  pamStage = c(TRUE)
  moduleMergeCutHeight = c(0.2)
  jackstrawnReplicate = 0
  TOMnReplicate = 0
  organism = "mmusculus"
  magma_gwas_dir = "/projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/"
  gwas_filter_traits <- "c('BMI', 't1d', 't2d')"
  n_cores = 6
}

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  
  make_option("--seurat_path", type="character",
              help = "Provide full path to Rdata input file with Seurat object with lognormalized expression data in the @data slot"),
  make_option("--project_dir", type="character", default=NULL,
              help = "Optional. Provide project directory. Must have subdirs RObjects, plots, tables. If not provided, assumed to be dir one level up from input data dir."),
  make_option("--data_prefix", type="character", default=paste0(substr(gsub("-","",as.character(Sys.Date())),3,1000), "_rWGCNA_run"),
              help = "Dataset prefix for output files"),
  make_option("--autosave", type="logical", default=T,
              help = "Autosave session images at regular intervals to resume later? Defaults to TRUE."),
  make_option("--resume", type="character", default=NULL,
              help = "Resume from a checkpoint? Must have same path and data_prefix as provided. Options are 'checkpoint_1' - 'checkpoint_4'"),
  make_option("--quit_session", type="character", default=NULL,
              help = "Quit after saving session image at a checkpoint? Options are 'checkpoint_1' - 'checkpoint_4'"),
  make_option("--metadata_subset_col", type="character", default=NULL,
              help = "Specify a seurat@meta.data$... column to use for subsetting the Seurat object. If NULL (default) uses the @ident slot."),
  make_option("--metadata_corr_col", type="character", default=NULL,
              help = "Specify seurat_obj@meta.data$... column(s) for which to compute correlations with gene modules. Takes a character with a vector of seurat_obj@meta.data column names e.g. 'nUMI' or 'c('nUMI', 'Age')'. For factor or character metadata, each levels is analysed as a dummy variable, so exercise caution. Defaults to NULL."),
  make_option("--metadata_corr_filter_vals", type="character", default=NULL,
              help = "Specify one or more values within the seurat@meta.data$... column(s). Retain only modules which are significantly (anti-) correlated (at present, there is no threshold value for the correlation). Takes a character with a vector of meta.data column names e.g. 'Female' or 'c('fasted', 'HFD')'. Case-insensitive. Defaults to NULL."),
  make_option("--use.imputed", type="logical", default=F,
              help="Use data in the obj@imputed slot for the computations to replace the @data slot? If the @imputed slot is empty, will revert to the default (FALSE)"),
  make_option("--min.cells", type="integer", default=5L,
              help="What is the minimum number of cells in each subset in the data in which a gene should be detected to not be filtered out? Integer, defaults to 5."),
  # make_option("--do.center", type="logical", default=T,
  #             help="Use centered data? In either case data is scaled and nUMI and mitochrondrial genes are regressed out."),
  make_option("--genes_use", type="character", default="PCA",
              help="One of 'all', 'var.genes' or 'PCA' for genes with significant loading on at least one significant PC. 'All' is not recommended. Defaults to 'PCA'"), 
  make_option("--pca_genes", type="character", default="var.genes",
              help="'all' or 'var.genes'. 'All' is computationally expensive but allows for selecting genes based on PC loading p-values rather than magnitudes as with the default (var.genes)"), 
  make_option("--corFnc", type="character", default="bicor",
              help="Use 'cor' for Pearson or 'bicor' for midweighted bicorrelation function (https://en.wikipedia.org/wiki/Biweight_midcorrelation). Defaults to 'bicor'"), 
  make_option("--networkType", type="character", default = "signed",
              help="'signed' scales correlations to [0:1]; 'unsigned' takes the absolute value (but the TOM can still be 'signed'); 'signed hybrid' sets negative correlations to zero."),
  make_option("--anti_cor_action", type="character", default=NULL, 
              help = "Optional. 'kME_reassign' reassigns genes with a negative kME more than 1.5 the kME w.r.t. their own (primary) module. Should be used only with networkType 'signed hybrid'."),
  make_option("--hclustMethod", type="character", default="average",
              help = "Hierarchical clustering agglomeration method. One of 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA), 'median' (= WPGMC) or 'centroid' (= UPGMC). See hclust() documentation for further information. Defaults to 'average'"),
  make_option("--minClusterSize", type="character", default="15L",
              help = "Minimum genes needed to form a module, or an initial cluster before the Partitioning Around Medoids-like step. WGCNA authors recommend decreasing the minimum cluster size when using higher settings of deepSplit. Takes a character with a vector of integer values to try 'c(15,20,25)', defaults to 15."),
  make_option("--deepSplit", type="character", default="3L",
              help = "Controls the sensitivity of the cutreeDynamic/cutreeHybrid algorithm. Takes a character with a vector of integer values between 0-4 to try, e.g. 'c(1,2,3)', defaults to '3L'"),
  make_option("--moduleMergeCutHeight", type="character", default="c(0.2)",
              help = "Cut-off level for the variable (1-correlation) for merging eigengenes. Takes a character with a vector of double values to try, e.g. 'c(0.1, 0.2)'. Defaults to '0.2'"),
  make_option("--pamStage", type="character", default="c(TRUE)",
              help = "cutreeHybrid. Perform additional Partition Around Medroids step? Users favoring speciﬁcity over sensitivity (that is, tight clusters with few mis-classiﬁcations) should set pamStage = FALSE, or at least specify the maximum allowable object–cluster dissimilarity to be zero (in which case objects are assigned to clusters only if the object–cluster dissimilarity is smaller than the radius of the cluster). Takes a character with a vector of logicals to try, e.g. 'c(TRUE,FALSE)', default 'c(FALSE)'"),
  make_option("--jackstrawnReplicate", type="integer", default=500L,
              help = "Number of times to re-run PCA after permuting a small proportion of genes to perform empirical significance tests, i.e. the `JackStraw` procedure (see `pca_genes` above). Integer, defaults to 1000."),
  make_option("--TOMnReplicate", type="integer", default=100L,
              help = "Number of times to resample the dataset when finding the consensus TOM, defaults to 100"),
  # make_option("--replace", type="logical", default=TRUE,
  #             help = "In TOM, sample with replacement? Defaults to TRUE. If TRUE, uses all samples, if FALSE, uses 66% each time."),
  make_option("--checkPPI", type="logical", default=T,
              help="Valiate gene modules using Protein-Protein Interactions?"),
  make_option("--organism", type="character", default="mmusculus",
              help = "'hsapiens' or 'mmusculus'"),
  make_option("--magma_gwas_dir", type="character", default = NULL,
              help = "MAGMA input GWAS data directory, defaults to NULL. E.g. '/projects/jonatan/tmp-epilepsy/data/magma/ilae-lancet-2014/', '/projects/jonatan/tmp-bmi-brain/data/magma/BMI-brain/'. NULL skips the magma GWAS step"),
  make_option("--gwas_filter_traits", type="character", default = NULL,
              help = "Filter out modules not significantly correlated with matching gwas studies within the magma_gwas_dir. Takes a character with a vector of character names to match within the filename of the GWAS , e.g. 'body_BMI_Locke2015' or 'c('BMI', 'T1D', 'T2D')'. Case-insensitive. Defaults to NULL"),
  make_option("--n_cores", type="integer", default=5,
              help = "Number of cores to use for parallelization [default %default]")
)

######################################################################
############################## PACKAGES #############################
######################################################################

#suppressMessages(library(profvis)) #  for processing time and memory profiling
#suppressMessages(library(pryr)) # for processing time and memory profiling
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(WGCNA)) 

message("Packages loaded")

######################################################################
################### GET COMMAND LINE OPTIONS #########################
######################################################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 

opt <- parse_args(OptionParser(option_list=option_list))

seurat_path <- opt$seurat_path 

project_dir <- opt$project_dir

data_prefix <- opt$data_prefix 

autosave <- opt$autosave

resume <- opt$resume

quit_session <- opt$quit_session

metadata_subset_col <- opt$metadata_subset_col

metadata_corr_col <- opt$metadata_corr_col

if (!is.null(metadata_corr_col)) {
  metadata_corr_col <- eval(parse(text=metadata_corr_col))
}

metadata_corr_filter_vals <- opt$metadata_corr_filter_vals

if (!is.null(metadata_corr_filter_vals)) {
  metadata_corr_filter_vals <- eval(parse(text=metadata_corr_filter_vals))
}

use.imputed <- opt$use.imputed

min.cells <- opt$min.cells

#do.center <- opt$do.center

genes_use <- opt$genes_use

pca_genes <- opt$pca_genes

corFnc <- opt$corFnc 

networkType <- opt$networkType

anti_cor_action <- opt$anti_cor_action

hclustMethod <- opt$hclustMethod

minClusterSize <- eval(parse(text=opt$minClusterSize))

deepSplit <- eval(parse(text=opt$deepSplit))

moduleMergeCutHeight <- eval(parse(text=opt$moduleMergeCutHeight))

pamStage <- eval(parse(text=opt$pamStage))

jackstrawnReplicate <- opt$jackstrawnReplicate

TOMnReplicate <- opt$TOMnReplicate

checkPPI <- opt$checkPPI

magma_gwas_dir <- opt$magma_gwas_dir 

gwas_filter_traits <- opt$gwas_filter_traits

if (!is.null(gwas_filter_traits)) {
  gwas_filter_traits <- eval(parse(text=gwas_filter_traits))
}
organism <- opt$organism

n_cores <- opt$n_cores  

######################################################################
############################## CONSTANTS #############################
######################################################################

if (!file.exists(seurat_path) & is.null(resume)) stop("Input data path not found and no previous session to resume given")
    
# if no output directory provided, use the dir one level above that of the input file.

if (is.null(project_dir)) {
  pos <- tail(gregexpr("/", seurat_path)[[1]],2)[1]
  project_dir = substr(seurat_path, 1, pos)
}

# if specified output directory doesn't exist, create it 
if (!file.exists(project_dir)) {
  dir.create(project_dir) 
  message("Project directory not found, new one created")
}

plots_dir = paste0(project_dir,"plots/")
if (!file.exists(plots_dir)) dir.create(plots_dir) 

tables_dir = paste0(project_dir,"tables/")
if (!file.exists(tables_dir)) dir.create(tables_dir)

RObjects_dir = paste0(project_dir,"RObjects/")
if (!file.exists(RObjects_dir)) dir.create(RObjects_dir)

log_dir = paste0(project_dir,"log/")
if (!file.exists(log_dir)) dir.create(log_dir)

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

# Load parameter values and utility functions
source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
source(file = "/projects/jonatan/functions-src/functions.R")

# Entrezgene to ensemb_gene_id mapping file path for MAGMA
mapping_hs_entrez2ensembl_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz"

if (organism == "mmusculus") {
  # For mapping symbol to ensembl
  mapping_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
  # Synonyms
  mapping_mm_synonyms_filepath = "/data/genetic-mapping/ncbi/Mus_musculus.gene_info_symbol2ensembl.gz"
  # Mouse to human ortholog mapping  
  mapping_hs_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
  
  # mendelian and coding variants paths
  coding_variants_path <- "/projects/jonatan/genesets/coding_variants_ensembl_mm.RData"
  mendelian_genes_path <- "/projects/jonatan/genesets/mendelian_ensembl_mm.RData"

} else if (organism == "hsapiens") {
  
  mapping_hs_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz" # columns: ensembl_gene_id, entrezgene, hgnc_symbol
 
  # mendelian and coding variants paths
  mendelian_genes_path <- "/projects/jonatan/genesets/mendelian_ensembl_hs.RData"
  coding_variants_path <- "/projects/jonatan/genesets/coding_variants_ensembl_hs.RData"
}

try(disableWGCNAThreads())

if (is.null(resume)) {
  
  ######################################################################
  ############################ CHECK INPUT #############################
  ######################################################################
  
  if (min.cells < 0) stop("min.cells must be a non-negative integer") 
  if (min.cells > 25 | min.cells < 5) warning("Recommended values for min.cells between 5 and 10")
  
  if (!(sapply(c("all", "var.genes", "PCA"), function(x) grepl(x, genes_use, ignore.case=T)) %>% any())) {
    stop("genes_use must be one of 'all', 'var.genes' or 'PCA'")
  }
  
  if (!pca_genes %in% c('all', 'var.genes')) {
    warning("Invalid pca_genes argument, reverted to var.genes")
    pca_genes <- var.genes
  }
  
  if (!corFnc %in% c("cor", "bicor")) stop("corFnc must be one of 'cor' for Pearson's or 'bicor' for biweighted midcorrelation")
  
  if (!networkType %in% c('signed', 'unsigned', 'signed hybrid')) stop("networkType must be one of 'signed', 'unsigned' or 'signed hybrid' (not 'signed_hybrid')")
  
  if (length(c(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight))>4) warning("Comparing different parameters increases processing time")
  
  if (min(minClusterSize) < 5) stop("minClusterSize must be a vector of integers over 5")
  
  if (max(deepSplit) > 4 | min(deepSplit) < 0) stop("deepSplit must be a vector of integers between 0 and 4")
  
  if (min(moduleMergeCutHeight) < 0 | max(moduleMergeCutHeight) > 1) stop("moduleMergeCutHeight must be a vector of doubles between 0 and 1, recommended range is between 0.1 and 0.2")
  
  if (! (jackstrawnReplicate >= 0)) stop("jackstrawnReplicate must be 0 or higher")
  
  if (! (TOMnReplicate >= 0)) stop("TOMnReplicate must be non-negative")
  
  if (! (n_cores >= 0 & n_cores <= 100)) stop("n_cores must be in the range 0-100")
  
  ######################################################################
  ################# LOAD AND SUBSET SEURAT OBJECT ######################
  ######################################################################

  message("Loading and subsetting seurat object..")
  
  seurat_obj <- load_obj(f=seurat_path)
  
  ######################################################################
  ################# IF SYMBOL, REMAP TO ENSEMBL ID #####################
  ######################################################################

  if (!any(grepl("ENSG|ENSMUSG",rownames(seurat_obj@data)))) {

    if (is.null(organism)) {
      stop("Error: seurat object gene names are not in ensembl ID format. To remap from hgnc to ensembl, please provide an organism 'mmusculus' or 'hsapiens'")
    
    } else if (!is.null(organism)) {

      if (organism == "mmusculus") { 
        
        # Step 1: direct mapping
        mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
        mapping = data.frame(symbol=row.names(seurat_obj@data), ensembl = mapping_direct$ensembl_gene_id[ match(toupper(row.names(seurat_obj@data)), toupper(mapping_direct$gene_name_optimal)) ])
        
        # Step 2: map remaing using synonyms
        mapping_synonyms = read.delim(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T, stringsAsFactors = F)
        mapping$ensembl[ which(is.na(mapping$ensembl)) ] = mapping_synonyms$ensembl[ match( toupper(mapping$symbol[which(is.na(mapping$ensembl)) ]) , toupper(mapping_synonyms$symbol)) ]
        
        # save mapping file for reference and later use
        write.csv(mapping, file=sprintf("%s%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, organism), quote = F, row.names = F)
        
      } else if (organism == "hsapiens") {
        
        mapping_direct = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T, stringsAsFactors = F) # columns: ensembl_gene_id, entrezgene, hgnc_symbol
        # Step 1: direct mapping
        mapping = data.frame(symbol=row.names(seurat_obj@data), ensembl = mapping_direct$ensembl_gene_id[ match(toupper(row.names(seurat_obj@data)), toupper(mapping_direct$hgnc_symbol)) ])
        
        # save mapping file for reference and later use
        write.csv(mapping, file=sprintf("%s%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, organism), quote = F, row.names = F)
        
      }

      # Make a log of unmapped genes 
      log_not_mapped_filepath = sprintf("%s%s_not_mapped_to_ensembl_%s", log_dir, data_prefix, flag_date,".tab")
      log_duplicate_filepath = sprintf("%s%s_duplicate_ensembl_id_genenames_%s", log_dir, data_prefix, flag_date,".tab")
      
      df_not_mapped = mapping[is.na(mapping$ensembl),]
      write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
      
      # Make a log of duplicate genes
      idx_duplicate_genes <- duplicated(mapping$ensembl)
      df_duplicate <- mapping[idx_duplicate_genes,]
      write.table(df_duplicate,log_duplicate_filepath,quote=F,sep="\t",row.names=F)
      
      # Filter out unmapped and duplicate genes from Seurat object
      seurat_obj@data <- seurat_obj@data[!is.na(mapping$ensembl) & !idx_duplicate_genes,]
      seurat_obj@raw.data <- seurat_obj@raw.data[!is.na(mapping$ensembl) & !idx_duplicate_genes,]
  
      # rename Seurat object rows where mapping was successful to ensembl ID
      rownames(seurat_obj@data) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
      rownames(seurat_obj@raw.data) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
      
      # Save ensembl seurat object for later use
      save(seurat_obj, file=sprintf("%sseurat_obj_ensembl.RData",RObjects_dir))
      
      }
    } else if (any(grepl("ENSG", rownames(seurat_obj@data)))) {
      
      message("Homo sapiens ensembl id gene names detected")
      organism <- "hsapiens"
      # Save ensembl seurat object for later use
      save(seurat_obj, file=sprintf("%sseurat_obj_ensembl.RData",RObjects_dir))
    } else if (any(grepl("ENSMUSG", rownames(seurat_obj@data)))) {
      message("mus musculus ensembl id gene names detected")
      organism <- "mmusculus"
      # Save ensembl seurat object for later use
      save(seurat_obj, file=sprintf("%sseurat_obj_ensembl.RData",RObjects_dir))
    }
  
  # Use imputed data if user opted to do so 
  if (is.null(seurat_obj@imputed) | any(dim(seurat_obj@imputed)) == FALSE) use.imputed <- F
  if (use.imputed==T) seurat_obj@data <- seurat_obj@imputed
  
  # Set seurat object ident
  if (!is.null(metadata_subset_col)) {
    seurat_obj <- SetAllIdent(object = seurat_obj,id = metadata_subset_col)
  }
  sNames <- names(table(seurat_obj@ident))
  
  # Convert any character or factor meta.data to numeric dummy variables each level with its own numeric column
  if (!is.null(metadata_corr_col)) {
    if (any(colnames(seurat_obj@meta.data) %in% metadata_corr_col)) {
      metadata <- matrix(NA, nrow=nrow(seurat_obj@meta.data), ncol=1)
      include <- seurat_obj@meta.data[,colnames(seurat_obj@meta.data) %in% metadata_corr_col, drop=F]
      for (i in 1:ncol(include)) {
        if (class(include[,i]) %in% c("factor", "character")) {
          metadata <- cbind(metadata, factorToIndicator(include[,i]))        
        } else {
          metadata <- cbind(metadata, include[,i, drop=F])
        }
      }
      rownames(metadata) <- rownames(seurat_obj@meta.data)
      metadata <- metadata[,-1, drop=F]
      if (!is.null(metadata_corr_filter_vals)) metadata <- metadata[, toupper(colnames(metadata)) %in% toupper(metadata_corr_filter_vals), drop = F]
      metadata <- as.data.frame(metadata)
      } else metadata <- NULL
    } else metadata <- NULL
     
    # meta <- seurat_obj@metadata[,colnames(seurat_obj@metadata) %in% meta.data_corr_col, drop=F]
    # model = paste0("~ ", paste0(metadata_corr_col, " + ", collapse=""), collapse="") 
    # model<- substr(x = model, start = 1, stop = nchar(model)-3)
    # metadata <- model.matrix(eval(parse(text=model)), data=meta)
    # metadata <- metadata[,-1, drop=F] # remove the intersect

    # #formula <- paste0("~", paste0(colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %in% metadata_corr_col], "+"))
    # meta.numeric <- metadata[,sapply(metadata, function(x) x[!class(x) %in% c("factor", "character")])]
    # #metadata <- model.matrix(metadata_corr_col[metadata_corr_col %in% colnames(seurat_obj@meta.data)], data=seurat_obj@meta.data)
    # if (any(sapply(metadata, function(x) class(x) %in% c("character", "factor"), simplify = T))) {
    #   for (j in colnames(metadata[sapply(metadata, function(x) class(x) %in% c("character", "factor"), simplify = T)])) {
    #     metadata <- cbind(metadata, factorToIndicator(metadata[[j]]))
    #     metadata[[j]] <- NULL
    #   }
    # }


  ### 180523 v1.8_dev4
  # Check metadata_corr_filter_vals
  # if (!is.null(metadata)) if(!is.null(metadata_corr_filter_vals)) if (!any(metadata_corr_filter_vals %in% colnames(metadata))) {
  #   message("metadata_corr_filter_vals not found in metadata_corr_col. Will not filter modules on meta data correlation")
  #   metadata_corr_filter_vals <- NULL
  # }
  ###  
  
  
  # Subset the seurat object
  subsets <- lapply(sNames, function(x) SubsetData(seurat_obj,
                                                   ident.use = x, 
                                                   do.scale = F, 
                                                   do.center = F,
                                                   subset.raw = if (!is.null(seurat_obj@raw.data)) T else F))
  
  names(subsets) <- sNames 
  
  message("Seurat object loaded and subsetted")
  
  ######################################################################
  ##################### SAVE PARAMETERS TO FILE ########################
  ######################################################################
  
  subsets_n_cells <- data.frame(subset = sNames, n_cells = as.numeric(table(seurat_obj@ident)))
  
  userParams <- cbind(list("seurat_path" = seurat_path,
                           "project_dir" = project_dir,
                           "data_prefix" = data_prefix,
                           "autosave" = autosave,
                           "resume" = resume,
                           "quit_session" = quit_session,
                           "metadata_subset_col" = metadata_subset_col,
                           "metadata_corr_col" = metadata_corr_col,
                           "use.imputed" = use.imputed,
                           "min.cells" = min.cells,
                           "do.center" = do.center,
                           "genes_use" = genes_use,
                           "corFnc" = corFnc,
                           "networkType" = networkType,
                           "anti_cor_action" = anti_cor_action,
                           "hclustMethod" = hclustMethod,
                           "minClusterSize" = minClusterSize,
                           "deepSplit" = deepSplit,
                           "pamStage" = pamStage,
                           "moduleMergeCutHeight" = moduleMergeCutHeight,
                           "jackstrawnReplicate" = jackstrawnReplicate,
                           "TOMnReplicate" = TOMnReplicate,
                           #"replace" = replace,
                           "organism" = organism,
                           "magma_gwas_dir" = magma_gwas_dir,
                           "n_cores" = n_cores))
  
  builtInParams <- cbind(list("minCoreKME" = minCoreKME,
                              "minKMEtoStay" = minKMEtoStay,
                              "consensusQuantile" = consensusQuantile,
                              "fraction" = fraction,
                              "replace" = replace,
                              "impute" = impute,
                              "nPC_seurat" = nPC_seurat,
                              "TOMType" = TOMType,
                              "p.val.threshold" = p.val.threshold))
  
  # Save run params to file
  write.csv(userParams, file=sprintf("%s%s_INFO_user_parameters_%s.csv", log_dir, data_prefix, flag_date), quote = F, row.names = T)
  write.csv(builtInParams, file=sprintf("%s%s_INFO_built_in_parameters_%s.csv", log_dir, data_prefix, flag_date), quote = F, row.names = T)
  write.csv(subsets_n_cells, file=sprintf("%s%s_INFO_subsets_n_cells_%s.csv", tables_dir, data_prefix, flag_date), quote = F, row.names=F)
  
  ######################################################################
  ##################### SEURAT PROCESSING ##############################
  ######################################################################
  
  ### 180502_v1.7_1
  subsets <- lapply(subsets, function(x) FilterGenes(x, min.cells = min.cells))

  vars.to.regress = c("nUMI", "percent.mito", "percent.ribo")[c("nUMI", "percent.mito", "percent.ribo") %in% names(seurat_obj@meta.data)]

  subsets <- lapply(subsets, function(x) ScaleData(object = x,
                                                  vars.to.regress = if (length(vars.to.regress) > 0) vars.to.regress else NULL,
                                                  model.use="linear",
                                                  do.par=T,
                                                  num.cores = min(n_cores, detectCores()-1),
                                                  do.scale=T,
                                                  do.center=do.center))

  subsets <- lapply(subsets, function(x) FindVariableGenes(object = x,
                                                           x.low.cutoff = 0.0125,
                                                           x.high.cutoff = 8,
                                                           y.cutoff = 1,
                                                           do.plot=F,
                                                           display.progress=F))
  ###
  invisible(gc())

  if (grepl("PCA", genes_use, ignore.case = T)) {

    message("Performing Principal Component Analysis..")
    
    start_time <- Sys.time()
  
    cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "log_parPCA.txt"))

    tryCatch({
      subsets <- parLapply(cl, subsets, function(x) RunPCA(object = x,
                                                           pc.genes = if (pca_genes == 'all') rownames(x@data) else x@var.genes,
                                                           pcs.compute = min(nPC_seurat, (if (pca_genes == 'all') nrow(x@data) else length(x@var.genes)) %/% 2, ncol(x@data) %/% 2),
                                                           use.imputed = F, # if use.imputed=T the @imputed slot has been copied to @data
                                                           weight.by.var = F,
                                                           do.print = F,
                                                           seed.use = randomSeed,
                                                           maxit = maxit, # set to 500 as default
                                                           fastpath = fastpath)) 
      
    }, warning = function(c) {
      subsets <- parLapply(cl, subsets, function(x) RunPCA(object = x,
                                                           pc.genes = if (pca_genes == 'all') rownames(x@data) else x@var.genes,
                                                           pcs.compute = min(nPC_seurat, (if (pca_genes == 'all') nrow(x@data) else length(x@var.genes)) %/% 2, ncol(x@data) %/% 2),
                                                           use.imputed = F, # if use.imputed=T the @imputed slot has been copied to @data
                                                           weight.by.var = F,
                                                           do.print = F,
                                                           seed.use = randomSeed,
                                                           maxit = maxit*2, # set to 500 as default
                                                           fastpath = F))
    })
    ### 
    stopCluster(cl)
    invisible(gc())
    end_time <- Sys.time()
   
    message(sprintf("PCA done, time elapsed: %s seconds", round(end_time - start_time,2)))
    
    # Select significant PCs using empirical p-value based on JackStraw resampling
    # Source: https://rdrr.io/cran/Seurat/src/R/plotting.R
    # score the PCs using JackStraw resampling to get an empirical null distribution to get p-values for the PCs based on the p-values of gene loadings
    if (jackstrawnReplicate > 0) message(sprintf("Performing JackStraw with %s replications to select genes that load on significant PCs", jackstrawnReplicate))
    list_datExpr <- lapply(subsets, function(x) wrapJackStraw(seurat_obj_sub = x, n_cores = n_cores, jackstrawnReplicate = jackstrawnReplicate, p.val.threshold = p.val.threshold))
    
    invisible(gc())
    
  } else if (genes_use == "all") {
    list_datExpr <- lapply(subsets, function(x) seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rep(TRUE, nrow(x@scale.data))))
  } else if (genes_use == "var.genes") {
    list_datExpr <- lapply(subsets, function(x) seurat_to_datExpr(seurat_obj_sub = x, idx_genes_use = rownames(x@scale.data) %in% x@var.genes))
  } 
  
  # Clear up
  rm(subsets, userParams, builtInParams, subsets_n_cells, seurat_obj)
  invisible(gc())
 
  # Save or load session image 
  resume = "checkpoint_1"
  if (autosave==T) save(list = ls(all.names = TRUE)[-grep("data_prefix$|corFnc$|networkType$|anti_cor_action$|hclustMethod$|minClusterSize$|deepSplit$|moduleMergeCutHeight$|pamStage$|TOMnReplicate$|checkPPI$|organism$|magma_gwas_dir$|gwas_filter_traits$|n_cores$|resume$" , ls(), ignore.case=F) ],  envir = .GlobalEnv, file=sprintf("%s%s_checkpoint_1_image.RData", RObjects_dir, data_prefix))
  if (!is.null(quit_session)) if (quit_session=="checkpoint_1") quit(save="no")
    
} else if (!is.null(resume)) {
  
  if (resume == "checkpoint_1") {
    
    load((file=sprintf("%s%s_checkpoint_1_image.RData", RObjects_dir, data_prefix)))
    source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
    source(file = "/projects/jonatan/functions-src/functions.R")
    options(stringsAsFactors = F)
    disableWGCNAThreads() 
  }
}

if (resume == "checkpoint_1") {

  ######################################################################
  ###################### (UN) LOAD PACKAGES ############################
  ######################################################################
  
  # Unload packages
  try(detach("package:Seurat", unload=TRUE))

  # Free up DLLs
  invisible(R.utils::gcDLLs())

  ######################################################################
  ####### PICK SOFT THRESHOLD POWER FOR ADJACENCY MATRIX ###############
  ######################################################################
  message("Computing soft threshold powers to maximise the fit of a scale free topology to the adjacency matrix")
  
    list_softPower <- mapply(function(x,y) sft_for_par(datExpr=x, subsetName=y), list_datExpr, sNames, SIMPLIFY=F)
  
    ###
    invisible(gc())
    names(list_softPower) <- sNames
    
    # TOM files will be saved here by consensusTOM / blockwiseConsensusModules
    setwd(RObjects_dir)
    
    if (TOMnReplicate > 0) {
      
      ######################################################################
      ############### RESAMPLE THE DATA FOR ROBUSTNESS #####################
      ######################################################################
      
      message(sprintf("Resampling the data %s times with replace = %s", TOMnReplicate, replace))
      
      invisible(gc())
      
      cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_bootStrap.txt"))
      
      list_multiExpr <- parLapply(cl, list_datExpr, function(x) bootstrap(datExpr=x, 
                                                                          ### 180508_v1.8_dev2
                                                                          nPermutations = TOMnReplicate,
                                                                          ###
                                                                          replace = replace,
                                                                          fraction = fraction,
                                                                          randomSeed = randomSeed))
      stopCluster(cl)
      
      invisible(gc())
      
      names(list_multiExpr) <- sNames 
    
      ######################################################################
      ######### RUN CONSENSUSTOM ON THE RESAMPLED DATASETS #################
      ######################################################################
      
      message(sprintf("Computing the consensus Topological Overlap Matrix with %s permutations", TOMnReplicate))
    
      invisible(gc())
      cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_consensusTOM.txt"))
      list_consensus <- clusterMap(cl, function(x,y)  consensusTOM(multiExpr = x, 
                                                                  checkMissingData = checkMissingData,
                                                                  maxBlockSize = maxBlockSize, 
                                                                  blockSizePenaltyPower = blockSizePenaltyPower, 
                                                                  randomSeed = randomSeed,
                                                                  corType = corType,
                                                                  maxPOutliers = maxPOutliers,
                                                                  quickCor = quickCor,
                                                                  pearsonFallback = pearsonFallback,
                                                                  cosineCorrelation = cosineCorrelation,
                                                                  replaceMissingAdjacencies = replaceMissingAdjacencies,
                                                                  power = list_softPower[[y]],
                                                                  networkType = networkType,
                                                                  TOMDenom = TOMDenom,
                                                                  saveIndividualTOMs = saveIndividualTOMs,
                                                                  individualTOMFileNames = paste0(y, "_individualTOM-Set%s-Block%b.RData"),
                                                                  networkCalibration = networkCalibration,
                                                                  sampleForCalibration = sampleForCalibration,
                                                                  sampleForCalibrationFactor = sampleForCalibrationFactor,
                                                                  getNetworkCalibrationSamples = getNetworkCalibrationSamples,
                                                                  consensusQuantile = consensusQuantile,
                                                                  useMean = useMean,
                                                                  saveConsensusTOMs = saveConsensusTOMs,
                                                                  consensusTOMFilePattern = paste0(y,"_consensusTOM-block.%b.RData"),
                                                                  returnTOMs = F,
                                                                  useDiskCache = T,
                                                                  cacheDir = RObjects_dir,
                                                                  cacheBase = ".blockConsModsCache",
                                                                  verbose = verbose,
                                                                  indent = indent), 
                                   x=list_multiExpr, 
                                   y=sNames, 
                                   SIMPLIFY=F)
      stopCluster(cl)
      invisible(gc())
      
      # For each consensusTOM, get a logical vector where 'good_genes' that were used are TRUE
      list_goodGenesTOM_idx <- lapply(list_consensus, function(x) as.logical(x$goodSamplesAndGenes$goodGenes))
      
      # Use the goodGenesTOM_idx to filter the datExpr matrices
      list_datExpr_gg <- mapply(FUN=function(x,y) x[,y], x=list_datExpr, y=list_goodGenesTOM_idx, SIMPLIFY=F)
      
    } else if (TOMnReplicate==0) {
      
      message("Computing the Topological Overlap Matrix")
      
      ######################################################################
      ##### COMPUTE THE ADJACENCY MATRIX AND TOM WITHOUT RESAMPLING ########
      ######################################################################
      
      invisible(gc())
      
      cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_TOM_for_par.txt"))
      list_goodGenesTOM_idx <- clusterMap(cl, function(x,y) TOM_for_par(datExpr=x, 
                                                                        subsetName=y, 
                                                                        softPower=list_softPower[[y]]), 
                                          x=list_datExpr, 
                                          y=sNames, SIMPLIFY=F)
      stopCluster(cl)
      
      invisible(gc())
      
      list_datExpr_gg <- list_datExpr
  }
  
  message("Done computing the consensus Topological Overlap Matrix")
    
  if (TOMnReplicate > 0) rm(list_multiExpr) #TODO: we'll need it for plotting permuted
  
  ######################################################################
  ####################### CLUSTER ON THE TOM ###########################
  ######################################################################
  # in:
  #   consTomDS (from disk!)
  #
  # out:
  #   list_geneTree
  
  # Convert TOM to distance matrix
  list_dissTOM <- lapply(sNames, dissTOM_for_par)
  
  # Cluster
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parHclust.txt"))
  list_geneTree <- parLapply(cl, list_dissTOM, function(x) hclust(d=x, method=hclustMethod))
  stopCluster(cl)
  
  names(list_geneTree) = sNames # used for PlotDendro
  
  ######################################################################
  ######################### CUT OUT MODULES ############################
  ######################################################################
  # in:
  #   list_dissTOM
  #   list_geneTree
  #
  # out:
  #   list_list_plot_label_ok 
  #   list_list_MEs_ok 
  #   list_list_colors_ok 
  #   sNames_ok 
  #   list_datExpr_ok 
  #   list_geneTree_ok
  #   list_dissTOM_ok 

  message("Computing modules on the Topological Overlap Matrix")
  
  dims = sapply(list(minClusterSize, deepSplit, pamStage, moduleMergeCutHeight), length) # list of number of values per parameter
  n_combs <- prod(dims) # Total number of combinations of parameter values
  comb_list = vector("list", length = n_combs)
    
  # Make a vector of different parameter combinations
  k=1
  for (minModSize in 1:dims[1]) {
    for (ds in 1:dims[2]) { # Gandal et al 2018 use c(50,100, 200)
      for (pam in 1:dims[3]) {
        for (dthresh in 1:dims[4]) {
          comb_list[[k]] <- c(minClusterSize[minModSize],deepSplit[ds], pamStage[pam], moduleMergeCutHeight[dthresh])
          k = k+1
        }
      }
    }
  }
  
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_cutreeHybrid_for_vec.txt"))
  # nb: comb_list is in the environment - no need to pass it through clusterMap
  list_list_cutree <- clusterMap(cl, function(x,y) lapply(comb_list, function(z) cutreeHybrid_for_vec(comb=z, 
                                                                                                      geneTree=x, 
                                                                                                      dissTOM = y,
                                                                                                      maxPamDist = maxPamDist, 
                                                                                                      useMedoids = useMedoids)), 
                                 x=list_geneTree,
                                 y=list_dissTOM, 
                                 SIMPLIFY=F)
  stopCluster(cl)
  
  names(list_list_cutree) = sNames
  
  # Merge close modules
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parvecMergeCloseModules.txt"))
  list_list_merged <- parLapply(cl, sNames, function(x) lapply(1:n_combs, function(y) mergeCloseModules_for_vec(cutree=list_list_cutree[[x]][[y]], comb=comb_list[[y]], datExpr=list_datExpr_gg[[x]], excludeGrey = excludeGrey)))
  stopCluster(cl)
  invisible(gc())
  
  names(list_list_merged) = sNames
  
  # Produce labels for plotting the modules found by different parameters
  list_plot_label <- lapply(comb_list, function(x) plotLabel_for_vec(comb=x)) # list of labels for plot
  # Make a copy for each subset
  list_list_plot_label <- list()
  list_list_plot_label <- lapply(sNames, function(x) list_list_plot_label$x = list_plot_label)
  names(list_list_plot_label) = sNames
  
  # Extract the Module Eigengenes from the list returned by mergeCloseModules
  list_list_MEs <- lapply(list_list_merged, parGetMEs) # nested lapply
  names(list_list_MEs) = sNames
  
  # Extract the colors from the list returned by mergeCloseModules
  list_list_colors <- mapply(parGetColors, list_list_merged, list_datExpr_gg, SIMPLIFY=F)
  names(list_list_colors) <- sNames
  
  # Test if for a given vector all the genes are grey
  # count the grey
  list_vec_n_grey <- lapply(list_list_colors, count_grey_in_list_of_vec)
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  list_logical_params_ok <- mapply(function(x,y) as.logical(mapply(function(a,b) a!=length(b), a=x, b=y, SIMPLIFY=F)), x=list_vec_n_grey, y=list_list_colors, SIMPLIFY=F)
  logical_subsets_ok <- sapply(list_logical_params_ok, any, simplify = T)
  
  # First, for each run, remove results of any parametrisation for which all the genes were assigned to grey
  # If for a subset all parameters gave only grey modules, take it out of the top level list.
  ### 180510_v1.8_dev2_1
  list_list_plot_label_ok <- mapply(function(x,y) x[y], x = list_list_plot_label, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_cutree_ok <- mapply(function(x,y) x[y], x = list_list_cutree, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_MEs_ok <- mapply(function(x,y) x[y], x = list_list_MEs, y = list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]
  list_list_colors_ok <- mapply(function(x,y) x[y], x=list_list_colors, y=list_logical_params_ok, SIMPLIFY = F)[logical_subsets_ok]

  sNames_ok <- sNames[logical_subsets_ok]
  list_datExpr_ok <- list_datExpr_gg[logical_subsets_ok] # If for a subset all parameters gave only grey modules, take it out of the top level list.
  list_geneTree_ok <- list_geneTree[logical_subsets_ok]
  #list_dissTOM_ok <- list_dissTOM[logical_subsets_ok]
  
  # Assign gene names to each color vector
  list_list_colors_ok <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), dimension = NULL)), x = list_list_colors_ok, y = list_datExpr_ok, SIMPLIFY=F)
  
  # Make a warning and report if any whole subsets produced no modules
  discarded <- setdiff(sNames, sNames_ok)
  if (length(discarded > 0)) {
    warning(paste0(discarded, " were dropped because the script didn't find any modules  "))
    fileConn <- file(description = sprintf("%s%s_ERROR_report_%s.txt", log_dir,  data_prefix, flag_date), open = 'a')
    writeLines(paste0(discarded, " were dropped because the script didn't find any modules"), fileConn)
    close(fileConn)
  }
  
  # clear up
  rm(list_datExpr, list_list_cutree, list_list_MEs, list_datExpr_gg, list_geneTree, list_list_plot_label, list_dissTOM)#, list_idx_mapped_gg)
  invisible(gc())
    
  ######################################################################
  ###################### compute kMEs & pkMEs ##########################
  ######################################################################
  # in:
  #   list_list_MEs_ok
  #   list_datExpr_ok
  #
  # out:
  #   list_list_kMEs
  #   list_list_pkMEs
  
  message("Computing kMEs: correlation between each gene and each module 'eigengene'")
 
  # Compute kMEs
  invisible(gc())
  cl <- makeCluster(spec=n_cores, type = "FORK", outfile = paste0(log_dir, "log_parkMEs_parPkMEs.txt"))
  
  list_list_kMEs <- clusterMap(cl, 
                               function(x,y) parkMEs(list_MEs=x, 
                                                     datExpr=y), 
                               x=list_list_MEs_ok, 
                               y=list_datExpr_ok, 
                               SIMPLIFY = F)
  
  # Extract the primary kMEs - i.e. those of each gene w.r.t. the module to which it belongs
  # We use these to filter genes that we submit to STRINGdb 
  list_list_pkMEs <- clusterMap(cl, function(x,y) parPkMEs(list_kMEs=x, list_colors=y), x=list_list_kMEs, y=list_list_colors_ok, SIMPLIFY = F)
  stopCluster(cl)
  invisible(gc())
  
  list_list_pkMEs <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_pkMEs, y=list_list_plot_label_ok, SIMPLIFY=F)
  names(list_list_pkMEs) <- sNames_ok

  ######################################################################
  ########### Filter out genes with poor kME p-values ##################
  ######################################################################
  # Note: this seems to 'punish' smaller gene modules disproportionately due to lower sample size
  if (FALSE) {
    invisible(gc())
    
    cl <- makeCluster(spec=n_cores, 
                      type = "FORK", 
                      outfile = paste0(log_dir, "log_kMEs_pvals_par.txt"))
    
    list_list_ppvals <- clusterMap(cl, function(a,b) lapply(a, 
                                                    function(x) WGCNA::corPvalueStudent(cor = x, 
                                                                                        nSamples = nrow(b))),
                                 a = list_list_pkMEs,
                                 b = list_datExpr_ok,
                                 SIMPLIFY = F)
    
      # Name the p-value vectors
      list_list_ppvals <- clusterMap(cl,
                                        function(a,b) lapply(a, function(x) name_for_vec(to_be_named=x,
                                                                                        given_names = colnames(b),
                                                                                        dimension=NULL)),
                                      a = list_list_ppvals,
                                      b = list_datExpr_ok,
                                      SIMPLIFY=F)
  
      
    list_list_colors_ok_sig <- clusterMap(cl,
                                      function(a,b) mapply(function(x,y) ifelse(y<p.val.threshold,
                                                                                yes=x,
                                                                                no="grey"),
                                                           x = a,
                                                           y = b,
                                                           SIMPLIFY=F),
                                      a = list_list_colors_ok,
                                      b = list_list_ppvals,
                                      SIMPLIFY=F)
    # name the new color vectors
    list_list_colors_ok_sig <- clusterMap(cl,
                                          function(a,b) mapply(function(x,y) name_for_vec(to_be_named = x,
                                                                                          given_names = names(y),
                                                                                          dimension=NULL),
                                                               x = a,
                                                               y = b,
                                                               SIMPLIFY=F),
                                          a = list_list_colors_ok_sig,
                                          b = list_list_colors_ok,
                                          SIMPLIFY=F)
    stopCluster(cl)
    invisible(gc())
  
    list_list_colors_ok <- list_list_colors_ok_sig
    
    # TODO: either remove or fix section below
    
    ######################################################################
    ################### Compute new MEs, kMEs, pkMEs #####################
    ######################################################################
  
   # in:
    #   list_list_colors_ok_sig
    #   list_datExpr_ok
    #
    # out:
    #   list_list_MEs
    #   list_list_kMEs
    #   list_list_pkMEs
  
    
    invisible(gc())
    cl <- makeCluster(spec=n_cores, type = "FORK", outfile = paste0(log_dir, "log_parkMEs_parPkMEs.txt"))
    # Compute MEs
    
    list_list_MEs_ok <- clusterMap(cl,
                                   function(a,b) lapply(b, function(x) moduleEigengenes(expr = a, 
                                                                                     colors = x, 
                                                                                     excludeGrey = T)),
                                   a = list_datExpr_ok,
                                   b = list_list_colors_ok,
                                   SIMPLIFY=F)
      
    # Compute kMEs
    list_list_kMEs <- clusterMap(cl, 
                                 function(x,y) parkMEs(list_MEs=x, 
                                                       datExpr=y), 
                                 x=list_list_MEs_ok, 
                                 y=list_datExpr_ok, 
                                 SIMPLIFY = F)
    
    # Extract the primary kMEs - i.e. those of each gene w.r.t. the module to which it belongs
    # We use these to filter genes that we submit to STRINGdb 
    list_list_pkMEs <- clusterMap(cl, function(x,y) parPkMEs(list_kMEs=x, list_colors=y), 
                                  x=list_list_kMEs, 
                                  y=list_list_colors_ok, 
                                  SIMPLIFY = F)
    stopCluster(cl)
    invisible(gc())
    
    list_list_pkMEs <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y[[1]], dimension = NULL), x=list_list_pkMEs, y=list_list_plot_label_ok, SIMPLIFY=F)
    names(list_list_pkMEs) <- sNames_ok
  }
  ######################################################################
  ############## Match colors between parameter settings ###############
  ######################################################################
  # in:
  #   list_list_colors_ok
  #
  # out:
  #   list_list_colors_ok_matched

  message("Matching module color labels between parameter settings")

  # Align / match the colors so similar modules found with different sets of parameters have the same name
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_parMatchColors.txt"))
  list_list_colors_ok_matched <- parLapply(cl, list_list_colors_ok, parMatchColors)
  stopCluster(cl)
  invisible(gc())
  
  # Rename highest level list entries (cell clusters)
  names(list_list_colors_ok_matched) = sNames_ok
  
  # Name each entry of the vectors of color assignments with the corresponding genes
  list_list_colors_ok_matched <- mapply(function(x,y) lapply(x, function(z) name_for_vec(to_be_named=z, given_names=colnames(y), dimension = NULL)), x=list_list_colors_ok_matched, y=list_datExpr_ok, SIMPLIFY=F)

  # clear up 
  rm(list_list_colors_ok)
  invisible(gc())

  resume="checkpoint_2"
  
  #save.image(file=sprintf("%s%s_checkpoint_2_image.RData", RObjects_dir, data_prefix))
  if (autosave==T) save(list = ls(all.names = TRUE)[-grep("data_prefix$|corFnc$|networkType$|anti_cor_action$|hclustMethod$|minClusterSize$|deepSplit$|moduleMergeCutHeight$|pamStage$|TOMnReplicate$|checkPPI$|organism$|magma_gwas_dir$|gwas_filter_traits$|n_cores$|resume$" , ls(), ignore.case=F) ], envir = .GlobalEnv, file=sprintf("%s%s_checkpoint_2_image.RData", RObjects_dir, data_prefix))
  if (!is.null(quit_session)) if (quit_session=="checkpoint_2") quit(save="no")
  
} else if (!is.null(resume)) {
  if (resume == "checkpoint_2") {
    load((file=sprintf("%s%s_checkpoint_2_image.RData", RObjects_dir, data_prefix)))
    source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
    source(file = "/projects/jonatan/functions-src/functions.R")
    options(stringsAsFactors = F)
  }
}

if (resume == "checkpoint_2") {
  
  ######################################################################
  #################### (UN) LOAD PACKAGES #############################
  ######################################################################

  # Free up DLLs
  invisible(R.utils::gcDLLs())

  if (checkPPI==T)  suppressPackageStartupMessages(library(STRINGdb))

  ######################################################################
  #################### CHECK MODULES FOR PPI ENRICHMENT ################
  ######################################################################
  # in:
  #   list_list_colors_ok_matched
  #   list_list_pkMEs
  #
  # out:
  #   list_list_colors_ok_matched_PPI
  if (checkPPI == T) {
    message("Checking modules for significant Protein-Protein Interactions through STRINGdb")
  
    if (organism == "hsapiens") STRINGdb_species <- 9606 else if (organism == "mmusculus") STRINGdb_species <- 10090
    
    invisible(gc())
    
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_PPI_outer_for_vec.txt"))
    
    list_list_colors_PPI <- clusterMap(cl, function(a,b) mapply(function(x,y) PPI_outer_for_vec(colors = x,
                                                                                                pkMEs = y,
                                                                                                STRINGdb_species = STRINGdb_species, 
                                                                                                PPI_pkME_threshold = PPI_pkME_threshold,
                                                                                                p.val.threshold = p.val.threshold, 
                                                                                                project_dir = project_dir, 
                                                                                                data_prefix = data_prefix, 
                                                                                                flag_date = flag_date), 
                                                                x = a, y = b, SIMPLIFY=F), 
                                       a = list_list_colors_ok_matched, 
                                       b = list_list_pkMEs, 
                                       SIMPLIFY=F)
    stopCluster(cl)
    invisible(gc())
    
  } else if (checkPPI==F) list_list_colors_PPI <- list_list_colors_ok_matched
  
  ######################################################################
  ########### ORDER PARAMETER SETS BY PPI ENRICHMENT AND PLOT ##########
  ######################################################################
  # in:
  #   list_list_colors_ok_matched
  #   list_list_colors_ok_matched_PPI
  #   list_list_plot_label_ok
  #   
  # out:
  #   list_list_colors_ok_matched_PPI_order
  #   list_list_plot_label_ok_order
  #   list_list_colors_ok_matched_order
  #   list_list_colors_PPI_order

  message("Selecting parameters with the highest number of genes assigned to modules significantly enriched for Protein-Protein Interactions")
    
  # Count how many genes were assigned to grey under each parametrisation and order the parametrisations
  list_PPI_vec_n_grey <- lapply(list_list_colors_PPI, count_grey_in_list_of_vec)

  # Order all the outputs by how many genes were assigned to a (non-grey) module
  list_list_plot_label_ok_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_plot_label_ok, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  list_list_colors_ok_matched_order <- mapply(function(x,y) x[order(y, decreasing=F)], x  =  list_list_colors_ok_matched , y = list_PPI_vec_n_grey, SIMPLIFY = F )
  list_list_colors_PPI_order <- mapply(function(x,y) x[order(y, decreasing=F)], x = list_list_colors_PPI, y = list_PPI_vec_n_grey, SIMPLIFY=F)
  
  ######################################################################
  ##### FOR EACH SUBSET SELECT PARAMETERS WITH BEST PPI ENRICHMENT #####
  ######################################################################
  
  # Eliminate a layer of nesting by selecting only the best parametrisation per celltype
  list_plot_label_final <- lapply(list_list_plot_label_ok_order, function(x) x[[1]])
  list_colors_PPI <- lapply(list_list_colors_PPI_order, function(x) x[[1]])
  list_colors <- lapply(list_list_colors_ok_matched_order, function(x) x[[1]])

  # Name by cell clusters
  names(list_plot_label_final) <- sNames_ok
  names(list_colors_PPI) <- sNames_ok
  names(list_colors) <- sNames_ok

  # Make list of list of final parameters
  param_names = c("minClusterSize", "deepSplit","pamStage", "moduleMergeCutHeight")
  list_list_cutree_params_final <- lapply(list_PPI_vec_n_grey, function(x) comb_list[order(x,decreasing = F)][[1]])
  list_list_cutree_params_final <- lapply(list_list_cutree_params_final, function(x) name_for_vec(to_be_named = x, given_names = param_names, dimension = NULL)) 

  ######################################################################
  ################### DELETE RUNS WITH ONLY GREY #######################
  ######################################################################
  # in:
  #   list_colors_PPI
  #   list_datExpr_ok
  #   sNames_ok
  #
  # out:
  #   list_colors_PPI_ok
  #   list_datExpr_PPI_ok
  #   sNames_PPI_ok
  
  PPI_vec_n_grey <- count_grey_in_list_of_vec(list_colors_PPI)
  
  # For any parameter setting, does the number of genes assigned to grey correspond to the length of the vector of assignments? If not, it's ok.
  logical_subsets_PPI_ok <- as.logical(mapply(function(x,y) x != length(y), x = PPI_vec_n_grey, y = list_colors_PPI, SIMPLIFY = F))
  
  sNames_PPI <- sNames_ok[logical_subsets_PPI_ok]
  list_colors_PPI <- list_colors_PPI[logical_subsets_PPI_ok]
  list_datExpr_PPI <- list_datExpr_ok[logical_subsets_PPI_ok]
  #list_dissTOM_PPI <- list_dissTOM_ok[logical_subsets_PPI_ok]
  
  list_colors <- list_colors[logical_subsets_PPI_ok]
  list_list_cutree_params_final_PPI <- list_list_cutree_params_final[logical_subsets_PPI_ok]
  list_plot_label_final_PPI  <- list_plot_label_final[logical_subsets_PPI_ok]
  list_geneTree_PPI <- list_geneTree_ok[logical_subsets_PPI_ok]
  
  # Checkpoint
  resume = "checkpoint_3"
  message("Reached checkpoint 3, saving session image")
  #save.image(file=sprintf("%s%s_checkpoint_3_image.RData", RObjects_dir, data_prefix))
  if (autosave==T) save(list = ls(all.names = TRUE)[-grep("data_prefix$|corFnc$|networkType$|anti_cor_action$|hclustMethod$|minClusterSize$|deepSplit$|moduleMergeCutHeight$|pamStage$|TOMnReplicate$|checkPPI$|organism$|magma_gwas_dir$|gwas_filter_traits$|n_cores$|resume$" , ls(), ignore.case=F) ], envir = .GlobalEnv, file=sprintf("%s%s_checkpoint_3_image.RData", RObjects_dir, data_prefix))
  if (!is.null(quit_session)) if (quit_session=="checkpoint_3") quit(save="no")
  
} else if (resume == "checkpoint_3") {
  load(file=sprintf("%s%s_checkpoint_3_image.RData", RObjects_dir, data_prefix))
  # Load parameter values and utility functions anew (in case a bug was fixed)
  source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
  source(file = "/projects/jonatan/functions-src/functions.R")
  options(stringsAsFactors = F)
}  

if (resume == "checkpoint_3") {
  
  ##########################################################################
  ############################ (UN)LOAD PACKAGES ############################
  ##########################################################################
  
  # Free up DLLs
  try(detach("package:STRINGdb", unload=TRUE))
  invisible(R.utils::gcDLLs())
  # Load packages
  suppressPackageStartupMessages(library(boot))
  
  ######################################################################
  ##################### MAKE MODULE COLORS UNIQUE ######################
  ######################################################################
  
  message("Making module colors unique across cell clusters")
  
  # Get nested list of modules
  list_mods <- lapply(list_colors_PPI, function(x) names(table(x))[-grep("^grey$", names(table(x)))])
  mods <- unlist(list_mods)
  names(mods) <- NULL
  
  all_cols_nogrey_uniq <- unique(gsub("\\d+", "", colors()[-grep("grey|gray", colors())]))
  all_cols_nogrey <- colors()[-grep("grey|gray", colors())]
  
  # Replace colors with new unique colors
  if (length(mods) <= length(all_cols_nogrey_uniq) ) { # if there are enough unique colors without adding numbers
    mods_uniq <- all_cols_nogrey_uniq[sample(x=1:length(all_cols_nogrey_uniq), size=length(mods), replace=F)]
  } else if (length(mods) > length(all_cols_nogrey_uniq) & length(mods) < length(all_cols_nogrey) ) { # if there aren't enough unique colors unless they have numbers added
    mods_uniq <- all_cols_nogrey[sample(x=1:length(all_cols_nogrey), size=length(mods), replace=F)]
  } else if (length(mods) > length(all_cols_nogrey)) { # if there aren't enough unique colors in R
    fakeCols <- paste0(all_cols_nogrey_uniq, "_", 1:(length(mods) - length(all_cols_nogrey)))
    mods_uniq <- mods
    mods_uniq[1:length(all_cols_nogrey)] <- all_cols_nogrey[sample(x=1:length(all_cols_nogrey), size=length(all_cols_nogrey), replace=F)]
    mods_uniq[(length(all_cols_nogrey)+1):length(mods_uniq)] <- fakeCols[sample(x=1:length(fakeCols), size=length(mods_uniq)-length(all_cols_nogrey), replace=F)]
  }

  # Nest the modules by celltype
  k = 0
  list_mods_uniq <- vector(mode="list", length=length(list_mods))
  for (i in 1:length(list_mods)) {
    list_mods_uniq[[i]] <- mods_uniq[(k+1):(k+length(list_mods[[i]]))]
    k <- k+length(list_mods[[i]])
  }
    
  names(list_mods_uniq) <- names(list_colors_PPI)
    
  # Update color vectors
  list_colors_PPI_uniq <- mapply(function(x,y,z) z[match(x, y)],
                             x = list_colors_PPI, 
                             y = list_mods, 
                             z = list_mods_uniq, SIMPLIFY = F)
  
  # Update also the old module assignment colors, but only change those that were validated by PPI (and therefore among those replaced with unique colors)
  list_colors_uniq <-  mapply(function(x,y,z) ifelse(x==y, yes=z, no = x),
                          x = list_colors, 
                          y = list_colors_PPI, 
                          z = list_colors_PPI_uniq,
                          SIMPLIFY = F)

  # Replace NAs (i.e. no match in non-grey modules) with "grey"
  list_colors_uniq <- lapply(list_colors_uniq, function(x) ifelse(is.na(x), yes="grey", no = x))
  list_colors_PPI_uniq <- lapply(list_colors_PPI_uniq, function(x) ifelse(is.na(x), yes="grey", no = x))
  
  # Give gene names to vectors
  list_colors_uniq <- mapply(function(x,y) name_for_vec(to_be_named= x, given_names = names(y), dimension=NULL), x = list_colors_uniq, y = list_colors, SIMPLIFY=F)
  list_colors_PPI_uniq <- mapply(function(x,y) name_for_vec(to_be_named= x, given_names = names(y), dimension=NULL), x = list_colors_PPI_uniq, y = list_colors_PPI, SIMPLIFY=F)
    
  ######################################################################
  ########## RECOMPUTE MEs, kMEs, kIMs, pkMEs AFTER PPI FILTER #########
  ######################################################################
  # in:
  #   list_colors_PPI_ok
  #   list_datExpr_PPI_ok
  #   list_dissTOM_PPI_ok
  #
  # out:
  #   list_MEs_PPI
  #   list_kMEs_PPI
  
  message("Computing kMEs after filtering modules for significant Protein-Protein Interactions")
  
  invisible(gc())
  cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_moduleEigengenes_kMEs_pkMEs_PPI.txt"))

  
  list_MEs_PPI <- clusterMap(cl, function(x,y) moduleEigengenes(expr = as.data.frame(x, col.names=col.names(x)),
                                                                      colors=y,
                                                                      excludeGrey=T), 
                                   x = list_datExpr_PPI, 
                                   y = list_colors_PPI_uniq, 
                                   SIMPLIFY = F)

  list_kMEs_PPI <- clusterMap(cl, function(x,y) signedKME(as.matrix(x),
                      y$eigengenes,
                      outputColumnName = "",
                      corFnc = corFnc), 
                      x=list_datExpr_PPI, 
                      y=list_MEs_PPI)
 
  stopCluster(cl)
  invisible(gc())
  
  # not needed
  # for (j in 1:length(list_kMEs_PPI)) {
  #   colnames(list_kMEs_PPI[[j]]) <- gsub("^kME", "", colnames(list_kMEs_PPI[[j]]), ignore.case = F)
  # }
  
  # Name rows in each dataframe - no, not necessary
  #list_kMEs_PPI <- mapply(function(x,y) name_for_vec(to_be_named=x, given_names = names(y), dimension = 1), x=list_kMEs_PPI, y=list_colors_PPI_uniq, SIMPLIFY=F)
  
  ###180509_v1.8_dev2
  # There should be no need since we have excludeGrey = T. 
  # But it may be necessary because there are definitely still grey modules
  ### 180531 we need MEs to make eigen matrix later
  #list_MEs_PPI <- deleteGrey(list_MEs_PPI) %>% Filter(f=length)
  #list_kMEs_PPI <- deleteGrey(list_kMEs_PPI) %>% Filter(f=length)

  # Delete consensus TOMs from disk
  for (subsetName in sNames) {
    if (file.exists(sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, subsetName))) {
      file.remove(sprintf("%s%s_consensusTOM-block.1.RData", RObjects_dir, subsetName))
    }
  }

  ##########################################################################
  ################## COMPUTE RARE VARIANTS ENRICHMENT ######################
  ##########################################################################

  message("Checking modules for enrichment with rare variants/mendelian genes")
  
  coding_variants <- load_obj(f=coding_variants_path)
  mendelian_genes <- load_obj(f=mendelian_genes_path)
  
  # 'Looping' over celltypes, compute logistic model for each module
  cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "mendelian_rare_variants_par.txt"))
  
  rare_variants_results <- clusterMap(cl, function(x,y) mendelianGenes(cellType = x,
                                                                       colors = y,
                                                                       coding_variants = coding_variants,
                                                                       mendelian_genes = mendelian_genes), 
                                  x = names(list_colors_PPI_uniq),
                                  y = list_colors_PPI_uniq,
                                  SIMPLIFY=F)
  stopCluster(cl)
  invisible(gc())
  # Put odd-ratios and p-values into separate celltype-length lists
  list_variants.or <- list_variants.p <- vector(mode = "list", length = length(rare_variants_results))
  
  for (i in 1:length(rare_variants_results)) { 
    list_variants.or[[i]] <- rare_variants_results[[i]][['or']]
    list_variants.p[[i]] <- rare_variants_results[[i]][['p.val']]
  }
  
  # Merge each list into a single table
  for (t in 1:length(list_variants.or)) {
    if (t == 1) {
      variants.or.all <- list_variants.or[[t]]
      variants.p.all <- list_variants.p[[t]]
    } else {
      variants.or.all <- rbind(variants.or.all, list_variants.or[[t]]) 
      variants.p.all <- rbind(variants.p.all, list_variants.p[[t]]) 
    }
  }
  
  # adjust analytical p-values for multiple testing
  variants.p.fdr.all = p.adjust(variants.p.all, method="fdr")
  dim(variants.p.fdr.all) = dim(variants.p.all);  dimnames(variants.p.fdr.all) = dimnames(variants.p.all)
  
  # Transform fdr-corrected p-values to -log10(p.val)
  variants.p.fdr.log <- -log10(variants.p.fdr.all)
  
  # Convert everything to dataframe
  variants.or.all <- variants.or.all %>% as.data.frame 
  variants.p.all <- variants.p.all  %>% as.data.frame 
  variants.p.fdr.log  <- variants.p.fdr.log  %>% as.data.frame 
  
  # Add module and celltype as columns
  module <- gsub(".+__","",rownames(variants.p.fdr.log))
  celltype <- gsub("__.+","",rownames(variants.p.fdr.log))
  
  variants.p.fdr.log <- cbind(celltype, module, variants.p.fdr.log) 
  variants.p.all <- cbind(celltype, module,variants.p.all)
  variants.or.all <- cbind(celltype,module,variants.or.all)
  
  ##########################################################################
  ##############COMPUTE MAGMA COMMON VARIANTS GWAS ENRICHMENT ##############
  ########  ##################################################################
  
  # MAGMA docs: http://ctg.cncr.nl/software/MAGMA/doc/manual_v1.06.pdf
  
  if (!is.null(magma_gwas_dir)) {
    
    message("Scoring PPI enriched modules for enrichment with genes linked by GWAS to phenotypes of interest")
      
    if (organism == "mmusculus") {
      
      # map mouse to human gene orthologs  
      mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
      
      cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "mapMMtoHs_par.txt"))
      list_kMEs_hs <- parLapply(cl, list_kMEs_PPI, function(x) mapMMtoHs(modulekME = x, 
                                                                    log_dir = log_dir, 
                                                                    flag_date = flag_date, 
                                                                    data_prefix = data_prefix, 
                                                                    mapping_orthology = mapping_orthology))
      stopCluster(cl)
      invisible(gc())
     
      names(list_kMEs_hs) <- names(list_kMEs_PPI)
      #rm(list_kMEs_PPI)
      
    } else if (organism == "hsapiens") {
      list_kMEs_hs <- list_kMEs_PPI
      names(list_kMEs_hs) <- names(list_kMEs_PPI)
      #rm(list_kMEs_PPI)
    }
    
    # Load gene mapping and annotation files
    
    # columns: ensembl_gene_id, entrezgene, hgnc_symbol. Use this for mapping entrezgene to ensembl
    mapping_hs_entrez2ensembl = read.csv(gzfile(mapping_hs_entrez2ensembl_filepath),sep="\t",header=T, stringsAsFactors = F)
    
    # Load MAGMA genes and remap to Ensembl gene IDs
    
    d = dir(path=magma_gwas_dir, pattern="[.]genes.out", recursive = T)
    gwas = vector(mode="list")
    for(i in 1:length(d)) {
      gwas[[i]] = read.table(paste(magma_gwas_dir, d[[i]],sep=""),head=T, check.names = FALSE, stringsAsFactors = F)
    }
    names(gwas) = gsub(".genes.out", "", d)
    
    # Remap from human Entrez to human Ensembl gene IDs
    for(i in 1:length(gwas)) {
      idx = match(gwas[[i]]$GENE, mapping_hs_entrez2ensembl$entrezgene)
      mapping = data.frame(entrez=gwas[[i]]$GENE, ensembl=mapping_hs_entrez2ensembl$ensembl_gene_id[idx])
      gwas[[i]]$gene_name = mapping$ensembl
    }

    cl <- makeCluster(n_cores, type="FORK", outfile = paste0(log_dir, "kME_magma_par.txt"))
    magma_results <- clusterMap(cl, function(x,y) kME_magma(cellType = x, 
                                                            modulekME = y,
                                                            gwas = gwas),
                                x = names(list_kMEs_hs),
                                y = list_kMEs_hs,
                                SIMPLIFY=F)
    stopCluster(cl)
    invisible(gc())
    
    # Prepare tables for results
    list_magma.r <- list_magma.p <- list_magma.emp.p <- vector(mode="list", length = length(list_kMEs_hs))
    
    # Extract the coefficient, analytical p-value and empirical p-value dataframes for each module, put them into lists
    
    for (i in 1:length(magma_results)) { 
      list_magma.r[[i]] <- magma_results[[i]][['corrCoef']]
      list_magma.p[[i]] <- magma_results[[i]][['p.val']]
      list_magma.emp.p[[i]] <- magma_results[[i]][['emp.p.val']]
    }
    
    # Merge each list into a single table
    for (t in 1:length(list_magma.p)) {
      if (t == 1) {
        magma.r.all <- list_magma.r[[t]]
        magma.p.all <- list_magma.p[[t]]
        magma.emp.p.all <- list_magma.emp.p[[t]]
      } else {
        magma.r.all <- rbind(magma.r.all, list_magma.r[[t]] )
        magma.p.all <- rbind(magma.p.all, list_magma.p[[t]] )
        magma.emp.p.all <- rbind(magma.emp.p.all, list_magma.emp.p[[t]] )
      }
    }
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
    
    # Add module and celltype as columns
    #magma.p.fdr.log$module <- magma.p.all$module <- magma.emp.p.all$module <- magma.r.all$module <- gsub("kME|kME\\.","",rownames(magma.p.fdr.log))
    
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
  }
  
  ##########################################################################
  ############## FILTER MODULES ON GENES WITH GWAS ENRICHMENT #############
  ##########################################################################
  
  # TODO: resume filtering once we've sorted MAGMA
  
  if (FALSE & !is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) {
    
    # If filtering out non-gwas enriched modules, also retain modules with rare variant / mendelian gene enrichment 
    idx_row_rareVariants <- rowMin(as.matrix(variants.p.fdr.all)) < p.val.threshold
    
    # Check if the gwas_filter_traits strings provided by the user can be found in the magma_gwas_dir files
    
    if (sapply(gwas_filter_traits, function(x) any(grepl(x, colnames(magma.p.fdr.log), ignore.case=T)), simplify = T) %>% all) {
    #if (sapply(gwas_filter_traits, function(x) any(grepl(x, colnames(magma.p.fdr.log[,-grep("module", colnames(magma.p.fdr.log))]), ignore.case=T)), simplify = T) %>% all) {
      
      # Identify columns to keep
      idx_col_keep <- sapply(gwas_filter_traits, function(x) grep(paste0(x,"|module|celltype"), colnames(magma.p.fdr.log), ignore.case=T), simplify = T) %>% Reduce(union,.)# %>% as.logical 
        
      # fdr-corrected log-transformed analytical p-values
      magma.p.fdr.log.sig <- magma.p.fdr.log[,idx_col_keep]
     
      # correlation coefficients
      magma.r.sig <- magma.r.all[,idx_col_keep]
      
      # Identify rows to keep
      idx_row_gwas <- apply(magma.p.fdr.log.sig[,-grep("module|celltype", colnames(magma.p.fdr.log.sig))], MARGIN = 1, max) > -log10(p.val.threshold)
      
      # rows enriched for gwas common variants and/or for rare variants
      idx_row_keep <- idx_row_gwas | idx_row_rareVariants
        
      if (sum(idx_row_keep)>0) {

        magma.p.fdr.log.sig <- magma.p.fdr.log.sig[idx_row_keep,]
        magma.r.sig <- magma.r.sig[idx_row_keep,]
        
      } else if (sum(idx_row_keep)==0) { # Warn if there are no significant enrichments
        
        gwas_filter_traits <- NULL
        warning("No modules enriched for gwas_filter_traits or rare variants - skipping gwas significance filtering step")
        
      }
    } else {
      
      gwas_filter_traits <- NULL 
      message("gwas_filter_traits not found in magma_gwas_dir files - skipping gwas significance filtering step")
      
    }
  }  
  
  if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
    magma.p.fdr.log.sig <- magma.p.fdr.log
    magma.r.sig <- magma.r.all
  }
  
  ######################################################################
  ####### FILTER MODULES ON GWAS SIGNIFICANCE ALSO IN OTHER FILES ######
  ######################################################################
    
    if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) { # NB: if no significant gwas correlations found in previous section, gwas_filter_traits <- NULL
      
      # Also filter out non-MAGMA/rare variant enriched modules in rare variants results (although we can still keep and plot the full results)
      # variants.or.all_gwas <- variants.or.all[idx_row_keep,] 
      # variants.p.fdr.log_gwas <- variants.p.fdr.log[idx_row_keep,]
  
      # get cell types with modules with gwas enrichment
      sNames_gwas <- names(table(magma.p.fdr.log.sig$celltype))
      
      list_module_gwas <- vector(mode="list", length=length(sNames_gwas))
      names(list_module_gwas) <- sNames_gwas
      
      # Get a list of vectors with modules per celltype
      for (name in sNames_gwas) {
        idx <- grep(pattern = name, x = magma.p.fdr.log.sig$celltype)
        list_module_gwas[[name]] <- sapply(idx, function(j) magma.p.fdr.log.sig$module[j], simplify = T)
      }
      
      # Filter out cell types with no significant modules
      list_colors_gwas <- list_colors_PPI_uniq[sNames_PPI  %in% sNames_gwas]
        
      # relabel non-GWAS modules 'grey'
      list_colors_gwas <- mapply(function(x,y) ifelse(x %in% y, yes = x, no = "grey"),
                                 ### 180611
                                 # x = list_colors_gwas,
                                 x = list_colors_gwas, 
                                 y = list_module_gwas,
                                 SIMPLIFY = F)
      
      # give gene names to color assignment vectors
      list_colors_gwas <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(y), dimension = NULL), 
                                 x = list_colors_gwas,
                                 y = list_colors_PPI_uniq[sNames_PPI %in% sNames_gwas, drop=F],
                                 SIMPLIFY = F)
      
      # remove any cell clusters without gwas enrichment: expression matrices
      list_datExpr_gwas <- list_datExpr_PPI[sNames_PPI %in% sNames_gwas, drop=F] 
      
    } else if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
      
      sNames_gwas <- sNames_PPI
      list_module_gwas <- lapply(list_kMEs_PPI, function(x) colnames(x)) 
      names(list_module_gwas) <- sNames_gwas
      list_colors_gwas <- list_colors_PPI_uniq
      list_datExpr_gwas <- list_datExpr_PPI
    }
    
  ######################################################################
  #################### RECOMPUTE MEs, kMEs AFTER MAGMA #################
  ######################################################################
  
  if (!is.null(magma_gwas_dir) & !is.null(gwas_filter_traits)) {
  
    message("Re-computing MEs and kME after filtering modules for MAGMA enrichment")
    
    invisible(gc())
    cl <- makeCluster(n_cores, type = "FORK", outfile = paste0(log_dir, "log_moduleEigengenes_kMEs_magma.txt"))
    
    list_MEs_gwas <- clusterMap(cl, function(x,y) moduleEigengenes(expr=as.data.frame(x,col.names=col.names(x)),
                                                                  colors=y,
                                                                  excludeGrey=T), 
                               x = list_datExpr_gwas, 
                               y = list_colors_gwas, 
                               SIMPLIFY = F)
  
    list_kMEs_gwas <- clusterMap(cl, function(x,y) signedKME(as.matrix(x),
                                                            y$eigengenes,
                                                            outputColumnName = "",
                                                            corFnc = corFnc), 
                                x=list_datExpr_gwas, 
                                y=list_MEs_gwas)
    
    ### 180524 v1.8_dev4
    # Temporarily removed. In AgRP run 21, list_kIMs_PPI[[1]] has an extra column which has name NA
    # list_kIMs_PPI <- clusterMap(cl, function(x,y) kIM_eachMod_norm(dissTOM=x, colors=y, genes=names(y)),
    #                             x = list_dissTOM_PPI_ok,
    #                             y = list_colors_PPI_ok,
    #                             SIMPLIFY=F)
    ###
    
    stopCluster(cl)
    invisible(gc())
    
    # Remove grey modules (if any and filter out runs with only grey (probably not possible but still))
    list_MEs_gwas <- deleteGrey(list_MEs_gwas) %>% Filter(f=length)
    list_kMEs_gwas<- deleteGrey(list_kMEs_gwas) %>% Filter(f=length)
    
  } else if (is.null(magma_gwas_dir) | is.null(gwas_filter_traits)) {
    list_MEs_gwas <- list_MEs_PPI
    list_kMEs_gwas <- list_kMEs_PPI
  }
  
  resume = "checkpoint_4"
  message("Reached checkpoint 4, saving session image")
  #save.image(file=sprintf("%s%s_checkpoint_4_image.RData", RObjects_dir, data_prefix))
  if (autosave==T) save(list = ls(all.names = TRUE)[-grep("data_data_prefix$|corFnc$|networkType$|anti_cor_action$|hclustMethod$|minClusterSize$|deepSplit$|moduleMergeCutHeight$|pamStage$|TOMnReplicate$|checkPPI$|organism$|magma_gwas_dir$|gwas_filter_traits$|n_cores$|resume$" , ls(), ignore.case=F) ], envir = .GlobalEnv, file=sprintf("%s%s_checkpoint_4_image.RData", RObjects_dir, data_prefix))
  if (!is.null(quit_session)) if (quit_session=="checkpoint_4") quit(save="no")
  
} else if (resume == "checkpoint_4") {
  load(file=sprintf("%s%s_checkpoint_4_image.RData", RObjects_dir, data_prefix))
  # Load parameter values and utility functions anew (in case a bug was fixed)
  source(file = "/projects/jonatan/wgcna-src/rwgcna-pipeline/rwgcna_params.R")
  source(file = "/projects/jonatan/functions-src/functions.R")
  options(stringsAsFactors = F)
}

if (resume == "checkpoint_4") {

  ######################################################################
  #### COMPUTE EIGENGENE - METADATA CORRELATION IN EACH CELL CLUSTER ###
  ######################################################################
  
  if (!is.null(metadata_corr_col)) if (!is.null(metadata)) {
    
    list_metadata <- lapply(list_datExpr_gwas, function(x) metadata[match(rownames(x), rownames(metadata)), , drop=F]) # get list of cell * metadata
    
    # Compute correlation between metadata (columns) and eigengenes (columns). 
    list_eigen_metadata_corr_rho <- mapply(function(x,y) WGCNA::cor(x=x, 
                                                                 y=y$eigengenes, 
                                                                 method = c("pearson"), 
                                                                 verbose = 0, 
                                                                 #use = 'all.obs'),
                                                                 use = 'pairwise.complete.obs'), 
                                        x=list_metadata, 
                                        y=list_MEs_gwas,  
                                        SIMPLIFY = F) 
      
    list_eigen_metadata_corr_rho <- lapply(list_eigen_metadata_corr_rho, function(x) name_for_vec(to_be_named = x, given_names = colnames(metadata), dimension = 1)) 
      
    # Remove 'ME' from eigengene names
    for (j in 1:length(list_eigen_metadata_corr_rho)) {
      colnames(list_eigen_metadata_corr_rho[[j]]) <- gsub("^ME", "", colnames(list_eigen_metadata_corr_rho[[j]]), ignore.case=F)
    }
    
    # Compute p values
    list_eigen_metadata_corr_pval <- mapply(function(x,y) WGCNA::corPvalueStudent(x, 
                                                                                   n = nrow(y)), 
                                             x=list_eigen_metadata_corr_rho, 
                                             y=list_metadata, 
                                             SIMPLIFY=F)
    
    # Convert p values into one big matrix in order to adjust p-values for the number of modules tested across *all* celltypes
    for (j in 1:length(list_eigen_metadata_corr_pval)) {
      if (j==1) {
        corr_pval <- as.data.frame(list_eigen_metadata_corr_pval[[j]])
        colnames(corr_pval) <- paste0(names(list_eigen_metadata_corr_pval)[j], "__", colnames(corr_pval)) 
      } else {
        tmp <- list_eigen_metadata_corr_pval[[j]]
        colnames(tmp) <- paste0(names(list_eigen_metadata_corr_pval)[j], "__", colnames(tmp)) 
        corr_pval <- cbind(corr_pval, tmp)
      }
    }
    
    # set NAs to 1
    corr_pval[is.na(corr_pval)] <- 1
      
    # Compute the false discovery rates
    corr_fdr <- t(apply(t(corr_pval), MARGIN=2, FUN = function(x) p.adjust(x, method = "fdr")))
    corr_fdr.log <- -log10(corr_fdr)     
    
    # Split single dataframe back into a list of celltypes
    list_eigen_metadata_corr_fdr <- list_eigen_metadata_corr_fdr.log <- vector(mode = "list", length = length(list_eigen_metadata_corr_pval))
    names(list_eigen_metadata_corr_fdr) <- names(list_eigen_metadata_corr_fdr.log) <- names(list_eigen_metadata_corr_rho)
    
    k = 0
    for (j in 1:length(list_eigen_metadata_corr_fdr.log)) {
      list_eigen_metadata_corr_fdr[[j]] <- as.data.frame(corr_fdr[,(k+1):(k+ncol(list_eigen_metadata_corr_pval[[j]])), drop=F])
      list_eigen_metadata_corr_fdr.log[[j]] <- as.data.frame(corr_fdr.log[,(k+1):(k+ncol(list_eigen_metadata_corr_pval[[j]])), drop=F])
      k <- k + ncol(list_eigen_metadata_corr_pval[[j]])
    }
    
    # Also make a single dataframe for saving 
    dim(corr_fdr.log) <- c(ncol(metadata), length(unlist(list_module_gwas)))
    rownames(corr_fdr.log) <- rownames(corr_pval)
    colnames(corr_fdr.log) <- colnames(corr_pval)
    corr_fdr.log <- as.data.frame(corr_fdr.log)
  } 
  
  ##########################################################################
  #### FILTER COLORS VECS, GENE AND KME LISTS FOR METADATA CORRELATIONS ####
  ##########################################################################
  
  if (!is.null(metadata_corr_col) & !is.null(metadata_corr_filter_vals)) {
    
    if (!is.null(metadata)) {
      # Get a list of logical vectors indicating significant correlations
      list_idx_module_sig <- lapply(list_eigen_metadata_corr_fdr.log, function(x) apply(x, 2, function(y) any(y>-log10(p.val.threshold)))) # logical
     
      # Only keep significantly correlated modules within each celltype
      list_module_meta <- mapply(function(x,y) colnames(x[,y, drop=F]),x=list_eigen_metadata_corr_fdr.log, y=list_idx_module_sig, SIMPLIFY=F) 
      list_module_meta <- lapply(list_module_meta, function(x) gsub("^.*__", "", x)) 
  
      sNames_meta <- sNames_gwas[sapply(list_idx_module_sig, any, simplify = T)]
      
      # Filter out celltypes with no modules significantly correlated with metadata
      list_module_meta <- list_module_meta[sNames_gwas %in% sNames_meta]
  
      # reassign genes of filtered out modules to grey and remove any empty cell clusters
      list_colors_meta <- mapply(function(x,y) ifelse(x %in% y, yes = x, no = "grey"),
                                   x = list_colors_gwas[names(list_colors_gwas) %in% sNames_meta, drop=F],
                                 y = list_module_meta,
                                 SIMPLIFY = F)
      
      # give gene names to color assignment vectors
      list_colors_meta <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = names(y), dimension = NULL), 
                                 x = list_colors_meta,
                                 y = list_colors_gwas[names(list_colors_gwas) %in% sNames_meta, drop=F],
                                 SIMPLIFY = F)
      
      # Filter kME lists
      list_kMEs_meta <- list_kMEs_gwas[names(list_kMEs_gwas) %in% sNames_meta, drop = F]
      list_kMEs_meta <- mapply(function(x,y) x[,colnames(x)%in%y],x=list_kMEs_meta, y = list_module_meta) %>% Filter(f=length)
  
      # if (!is.null(magma_gwas_dir)) {
      #   idx <- unlist(list_module_meta) %in% unlist(list_module_gwas)
      #   magma.p.fdr.log.sig_meta <- magma.p.fdr.log.sig[idx,]
      #   magma.r.sig_meta <- magma.r.sig[idx,]
      # }
    } else 
    
    # TODO: remove, just a backstop
    sNames_meta <- sNames_gwas
    list_colors_meta <- list_colors_gwas
    list_module_meta <- list_module_gwas
    list_kMEs_meta <- list_kMEs_gwas 
    
  } else if (is.null(metadata_corr_col) | is.null(metadata_corr_filter_vals)) {
    
    sNames_meta <- sNames_gwas
    list_colors_meta <- list_colors_gwas
    list_module_meta <- list_module_gwas
    list_kMEs_meta <- list_kMEs_gwas 
    # if (!is.null(magma_gwas_dir)) {
    #   magma.p.fdr.log.sig_meta <- magma.p.fdr.log.sig
    #   magma.r.sig_meta <- magma.r.sig
    # }
  }
  
  ######################################################################
  ##################### PREPARE ORDERED GENE LISTS #####################
  ######################################################################
  
  # prepare nested lists of module genes
  list_list_module_meta_genes <- mapply(function(a,b) lapply(b, function(x) names(a)[a==x]), a=list_colors_meta, b=list_module_meta, SIMPLIFY=F)
  list_list_module_meta_genes <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = NULL), x=list_list_module_meta_genes, y=list_module_meta, SIMPLIFY=F)
  
  # order the gene lists by kME
  tmp <- list_list_module_meta_genes
  
  # iterate over celltypes
  for (i in 1:length(list_list_module_meta_genes)) {
    # iterate over modules
    for (j in 1:length(list_list_module_meta_genes[[i]])) {
      kMEs <- list_kMEs_meta[[i]] 
      genes <- list_list_module_meta_genes[[i]][[j]]
      tmp[[i]][[j]] <- genes[order(kMEs[rownames(kMEs) %in% genes, j], decreasing=T)]
    }
  }
  
  list_list_module_meta_genes <- tmp
  
  ######################################################################
  ##################### COMPUTE CELL x EIGENGENE MATRIX ################
  ######################################################################

  message("Computing all cell embeddings on all gene module eigengenes")

  eigen_mat <- make_eigen_mat(RObjects_dir = RObjects_dir,
                            list_list_module_genes = list_list_module_meta_genes,
                            n_cores = n_cores,
                            log_dir = log_dir)
  
  ##########################################################################
  ######### PREPARE GENES LISTS AND DATAFRAME WITH MODULES, GENES ##########
  ##########################################################################
  
  # Prepare module genes dataframe
  cell_cluster <- rep(sNames_meta, times=unlist(sapply(list_list_module_meta_genes, FUN=function(x) sum(sapply(x, function(y) length(y), simplify=T)), simplify=T)))
  module <- unlist(sapply(list_list_module_meta_genes, function(x) rep(names(x), sapply(x, function(y) length(y), simplify = T)), simplify=T), use.names = F)
  ensembl <- unlist(list_list_module_meta_genes, recursive = T, use.names = F)
  
  # if we mapped from hgnc to ensembl, retrieve the hgnc symbols
  if (file.exists(sprintf("%s%s_%s_hgnc_to_ensembl_mapping_df.RData", tables_dir, data_prefix, organism))) {
    mapping <- read.csv(file=sprintf("%s%s_%s_hgnc_to_ensembl_mapping_df.csv", tables_dir, data_prefix, organism), stringsAsFactors = F)
    list_list_module_meta_genes_hgnc <- lapply(list_list_module_meta_genes, function(x) lapply(x, function(y) mapping$symbol[match(y, mapping$ensembl)]))
    list_list_module_meta_genes_hgnc <- mapply(function(x,y) name_for_vec(to_be_named = x, given_names = y, dimension = NULL), x=list_list_module_meta_genes_hgnc, y=list_module_meta)
    hgnc <- unlist(list_list_module_meta_genes_hgnc, recursive = T, use.names=F)
    df_meta_module_genes <- data.frame(cell_cluster, module, ensembl, hgnc, row.names = NULL)
  } else {
    df_meta_module_genes <- data.frame(cell_cluster, module, ensembl, row.names = NULL)
  }
  
  ##########################################################################
  ############################# OUTPUT TABLES ##############################
  ##########################################################################
  
  invisible(write.csv(df_meta_module_genes, file = sprintf("%s%s_cell_cluster_module_genes_%s.csv",tables_dir, data_prefix, flag_date), row.names = F, quote = F))
  
  # also output ordered gene lists
  for (i in names(list_list_module_meta_genes)) {
   for (j in names(list_list_module_meta_genes[[i]])) {
     invisible(write.csv(list_list_module_meta_genes[[i]][[j]], file = sprintf("%s%s_%s_%s_module_genes.csv",tables_dir, data_prefix, i, j), row.names = F, quote = F))
   }
  }
  
  save(list_list_module_meta_genes, file=sprintf("%s%s_list_list_module_genes.RData", RObjects_dir, data_prefix))
  
  ################## SAVE KMEs FOR BMI-ENRICHED MODULES ####################
  ##########################################################################

  # Prepare dfs with a gene column followed by kMEs 
  list_kMEs_meta_out <- lapply(list_kMEs_meta, function(x) cbind(genes=rownames(x), x))
  rownames(list_kMEs_meta) <- NULL
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_kMEs_meta_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=F, quote = F), list_kMEs_meta_out, sNames_meta, SIMPLIFY = F))

  ############### OUTPUT MENDELIAN/RARE VARIANT RESULTS #####################
  ###########################################################################
  
  # Full set of results
  invisible(write.csv(variants.p.fdr.log, file=sprintf("%s%s_mendelian_rareVariants_fdr_log_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
  invisible(write.csv(variants.or.all, file=sprintf("%s%s_mendelian_rareVariants_OR_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))

  ######################## OUTPUT MAGMA RESULTS #############################
  ###########################################################################
  
  invisible(write.csv(magma.p.fdr.log.sig, file=sprintf("%s%s_magma.fdr.log.sig_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
  invisible(write.csv(magma.r.sig, file=sprintf("%s%s_magma.r.all_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
  invisible(write.csv(magma.p.all, file=sprintf("%s%s_magma.p.all_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
  invisible(write.csv(magma.emp.p.all, file=sprintf("%s%s_magma.emp.p.all_%s.csv", tables_dir, data_prefix, flag_date), row.names=F, quote = F))
  
  ####################### OUTPUT METADATA CORRELATION RESULTS ###############
  ###########################################################################
  
  invisible(write.csv(corr_fdr.log, file=sprintf("%s%s_all_metadata_corr_logfdr_%s.csv", tables_dir, data_prefix, flag_date), row.names=T, quote = F))
  
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_metadata_corr_rho_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=T, quote = F), list_eigen_metadata_corr_rho, sNames_meta, SIMPLIFY = F))
  invisible(mapply(function(x,y) write.csv(x, file=sprintf("%s%s_%s_metadata_corr_logfdr_%s.csv", tables_dir, data_prefix, y, flag_date), row.names=T, quote = F), list_eigen_metadata_corr_fdr.log, sNames_meta, SIMPLIFY = F))

  ############################ OUTPUT EIGEN MAT #############################
  ###########################################################################
  
  invisible(write.csv(eigen_mat, file=sprintf("%s%s_eigen_mat_%s.csv", tables_dir, data_prefix, flag_date), row.names=T, quote = F))

}
######################################################################
########################### CLEAR UP #################################
######################################################################

message("Saving final session image")

save.image(file=sprintf("%s%s_final_session_image.RData", RObjects_dir, data_prefix, flag_date))

rm(list=ls())
invisible(gc())

##########################################################################
################################ FINISH ##################################
##########################################################################

message("Script DONE!")
