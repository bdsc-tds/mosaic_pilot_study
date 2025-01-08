# -------------------------------------------------------------------------
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5/B3"

# noSpotClean + logNorm + BayesSpace
sample_name = "B3_2"
disease = "breast"
sce <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered.rds") # noSpotClean + logNorm + BayesSpace

source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Visium/01_params.R")
set.seed(100)
#sce <- runUMAP(sce, assay.type = "logcounts")
sce <- runUMAP(sce, assay.type = "logcounts", n_neighbors = 5) 
set.seed(100)
# sce <- runTSNE(sce, assay.type = "logcounts")
sce <- runUMAP(sce, assay.type = "logcounts", n_neighbors = 5) 

sce$spatial.cluster <- as.factor(sce$spatial.cluster)
plotDimRed(sce, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster") |
  plotDimRed(sce, type = "TSNE", annotate = "spatial.cluster", text_by = "spatial.cluster")

sce$pxl_col_in_fullres <- NULL; sce$pxl_row_in_fullres <- NULL

p <- plotSpots(sce, annotate = "spatial.cluster", text_by = "spatial.cluster", y_reverse = TRUE,
               pt.size = 0.5)

plot_title = "Vis_B3_2_Spatial_clus.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 2,
    height = 3.5)
print(p)
dev.off()

# Tumor pure percentage
round(table(sce$Region, sce$spatial.cluster)["Tumor_pure",]/colSums(table(sce$Region, sce$spatial.cluster)), 1) # 1, 5, 7
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14 
# 0.8 0.2 0.0 0.3 0.7 0.0 0.7 0.1 0.6 0.0 0.0 0.2 0.1 0.6 

# round(table(sce$Region, sce$spatial.cluster)["Tumor_Stroma_mix", ]/colSums(table(sce$Region, sce$spatial.cluster)), 2) # 1, 5, 7
# # 1   2   3   4   5   6   7   8   9  10  11  12  13  14 
# # 0.2 0.2 0.0 0.5 0.2 0.2 0.2 0.5 0.1 0.5 0.1 0.3 0.2 0.4 
# 
# round(table(sce$Region, sce$spatial.cluster)["Intratumoral_Stroma", ]/colSums(table(sce$Region, sce$spatial.cluster)), 2) # 1, 5, 7
# # 1   2   3   4   5   6   7   8   9  10  11  12  13  14 
# # 0.0 0.7 1.0 0.1 0.1 0.8 0.2 0.4 0.3 0.4 0.9 0.6 0.7 0.0 


round(table(sce$Level4_decon_max, sce$spatial.cluster)["Tu_B3_CYP4F8",]/colSums(table(sce$Level4_decon_max, sce$spatial.cluster)), 1) # 1, 4, 5, 7, 9, 14
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 0.96 0.28 0.00 0.65 0.78 0.05 0.71 0.64 0.67 0.32 0.00 0.19 0.08 1.00 

round(table(sce$Level4_decon_max, sce$spatial.cluster)["Fibroblast_B3",]/colSums(table(sce$Level4_decon_max, sce$spatial.cluster)), 2) # 1, 4, 5, 7, 9, 14

# decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/"
# decon <- read.csv(file.path(decon_path, "B3_2_spot_level4_RCTD.csv")) %>% # RCTD
# decon <- read.csv(file.path(decon_path, "B3_2_spot_Level4_decon.csv")) %>% # CARD
#   dplyr::rename(Barcode = X)

# # 1st dominant cell type
# decon_ <- decon %>% 
#   column_to_rownames('Barcode')
# 
# decon_$max <- colnames(decon_)[max.col(decon_, ties.method="first")]
# decon$max <- decon_$max
# 
# 
# CD <- as.data.frame(colData(sce)) %>%
#   select(Barcode, spatial.cluster)
# 
# decon_sub <- decon %>%
#   left_join(CD) %>%
#   filter(spatial.cluster %in% c("1", "4", "5", "7", "9", "14")) %>%
#   column_to_rownames("Barcode")
# 
# ## Dumb down
# # > round(table(decon_sub$max, decon_sub$spatial.cluster)["Tu_B3_CYP4F8",]/colSums(table(decon_sub$max, decon_sub$spatial.cluster)), 4) # 1, 4, 5, 7, 9, 14
# # 1      2      3      4      5      6      7      8      9     10     11     12     13     14 
# # 1.0000    NaN    NaN 0.9806 0.9874    NaN 0.9855    NaN 0.9552    NaN    NaN    NaN    NaN 1.0000 
# 
# decon_sub_stats <- decon_sub %>%
#   select(Tu_B3_CYP4F8, Fibroblast_B3, spatial.cluster) %>%
#   group_by(spatial.cluster) %>%
#   summarise(mean_tu = mean(Tu_B3_CYP4F8),
#             mean_fi = mean(Fibroblast_B3))
#               #quantile(Tu_B3_CYP4F8, 0.80))

# spatial.cluster mean_val # Mean > 60%
# <fct>              <dbl>
# 1 1                  0.788
# 2 4                  0.580
# 3 5                  0.670
# 4 7                  0.656
# 5 9                  0.620
# 6 14                 0.796


# Save plots --------------------------------------------------------------
p <- plotDimRed(sce, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster")
plot_title = "Vis_B3_2_UMAP.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 11,
    height = 5)
print(p)
dev.off()


# Merge DE ----------------------------------------------------------------
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Manuscript_Figure/Fig5_Visium_clustering_biology/B3/00_Seurat_DE_heatmap_ClusPlot_by_Cluster_helper.R")
# p <- clusterPlot_by_cluster(sce, type = "spot", sample_name = sample_name)

# # 1,5,7,9 vs 4, 14 ----------------------------------------------------- # 14 is stroma mix
# sce$spatial.cluster_merge_final <- ifelse(sce$spatial.cluster %in% c(1, 5, 7, 9), "1_5_7_9", 
#                                           ifelse(sce$spatial.cluster %in% c(4, 14), "4_14", sce$spatial.cluster))
# 
# sce$spatial.cluster_merge_final_consensus_TME <- 
#   case_when(sce$spatial.cluster_merge_final == "1_5_7_9" & sce$Region == "Tumor_pure" & sce$Level4_decon_max == "Tu_B3_CYP4F8" ~ "1_5_7_9_Consensus",
#             sce$spatial.cluster_merge_final == "4_14" & sce$Region == "Tumor_pure" & sce$Level4_decon_max == "Tu_B3_CYP4F8" ~ "4_14_Consensus",
#             TRUE ~ sce$spatial.cluster_merge_final)
# sce$spatial.cluster_merge_final_consensus_TME <- 
#   case_when(sce$spatial.cluster_merge_final_consensus_TME == "1_5_7_9" ~ "1_5_7_9_TME",
#             sce$spatial.cluster_merge_final_consensus_TME == "4_14" ~ "4_14_TME",
#             TRUE ~ sce$spatial.cluster_merge_final_consensus_TME)  

# # 1,5,7,9 vs 14 -----------------------------------------------------
# sce$spatial.cluster_merge_final <- ifelse(sce$spatial.cluster %in% c(1, 5, 7, 9), "1_5_7_9", 
#                                           ifelse(sce$spatial.cluster %in% c(14), "4_14", sce$spatial.cluster))
# 
# sce$spatial.cluster_merge_final_consensus_TME <- 
#   case_when(sce$spatial.cluster_merge_final == "1_5_7_9" & sce$Region == "Tumor_pure" & sce$Level4_decon_max == "Tu_B3_CYP4F8" ~ "1_5_7_9_Consensus",
#             sce$spatial.cluster_merge_final == "4_14" & sce$Region == "Tumor_pure" & sce$Level4_decon_max == "Tu_B3_CYP4F8" ~ "4_14_Consensus",
#             TRUE ~ sce$spatial.cluster_merge_final)
# sce$spatial.cluster_merge_final_consensus_TME <- 
#   case_when(sce$spatial.cluster_merge_final_consensus_TME == "1_5_7_9" ~ "1_5_7_9_TME",
#             sce$spatial.cluster_merge_final_consensus_TME == "4_14" ~ "4_14_TME",
#             TRUE ~ sce$spatial.cluster_merge_final_consensus_TME) 

# 1,5,9 vs 14 ----------------------------------------------------- # EMT order low
sce$spatial.cluster_merge_final <- ifelse(sce$spatial.cluster %in% c(1, 5, 9), "1_5_7_9",
                                          ifelse(sce$spatial.cluster %in% c(14), "4_14", sce$spatial.cluster))

sce$spatial.cluster_merge_final_consensus_TME <-
  case_when(sce$spatial.cluster_merge_final == "1_5_7_9" & sce$Region == "Tumor_pure" & sce$Level4_decon_max == "Tu_B3_CYP4F8" ~ "1_5_7_9_Consensus",
            sce$spatial.cluster_merge_final == "4_14" & sce$Region == "Tumor_pure" & sce$Level4_decon_max == "Tu_B3_CYP4F8" ~ "4_14_Consensus",
            TRUE ~ sce$spatial.cluster_merge_final)
sce$spatial.cluster_merge_final_consensus_TME <-
  case_when(sce$spatial.cluster_merge_final_consensus_TME == "1_5_7_9" ~ "1_5_7_9_TME",
            sce$spatial.cluster_merge_final_consensus_TME == "4_14" ~ "4_14_TME",
            TRUE ~ sce$spatial.cluster_merge_final_consensus_TME)

# # 1 vs 14 ----------------------------------------------------- # Volcano with NAMPT, no ADRADA2; Pathways no EMT
# sce$spatial.cluster_merge_final <- ifelse(sce$spatial.cluster %in% c(1), "1_5_7_9", 
#                                           ifelse(sce$spatial.cluster %in% c(14), "4_14", sce$spatial.cluster))
# 
# sce$spatial.cluster_merge_final_consensus_TME <- 
#   case_when(sce$spatial.cluster_merge_final == "1_5_7_9" & sce$Region == "Tumor_pure" & sce$Level4_decon_max == "Tu_B3_CYP4F8" ~ "1_5_7_9_Consensus",
#             sce$spatial.cluster_merge_final == "4_14" & sce$Region == "Tumor_pure" & sce$Level4_decon_max == "Tu_B3_CYP4F8" ~ "4_14_Consensus",
#             TRUE ~ sce$spatial.cluster_merge_final)
# sce$spatial.cluster_merge_final_consensus_TME <- 
#   case_when(sce$spatial.cluster_merge_final_consensus_TME == "1_5_7_9" ~ "1_5_7_9_TME",
#             sce$spatial.cluster_merge_final_consensus_TME == "4_14" ~ "4_14_TME",
#             TRUE ~ sce$spatial.cluster_merge_final_consensus_TME) 

saveRDS(sce, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered.rds")
# > table(sce$spatial.cluster_merge_final_consensus_TME)
# 1_5_7_9_Consensus       1_5_7_9_TME                10                11                12                13                 2                 3    4_14_Consensus          4_14_TME 
# 674               431               107               127                63                76                36                 1               197               234 
# 6                 8 
# 140                70 

p <- clusterPlot_by_cluster(sce, type = "spot", sample_name = sample_name, annotate = "spatial.cluster_merge_final_consensus_TME")

plot_title = "Vis_B3_2_Spatial_clus_by_clus_final.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 7,
    height = 10)
print(p)
dev.off()


# Consensus DE -----------------------------------------------------------
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2"
DE_result <- seurat_cluster_DE(sce, clusters = c("1_5_7_9_Consensus", "4_14_Consensus"), cluster_col="spatial.cluster_merge_final_consensus_TME", n_markers = 40)
contrast <- "ConsensusA_B"

write.csv(DE_result, file.path(save_path_DE, paste0(contrast, ".csv")))
head(DE_result)


# TME DE -----------------------------------------------------------------
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2"
DE_result <- seurat_cluster_DE(sce, clusters = c("1_5_7_9_TME", "4_14_TME"), cluster_col="spatial.cluster_merge_final_consensus_TME", n_markers = 40)
contrast <- "TMEA_B"

write.csv(DE_result, file.path(save_path_DE, paste0(contrast, ".csv")))
head(DE_result)
 
# plotDimRed(sce, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster") + facet_wrap(~ spatial.cluster)
# 
# sce$spatial.cluster_merge <- ifelse(sce$spatial.cluster %in% c(5, 9), "5_9", sce$spatial.cluster)
sce$spatial.cluster_merge2 <- ifelse(sce$spatial.cluster %in% c(1, 5, 9), "1_5_9", sce$spatial.cluster)
# 
# # > table(sce$Tu_Consensus, sce$spatial.cluster159reclus)
# # 
# #        10  11  12  13  14 159_1 159_2   2   3   4   6   7   8
# # FALSE 104 127  59  75  82   258   142  33   1 152 140  31  65
# # TRUE    3   0   4   1 143   392   244   3   0  54   0  38   5
# 
# sce$spatial.cluster159reclus_tu_consensus <- case_when(sce$spatial.cluster159reclus == "14" & sce$Tu_Consensus ~ "14", 
#                                                        sce$spatial.cluster159reclus == "14" & !sce$Tu_Consensus ~ "14_TME", 
#                                                        sce$spatial.cluster159reclus == "4" ~ "14_TME",
#                                                        sce$spatial.cluster159reclus == "159_1" & sce$Tu_Consensus ~ "159_1", 
#                                                        sce$spatial.cluster159reclus == "159_1" & !sce$Tu_Consensus ~ "159_1_TME",
#                                                        sce$spatial.cluster159reclus == "159_2" & sce$Tu_Consensus ~ "159_2",
#                                                        sce$spatial.cluster159reclus == "159_2" & !sce$Tu_Consensus ~ "159_2_TME",
#                                                        .default = sce$spatial.cluster159reclus)
# 
# sce$spatial.cluster_tu_consensus <- case_when(sce$spatial.cluster159reclus == "14" & sce$Tu_Consensus ~ "14", 
#                                               sce$spatial.cluster159reclus == "14" & !sce$Tu_Consensus ~ "14_TME", 
#                                               sce$spatial.cluster159reclus == "4" ~ "14_TME",
#                                               sce$spatial.cluster159reclus == "159_1" & sce$Tu_Consensus ~ "159", 
#                                               sce$spatial.cluster159reclus == "159_1" & !sce$Tu_Consensus ~ "159_TME",
#                                               sce$spatial.cluster159reclus == "159_2" & sce$Tu_Consensus ~ "159",
#                                               sce$spatial.cluster159reclus == "159_2" & !sce$Tu_Consensus ~ "159_TME",
#                                               .default = sce$spatial.cluster159reclus)
# 
# # TME clus 14: Tu consensus FALSE in clus 14 + entire clus 4
# # TME clus 159: Tu consensus FALSE in clus 159_1 and clus 159_2
# 
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2"
DE_result <- seurat_cluster_DE(sce, clusters = c("1_5_9", "14"), cluster_col="spatial.cluster_merge2", n_markers = 40)
contrast <- "1&5&9_14"

# DE_result <- seurat_cluster_DE(sce, cluster_col="spatial.cluster159reclus_tu_consensus", clusters = c("14", "159_1", "159_2"), n_markers = 100)
DE_result <- seurat_cluster_DE(sce, cluster_col="spatial.cluster_tu_consensus", clusters = c("14", "159"), n_markers = 100)
DE_result2 <- seurat_cluster_DE(sce, cluster_col="spatial.cluster_tu_consensus", clusters = c("14_TME", "159_TME"), n_markers = 100)


write.csv(DE_result, file.path(save_path_DE, paste0(contrast, ".csv")))
head(DE_result)
