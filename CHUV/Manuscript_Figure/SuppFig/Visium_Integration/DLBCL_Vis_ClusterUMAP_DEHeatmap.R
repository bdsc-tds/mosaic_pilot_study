library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)

disease = "dlbcl"

# After spotclean -----------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
UMAP_name = paste0("/", str_to_title(disease), "-UMAP-")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_UMAP_Heatmap"


# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(savepath.spotclean, save_rds_name))

pdf(paste0(figpath, paste0(UMAP_name, "Cluster", ".pdf")), width = 7, height = 6)
print(DimPlot(seu, label=TRUE) + labs(y= "UMAP_2", x = "UMAP_1"))
dev.off()


# DLBCL post spotclean DE ----------------------------------------------------------------
seu_de <- PrepSCTFindMarkers(seu)
seu_de <- ScaleData(seu_de, assay = "SCT")

all.markers <- FindAllMarkers(object = seu_de, only.pos = TRUE)
all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15


p <- DoHeatmap(seu_de, features = top15$gene, size = 20, angle = 45, hjust = 0.8) + # color bar text size
  theme(axis.text = element_text(size = 30),                                        # gene text size
        legend.position = "bottom",
        legend.title = element_text(size = 40, face = "bold", hjust = 100),
        legend.text = element_text(size = 40),
        plot.margin = unit(c(1,2,1,1), "cm")) +
  guides(colour="none",
         fill = guide_colorbar(title = "Avg. Expression", title.position = "top", title.hjust = 0.5,
                               barwidth = unit(140, "cm"), barheight = unit(1, "cm")))

# p2 <- p + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white")
p2 <- p + scale_fill_gradientn(colors = c("purple", "royalblue","brown", "orange"), na.value = "white")



pdf(paste0(figpath, paste0("/DLBCL_Vis_DE_Heatmap_15_prepsct_scale.pdf")), width = 60, height = 50)
print(p)
dev.off()

# pdf(paste0(figpath, paste0("/DLBCL_Vis_DE_Heatmap_15_prepsct_scale_bluewhite.pdf")), width = 60, height = 50)
# print(p2)
# dev.off()

pdf(paste0(figpath, paste0("/DLBCL_Vis_DE_Heatmap_15_prepsct_scale_purpleyellow.pdf")), width = 60, height = 50)
print(p2)
dev.off()

# write.csv(all.markers, paste0(savepath, paste0("/DLBCL_Integrated_Visium_DE.csv")))


# -------------------------------------------------------------------------
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"

# Heatmap of tumor clus versus sample id ----------------------------------
seu_sub <- seu 
seu_sub$sample_id <- paste0(substr(seu$sample_id, 1, 1), substr(seu$sample_id, 7, 7))
seu_sub$seurat_clusters <- droplevels(seu_sub$seurat_clusters)
Idents(seu_sub) <- seu_sub$seurat_clusters
df <- as.data.frame(table(seu_sub$seurat_clusters, seu_sub$sample_id)) %>%
  dplyr::rename(cluster = Var1,
                patient = Var2,
                Count = Freq)

# df2 <- df %>%
#   group_by(cluster) %>%
#   mutate(countT= sum(Count)) %>%
#   mutate(per=round(100*Count/countT,0))

# Weight by sample size to compare D5 and D6 fairly
df2_ <- df %>%
  group_by(patient) %>%
  mutate(SampleSize = sum(Count)) %>%
  ungroup() %>%
  mutate(Weight = 1/SampleSize, 
         WeightedCount = Count * Weight) %>%
  group_by(cluster) %>%
  mutate(per = round(100 * WeightedCount/sum(WeightedCount), 0))
df2_$patient <- factor(df2_$patient, levels = rev(levels(factor(df2_$patient))))

# Create the heatmap pt per cluster
p <- ggplot(df2_, aes(x = cluster, y = patient, fill = per)) + 
  geom_tile() + 
  scale_x_discrete(position = "top") +
  #scale_fill_viridis_c(option = "mako", direction = -1) + 
  scale_fill_gradientn(colors = c("white", "royalblue", "purple")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15)) + 
  labs(fill = "Per-cluster Percentage", title = "", 
       x = "Integrated DLBCL Visium Seurat Clusters", y = "Sample") + 
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                               barwidth = unit(23, "cm"), barheight = unit(0.5, "cm")))
p

write.csv(df2_, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS13b.csv")

pdf(file.path(figpath, "/DLBCL_Vis_Pt_Clus_Heatmap.pdf"), width = 13, height = 6)
print(p)
dev.off()


# Heatmap of tumor clus versus Region ----------------------------------
seu_sub <- seu 
df <- as.data.frame(table(seu_sub$seurat_clusters, seu_sub$Region)) %>%
  dplyr::rename(cluster = Var1,
                Region = Var2,
                Count = Freq)
df2 <- df %>%
  group_by(cluster) %>%
  mutate(countT= sum(Count)) %>%
  mutate(per=round(100*Count/countT,0))
df2$Region <- factor(df2$Region, levels = rev(levels(factor(df2$Region))))


# Create the heatmap pt per cluster
p <- ggplot(df2, aes(x = cluster, y = Region, fill = per)) + 
  geom_tile() + 
  scale_x_discrete(position = "top") +
  # scale_fill_viridis_c(option = "rocket", direction = -1) + 
  scale_fill_gradientn(colors = c("white", "orange", "brown")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15)) + 
  labs(fill = "Per-cluster Percentage", title = "", 
       x = "Integrated DLBCL Visium Seurat Clusters", y = "Pathology Annotation") + 
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                               barwidth = unit(23, "cm"), barheight = unit(0.5, "cm")))

write.csv(df2, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS13c.csv")

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "/DLBCL_Vis_Patho_Clus_Heatmap.pdf"), width = 15, height = 6.8)
print(p)
dev.off()


