library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SCpubr)
library(patchwork)

disease = "dlbcl"

# After spotclean -----------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
UMAP_name = paste0("/", str_to_title(disease), "-UMAP-")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
savepath = savepath.spotclean

# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(savepath, save_rds_name))

DimPlot(seu, label = TRUE, group.by = "Level4_decon_max") + facet_wrap(~Level4_decon_max) + 
  theme(legend.position = "none")
DimPlot(seu, label = TRUE, group.by = "seurat_clusters") + facet_wrap(~seurat_clusters) + 
  theme(legend.position = "none")
DimPlot(seu, label = TRUE, group.by = "Region") + facet_wrap(~Region) + 
  theme(legend.position = "none")




# Plasma cells
FeaturePlot(seu, features = c("IGKC", "IGHG3", "IGHA", "JCHAIN", "IGHG1"))
# Epithelium
FeaturePlot(seu, features = c("MUCL3", "LCN2", "TFF1", "MUC5AC"))
# Stroma
FeaturePlot(seu, features = c("THBS1", "POSTN", "COL3A1", "COL1A2", "COL1A1"), ncol = 3)
# Vessels
FeaturePlot(seu, features = c("ACKR1", "MEOX2", "CCL14", "TSPAN7", "PLA1A", "VWF", "GIMAP4"), ncol = 4)

# Marcophage 
FeaturePlot(seu, features = c("CD68", "CD163", "CD14", "ITGAM", "MRC1", "MSR1", "SIRPA", "ARG1"))
FeaturePlot(seu, features = rownames(seu)[grepl("^MMP", rownames(seu))])
FeaturePlot(seu, features = rownames(seu)[grepl("^TGFB", rownames(seu))])


table# Tu_D1_LMO2
FeaturePlot(seu, features = c("LMO2"))
# Tu_D1_SMIM14
FeaturePlot(seu, features = c("SMIM14"))

# # Tu_D2_mito
# FeaturePlot(seu, features = c("mito?"))
FeaturePlot(seu, features = rownames(seu)[grepl("^MT-", rownames(seu))])

# check marker gene from heatmap for clus 14 and 16
FeaturePlot(seu, features = c("SLC40A1"))

FeaturePlot(seu, features = c("MMP12"))

# Tu_D3_FAM3C
FeaturePlot(seu, features = c("FAM3C"))
# Tu_D4_BCL7A
FeaturePlot(seu, features = c("BCL7A"))
# Tu_D5_CCL22
FeaturePlot(seu, features = c("CCL22"))
# Tu_D6_BCL2
FeaturePlot(seu, features = c("BCL2"))

seu_sub <- seu[, seu$seurat_clusters == "18"]
table(seu_sub$sample_id, seu_sub$Level4_decon_max)

# SuppFig (cluster vs marker genes dotplot) -------------------------------
# Key tumor and annotation markers DLBCL
tu_dlbcl <- c(# "CD19", "CD79A", "BCL6", "IRF4", "MYC", "BCL2", "MKI67", "TP53", # Tumor
  "IGKC", "IGHG3", "JCHAIN", "IGHG1",                      # Plasma
  "MUCL3", "LCN2", "TFF1", "MUC5AC",                               # Epithelium
  "THBS1", "POSTN", "COL3A1", "COL1A2", "COL1A1",                  # Stroma
  "ACKR1", "MEOX2", "CCL14", "TSPAN7", "PLA1A", "VWF", "GIMAP4",    # Vessels
  "CD4", "CD28", "CD68", "CD163", "CD14", "SIRPA"            # Immune cells
)

p <- SCpubr::do_DotPlot(sample = seu, features = tu_dlbcl, 
                        axis.text.x.angle = 0,
                        flip = TRUE# , legend.position = "right"
                        ) + 
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8), 
        axis.line = element_blank()) + 
  scale_fill_viridis_c(option = "magma", direction = -1) + 
  scale_x_discrete(position = "top") +
  guides(fill = guide_colorbar(title = "Avg. Expression", title.position = "top", title.hjust = 0.5, 
                               barwidth = unit(5, "cm"), barheight = unit(0.5, "cm")))
p

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "/DLBCL_Vis_Marker_DotPlot.pdf"), width = 8, height = 7)
print(p)
dev.off()


# Chrom pt spec marker in Visium  -------------------------------------------
plot_UMAP_agg_markers <- function(seu, genelist, genelistname = "Tu_D1"){
  seu <- AddModuleScore(object = seu, features = genelist, name = genelistname)
  p <- FeaturePlot(object = seu, features = paste0(genelistname, "1"), label = TRUE) + 
    scale_color_viridis_c(option = "mako") + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank()) + 
    ggtitle(paste0(genelistname))
  return(p)
}

Tu_D1 <- c("NME2", "FCRL1", "LMO2", "IGKC", "MT-CO3", "TMSB4X", "GRHPR",
           "MPEG1", "CD1C", "CD52", "LTB", "SMIM14", "TNFRSF13C", "CD79A",
           "CD22", "FCRL2", "FCRL5")
Tu_D2 <- c("MT-CYB", "NIBAN3", "SPIB", "MT-ATP6", "HIST1H1C", "IGHM",
           "MT-ND4", "TCL1A", "MT-ND4L", "HIST1H1R") 
Tu_D3 <- c("CD83", "SWAP70", "SEL1L3", "FCRL3", "ACTB", "IL4I1", "FAM3C",
           "ACTG2", "CD74", "MS4A1", "MT-CYB", "CNN2", "LRMP", "LCP1")
Tu_D4 <- c("GGA2", "HELLS", "RASGRP2", "NKX6-3", "DTX1", "BCL7A", "MAT2A",
           "ARGLU1", "ALOX5", "AKNA", "ARHGEF1", "PNN", "TMC8")
Tu_D5 <- c("TSPAN33", "SPATC1", "NFKB2", "BCL2L1", "OGI", "CCL17",
           "CD40", "DUSP2", "ITPKB", "CCL22")
Tu_D6 <- c("POU2F2", "CCDC88A", "LENG8", "KLHL6", "TNFRSF13B", "BCL2",
           "NFATC1", "MT-ND6", "NIBAN3", "FCRL2")

markers <- lapply(list(Tu_D1, Tu_D2, Tu_D3, Tu_D4, Tu_D5, Tu_D6), function(x) x[x %in% rownames(seu)])
names(markers) <- c("Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6")

for(i in 1:6){
  seu <- AddModuleScore(object = seu, 
                        features = markers[i], 
                        name = names(markers)[i])
}

df_save <- data.frame(Tu_D1 = seu$Tu_D11,
                      Tu_D2 = seu$Tu_D21,
                      Tu_D3 = seu$Tu_D31,
                      Tu_D4 = seu$Tu_D41,
                      Tu_D5 = seu$Tu_D51,
                      Tu_D6 = seu$Tu_D61,
                      cluster = seu$seurat_clusters,
                      seu@reductions$umap@cell.embeddings
                      )
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS13ae.csv")

UMAP_plots <- purrr::map2(markers, names(markers), function(genelist, genelistname) plot_UMAP_agg_markers(seu, genelist, genelistname))
plot <- patchwork::wrap_plots(UMAP_plots, ncol=3) + plot_annotation(title = "Aggregated Patient-specific Tumor Marker Expression from DLBCL Chromium", 
                                                                    theme = theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold")))

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "DLBCL_Vis_AggMarker_UMAP.pdf"), width = 15, height = 10)
print(plot)
dev.off()


# -------------------------------------------------------------------------
## Single cell DLBCL markers 
# tu_dlbcl <- c("NME2", "FCRL1", "LMO2", "IGKC", "MT-CO3", "TMSB4X", "GRHPR",
#               "MPEG1", "CD1C", "CD52", "LTB", "SMIM14", "TNFRSF13C", "CD79A",
#               "CD22", "FCRL2", "FCRL5",                                    # D1 clus 2, 5, 13, 19:
#               "MT-CYB", "NIBAN3", "SPIB", "MT-ATP6", "HIST1H1C", "IGHM",
#               "MT-ND4", "TCL1A", "MT-ND4L", "HIST1H1R",                    # D2 clus 0:
#               "CD83", "SWAP70", "SEL1L3", "FCRL3", "ACTB", "IL4I1", "FAM3C",
#               "ACTG2", "CD74", "MS4A1", "MT-CYB", "CNN2", "LRMP", "LCP1",  # D3 clus 1, 6, 26:
#               "GGA2", "HELLS", "RASGRP2", "NKX6-3", "DTX1", "BCL7A", "MAT2A",
#               "ARGLU1", "ALOX5", "AKNA", "ARHGEF1", "PNN", "TMC8",         # D4 clus 9&22:
#               "TSPAN33", "SPATC1", "NFKB2", "BCL2L1", "OGI", "CCL17",
#               "CD40", "DUSP2", "ITPKB", "CCL22",                           # D5 clus 17:
#               "POU2F2", "CCDC88A", "LENG8", "KLHL6", "TNFRSF13B", "BCL2",
#               "NFATC1", "MT-ND6", "NIBAN3", "FCRL2"                        # D6 clus 12:
# )
# 
# 
# # p <- SCpubr::do_DotPlot(sample = seu, features = tu_dlbcl, 
# #                         axis.text.x.angle = 45,
# #                         flip = FALSE, legend.position = "right"
# # ) + 
# #   theme(axis.ticks = element_blank(),
# #         axis.line = element_blank()) + 
# #   scale_fill_viridis_c(option = "magma", direction = -1) + 
# #   scale_x_discrete(position = "top")
# # p
# 
# ## C'est impossible


# -------------------------------------------------------------------------
# long format with percentage per-clus
plt_df <- as.data.frame(table(seu$Level4_decon_max, seu$seurat_clusters)) %>%
  group_by(Var2) %>%
  mutate(countT= sum(Freq)) %>%
  mutate(per=round(100*Freq/countT,0))

data_frame <- data.frame(cluster = plt_df$Var2, decon = plt_df$Var1, per = plt_df$per)
data_frame$decon <- factor(data_frame$decon, levels = rev(levels(factor(data_frame$decon))))

# Create the heatmap
p <- ggplot(data_frame, aes(x = cluster, y = decon, fill = per)) + 
  geom_tile() + 
  scale_x_discrete(position = "top") +
  # scale_fill_viridis_c(option = "plasma") + 
  scale_fill_gradientn(colors = c("white", "orchid", "#FC0FC0")) +
  theme_minimal() +
  labs(fill = "Per-cluster\nPercentage", title = "", 
       x = "Integrated DLBCL Visium Seurat Clusters", y = "Level 4 Deconvolution Majority Vote") 

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "DLBCL_Vis_Decon_Clus_Heatmap_pink.pdf"), width = 12, height = 7)
print(p)
dev.off()




