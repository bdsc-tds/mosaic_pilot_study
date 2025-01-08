library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig6/"

# DLBCLnChromium_Marker_Gene_List <- c(
#   "MS4A1", "TNFRSF13C", "CD79B", "CD37", "PSMB8", "CD19", "TYMS",
#   "TUBB", "TOP2A", "POLD4", "CD47", "CD52", "BLK", "CD38", "MAP2K1",
#   "CD40", "BCL2L1", "TNFRSF8", "SMO", "RARA", "TYK2", "TNFRSF10B"
# )

DLBCLnChromium_Marker_Gene_List <- c(
  "MS4A1", "TNFRSF13C", "CD79B", "CD37", "CD19", "BTK",
  "TUBB",  "CD47", "CD52", "CD38", "MAP2K1",
  "CD40", "BCL2L1", "TNFRSF8", "RARA", "TYK2", "TNFRSF10B",
  "MEF2B", "NOP56", "MYO1G", "SEMA7A", "FCRL5",
  "FCRL3", "BCL6", "BCL11A", "BCL2L11", "BCL7A")

# Chromium ----------------------------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

## Separate tumor subtypes
chrom$level1_5 <- chrom$Harmonised_Level4
chrom$level1_5 <- ifelse(chrom$Level1 == "B", "B cells", 
                         ifelse(chrom$Level1 == "Epithelia", "Epithelia", 
                                ifelse(chrom$Level1 == "Myeloid", "Myeloid", 
                                       ifelse(chrom$Level1 == "Stroma", "Stroma", 
                                              ifelse(chrom$Level1 == "T_NK", "T/NK", 
                                                     chrom$Harmonised_Level4)))))
Idents(chrom) <- factor(chrom$level1_5, levels = c("Epithelia", "Stroma", "B cells", "Myeloid", "T/NK", 
                                                   names(table(chrom$level1_5))[grepl("Tu_D", names(table(chrom$level1_5)))]))
                      
# chrom <- chrom[rownames(chrom) %in% DLBCLnChromium_Marker_Gene_List, ]

p_chrom <- DotPlot(
  # assay = "RNA",
  chrom,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 10.5, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank()
        ,legend.position = "none"
        )

# For legend
p_chrom_legen <- DotPlot(
  chrom,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 10.5, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank()
        ,legend.position = "bottom"
  )

# Visium ------------------------------------------------------------------
vispath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/dlbcl/spotclean/Results/"
vis <- readRDS(file.path(vispath, "Dlbcl-merge-SCTpostSpotClean.rds"))

Idents(vis) <- factor(vis$annot, levels = c("Epithelium", "Stroma", "Necrosis", "Plasma", "Vessels/Immune", 
                                            "Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6"))

# vis_small <- subset(vis, features = DLBCLnChromium_Marker_Gene_List)
# vis[rownames(vis) %in% DLBCLnChromium_Marker_Gene_List, ]

p_vis <- DotPlot(
  vis,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 10.5, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank()
        ,legend.position = "none"
        )

# For legend
p_vis_legen <- DotPlot(
  vis,
  features = rev(DLBCLnChromium_Marker_Gene_List)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank()
        ,legend.position = "bottom"
  )


# GeoMx ----------------------------------------------------------------
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/dlbcl_seu_ruv.rds")
# head(Idents(geo)) # Levels: Other, Macrophage, T cells, Tu_D1 Tu_D2 Tu_D3 Tu_D4 Tu_D5 Tu_D6
Idents(geo) <- factor(geo$clusters, levels = c("Other", "Macrophage", "T cells", "Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6"))
# geo <- geo[rownames(geo) %in% DLBCLnChromium_Marker_Gene_List, ]

p_geo <- DotPlot(
  geo,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 10.5, hjust = 1, vjust = 0.5),
        # axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank()
        ,legend.position = "none"
        ) + scale_size(range = c(1, 9))

# For legend
p_geo_legen <- DotPlot(
  geo,
  features = rev(DLBCLnChromium_Marker_Gene_List)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, size = 10.5, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank()
        ,legend.position = "bottom"
  ) + scale_size(range = c(1, 9))

# Stitch plots ------------------------------------------------------------
p_final <- p_chrom + theme(plot.margin = margin(0, 4, 0, 0.5)) + 
  p_vis + theme(plot.margin = margin(0, 4, 0, 4)) +
  p_geo + theme(plot.margin = margin(0, 0.5, 0, 4)) +
  plot_layout(widths = c(1.65, 1.05, 0.93))

# -------------------------------------------------------------------------
# plot_title = "d_Three_Tech_Cluster_Dotplot_level2.pdf"
plot_title = "d_Three_Tech_Cluster_Dotplot_level2_new_gene.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10.2,
    height = 7.7)
print(p_final)
dev.off()

# Legend chrom
plot_title = "d_Three_Tech_Cluster_Dotplot_level2_new_gene_chrom_legend.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 8)
print(p_chrom_legen)
dev.off()

# Legend vis
plot_title = "d_Three_Tech_Cluster_Dotplot_level2_new_gene_vis_legend.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 8)
print(p_vis_legen)
dev.off()

# Legend geo
plot_title = "d_Three_Tech_Cluster_Dotplot_level2_new_gene_geo_legend.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 8)
print(p_geo_legen)
dev.off()



