chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

Idents(chrom) <- chrom$Harmonised_Level4

df_save <- data.frame(
  chrom@reductions$umap@cell.embeddings,
  clus = chrom$Harmonised_Level4
)

write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig6a_chrom.csv")

p <- DimPlot(chrom) + labs(y= "UMAP_2", x = "UMAP_1") + 
  scale_color_manual(values = level4_cellcolors)

plot <- LabelClusters(plot = p, id = "ident", box = TRUE, color = "white") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        axis.line = element_blank()) + 
  labs(x = "", y = "")


figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig6"
pdf(file.path(figpath, "chrom_cluster_by_annot.pdf"), width = 8, height = 8)
print(plot)
dev.off()