library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
disease = "breast"

# Before spotclean ------------------------------------------------------
datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
datapath = paste0(datapath.filtered, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTnoSpotClean_filtered.rds")

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig1/Fig1d"
# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(datapath, save_rds_name))
p <- DimPlot(seu, label=FALSE, group.by="sample_id") + labs(y= "UMAP_2", x = "UMAP_1") + ggtitle("")

red_df <- data.frame(seu@reductions$umap@cell.embeddings)
red_df$section <- seu$sample_id
write.csv(red_df, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig1d_vis.csv")

pdf(file.path(figpath, "Vis_UMAP_Breast_Sample.pdf"), width = 7, height = 6)
print(p)
dev.off()