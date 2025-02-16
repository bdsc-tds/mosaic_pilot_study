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

# Better plots for SCT before and after spotclean -----------------------
vis <- readRDS(paste0(savepath.spotclean, save_rds_name))

vis$annot <- case_when(vis$seurat_clusters == "0" ~ "Tu_D1",
                       vis$seurat_clusters %in% c("15", "18") ~ "Tu_D2",
                       vis$seurat_clusters %in% c("3", "8", "11") ~ "Tu_D3",
                       vis$seurat_clusters == "17" ~ "Tu_D4",
                       vis$seurat_clusters == "4" ~ "Tu_D5",
                       vis$seurat_clusters == "5" ~ "Tu_D6",
                       vis$seurat_clusters == "1" ~ "Stroma",
                       vis$seurat_clusters == "6" ~ "Epithelium",
                       vis$seurat_clusters %in% c("2", "12", "14", "16") ~ "Necrosis",
                       vis$seurat_clusters %in% c("9", "10", "13") ~ "Plasma",
                       vis$seurat_clusters == "7" ~ "Vessels/Immune")


df_save <- data.frame(
  vis@reductions$umap@cell.embeddings,
  clus = vis$annot
)

write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig6a_vis.csv")

Idents(vis) <- vis$annot
p <- DimPlot(vis) + labs(y= "UMAP_2", x = "UMAP_1") + 
  scale_color_manual(values = c(
    Tu_D1 = "#FF8C00",
    Tu_D2 = "#EEEE00",
    Tu_D3 = "#FFD700",
    Tu_D4 = "#A2CD5A",
    Tu_D5 = "#00EE76",
    Tu_D6 = "#ADFF2F",
    Stroma = "#388E8E",
    `Vessels/Immune` =   "#FF7256",
    Plasma = "#CD00CD",
    Necrosis = "#666666",
    Epithelium = "#BC8F8F"
  ))

plot <- LabelClusters(plot = p, id = "ident", box = TRUE, color = "white") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        axis.line = element_blank()) + 
  labs(x = "", y = "")


figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig6"
pdf(file.path(figpath, "vis_cluster_by_annot.pdf"), width = 6.5, height = 6)
print(plot)
dev.off()
