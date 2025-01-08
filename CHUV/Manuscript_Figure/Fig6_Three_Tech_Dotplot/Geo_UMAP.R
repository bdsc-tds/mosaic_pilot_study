library(ggspavis)
library(patchwork)

absolute_path_cur <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/"
source(paste0(absolute_path_cur, "env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/GeoMx_init.R"))
read_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_test_norm/")


# -------------------------------------------------------------------------
disease = "dlbcl"
spe_ruv <- readRDS(file.path("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", paste0(disease, "_spe_ruv.rds")))

spe_ruv$clusters <- case_when(spe_ruv$cell_fraction == "Macro" ~ "Macrophage",
                              spe_ruv$cell_fraction == "B cells" & spe_ruv$patient == "D1" ~ "Tu_D1", 
                              spe_ruv$cell_fraction == "B cells" & spe_ruv$patient == "D2" ~ "Tu_D2", 
                              spe_ruv$cell_fraction == "B cells" & spe_ruv$patient == "D3" ~ "Tu_D3", 
                              spe_ruv$cell_fraction == "B cells" & spe_ruv$patient == "D4" ~ "Tu_D4", 
                              spe_ruv$cell_fraction == "B cells" & spe_ruv$patient == "D5" ~ "Tu_D5", 
                              spe_ruv$cell_fraction == "B cells" & spe_ruv$patient == "D6" ~ "Tu_D6",
                              .default = spe_ruv$cell_fraction
)

dlbcl_cluster_type <- c("T cells", "Macrophage", "Other", "Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6")
dlbcl_cluster_color <- c("#4169E1", "#9A32CD", "#388E8E", "#FF8C00", "#EEEE00", "#FFD700", "#A2CD5A", "#00EE76", "#ADFF2F")

names(dlbcl_cluster_color) <- dlbcl_cluster_type

spe_ruv$clusters <- factor(spe_ruv$clusters, levels = dlbcl_cluster_type)


library(Seurat)
seu_ruv <- as.Seurat(spe_ruv)
Idents(seu_ruv) <- seu_ruv$clusters
datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched"
saveRDS(seu_ruv, file.path(datapath, "dlbcl_seu_ruv.rds"))

plot <- DimPlot(seu_ruv, reduction = "UMAP", cols = dlbcl_cluster_color) + NoLegend()
p <- LabelClusters(plot = plot, id = "ident", box = TRUE, color = "white") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  labs(x = "", y = "")

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig6"
pdf(file.path(figpath, "geo_cluster_by_annot.pdf"), width = 6, height = 6)
print(p)
dev.off()
