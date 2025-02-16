library(BayesSpace)
library(patchwork)
library(ggplot2)

# disease = "DLBCL"
# sample = "DLBCL_6"
disease = "Lung"
sample = "L1_4"
# disease = "Breast"
# sample = "B4_2"
# resolution = "spots"
# resolution = "subspots"
read_path = paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/", disease, "/", sample)

# if(resolution == "spot"){ # spot
sce <- readRDS(file.path(read_path, paste0(sample, "_baye_clustered.rds"))) # spot
assay(sce, "log1p") <- log1p(counts(sce))
# raw_assay_name = "counts"

# }else{ # subspot
sce_enhanced <- readRDS(file.path(read_path, paste0(sample, "_baye_clustered_all_enhanced_expr.rds"))) # subspot
# raw_assay_name = "decont_subspot"
# }

dim(rowData(sce_enhanced))

# -------------------------------------------------------------------------
markers <- c("MS4A1",
             "CXCL13",
             "CXCR5",
             "CCL19")

# markers <- c("CD14", "CD68", "CSF1R", "ITGAM",# Marcrophage
#              "CD80", "CD86", "IL1B", # M1
#              "CD163", "MRC1", "IL10") # M2

markers <- c("MS4A1", "CXCL13", # B cells
             "CD3D", "CD4",     # T cells
             "CD14", "CD68", "CD163") # Marcrophage

df_save <- data.frame(as.data.frame(t(assay(sce, "log1p")[markers, ])),
                      x = sce$array_col,
                      y = sce$array_row)
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig5a_vis.csv")
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS9b_vis.csv")

df_save <- data.frame(as.data.frame(t(assay(sce_enhanced, "log1p")[markers, ])),
                      x = sce_enhanced$array_col,
                      y = sce_enhanced$array_row)
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig5a_vis_enhanced.csv")
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS9b_vis_enhanced.csv")

plot_expression <- function(sce_obj = sce, marker = "CD4"){
  featurePlot(sce_obj, marker, assay.type = "log1p", color=NA) +
    scale_x_reverse() + 
    viridis::scale_fill_viridis(option="A") +
    labs(title=marker, fill="Log\nexpression") +
    theme(plot.title = element_blank())
  # theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
}

# With legend
feature.plots1 <- purrr::map(markers, function(x) plot_expression(sce, x)) 
enhanced.feature.plots1 <- purrr::map(markers, function(x) plot_expression(sce_enhanced, x))

# Without legend
feature.plots2 <- purrr::map(markers, function(x) plot_expression(sce, x) + theme(legend.position = "none")) 
enhanced.feature.plots2 <- purrr::map(markers, function(x) plot_expression(sce_enhanced, x) + theme(legend.position = "none"))

# -------------------------------------------------------------------------
saveFeaturePlots <- function(feature.plots, enhanced.feature.plots, legend, width, height){
  p1 <- patchwork::wrap_plots(feature.plots, ncol = 7) + plot_layout(guides = "collect") 
  p2 <- patchwork::wrap_plots(enhanced.feature.plots, ncol = 7) + plot_layout(guides = "collect") 
  
  fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/All_samples_TLS"
  # pdf(file = file.path(file.path(fig_path, legend), paste0(sample, "_TLS_subspot_spot.pdf")),
  pdf(file = file.path(file.path(fig_path, legend), paste0(sample, "_TLS_subspot_spot_BTMacro.pdf")),
      width = width,
      height = height)
  print(p1 / p2)
  dev.off()
}

saveFeaturePlots(feature.plots1, enhanced.feature.plots1, "legend", width = 17, height = 25)
saveFeaturePlots(feature.plots2, enhanced.feature.plots2, "nolegend", width = 17, height = 4.5)



