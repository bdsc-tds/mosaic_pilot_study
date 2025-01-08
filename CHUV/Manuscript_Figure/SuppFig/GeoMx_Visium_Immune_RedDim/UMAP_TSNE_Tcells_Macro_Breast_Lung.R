library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")

# helper ------------------------------------------------------------------
plotCTFracDimRed <- function(seu, feat_name = "Macrophage_frac", 
                             reduction = "umap", pt.size = 4,
                             col_vec = c("#F0F0F0","purple", "#301934"),
                             shape.by = NULL){
  p <- FeaturePlot(seu, feat_name, reduction = reduction, pt.size = pt.size, 
                   shape.by = shape.by) +
    scale_color_gradientn(colors = col_vec) + 
    theme(axis.line = element_blank(), 
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank()) +
    ggtitle(NULL) 
  
  return(p)
}

# Read Data ----------------------------------------------------------------

disease_list <- c("breast", "lung")

for(disease in disease_list){
  ## Geo
  datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
  spe_ruv <- readRDS(file.path(datapath, paste0(disease, "_spe_ruv.rds")))
  if(class(assay(spe_ruv, "logcounts"))[1] != "dgCMatrix"){assay(spe_ruv, "logcounts") <- as(as.matrix(assay(spe_ruv, "logcounts")), "dgCMatrix")}
  seu_ruv <- as.Seurat(spe_ruv)
  seu_ruv$cell_fraction_Tcells <- factor(ifelse(seu_ruv$cell_fraction == "T cells", "T cells", "Not T cells"), levels = c( "Not T cells", "T cells"))
  seu_ruv$cell_fraction_Macrophage <- factor(ifelse(seu_ruv$cell_fraction == "Macro", "Macrophage", "Not Macrophage"), levels = c( "Not Macrophage", "Macrophage"))
  
  ## Vis
  datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
  savepath.spotclean = paste0(datapath.spotclean, "/Results")
  save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
  seu <- readRDS(paste0(savepath.spotclean, save_rds_name))
  
  ## Geo
  p1 <- plotCTFracDimRed(seu_ruv, feat_name = "Tcells_frac", pt.size = 2, reduction = "TSNE", 
                         col_vec = c("grey", "blue", "darkblue"), shape.by = "cell_fraction_Tcells") + 
    guides(shape = "none", color = guide_colorbar())
  p2 <- plotCTFracDimRed(seu_ruv, feat_name = "Macrophage_frac",  pt.size = 2, reduction = "TSNE", 
                         col_vec = c("grey","purple", "#301934"), shape.by = "cell_fraction_Macrophage")+
    guides(shape = "none", color = guide_colorbar())
  
  ## Vis
  p3 <- plotCTFracDimRed(seu, feat_name = "Tcells_frac",  pt.size = 0.3, col_vec = c("#F0F0F0", "blue", "darkblue")) 
  p4 <- plotCTFracDimRed(seu, feat_name = "Macrophage_frac",  pt.size = 0.3, col_vec = c("#F0F0F0","purple", "#301934"))
  
  
  p_indication_geo <- p1 + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) | 
    p2 + theme(plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
  
  p_indication_vis <- p3 + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) | 
    p4 + theme(plot.margin = unit(c(0, 0, 0, 0.5), "cm")) 
  
  assign(paste0("p_", disease, "_geo"), p_indication_geo)
  assign(paste0("p_", disease, "_vis"), p_indication_vis)
}

p_breast_geo
p_lung_geo
p_breast_vis
p_lung_vis


# -------------------------------------------------------------------------
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig3"

plot_title = "Geo_breast_TSNE.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 4)
print(p_breast_geo)
dev.off()

plot_title = "Geo_lung_TSNE.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 4)
print(p_lung_geo)
dev.off()

plot_title = "Vis_breast_UMAP.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 4)
print(p_breast_vis)
dev.off()

plot_title = "Vis_lung_UMAP.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 4)
print(p_lung_vis)
dev.off()


