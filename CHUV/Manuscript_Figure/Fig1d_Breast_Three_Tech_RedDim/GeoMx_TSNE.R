library(ggspavis)
library(patchwork)
library(SpatialExperiment)

disease = "breast"
geopath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig1/Fig1d"

spe <- readRDS(file.path(geopath, paste0(disease, "_spe.rds")))
spe_ruv <- readRDS(file.path(geopath, paste0(disease, "_spe_ruv.rds")))


red_df <- as.data.frame(reducedDim(spe, "TSNE"))
red_df$section = spe$section_id
write.csv(red_df, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig1d_geo.csv")

# TSNE --------------------------------------------------------------------
plotGeoBreastFig1d <- function(spe, type = "TSNE", legend_pos = c(0.65, 0.3), 
                               annotx = 5, annoty = 2.5, 
                               plot.title = "geo_breast_nobatch_PCA_fig1d.pdf"){
  p <- plotDimRed(spe, type = type, annotate = "section_id", pt.size = 3) + 
    theme(legend.position = legend_pos,
          legend.title=element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.5, "lines")
    ) 
  p <- p + annotate("text", x = annotx, y = annoty, label = "Breast", 
                    size = 10, fontface = "bold")
  
  pdf(file = file.path(figpath, plot.title = plot.title),
      width = 5.5,
      height = 5)
  print(p)
  dev.off()
  
  return(p)
}

plotGeoBreastFig1d(spe, type = "PCA", legend_pos = c(0.9, 0.2), annotx = 22, annoty = 11, plot.title = "geo_breast_nobatch_PCA_fig1d.pdf")
plotGeoBreastFig1d(spe, type = "UMAP", legend_pos = c(0.9, 0.2), annotx = -5, annoty = 7, plot.title = "geo_breast_nobatch_UMAP_fig1d.pdf")
plotGeoBreastFig1d(spe, type = "TSNE", legend_pos = c(0.65, 0.3), annotx = 5, annoty = 2.5, plot.title = "geo_breast_nobatch_TSNE_fig1d.pdf")
