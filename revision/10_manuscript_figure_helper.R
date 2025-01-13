#######################################################################
#                                 Breast                              #
#######################################################################

# PCA ---------------------------------------------------------------------
p <- plotDimRed(spe, type = "PCA", annotate = "section_id", pt.size = 3) + 
  theme(legend.position = c(0.9, 0.2),
        legend.title=element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.5, "lines")
        ) 
p <- p + annotate("text", x = 22, y = 11, label = "Breast", 
                  size = 10, fontface = "bold")

pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1d/", "geo_breast_nobatch_PCA_fig1c.pdf"),
    width = 5.5,
    height = 5)
print(p)
dev.off()



# UMAP --------------------------------------------------------------------
p <- plotDimRed(spe, type = "UMAP", annotate = "section_id", pt.size = 3) + 
  theme(legend.position = c(0.9, 0.2),
        legend.title=element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.5, "lines")
  ) 
p <- p + annotate("text", x = -5, y = 7, label = "Breast", 
                  size = 10, fontface = "bold")

pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1d/", "geo_breast_nobatch_UMAP_fig1c.pdf"),
    width = 5.5,
    height = 5)
print(p)
dev.off()

# TSNE --------------------------------------------------------------------
p <- plotDimRed(spe, type = "TSNE", annotate = "section_id", pt.size = 3) + 
  theme(legend.position = c(0.65, 0.3),
        legend.title=element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.5, "lines")
  ) 
p <- p + annotate("text", x = 5, y = 2.5, label = "Breast", 
                  size = 10, fontface = "bold")

pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1d/", "geo_breast_nobatch_TSNE_fig1c.pdf"),
    width = 5.5,
    height = 5)
print(p)
dev.off()

saveRDS(spe, file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", paste0(disease, "_spe.rds")))
saveRDS(spe_ruv, file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", paste0(disease, "_spe_ruv.rds")))


#######################################################################
#                             Lung & DLBCL                            #
#######################################################################
disease = "dlbcl"
spe <- readRDS(file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", paste0(disease, "_spe.rds")))
spe_ruv <- readRDS(file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", paste0(disease, "_spe_ruv.rds")))

plotDimRed_publication <- function(spe, reddim = "PCA", annotate = "section_id", legen.lab = "Section"){
  p <- plotDimRed(spe, type = reddim, annotate = annotate, pt.size = 1.5) + 
    theme(legend.position = "bottom",
          legend.title= element_text(size = 10),
          legend.key.size = unit(0.8, "lines")) +
    labs(color = legen.lab)
  p
}

# Lung/DLBCL - PCA/UMAP/TSNE (SuppFig GeoMx b)
plot_save_reddim <- function(reddim = "PCA", plot_title = "geo_lung_nobatch_PCA.pdf"){
  # p1 <- plotDimRed_publication(spe, reddim = reddim, annotate = "section_id", legen.lab = "Section") + # lung
  p1 <- plotDimRed_publication(spe, reddim = reddim, annotate = "patient", legen.lab = "Patient") + # dlbcl
    theme(plot.margin = margin(5.5, 35, 5.5, 5.5)) |
    plotDimRed_publication(spe, reddim = reddim, annotate = "cell_fraction", legen.lab = "Cell fraction") +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 35))
  
  pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1d/", plot_title),
      width = 10.5,
      height = 5)
  print(p1)
  dev.off()
}

plot_save_reddim(reddim = "PCA", plot_title = paste0("geo_", disease, "_nobatch_PCA.pdf"))
plot_save_reddim(reddim = "UMAP", plot_title = paste0("geo_", disease, "_nobatch_UMAP.pdf"))
plot_save_reddim(reddim = "TSNE", plot_title = paste0("geo_", disease, "_nobatch_TSNE.pdf"))


# Do batch correction -----------------------------------------------------
# saveRDS(spe, file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", paste0(disease, "_spe.rds")))
# saveRDS(spe_ruv, file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", paste0(disease, "_spe_ruv.rds")))

# -------------------------------------------------------------------------
# Lung - slide before & after batch + cell_fraction before & after batch (SuppFig GeoMx c)

plot_save_reddim_batch <- function(reddim = "PCA", plot_title = "geo_lung_before_after_batch_PCA.pdf"){
  p1 <- plotDimRed_publication(spe, reddim = reddim, annotate = "slide_name", legen.lab = "Slide") +
    theme(plot.margin = margin(5.5, 35, 5.5, 5.5)) |
    plotDimRed_publication(spe_ruv, reddim = reddim, annotate = "slide_name", legen.lab = "Slide") +
    theme(plot.margin = margin(5.5, 35, 5.5, 35)) |
    plotDimRed_publication(spe_ruv, reddim = reddim, annotate = "cell_fraction", legen.lab = "Cell fraction") +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 35))
  
  pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1d/", plot_title),
      width = 16,
      height = 5)
  print(p1)
  dev.off()
}

plot_save_reddim_batch(reddim = "PCA", plot_title = paste0("geo_", disease, "_before_after_batch_PCA.pdf"))
plot_save_reddim_batch(reddim = "UMAP", plot_title = paste0("geo_", disease, "_before_after_batch_UMAP.pdf"))
plot_save_reddim_batch(reddim = "TSNE", plot_title = paste0("geo_", disease, "_before_after_batch_TSNE.pdf"))


##########################################################################
#                            Fig 2b                                      #
##########################################################################
geo_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results/"
breast_decon <- read.csv(file.path(geo_decon_path, "breast_batched_decon_long.csv"))
lung_decon <- read.csv(file.path(geo_decon_path, "lung_batched_decon_long.csv"))
dlbcl_decon <- read.csv(file.path(geo_decon_path, "dlbcl_batched_decon_long.csv"))

# Individuals
# Plot
# p <- ggplot(breast_decon, aes(x=CellType, y=Fraction)) +
# p <- ggplot(lung_decon, aes(x=CellType, y=Fraction)) +
p <- ggplot(dlbcl_decon, aes(x=CellType, y=Fraction)) +
  geom_boxplot() +
  # geom_jitter(aes(color=section_id), size=1, alpha=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~cell_fraction, ncol = 5) # + 
  # ylim(0.0, 1.0) + 
  # labs(color = "Section") +
  # ggtitle("Quantile Normed Batch Corrected Expr Decon Result - full genes")
  # ggtitle("Quantile Normed Batch Corrected Expr Decon Result - Breast & Lung")
p

# plot_title = "Geo_breast_decon.pdf"
# plot_title = "Geo_lung_decon.pdf"
plot_title = "Geo_dlbcl_decon.pdf"
pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig2b/", plot_title),
    width = 10.5,
    height = 5)
print(p)
dev.off()



# # Breast and Lung 
breast_lung_decon <- rbind(breast_decon, lung_decon)

# Plot
p <- ggplot(breast_lung_decon, aes(x=CellType, y=Fraction)) +
  geom_boxplot() +
  # geom_jitter(aes(color=section_id), size=1, alpha=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~cell_fraction, ncol = 5) # + 
  # ylim(0.0, 1.0) + 
  # labs(color = "Section") +
  # ggtitle("Quantile Normed Batch Corrected Expr Decon Result - full genes")
  # ggtitle("Quantile Normed Batch Corrected Expr Decon Result - Breast & Lung")
p

plot_title = "Geo_breast_lung_decon.pdf"
pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig2b/", plot_title),
    width = 10.5,
    height = 5)
print(p)
dev.off()



# Breast, Lung, and DLBCL 
dlbcl_decon <- dlbcl_decon %>% rename(section_id = patient)
breast_lung_dlbcl_decon <- rbind(breast_decon, lung_decon, dlbcl_decon)

# Plot
p <- ggplot(breast_lung_dlbcl_decon, aes(x=CellType, y=Fraction)) +
  geom_boxplot() +
  # geom_jitter(aes(color=section_id), size=1, alpha=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~cell_fraction, ncol = 5) # + 
  # ylim(0.0, 1.0) + 
  # labs(color = "Section") +
  # ggtitle("Quantile Normed Batch Corrected Expr Decon Result - full genes")
  # ggtitle("Quantile Normed Batch Corrected Expr Decon Result - Breast & Lung & DLBCL")
p

plot_title = "Geo_breast_lung_dlbcl_decon.pdf"
pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig2b/", plot_title),
    width = 10.5,
    height = 5)
print(p)
dev.off()








