library(magick)
library(ggplot2)
library(dplyr)
library(scatterpie)
library(stringr)

disease = "breast"
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")
sample = "B1_4"
if(disease == "lung"){patho_color <- lung_patho_color}else{patho_color <- breast_patho_color}

sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
sce <- readRDS(file.path(sce_path, paste0(sample, "_qcd.rds")))

datapath_rerun <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium/VisiumRerun/", sample)
coord_new <- read.csv(file.path(datapath_rerun, "tissue_positions.csv"))

CD <- data.frame(colData(sce))

CD_new <- coord_new %>% 
  mutate(Barcode = barcode) %>% 
  select(Barcode, pxl_col_in_fullres, pxl_row_in_fullres)

CD <- CD %>%
  left_join(CD_new, by = "Barcode")
rownames(CD) <- CD$Barcode

colData(sce) <- as(CD, "DFrame")

scalef <- rjson::fromJSON(file = file.path(datapath_rerun, "scalefactors_json.json"))
image <- image_read(file.path(datapath_rerun, "tissue_hires_image.png"))
info <- image_info(image)

sce$pxl_col_in_hires <- sce$pxl_col_in_fullres * scalef$tissue_hires_scalef
sce$pxl_row_in_hires <- sce$pxl_row_in_fullres * scalef$tissue_hires_scalef

sce$Region <- as.factor(sce$Region)

# C2L
# -------------------------------------------------------------------------
decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Visium_Decon/Level4/C2L"
res_C2L <- read.csv(file.path(decon_path, paste0(sample, ".csv")), row.names = 1)
res_C2L <- res_C2L %>%
  mutate(Fibroblast = Fibroblast + Fibroblast_B3) %>%
  select(-Fibroblast_B3)

location <- data.frame(x = sce$pxl_col_in_hires,
                       y = sce$pxl_row_in_hires)

rownames(location) <- colnames(sce)
colnames(location) <- c("x", "y")

overlap_barcodes <- intersect(rownames(res_C2L), rownames(location))
location <- location[overlap_barcodes, ]
res_C2L <- res_C2L[overlap_barcodes, ]
data = cbind(res_C2L, location)
ct.select = colnames(res_C2L)

colors = level4_cellcolors[names(level4_cellcolors) %in% sort(colnames(res_C2L))]
radius <- scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 1.7


########################### Decon HE full slide #########################
p0 <- ggplot(
  data.frame(x = 0, y = 0),
  aes(x, y)
) +
  geom_blank() +
  coord_fixed(
    expand = FALSE,
    xlim = c(min(sce$pxl_col_in_hires) - info$width * 0.05, max(sce$pxl_col_in_hires) + info$width * 0.05),
    ylim = c(min(sce$pxl_row_in_hires) - info$height * 0.05, max(sce$pxl_row_in_hires) + info$height * 0.05)
  ) +
  annotation_raster(
    image_flip(image),
    0,
    info$width,
    0,
    info$height
  ) +
  theme_bw()

p_decon_full <- p0 + 
  geom_scatterpie(aes(x = x, y = y, r = radius), 
                  data = data, cols = ct.select, color = NA) + 
  scale_fill_manual(values =  colors) +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 15),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm'),
        strip.text = element_text(size = 16,face="bold"),
        legend.position="right")+
  guides(fill=guide_legend(title="Cell Type", ncol = 1)) 

# Image size and save ---------
if(sample == "B1_4"){
  image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
  image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
}else{
  image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
  image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
}

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) - 2

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_full_slide_HE"
pdf(file = file.path(fig_path, paste0(sample, "_decon_HE_C2L.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_full)
dev.off()


########################################################################
# Decon_HE - B cell -------------------------------------------
p_decon_Bcells <- p_decon_full + xlim(c(480, 606)) + ylim(c(458, 626))

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_cropped"
pdf(file = file.path(fig_path, paste0(sample, "_decon_HE_Bcell_C2L.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_Bcells)
dev.off()


# Decon_HE - T cell -----------------
p_decon_Tcells <- p_decon_full + xlim(c(755, 920)) + ylim(c(100, 263))

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_cropped"
pdf(file = file.path(fig_path, paste0(sample, "_decon_HE_Tcell_C2L.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_Tcells)
dev.off()


################################ Legend ################################
data_legend <- data %>%
  select(B_plasma_IGHA1, B_plasma_IGHG1, B_plasma_IGHG3, B_plasma_IGHM, B_plasma_IGKC, B_plasma_IGLC1, 
         T_CD4, T_CD8_exhausted, T_CTL, T_CXCL13, T_reg, TNK_dividing, x, y)
  
p_decon_legend_sub <- p0 + 
  geom_scatterpie(aes(x = x, y = y, r = radius), 
                  data = data_legend, cols = colnames(data_legend)[!colnames(data_legend) %in% c("x", "y")], color = NA) + 
  scale_fill_manual(values =  colors) +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 15),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm'),
        strip.text = element_text(size = 16,face="bold"),
        legend.position="right")+
  guides(fill=guide_legend(title="Cell Type", ncol = 1)) 
  
legend_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_legend"
pdf(file = file.path(legend_path, paste0(sample, "_decon_legend_sub.pdf")),
    width = 8,
    height = 8)
print(p_decon_legend_sub)
dev.off()
  
############################### Decon ind ##############################
data_ind <- data %>%
  select(B_cell, B_plasma_IGHA1, B_plasma_IGHG1, B_plasma_IGHG3, B_plasma_IGHM, B_plasma_IGKC, B_plasma_IGLC1, 
         T_CD4, T_CD8_exhausted, T_CTL, T_CXCL13, T_reg, TNK_dividing,
         x, y) %>%
  mutate(B_cells = B_cell + B_plasma_IGHA1 + B_plasma_IGHG1 + B_plasma_IGHG3 + B_plasma_IGHM + B_plasma_IGKC + B_plasma_IGLC1,
         T_cells = T_CD4 + T_CD8_exhausted + T_CTL + T_CXCL13 + T_reg + TNK_dividing)

# Decon_HE ind - B cell ---------------------------------------------------
p_decon_ind_Bcells <- p0 + # xlim(c(480, 606)) + ylim(c(458, 626)) + 
  ggforce::geom_circle(
    data = data_ind,
    aes(
      x0 = x,
      y0 = y,
      r = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2,
      fill = B_cells, 
      colour = B_cells 
    ),
    inherit.aes = FALSE
  ) +
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right",
        legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size = 8)) +
  viridis::scale_fill_viridis(option="inferno") + 
  viridis::scale_color_viridis(option="inferno") +
  guides(fill = guide_colorbar(title.position = "left", title.hjust = 0.5, direction = "horizontal",
                               barwidth = unit(2.5, "cm"), barheight = unit(0.3, "cm")))

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_cropped"
# pdf(file = file.path(fig_path, paste0(sample, "_decon_ind_HE_Bcell_C2L.pdf")),
pdf(file = file.path(fig_path, paste0(sample, "_decon_ind_HE_Bcell_C2L_fullslide.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_ind_Bcells)
dev.off()



# Decon_HE ind - T cell ---------------
p_decon_ind_Tcells <- p0 + # xlim(c(755, 920)) + ylim(c(100, 263)) + 
  ggforce::geom_circle(
    data = data_ind,
    aes(
      x0 = x,
      y0 = y,
      r = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2,
      fill = T_cells, 
      colour = T_cells 
    ),
    inherit.aes = FALSE
  ) +
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right",
        legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size = 8)) +
  viridis::scale_fill_viridis(option="inferno") + 
  viridis::scale_color_viridis(option="inferno") +
  guides(fill = guide_colorbar(title.position = "left", title.hjust = 0.5, direction = "horizontal",
                               barwidth = unit(2.5, "cm"), barheight = unit(0.3, "cm")))

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_cropped"
# pdf(file = file.path(fig_path, paste0(sample, "_decon_ind_HE_Tcell_C2L.pdf")),
pdf(file = file.path(fig_path, paste0(sample, "_decon_ind_HE_Tcell_C2L_fullslide.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_ind_Tcells)
dev.off()






