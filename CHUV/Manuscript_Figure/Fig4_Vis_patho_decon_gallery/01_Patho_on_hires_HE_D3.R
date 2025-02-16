library(magick)
library(ggplot2)
library(dplyr)
library(ggspavis)

disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/Visium/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/color_palette.R")
sample = "DLBCL_3"
if(disease == "lung"){
  patho_color <- lung_patho_color
}else if(disease == "breast"){
  patho_color <- breast_patho_color
}else{
  patho_color <- dlbcl_patho_color
}

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

write.csv(data.frame(colData(sce)), "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig4h_patho_vis.csv")

# -------------------------------------------------------------------------
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


# -------------------------------------------------------------------------
p <- p0 + 
  ggforce::geom_circle(
    data = data.frame(colData(sce)),
    aes(
      x0 = pxl_col_in_hires,
      y0 = pxl_row_in_hires,
      r = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2,
      fill = Region,
      colour = Region,
    ),
    inherit.aes = FALSE
  ) +
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right",
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 15)) +
  scale_fill_manual(name = "Pathology", values = patho_color[names(table(sce$Region))],
                    guide = guide_legend(override.aes = list(shape = 19)),
                    na.value = "#d3d3d3") +
  scale_color_manual(values = patho_color[names(table(sce$Region))], guide = FALSE)+
  guides(fill = guide_legend(override.aes = list(color = NULL)))
p

if(sample == "B1_4"){
  image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
  image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
}else{
  image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
  image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
}

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) - 2

# figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/D3_full_slide_HE"
# pdf(file = file.path(figpath, paste0(sample, "_patho_HE.pdf")),
#     width = plot_width,
#     height = plot_height)
# print(p)
# dev.off()


# Check coord -------------------------------------------------------------
p0


# Patho ----------------------------------------------------------
p_patho_HE_sub <- p + xlim(c(780, 1060)) + ylim(c(975, 1160))

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/D3_cropped"
pdf(file = file.path(fig_path, paste0(sample, "_patho_HE_sub.pdf")),
    width = plot_width,
    height = plot_height)
print(p_patho_HE_sub)
dev.off()


# H&E ----------------------------------------------------------------------
p_HE_sub <- ggplot(
  data.frame(x = 0, y = 0),
  aes(x, y)
) +
  geom_blank() +
  coord_fixed(
    expand = FALSE,
    xlim = c(780, 1060),
    ylim = c(975, 1160)
  ) +
  annotation_raster(
    image_flip(image),
    0,
    info$width,
    0,
    info$height
  ) +
  theme_bw() 


fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/D3_cropped/archive"
pdf(file = file.path(fig_path, paste0(sample, "_HE_sub.pdf")),
    width = plot_width,
    height = plot_height)
print(p_HE_sub)
dev.off()



################################ Legend ################################
legend_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/D3_legend"
plot_width <- 8
plot_height <- 8


# Patho sub classes -------------------------------------------------------
sce_sub <- sce[, sce$Region %in% c("Small lymphocytes", "Tumor", "Vessels")]
sce_sub$pxl_col_in_fullres <- NULL
sce_sub$pxl_row_in_fullres <- NULL
p_patho_sub <- plotSpots(sce_sub, annotate = "Region", palette = patho_color[names(table(sce_sub$Region))]) + 
  theme(legend.key.size = unit(0.8, 'lines')) 

pdf(file = file.path(legend_path, paste0(sample, "_patho_sub.pdf")),
    width = plot_width,
    height = plot_height)
print(p_patho_sub)
dev.off()

