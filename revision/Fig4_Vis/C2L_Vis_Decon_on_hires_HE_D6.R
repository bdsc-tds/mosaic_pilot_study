library(magick)
library(ggplot2)
library(dplyr)
library(scatterpie)

disease = "dlbcl"
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Visium/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/color_palette.R")
sample = "DLBCL_6"
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

scaled_width <- (sce$pxl_col_in_fullres + abs(min(sce$pxl_col_in_fullres))) * scalef$tissue_hires_scalef
shrink_factor <- (range(scaled_width)[2] - range(scaled_width)[1])/(1925-15) # 3.036904 # 3.170465
sce$pxl_col_in_hires <- scaled_width/shrink_factor/1.15 + 303

scaled_height <- (sce$pxl_row_in_fullres + abs(min(sce$pxl_row_in_fullres))) * scalef$tissue_hires_scalef
shrink_factor <- (range(scaled_height)[2] - range(scaled_height)[1])/(1880-105) # 3.265216 # 3.679116
sce$pxl_row_in_hires <- scaled_height/shrink_factor - 10


# Rotate coordinates ------------------------------------------------------
# Define the coordinates of the point
x <- sce$pxl_col_in_hires
y <- sce$pxl_row_in_hires

# Convert angle from degrees to radians
angle <- -10 * pi / 180

# Define the rotation matrix
rotation_matrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2)

# Create a vector for the coordinates
coordinates <- cbind(x, y)

# Perform the rotation
rotated_coordinates <- rotation_matrix %*% t(coordinates)

# # Print the rotated coordinates
# print(rotated_coordinates)
sce$pxl_col_in_hires <- rotated_coordinates[1, ]
sce$pxl_row_in_hires <- rotated_coordinates[2, ]

sce$Region <- as.factor(sce$Region)


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
    image,
    0,
    info$width,
    0,
    info$height
  ) +
  theme_bw()

decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Visium_Decon/Level4/C2L"
res_C2L <- read.csv(file.path(decon_path, paste0(sample, ".csv")), row.names = 1)
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
radius <- scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 1.7 / 4

p_decon_full <- p0 + 
  geom_scatterpie(aes(x = x, y = y, r = radius), 
                  data = data, cols = ct.select, color = NA) + # coord_fixed(ratio = 1*max(data$x)/max(data$y)) + 
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

if(sample == "B1_4"){
  image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
  image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
}else{
  image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
  image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
}

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) - 2

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Vis_D6"
pdf(file = file.path(fig_path, paste0(sample, "_decon_HE_C2L.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_full)
dev.off()


########################################################################
# Decon_HE - sub -------------------------------------------
p_decon_sub <- p_decon_full + xlim(c(890, 1540)) + ylim(c(1190, 1840))

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Vis_D6"
pdf(file = file.path(fig_path, paste0(sample, "_decon_HE_sub_C2L.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_sub)
dev.off()


# ################################ Legend ################################
# data_legend <- data %>%
#   select(Epi_mucous, Epi_bronchus, Epi_Mucous_surface_gastric, Epi_parietal_gastric, Tu_D6_BCL2, x, y)
# 
# p_decon_legend_sub <- p0 + 
#   geom_scatterpie(aes(x = x, y = y, r = radius), 
#                   data = data_legend, cols = colnames(data_legend)[!colnames(data_legend) %in% c("x", "y")], color = NA) + 
#   scale_fill_manual(values =  colors) +
#   theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
#         panel.background = element_blank(),
#         plot.background = element_blank(),
#         panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
#         axis.text =element_blank(),
#         axis.ticks =element_blank(),
#         axis.title =element_blank(),
#         legend.title=element_text(size = 16,face="bold"),
#         legend.text=element_text(size = 15),
#         legend.key = element_rect(colour = "transparent", fill = "white"),
#         legend.key.size = unit(0.45, 'cm'),
#         strip.text = element_text(size = 16,face="bold"),
#         legend.position="right")+
#   guides(fill=guide_legend(title="Cell Type", ncol = 1)) 
# 
# fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Vis_D6"
# pdf(file = file.path(legend_path, paste0(sample, "_decon_legend_sub.pdf")),
#     width = 8,
#     height = 8)
# print(p_decon_legend_sub)
# dev.off()


############################### Decon ind ##############################
# Decon_HE ind - Tu_D6_BCL2 ---------------------------------------------------
p_decon_ind_Tu_D6_BCL2 <- p0 + xlim(c(890, 1540)) + ylim(c(1190, 1840)) +
  ggforce::geom_circle(
    data = data,
    aes(
      x0 = x,
      y0 = y,
      r = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2 / 4,
      fill = Tu_D6_BCL2, 
      colour = Tu_D6_BCL2 
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

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Vis_D6"
pdf(file = file.path(fig_path, paste0(sample, "_decon_ind_HE_Tu_D6_BCL2_C2L.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_ind_Tu_D6_BCL2)
dev.off()


# Decon_HE ind - Epi_mucous ---------------------------------------------------
p_decon_ind_Epi_mucous <- p0 + xlim(c(890, 1540)) + ylim(c(1190, 1840)) +
  ggforce::geom_circle(
    data = data,
    aes(
      x0 = x,
      y0 = y,
      r = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2 / 4,
      fill = Epi_mucous, 
      colour = Epi_mucous 
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

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Vis_D6"
pdf(file = file.path(fig_path, paste0(sample, "_decon_ind_HE_Epi_mucous_C2L.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_ind_Epi_mucous)
dev.off()


# Decon_HE ind - Epi_Mucous_surface_gastric ---------------------------------------------------
p_decon_ind_Epi_Mucous_surface_gastric <- p0 + xlim(c(890, 1540)) + ylim(c(1190, 1840)) +
  ggforce::geom_circle(
    data = data,
    aes(
      x0 = x,
      y0 = y,
      r = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2 / 4,
      fill = Epi_Mucous_surface_gastric, 
      colour = Epi_Mucous_surface_gastric 
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

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Vis_D6"
pdf(file = file.path(fig_path, paste0(sample, "_decon_ind_HE_Epi_Mucous_surface_gastric_C2L.pdf")),
    width = plot_width,
    height = plot_height)
print(p_decon_ind_Epi_Mucous_surface_gastric)
dev.off()
