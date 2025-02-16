options(java.parameters = "-Xmx12g")
library(RBioFormats)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(SpatialOmicsOverlay)
library(GeomxTools)
library(readxl)
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/lung_spe.rds")

# names(colData(geo))
tifFile <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/mosaic_pilot/lung/__ome_tiff_geomx/%CH_L_p003_4_WTA.ome.tiff"
LS <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/L4_new_LabWorksheet.txt"

# L4 <- readSpatialOverlay(ometiff = tifFile, annots = LS,
#                          slideName = "%CH_L_p003/4", image = FALSE,
#                          saveFile = FALSE, outline = FALSE)

annots <- as.data.frame(data.table::fread(file = LS))
xml <- xmlExtraction(ometiff = tifFile, saveFile = FALSE)
# AOIattrs <- parseOverlayAttrs(omexml = xml, annots = annots, labworksheet = TRUE) # manual step in to get AOIattrs, 
# modified annotMatching() to SpatialOmicsOverlay:::annotMatching()

omexml = xml;labworksheet = TRUE
ROIs <- omexml[which(names(omexml) == "ROI")]
names(ROIs) <- paste0(names(ROIs), 0:(length(ROIs) - 1))
AOIattrs <- NULL
for (ROI in names(ROIs)) {
  ROInum <- ROIs[[ROI]]$AnnotationRef
  ROInum <- as.numeric(gsub("Annotation:", "", ROInum))
  ROInum <- omexml$StructuredAnnotations[ROInum]$XMLAnnotation$Value$ChannelThresholds$RoiName
  ROInum <- gsub("\\W", "", ROInum)
  ROI <- ROIs[[ROI]]$Union
  masks <- which(names(ROI) == "Mask")
  for (mask in masks) {
    maskNum <- which(masks == mask)
    segmentation <- ifelse(length(masks) == 1, "Geometric", 
                           "Segmented")
    mask.attrs <- ROI[[mask]]$.attrs
    if (labworksheet == TRUE & !"Text" %in% names(mask.attrs)) {
      stop("Scan was not exported on version 2.4+, please use DA annotation instead of Lab Worksheet")
    }
    else if (labworksheet == TRUE) {
      maskText <- mask.attrs[["Text"]]
    }
    else {
      maskText <- NULL
    }
    ROIannot <- SpatialOmicsOverlay:::annotMatching(annots, ROInum, maskNum, 
                                                    maskText)
    if (is.null(ROIannot)) {
      next
    }
    else if (nrow(ROIannot) == 0) {
      next
    }
    AOIattr <- as.data.frame(c(ROILabel = ROInum, ROIannot, 
                               mask.attrs[c("Height", "Width", "X", "Y")], 
                               Segmentation = segmentation))
    AOIattr$Height <- as.numeric(AOIattr$Height)
    AOIattr$Width <- as.numeric(AOIattr$Width)
    AOIattr$X <- as.numeric(AOIattr$X)
    AOIattr$Y <- as.numeric(AOIattr$Y)
    AOIattrs <- rbind(AOIattrs, cbind(AOIattr, Position = ROI[[mask]]$BinData$text))
  }
}

AOIattrs <- SpatialPosition(position = AOIattrs)
scan_metadata <- parseScanMetadata(omexml = xml)
scan_metadata[["Segmentation"]] <- "Segmented"
labWorksheet = TRUE
slideName = "%CH_L_p003/4"

L4 <- SpatialOverlay(slideName = slideName, scanMetadata = scan_metadata, 
                     overlayData = AOIattrs, 
                     workflow = list(labWorksheet = labWorksheet, outline = TRUE, scaled = FALSE), 
                     image = list(filePath = NULL, imagePointer = NULL, resolution = NULL))
L4 <- createCoordFile(overlay = L4, outline = TRUE)


L4
head(sampNames(L4))
slideName(L4)
head(meta(overlay(L4)))
head(coords(L4))

# LS_read <- readLabWorksheet(lw = LS, slideName = "%CH_L_p003/4")
L4 <- addPlottingFactor(overlay = L4, annots = annots,
                        plottingFactor = "segment")

## Gene expr --------
geo_L4_3 <- geo[, geo$section_id == "L4_3"] #qcd


countmat <- assay(geo_L4_3, "quantile")
expr_df <- data.frame(MS4A1 = countmat["MS4A1", ], # TLS B cells 
                      CXCL13 = countmat["CXCL13", ],
                      CD3D = countmat["CD3D", ],   # T cells
                      CD4 = countmat["CD4", ],
                      CD14 = countmat["CD14", ], # Macrophage
                      CD68 = countmat["CD68", ], 
                      CD163 = countmat["CD163", ], # M2
                      Sample_ID = geo_L4_3$sample_id2)

write.csv(expr_df, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS9b_geo.csv")

LS_read_new <- annots %>%
  left_join(expr_df)


## Decon
library(tidyr)
geo_level4lung <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/lung_batched_decon_long.csv")
data <- geo_level4lung %>%
  filter(section_id == "L4_3")

data_ind <- data %>%
  pivot_wider(names_from = CellType, values_from = Fraction) %>%
  select(sample, section_id, cell_fraction,
         B_cell, B_plasma_IGHA1, B_plasma_IGHG1, B_plasma_IGHG3, B_plasma_IGHM, B_plasma_IGKC, B_plasma_IGLC1, 
         T_CD4, T_CD8_exhausted, T_CTL, T_CXCL13, T_reg, TNK_dividing,
         NK,
         Macrophage, Monocyte,
         DC_1, DC_2, DC_activated, DC_pc, Granulocyte, Mast_cell) %>%
  mutate(B_cells = B_cell + B_plasma_IGHA1 + B_plasma_IGHG1 + B_plasma_IGHG3 + B_plasma_IGHM + B_plasma_IGKC + B_plasma_IGLC1,
         T_cells = T_CD4 + T_CD8_exhausted + T_CTL + T_CXCL13 + T_reg + TNK_dividing,
         Macrophage = Macrophage + Monocyte,
         Myeloid_else = DC_1 + DC_2 + DC_activated + DC_pc + Granulocyte + Mast_cell) %>%
  select(sample, section_id, cell_fraction, B_cells, T_cells, NK, Macrophage, Myeloid_else) %>%
  dplyr::rename(Sample_ID = sample)

df_save <- data_ind %>% 
  select(c(B_cells, T_cells, Macrophage, cell_fraction))

write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS9c.csv")

# viridis::scale_fill_viridis(option="inferno") + 

LS_read_new <- LS_read_new %>%
  left_join(data_ind)

# TLS B cells 
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "MS4A1")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CXCL13")
# T cells 
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD3D")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD4")
# Macrophage
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD14")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD68")
# M2
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD163")

# Decon proportion B cell, T cell, Macrophage
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "B_cells")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "T_cells")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "Macrophage")

head(plotFactors(L4))

# TLS B cells
pMS4A1 <- plotSpatialOverlay(overlay = L4, colorBy = "MS4A1", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
pCXCL13 <- plotSpatialOverlay(overlay = L4, colorBy = "CXCL13", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

# T cells
pCD3D <- plotSpatialOverlay(overlay = L4, colorBy = "CD3D", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
pCD4 <- plotSpatialOverlay(overlay = L4, colorBy = "CD4", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


# Macrophage
pCD14 <- plotSpatialOverlay(overlay = L4, colorBy = "CD14", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
pCD68 <- plotSpatialOverlay(overlay = L4, colorBy = "CD68", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
pCD163 <- plotSpatialOverlay(overlay = L4, colorBy = "CD163", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


# Decon -------------------------------------------------------------------
pBcell <- plotSpatialOverlay(overlay = L4, colorBy = "B_cells", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
pTcell <- plotSpatialOverlay(overlay = L4, colorBy = "T_cells", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
pMacrophage <- plotSpatialOverlay(overlay = L4, colorBy = "Macrophage", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


# Saved from R plot panel 
# ## Save plots ------------------------------------------------
# fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_expr/"
# # pdf(file = file.path(fig_path, "L4_MS4A1_.pdf"),
# #     width = 6,
# #     height = 20)
# # print(pMS4A1)
# # dev.off()
# 
# ggsave(
#   plot = pMS4A1,
#   width = 5000,
#   height = 3457,
#   units = "px",
#   filename = file.path(fig_path, "L4_MS4A1.png")
# )

# -------------------------------------------------------------------------
CD <- as.data.frame(colData(geo_L4_3))

p <- ggplot(data = CD, aes(x = roi_coordinate_y, y = -roi_coordinate_x,
                           # label = sample_id2)) + 
                           label = roi)) + 
  # label = cell_fraction)) +
  geom_point(size = 2) + 
  ggrepel::geom_text_repel(size = 8, color = "blue", point.padding = 0.15, #lwd = 2,
                           min.segment.length = .1, box.padding = .2, max.overlaps = 100) + 
  theme_bw()

p
# D_Islet_4, D_TME_4
# D_Islet_1, D_TME_1 (to bar plot)

library(RColorBrewer)
LS_read_new_sub <- LS_read_new %>%
  filter(ROILabel %in% c("D_TME_1", "D_Islet_1")) %>%
  select(ROILabel, cell_fraction, 
         MS4A1, CXCL13, 
         CD3D, CD4,
         CD14, CD68, CD163) %>%
  mutate(cell_fraction = factor(cell_fraction, levels = c("Malignant", "Other", "T cells", "Macrophage")),
         ROILabel = factor(ROILabel, levels = rev(c("D_TME_1", "D_Islet_1"))))

range(na.omit(c(LS_read_new$MS4A1, LS_read_new$CXCL13, 
                LS_read_new$CD3D, LS_read_new$CD4,
                LS_read_new$CD14, LS_read_new$CD68, LS_read_new$CD163
))) # 1.687300 3.446389

range(na.omit(c(LS_read_new_sub$MS4A1, LS_read_new_sub$CXCL13, 
                LS_read_new_sub$CD3D, LS_read_new_sub$CD4,
                LS_read_new_sub$CD14, LS_read_new_sub$CD68, LS_read_new_sub$CD163
                ))) # 1.687300 3.446389

# Expression heatmap -------------------------------------------------
# B cells
pMS4A1_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = MS4A1)) +
  geom_tile() +
  theme_minimal() +
  viridis::scale_fill_viridis(option="magma", limits = c(1.65, 3.5)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pMS4A1_heat

pCXCL13_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CXCL13)) +
  geom_tile() +
  theme_minimal() +
  viridis::scale_fill_viridis(option="magma", limits = c(1.65, 3.5)) +
  theme(panel.grid = element_blank(),
       axis.text = element_blank(),
        axis.title = element_blank()); pCXCL13_heat

## T cells 
pCD3D_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CD3D)) +
  geom_tile() +
  theme_minimal() +
  viridis::scale_fill_viridis(option="magma", limits = c(1.65, 3.5)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pCD3D_heat

pCD4_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CD4)) +
  geom_tile() +
  theme_minimal() +
  viridis::scale_fill_viridis(option="magma", limits = c(1.65, 3.5)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pCD4_heat

## Macro
pCD14_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CD14)) +
  geom_tile() +
  theme_minimal() +
  viridis::scale_fill_viridis(option="magma", limits = c(1.65, 3.5)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pCD14_heat

pCD68_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CD68)) +
  geom_tile() +
  theme_minimal() +
  viridis::scale_fill_viridis(option="magma", limits = c(1.65, 3.5)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pCD68_heat

pCD163_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CD163)) +
  geom_tile() +
  theme_minimal() +
  viridis::scale_fill_viridis(option="magma", limits = c(1.65, 3.5)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pCD163_heat

######
save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_expr/final"

## B cells
plot_title = "L4_MS4A1_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pMS4A1_heat)
dev.off()

plot_title = "L4_CXCL13_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCXCL13_heat)
dev.off()

# T cells
plot_title = "L4_CD3D_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCD3D_heat)
dev.off()

plot_title = "L4_CD4_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCD4_heat)
dev.off()

# Macro
plot_title = "L4_CD14_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCD14_heat)
dev.off()

plot_title = "L4_CD68_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCD68_heat)
dev.off()

plot_title = "L4_CD163_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCD163_heat)
dev.off()

# Expression grouped bar plot -------------------------------------------
LS_read_new_sub_plt <- rbind(LS_read_new_sub,
                             c("D_Islet_1", "Other",   NA, NA, NA, NA, NA, NA, NA),
                             c("D_Islet_1", "Macrophage",   NA, NA, NA, NA, NA, NA, NA),
                             c("D_TME_1",   "Malignant", NA, NA, NA, NA, NA, NA, NA)) %>%
  mutate(MS4A1 = as.numeric(MS4A1),
         CXCL13 = as.numeric(CXCL13),
         CD3D = as.numeric(CD3D),
         CD4 = as.numeric(CD4),
         CD14 = as.numeric(CD14),
         CD68 = as.numeric(CD68),
         CD163 = as.numeric(CD163)) %>%
  mutate(cell_fraction = factor(cell_fraction, levels = c("Malignant", "Other", "T cells", "Macrophage"))) %>%
  dplyr::rename(`AOI label` = cell_fraction) %>%
  mutate(ROILabel = ifelse(ROILabel == "D_Islet_1", "Islet", "TME")) %>%
  mutate(ROILabel = factor(ROILabel, levels = c("TME", "Islet"))) %>%
  dplyr::rename(`ROI label` = ROILabel)

#### 
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

plot_bar_TLS_expression <- function(LS_read_new_sub_plt, gene){
  p <- ggplot(LS_read_new_sub_plt, aes(x=`AOI label`, y=get(gene), fill=`AOI label`, pattern = `ROI label`)) +
    geom_bar_pattern(stat='identity', 
                     position='dodge2',
                     width = 0.8,
                     color = "black", 
                     pattern_fill = "black",
                     pattern_angle = 45,
                     pattern_density = 0.1,
                     pattern_spacing = 0.025,
                     pattern_key_scale_factor = 0.6) + 
    theme_minimal() + 
    theme(axis.title = element_blank()) + 
    scale_fill_manual(values = c("Malignant" = "#8185D5", "Other" = "#AA93C0", "T cells" = "#30FF00", "Macrophage" = "#FD3B40")) +
    scale_pattern_manual(values = c("TME" = "none", "Islet" = "stripe")) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))
  
  p
}
pMS4A1_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "MS4A1")
pCXCL13_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CXCL13")
pCD3D_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CD3D")
pCD4_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CD4")
pCD14_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CD14")
pCD68_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CD68")
pCD163_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CD163")

save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_expr/"
plot_title = "L4_MS4A1_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 4)
print(pMS4A1_bar)
dev.off()

plot_title = "L4_CXCL13_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 4)
print(pCXCL13_bar)
dev.off()

plot_title = "L4_CD3D_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 4)
print(pCD3D_bar)
dev.off()

plot_title = "L4_CD4_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 4)
print(pCD4_bar)
dev.off()

plot_title = "L4_CD14_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 4)
print(pCD14_bar)
dev.off()

plot_title = "L4_CD68_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 4)
print(pCD68_bar)
dev.off()

plot_title = "L4_CD163_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 4)
print(pCD163_bar)
dev.off()


