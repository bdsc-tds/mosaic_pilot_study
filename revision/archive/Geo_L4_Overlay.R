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
                      CXCR5 = countmat["CXCR5", ],
                      CCL19 = countmat["CCL19", ],
                      CD79A = countmat["CD79A", ], # B cells
                      CD79B = countmat["CD79B", ], 
                      CD3D = countmat["CD3D", ],   # T cells
                      CD3E = countmat["CD3E", ],
                      CD4 = countmat["CD4", ],
                      CD8A = countmat["CD8A", ],
                      KLRK1 = countmat["KLRK1", ], # NK
                      CD14 = countmat["CD14", ], # Macrophage
                      CD68 = countmat["CD68", ], 
                      CSF1R = countmat["CSF1R", ], 
                      ITGAM = countmat["ITGAM", ], 
                      CD80 = countmat["CD80", ],  # M1
                      CD86 = countmat["CD86", ],
                      IL1B = countmat["IL1B", ],
                      CD163 = countmat["CD163", ], # M2
                      MRC1 = countmat["MRC1", ],
                      IL10 = countmat["IL10", ],
                      Sample_ID = geo_L4_3$sample_id2)

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

# viridis::scale_fill_viridis(option="inferno") + 

LS_read_new <- LS_read_new %>%
  left_join(data_ind)

# TLS B cells 
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "MS4A1")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CXCL13")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CXCR5")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CCL19")

# Macrophage
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD14")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD68")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CSF1R")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "ITGAM")

# M1
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD80")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD86")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "IL1B")

# M2
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "CD163")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "MRC1")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "IL10")

# Decon
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "B_cells")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "T_cells")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "NK")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "Macrophage")
L4 <- addPlottingFactor(overlay = L4, annots = LS_read_new,
                        plottingFactor = "Myeloid_else")

head(plotFactors(L4))
plotSpatialOverlay(overlay = L4, colorBy = "segment", scaleBar = FALSE) 
pB_cells <- plotSpatialOverlay(overlay = L4, colorBy = "B_cells", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") 
pT_cells <- plotSpatialOverlay(overlay = L4, colorBy = "T_cells", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pNK <- plotSpatialOverlay(overlay = L4, colorBy = "NK", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pMacrophage <- plotSpatialOverlay(overlay = L4, colorBy = "Macrophage", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pMyeloid_else <- plotSpatialOverlay(overlay = L4, colorBy = "Myeloid_else", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 

# TLS B cells
pMS4A1 <- plotSpatialOverlay(overlay = L4, colorBy = "MS4A1", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") 
pCXCL13 <- plotSpatialOverlay(overlay = L4, colorBy = "CXCL13", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pCXCR5 <- plotSpatialOverlay(overlay = L4, colorBy = "CXCR5", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pCCL19 <- plotSpatialOverlay(overlay = L4, colorBy = "CCL19", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 

# Macrophage
pCD14 <- plotSpatialOverlay(overlay = L4, colorBy = "CD14", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") 
pCD68 <- plotSpatialOverlay(overlay = L4, colorBy = "CD68", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pCSF1R <- plotSpatialOverlay(overlay = L4, colorBy = "CSF1R", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pITGAM <- plotSpatialOverlay(overlay = L4, colorBy = "ITGAM", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 

# M1
pCD80 <- plotSpatialOverlay(overlay = L4, colorBy = "CD80", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") 
pCD86 <- plotSpatialOverlay(overlay = L4, colorBy = "CD86", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pIL1B <- plotSpatialOverlay(overlay = L4, colorBy = "IL1B", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 

# M2
pCD163 <- plotSpatialOverlay(overlay = L4, colorBy = "CD163", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") 
pMRC1 <- plotSpatialOverlay(overlay = L4, colorBy = "MRC1", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 
pIL10 <- plotSpatialOverlay(overlay = L4, colorBy = "IL10", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") 

# # B cells
# pCD79A <- plotSpatialOverlay(overlay = L4, colorBy = "CD79A", scaleBar = FALSE, image = FALSE) + 
#   viridis::scale_fill_viridis(option="inferno") 
# pCD79B <- plotSpatialOverlay(overlay = L4, colorBy = "CD79B", scaleBar = FALSE, image = FALSE) + 
#   viridis::scale_fill_viridis(option="inferno") 
# 
# # T cells
# pCD3D <- plotSpatialOverlay(overlay = L4, colorBy = "CD3D", scaleBar = FALSE, image = FALSE) + 
#   viridis::scale_fill_viridis(option="inferno") 
# pCD3E <- plotSpatialOverlay(overlay = L4, colorBy = "CD3E", scaleBar = FALSE, image = FALSE) + 
#   viridis::scale_fill_viridis(option="inferno") 
# pCD4 <- plotSpatialOverlay(overlay = L4, colorBy = "CD4", scaleBar = FALSE, image = FALSE) + 
#   viridis::scale_fill_viridis(option="inferno") 
# pCD8A <- plotSpatialOverlay(overlay = L4, colorBy = "CD8A", scaleBar = FALSE, image = FALSE) + 
#   viridis::scale_fill_viridis(option="inferno") 



## Save plots with transparent background
save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_gallery/"
ggsave(
  plot = pTu_L4_BCL2,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "L4_Tu.png")
)

ggsave(
  plot = pEpi_mucous,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "L4_epi_mucuous.png")
)


ggsave(
  plot = pEpi_Mucous_surface_gastric,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "L4_epi_gas_.png")
)




# Indep validating --------------------------------------------------------
geo_L4 <- geo[, geo$patient == "L4"]
CD <- as.data.frame(colData(geo_L4))

p <- ggplot(data = df, aes(x = Y, y = X,
                           #  label = Sample_ID)) + 
                           # label = cell_fraction)) +
                           label = ROILabel)) +
  geom_point(size = 2) + 
  ggrepel::geom_text_repel(size = 8, color = "blue", point.padding = 0.15, #lwd = 2,
                           min.segment.length = .1, box.padding = .2, max.overlaps = 100) + 
  theme_bw()

p

df <- meta(overlay(L4))






