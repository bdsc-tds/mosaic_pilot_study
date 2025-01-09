options(java.parameters = "-Xmx12g")
library(RBioFormats)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(SpatialOmicsOverlay)
library(GeomxTools)
library(readxl)
# geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Final_Qcd/lung_qcd.rds")
# Do not QC
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/GeoMx/lung_raw/L1_1_0PSV.rds")


# names(colData(geo))
tifFile <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/mosaic_pilot/lung/__ome_tiff_geomx/_CH_L_p001_WTA_110523.ome.tiff"
# LSxlsx <- read_xlsx("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/lung_LabWorkSheet_owkin_ED.xlsx")
# table(LSxlsx$Patient, LSxlsx$`slide name`)
#    CH_D_P001 CH_D_P001,2 CH_D_P002
# D1         0          24         0
# D2        24           0         0
# D3        24           0         0
# D4         0           0        24
# D5         0           0        28
# L1         0           0        20
# "CH_D_P001"
LS <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/L1_1_new_LabWorksheet.txt"


# L1 <- readSpatialOverlay(ometiff = tifFile, annots = LS, 
#                          slideName = "CH_L_p001_WTA", image = FALSE,
#                          saveFile = FALSE, outline = FALSE)

annots <- as.data.frame(data.table::fread(file = LS))
xml <- xmlExtraction(ometiff = tifFile, saveFile = FALSE)
AOIattrs <- parseOverlayAttrs(omexml = xml, annots = annots, labworksheet = TRUE)
scan_metadata <- parseScanMetadata(omexml = xml)
scan_metadata[["Segmentation"]] <- "Segmented"
labWorksheet = TRUE
slideName = "CH_L_p001_WTA"

L1 <- SpatialOverlay(slideName = slideName, scanMetadata = scan_metadata, 
                     overlayData = AOIattrs, 
                     workflow = list(labWorksheet = labWorksheet, outline = TRUE, scaled = FALSE), 
                     image = list(filePath = NULL, imagePointer = NULL, resolution = NULL))
L1 <- createCoordFile(overlay = L1, outline = TRUE)


head(sampNames(L1))
slideName(L1)
head(meta(overlay(L1)))
head(coords(L1))

# LS_read <- readLabWorksheet(lw = LS, slideName = "CH_L_p001_WTA")
L1 <- addPlottingFactor(overlay = L1, annots = annots,
                        plottingFactor = "segment")

## Gene expr --------
# geo_L1_1 <- geo[, geo$section_id == "L1_1"] #qcd
geo_L1_1 <- geo

countmat <- assay(geo_L1_1, "counts")
# logcountmat <- assay(geo_L1_1, "logcounts") # logCPM
expr_df <- data.frame(MS4A1 = countmat["MS4A1", ],
                      logMS4A1 = log(countmat["MS4A1", ] + 1),
                      CXCL13 = countmat["CXCL13", ],
                      logCXCL13 = log(countmat["CXCL13", ] + 1),
                      CXCR5 = countmat["CXCR5", ],
                      logCXCR5 = log(countmat["CXCR5", ] + 1),
                      CCL19 = countmat["CCL19", ],
                      logCCL19 = log(countmat["CCL19", ] + 1),
                      Sample_ID = geo_L1_1$SampleID)

LS_read_new <- annots %>%
  left_join(expr_df)

# ## Decon
# library(tidyr)
# library(dplyr)
# geo_level4lung <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/lung_batched_decon_long.csv")
# data <- geo_level4lung %>%
#   filter(patient == "L1")
# 
# data_ind <- data %>%
#   pivot_wider(names_from = CellType, values_from = Fraction) %>%
#   select(sample, patient, cell_fraction,
#          Tu_L1_BCL2, Epi_mucous, Epi_Mucous_surface_gastric) %>%
#   dplyr::rename(Sample_ID = sample)
# 
# # viridis::scale_fill_viridis(option="inferno") + 
# 
# LS_read_new <- LS_read %>%
#   left_join(data_ind)

L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "logMS4A1")
L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "logCXCL13")
L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "logCXCR5")
L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "logCCL19")

# head(plotFactors(L1))
# df <- plotFactors(L1) # Even more dropped at plotting factor
# table(df$segment)
# CD3 Other PanCK 
# 2     7     7 
# table(geo_L1_1$segment)  # why so much 6 aois QCed out?
# CD3 Other PanCK 
# 3     7     8 

# table(annots$segment)
# CD3 Other PanCK 
# 7    10     7 

# options(java.parameters = "-Xmx12g")
# library(RBioFormats)
# checkValidRes(ometiff = tifFile)
# 
# res <- 10
# B1 <- addImageOmeTiff(overlay = B1, ometiff = tifFile, res = res)
#
# B1 <- flipY(B1)
# B1 <- flipY(B1)
plotSpatialOverlay(overlay = L1, colorBy = "segment", scaleBar = FALSE)
plogMS4A1 <- plotSpatialOverlay(overlay = L1, colorBy = "logMS4A1", scaleBar = FALSE) +
  viridis::scale_fill_viridis(option="magma") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


plogCXCL13 <- plotSpatialOverlay(overlay = L1, colorBy = "logCXCL13", scaleBar = FALSE) + 
  viridis::scale_fill_viridis(option="magma") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

plogCXCR5 <- plotSpatialOverlay(overlay = L1, colorBy = "logCXCR5", scaleBar = FALSE) + 
  viridis::scale_fill_viridis(option="magma") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

plogCCL19 <- plotSpatialOverlay(overlay = L1, colorBy = "logCCL19", scaleBar = FALSE) + 
  viridis::scale_fill_viridis(option="magma") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

## Save plots with transparent background
save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_expr/"
ggsave(
  plot = plogMS4A1,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "L1_logMS4A1.png")
)

ggsave(
  plot = plogCXCL13,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "L1_logCXCL13.png")
)


ggsave(
  plot = plogCXCR5,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "L1_logCXCR5.png")
)

ggsave(
  plot = plogCCL19,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "L1_logCCL19.png")
)

