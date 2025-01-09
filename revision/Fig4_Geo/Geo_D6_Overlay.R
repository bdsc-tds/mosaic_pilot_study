options(java.parameters = "-Xmx12g")
library(RBioFormats)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(SpatialOmicsOverlay)
library(GeomxTools)
library(readxl)
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Final_Qcd/dlbcl_qcd.rds")

# names(colData(geo))
tifFile <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/mosaic_pilot_2/geomx/__ome_tiff/CH_D_P002_19JUL23.ome.tiff"
# LSxlsx <- read_xlsx("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DLBCL_LabWorkSheet_owkin_ED.xlsx")
# table(LSxlsx$Patient, LSxlsx$`slide name`)
#    CH_D_P001 CH_D_P001,2 CH_D_P002
# D1         0          24         0
# D2        24           0         0
# D3        24           0         0
# D4         0           0        24
# D5         0           0        28
# D6         0           0        20
# "CH_D_P001"
LS <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DLBCL_owkin_ED_D3_D6_LabWorksheet.txt"


D6 <- readSpatialOverlay(ometiff = tifFile, annots = LS, 
                         slideName = "CH_D_P002", image = FALSE,
                         saveFile = FALSE, outline = FALSE)

D6
head(sampNames(D6))
slideName(D6)
head(meta(overlay(D6)))
head(coords(D6))

LS_read <- readLabWorksheet(lw = LS, slideName = "CH_D_P002")
D6 <- addPlottingFactor(overlay = D6, annots = LS_read,
                        plottingFactor = "segment")

# ## Gene expr --------
# geo_B1_3 <- geo[, geo$section_id == "B1_3"] 
# 
# countmat <- assay(geo_B1_3, "counts")
# expr_df <- data.frame(ACADM = countmat["ACADM", ],
#                       Sample_ID = colnames(geo_B1_3))
# 
# LS_read_new <- LS_read %>%
#   inner_join(expr_df)

## Decon
library(tidyr)
library(dplyr)
geo_level4dlbcl <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/dlbcl_batched_decon_long.csv")
data <- geo_level4dlbcl %>%
  filter(patient == "D6")

data_ind <- data %>%
  pivot_wider(names_from = CellType, values_from = Fraction) %>%
  select(sample, patient, cell_fraction,
         Tu_D6_BCL2, Epi_mucous, Epi_Mucous_surface_gastric) %>%
  dplyr::rename(Sample_ID = sample)

# viridis::scale_fill_viridis(option="inferno") + 

LS_read_new <- LS_read %>%
  left_join(data_ind)

D6 <- addPlottingFactor(overlay = D6, annots = LS_read_new,
                        plottingFactor = "Tu_D6_BCL2")
D6 <- addPlottingFactor(overlay = D6, annots = LS_read_new,
                        plottingFactor = "Epi_mucous")
D6 <- addPlottingFactor(overlay = D6, annots = LS_read_new,
                        plottingFactor = "Epi_Mucous_surface_gastric")
D6 <- addPlottingFactor(overlay = D6, annots = LS_read_new,
                        plottingFactor = "cell_fraction")

head(plotFactors(D6))

# options(java.parameters = "-Xmx12g")
# library(RBioFormats)
# checkValidRes(ometiff = tifFile)
# 
# res <- 10
# B1 <- addImageOmeTiff(overlay = B1, ometiff = tifFile, res = res)
#
# B1 <- flipY(B1)
# B1 <- flipY(B1)
plotSpatialOverlay(overlay = D6, colorBy = "segment", scaleBar = FALSE)
pTu_D6_BCL2 <- plotSpatialOverlay(overlay = D6, colorBy = "Tu_D6_BCL2", scaleBar = FALSE, image = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
pEpi_mucous <- plotSpatialOverlay(overlay = D6, colorBy = "Epi_mucous", scaleBar = FALSE, image = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
pEpi_Mucous_surface_gastric <- plotSpatialOverlay(overlay = D6, colorBy = "Epi_Mucous_surface_gastric", scaleBar = FALSE, image = FALSE) +  # , alpha = 1.5
  viridis::scale_fill_viridis(option="inferno") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

## Save plots with transparent background
save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_gallery/"
ggsave(
  plot = pTu_D6_BCL2,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "D6_Tu.png")
)

ggsave(
  plot = pEpi_mucous,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "D6_epi_mucuous.png")
)


ggsave(
  plot = pEpi_Mucous_surface_gastric,
  width = 2000,
  height = 1739,
  units = "px",
  filename = file.path(save_path, "D6_epi_gas_.png")
)




# Indep validating --------------------------------------------------------
geo_D6 <- geo[, geo$patient == "D6"]
CD <- as.data.frame(colData(geo_D6))

p <- ggplot(data = df, aes(x = Y, y = X,
                          #  label = Sample_ID)) + 
  # label = cell_fraction)) +
  label = ROILabel)) +
  geom_point(size = 2) + 
  ggrepel::geom_text_repel(size = 8, color = "blue", point.padding = 0.15, #lwd = 2,
                           min.segment.length = .1, box.padding = .2, max.overlaps = 100) + 
  theme_bw()

p

df <- meta(overlay(D6))






