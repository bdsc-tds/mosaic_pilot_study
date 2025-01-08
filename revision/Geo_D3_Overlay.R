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
tifFile <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/mosaic_pilot_2/geomx/__ome_tiff/CH_D_P001_18JUL23_bis.ome.tiff"
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

B1_3LS <- readLabWorksheet("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/B1_3_LabWorksheet.txt",
                           slideName = "%CH_BC_p002_WTA")


D3 <- readSpatialOverlay(ometiff = tifFile, annots = LS, 
                         slideName = "CH_D_P001", image = FALSE,
                         saveFile = FALSE, outline = FALSE)

D3
head(sampNames(D3))
slideName(D3)
head(meta(overlay(D3)))
head(coords(D3))

LS_read <- readLabWorksheet(lw = LS, slideName = "CH_D_P001")
D3 <- addPlottingFactor(overlay = D3, annots = LS_read,
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
geo_level4dlbcl <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/dlbcl_batched_decon_long.csv")
data <- geo_level4dlbcl %>%
  filter(patient == "D3")

data_ind <- data %>%
  pivot_wider(names_from = CellType, values_from = Fraction) %>%
  select(sample, patient, cell_fraction,
         Tu_D3_FAM3C, Tu_D3_dividing) %>%
  dplyr::rename(Sample_ID = sample)

# viridis::scale_fill_viridis(option="inferno") + 

LS_read_new <- LS_read %>%
  inner_join(data_ind)

D3 <- addPlottingFactor(overlay = D3, annots = LS_read_new,
                        plottingFactor = "Tu_D3_FAM3C")
D3 <- addPlottingFactor(overlay = D3, annots = LS_read_new,
                        plottingFactor = "Tu_D3_dividing")
D3 <- addPlottingFactor(overlay = D3, annots = LS_read_new,
                        plottingFactor = "cell_fraction")

head(plotFactors(D3))

# options(java.parameters = "-Xmx12g")
# library(RBioFormats)
# checkValidRes(ometiff = tifFile)
# 
# res <- 10
# B1 <- addImageOmeTiff(overlay = B1, ometiff = tifFile, res = res)
#
# B1 <- flipY(B1)
# B1 <- flipY(B1)
plotSpatialOverlay(overlay = D3, colorBy = "cell_fraction", scaleBar = FALSE)
pTu_D3_FAM3C <- plotSpatialOverlay(overlay = D3, colorBy = "Tu_D3_FAM3C", scaleBar = FALSE) +
  # viridis::scale_fill_viridis(option="inferno") + # all zero anyway
  scale_fill_gradient(low = "#02020C", high = "#08051D") +
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


pTu_D3_dividing <- plotSpatialOverlay(overlay = D3, colorBy = "Tu_D3_dividing", scaleBar = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_gallery/"
ggsave(
  plot = pTu_D3_FAM3C,
  width = 8000,
  height = 4069,
  units = "px",
  filename = file.path(save_path, "D3_Tu_D3_FAM3C.png")
)

ggsave(
  plot = pTu_D3_dividing,
  width = 8000,
  height = 4069,
  units = "px",
  filename = file.path(save_path, "D3_Tu_D3_dividing.png")
)





