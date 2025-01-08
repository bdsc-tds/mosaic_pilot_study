options(java.parameters = "-Xmx12g")
library(RBioFormats)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(SpatialOmicsOverlay)
library(GeomxTools)
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Final_Qcd/breast_qcd.rds")

# names(colData(geo))
tifFile <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/mosaic_pilot/breast/__ome_tiff_geomx/%CH_BC_p001_WTA_150523.ome.tiff"
LS <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/B1_3_LabWorksheet.txt"


B1 <- readSpatialOverlay(ometiff = tifFile, annots = LS, 
                         slideName = "%CH_BC_p002_WTA", image = FALSE,  
                         saveFile = FALSE, outline = FALSE)

B1
head(sampNames(B1))
slideName(B1)
head(meta(overlay(B1)))
head(coords(B1))

LS_read <- readLabWorksheet(lw = LS, slideName = "%CH_BC_p002_WTA")
LS_read$roi_sampleid <- paste0(LS_read$roi, LS_read$Sample_ID)
B1 <- addPlottingFactor(overlay = B1, annots = LS_read,
                        plottingFactor = "segment")
B1 <- addPlottingFactor(overlay = B1, annots = LS_read,
                        plottingFactor = "Sample_ID")
B1 <- addPlottingFactor(overlay = B1, annots = LS_read,
                        plottingFactor = "roi_sampleid")

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
geo_level4breast <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/breast_batched_decon_long.csv")
data <- geo_level4breast %>%
  filter(section_id == "B1_3")

data_ind <- data %>%
  pivot_wider(names_from = CellType, values_from = Fraction) %>%
  select(sample, section_id, cell_fraction,
         B_cell, B_plasma_IGHA1, B_plasma_IGHG1, B_plasma_IGHG3, B_plasma_IGHM, B_plasma_IGKC, B_plasma_IGLC1, 
         T_CD4, T_CD8_exhausted, T_CTL, T_CXCL13, T_reg, TNK_dividing) %>%
  mutate(B_cells = B_cell + B_plasma_IGHA1 + B_plasma_IGHG1 + B_plasma_IGHG3 + B_plasma_IGHM + B_plasma_IGKC + B_plasma_IGLC1,
         T_cells = T_CD4 + T_CD8_exhausted + T_CTL + T_CXCL13 + T_reg + TNK_dividing) %>%
  select(sample, section_id, cell_fraction, B_cells, T_cells) %>%
  dplyr::rename(Sample_ID = sample)

# viridis::scale_fill_viridis(option="inferno") + 

LS_read_new <- LS_read %>%
   inner_join(data_ind)

B1 <- addPlottingFactor(overlay = B1, annots = LS_read_new,
                        plottingFactor = "B_cells")
B1 <- addPlottingFactor(overlay = B1, annots = LS_read_new,
                        plottingFactor = "T_cells")
B1 <- addPlottingFactor(overlay = B1, annots = LS_read_new,
                        plottingFactor = "cell_fraction")

head(plotFactors(B1))

# Sanity ------------------------------------------------------------------
df <- plotFactors(B1)
dim(df)
dim(LS_read)
setdiff(LS_read$Sample_ID, df$Sample_ID)
# [1] "DSP-1001660018473-B-A02" "DSP-1001660018473-B-A05" "DSP-1001660018473-B-B05" "DSP-1001660018473-B-B08"
# [5] "DSP-1001660018473-B-B09" "DSP-1001660018473-B-B10"
LS_read_sub <- LS_read %>% filter(Sample_ID %in% setdiff(LS_read$Sample_ID, df$Sample_ID))


# -------------------------------------------------------------------------



# options(java.parameters = "-Xmx12g")
# library(RBioFormats)
# checkValidRes(ometiff = tifFile)
# 
# res <- 10
# B1 <- addImageOmeTiff(overlay = B1, ometiff = tifFile, res = res)
#
# B1 <- flipY(B1)
# B1 <- flipY(B1)
plotSpatialOverlay(overlay = B1, colorBy = "segment", scaleBar = FALSE)
p <- plotSpatialOverlay(overlay = B1, colorBy = "Sample_ID", scaleBar = FALSE)
p <- plotSpatialOverlay(overlay = B1, colorBy = "roi_sampleid", scaleBar = FALSE)
plotSpatialOverlay(overlay = B1, colorBy = "cell_fraction", scaleBar = FALSE)
pBcells <- plotSpatialOverlay(overlay = B1, colorBy = "B_cells", scaleBar = FALSE) + 
  viridis::scale_fill_viridis(option="inferno") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_gallery/"
ggsave(
  plot = pBcells,
  width = 8000,
  height = 5394,
  units = "px",
  filename = file.path(save_path, "B1_3_Bcells.png")
)



# # -------------------------------------------------------------------------
# geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Final_Qcd/breast_qcd.rds")
# geo_B1_3 <- geo[, geo$section_id == "B1_3"]
# CD <- as.data.frame(colData(geo_B1_3))
# table(CD$roi, CD$cell_fraction)
# 
# #               Macro Malignant Other PanCK- T_cells
# # A_Islet_2         0         1     1      0       1
# # A_Islet_4         1         1     1      0       1
# # A_Islet_TME_1     0         1     0      1       0
# # A_Islet_TME_2     0         1     0      1       0
# # A_Islet_TME_3     0         1     0      1       0
# # A_Islet_TME_4     0         1     0      1       0
# # A_TME_2           1         0     1      0       1
# # A_TME_4           1         0     1      0       1
# # A_TME_7           1         0     1      0       1
# 

# -------------------------------------------------------------------------

p <- ggplot(data = CD, aes(x = roi_coordinate_y, y = -roi_coordinate_x,
                           label = sample_id2)) + 
                           # label = cell_fraction)) +
                           # label = roi)) +
  geom_point(size = 2) + 
  ggrepel::geom_text_repel(size = 8, color = "blue", point.padding = 0.15, #lwd = 2,
                           min.segment.length = .1, box.padding = .2, max.overlaps = 100) + 
  theme_bw()

p

## B cell ROI - AOIs
# CD$sample_id2[grepl("18473-B-B11", CD$sample_id2)]
# CD$sample_id2[grepl("18473-B-C01", CD$sample_id2)]
# CD$sample_id2[grepl("18473-B-B12", CD$sample_id2)]

CD_Bcell <- CD %>%
  filter(grepl("18473-B-B11|18473-B-C01|18473-B-B12", sample_id2) )

CD_new <- CD %>%
  mutate(roi_coordinate_x = case_when(cell_fraction == "Malignant" ~ roi_coordinate_x - 100,
                                      cell_fraction == "Other" ~ roi_coordinate_x + 100,
                                      .default = roi_coordinate_x),
         roi_coordinate_y = case_when(cell_fraction == "T_cells" ~ roi_coordinate_y - 100,
                                      cell_fraction == "Macro" ~ roi_coordinate_y + 100,
                                      .default = roi_coordinate_y))

p <- ggplot(data = CD_new, aes(x = roi_coordinate_y, y = -roi_coordinate_x,
                           # label = sample_id2)) + 
  label = cell_fraction)) +
  # label = roi)) +
  geom_point(size = 2) + 
  ggrepel::geom_text_repel(size = 8, color = "blue", point.padding = 0.15, #lwd = 2,
                           min.segment.length = .1, box.padding = .2, max.overlaps = 100) + 
  theme_bw()

p


CD_B1 <- CD %>%
  filter(section_id == "B1_3")
CD_B1_atme4 <- CD_B1 %>% filter(roi == "A_TME_4") # not shown in QuPath




