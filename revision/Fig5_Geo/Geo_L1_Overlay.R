options(java.parameters = "-Xmx12g")
library(RBioFormats)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(SpatialOmicsOverlay)
library(GeomxTools)
library(readxl)
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/lung_spe_ruv.rds")
# Do not QC
# geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/GeoMx/lung_raw/L1_1_0PSV.rds")


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
geo_L1_1 <- geo[, geo$section_id == "L1_1"] #qcd
# geo_L1_1 <- geo

# countmat <- assay(geo_L1_1, "counts")
countmat <- assay(geo_L1_1, "quantile")
# logcountmat <- assay(geo_L1_1, "logcounts") # logCPM
# expr_df <- data.frame(MS4A1 = countmat["MS4A1", ],
#                       logMS4A1 = log(countmat["MS4A1", ] + 1),
#                       CXCL13 = countmat["CXCL13", ],
#                       logCXCL13 = log(countmat["CXCL13", ] + 1),
#                       CXCR5 = countmat["CXCR5", ],
#                       logCXCR5 = log(countmat["CXCR5", ] + 1),
#                       CCL19 = countmat["CCL19", ],
#                       logCCL19 = log(countmat["CCL19", ] + 1),
#                       Sample_ID = geo_L1_1$SampleID)
expr_df <- data.frame(MS4A1 = countmat["MS4A1", ],
                      CXCL13 = countmat["CXCL13", ],
                      CXCR5 = countmat["CXCR5", ],
                      CCL19 = countmat["CCL19", ],
                      Sample_ID = geo_L1_1$sample_id2)

LS_read_new <- annots %>%
  left_join(expr_df)


## Decon
library(tidyr)
geo_level4lung <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/lung_batched_decon_long.csv")
data <- geo_level4lung %>%
  filter(section_id == "L1_1")

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

LS_read_new <- LS_read_new %>%
  left_join(data_ind)

L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "B_cells")
L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "T_cells")
L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "cell_fraction")

# L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
#                         plottingFactor = "logMS4A1")
# L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
#                         plottingFactor = "logCXCL13")
# L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
#                         plottingFactor = "logCXCR5")
# L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
#                         plottingFactor = "logCCL19")

L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "MS4A1")
L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "CXCL13")
L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "CXCR5")
L1 <- addPlottingFactor(overlay = L1, annots = LS_read_new,
                        plottingFactor = "CCL19")

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
p_seg_L1 <- plotSpatialOverlay(overlay = L1, colorBy = "segment", scaleBar = FALSE) 
# plogMS4A1 <- plotSpatialOverlay(overlay = L1, colorBy = "logMS4A1", scaleBar = FALSE) +
#   viridis::scale_fill_viridis(option="magma") + 
#   theme_minimal() + 
#   theme(rect = element_rect(fill = "transparent"),
#         panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank())
# 
# 
# plogCXCL13 <- plotSpatialOverlay(overlay = L1, colorBy = "logCXCL13", scaleBar = FALSE) + 
#   viridis::scale_fill_viridis(option="magma") + 
#   theme_minimal() + 
#   theme(rect = element_rect(fill = "transparent"),
#         panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank())
#   
# plogCXCR5 <- plotSpatialOverlay(overlay = L1, colorBy = "logCXCR5", scaleBar = FALSE) + 
#   viridis::scale_fill_viridis(option="magma") + 
#   theme_minimal() + 
#   theme(rect = element_rect(fill = "transparent"),
#         panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank())
# 
# plogCCL19 <- plotSpatialOverlay(overlay = L1, colorBy = "logCCL19", scaleBar = FALSE) + 
#   viridis::scale_fill_viridis(option="magma") + 
#   theme_minimal() + 
#   theme(rect = element_rect(fill = "transparent"),
#         panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank())

range(na.omit(c(LS_read_new$MS4A1, LS_read_new$CXCL13, LS_read_new$CXCR5, LS_read_new$CCL19))) # , limits = c(1.326667, 3.808864)

pMS4A1 <- plotSpatialOverlay(overlay = L1, colorBy = "MS4A1", scaleBar = FALSE) +
  # viridis::scale_fill_viridis(option="magma", limits = c(1.31, 3.81)) + 
  viridis::scale_fill_viridis(option="magma") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

pCXCL13 <- plotSpatialOverlay(overlay = L1, colorBy = "CXCL13", scaleBar = FALSE) + 
  viridis::scale_fill_viridis(option="magma", limits = c(1.31, 3.81)) + 
  # viridis::scale_fill_viridis(option="magma") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

pCXCR5 <- plotSpatialOverlay(overlay = L1, colorBy = "CXCR5", scaleBar = FALSE) + 
  viridis::scale_fill_viridis(option="magma", limits = c(1.31, 3.81)) + 
  # viridis::scale_fill_viridis(option="magma") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

pCCL19 <- plotSpatialOverlay(overlay = L1, colorBy = "CCL19", scaleBar = FALSE) + 
  viridis::scale_fill_viridis(option="magma", limits = c(1.31, 3.81)) + 
  # viridis::scale_fill_viridis(option="magma") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


# Decon plot B cells T cells  ---------------------------------------------
pBcells <- plotSpatialOverlay(overlay = L1, colorBy = "B_cells", scaleBar = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

pTcells <- plotSpatialOverlay(overlay = L1, colorBy = "T_cells", scaleBar = FALSE) +
  viridis::scale_fill_viridis(option="inferno") + 
  theme_minimal() + 
  theme(rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


## Save plots with transparent background
save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_expr/"

# -------------------------------------------------------------------------
CD <- as.data.frame(colData(geo_L1_1))
CD <- CD[CD$roi != "A_Islet_1_", ]
CD$sample_id2 <- rownames(CD)

df_plotfactor <- plotFactors(L1)
df_plotfactor <- df_plotfactor %>%
  mutate(sample_id2 = rownames(df_plotfactor)) %>%
  left_join(CD %>% select(sample_id2, roi_coordinate_x, roi_coordinate_y))
df_plotfactor_new <- df_plotfactor %>%
  mutate(roi_coordinate_x = case_when(segment == "PanCK" ~ roi_coordinate_x - 100,
                                      segment == "Other" ~ roi_coordinate_x + 100,
                                      .default = roi_coordinate_x),
         roi_coordinate_y = case_when(segment == "CD3" ~ roi_coordinate_y - 100 ,
                                      # segment == "Macro" ~ roi_coordinate_y + 100,
                                      .default = roi_coordinate_y))


p <- ggplot(data = df_plotfactor_new, # CD, 
            aes(x = roi_coordinate_y, y = -roi_coordinate_x,
                color = CXCR5
                # , label = sample_id2
                )) + 
                # label = cell_fraction)) +
                # label = roi)) +
  geom_point(size = 5) + 
  viridis::scale_color_viridis(option="magma") + 
  # ggrepel::geom_text_repel(size = 8, color = "blue", point.padding = 0.15, #lwd = 2,
  #                          min.segment.length = .1, box.padding = .2, max.overlaps = 100) + 
  theme_bw()

p

exprdf <- meta(overlay(L1))

# "DSP-1001660018472-A-C01" %in% CD$sample_id2 # FALSE # A_TME_5


# -------------------------------------------------------------------------
CD <- as.data.frame(colData(geo_L1_1[, geo_L1_1$roi != "A_Islet_1_"]))

p <- ggplot(data = CD, aes(x = roi_coordinate_y, y = -roi_coordinate_x,
                               # label = sample_id2)) + 
                               label = roi)) + 
                               # label = cell_fraction)) +
  geom_point(size = 2) + 
  ggrepel::geom_text_repel(size = 8, color = "blue", point.padding = 0.15, #lwd = 2,
                           min.segment.length = .1, box.padding = .2, max.overlaps = 100) + 
  theme_bw()

p

# A_TME_5, 
# A_TME_2, 
# A_Islet_2

# Expression heatmap or boxplot --------------------------------------------
library(RColorBrewer)
LS_read_new_sub <- LS_read_new %>%
  filter(ROILabel %in% c("A_TME_2", "A_Islet_2")) %>% # "A_TME_5",
  select(ROILabel, cell_fraction, MS4A1, CXCL13, CXCR5, CCL19)

range(na.omit(c(LS_read_new_sub$MS4A1, LS_read_new_sub$CXCL13, LS_read_new_sub$CXCR5, LS_read_new_sub$CCL19))) # 1.408960 2.981307

pMS4A1_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = MS4A1)) +
  geom_tile() +
  theme_minimal() +
  # scale_fill_gradientn(colors = rev(brewer.pal(9, "RdYlBu")), limits = c(1.40, 3)) +
  viridis::scale_fill_viridis(option="magma", limits = c(1.40, 3)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pMS4A1_heat

pCXCL13_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CXCL13)) +
  geom_tile() +
  theme_minimal() +
  # scale_fill_gradientn(colors = rev(brewer.pal(9, "RdYlBu")), limits = c(1.40, 3)) +
  viridis::scale_fill_viridis(option="magma", limits = c(1.40, 3)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pCXCL13_heat

pCXCR5_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CXCR5)) +
  geom_tile() +
  theme_minimal() +
  # scale_fill_gradientn(colors = rev(brewer.pal(9, "RdYlBu")), limits = c(1.40, 3)) +
  viridis::scale_fill_viridis(option="magma", limits = c(1.40, 3)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pCXCR5_heat

pCCL19_heat <- ggplot(LS_read_new_sub, aes(x = cell_fraction, y = ROILabel, fill = CCL19)) +
  geom_tile() +
  theme_minimal() +
  # scale_fill_gradientn(colors = rev(brewer.pal(9, "RdYlBu")), limits = c(1.40, 3)) +
  viridis::scale_fill_viridis(option="magma", limits = c(1.40, 3)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()); pCCL19_heat


###########################################
save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_expr/final"
plot_title = "L1_MS4A1_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pMS4A1_heat)
dev.off()

plot_title = "L1_CXCL13_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCXCL13_heat)
dev.off()

plot_title = "L1_CXCR5_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCXCR5_heat)
dev.off()

plot_title = "L1_CCL19_heat.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 3.5,
    height = 1)
print(pCCL19_heat)
dev.off()


# Grouped bar plot --------------------------------------------------------
LS_read_new_sub_plt <- rbind(LS_read_new_sub,
                         c("A_Islet_2", "T cells",   NA, NA, NA, NA),
                         c("A_TME_2",   "Malignant", NA, NA, NA, NA)) %>%
  mutate(MS4A1 = as.numeric(MS4A1),
         CXCL13 = as.numeric(CXCL13),
         CXCR5 = as.numeric(CXCR5),
         CCL19 = as.numeric(CCL19)) %>%
  dplyr::rename(`AOI label` = cell_fraction) %>%
  mutate(ROILabel = ifelse(ROILabel == "A_Islet_2", "Islet", "TME")) %>%
  mutate(ROILabel = factor(ROILabel, levels = c("TME", "Islet"))) %>%
  dplyr::rename(`ROI label` = ROILabel)


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
    scale_fill_manual(values = c("Malignant" = "#36E7E4", "Other" = "#C2E2B9", "T cells" = "#30FF00")) +
    scale_pattern_manual(values = c("TME" = "none", "Islet" = "stripe")) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))
  
  p
}
pMS4A1_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "MS4A1")
pCXCL13_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CXCL13")
pCXCR5_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CXCR5")
pCCL19_bar <- plot_bar_TLS_expression(LS_read_new_sub_plt, gene = "CCL19")

save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_expr/"
plot_title = "L1_MS4A1_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 3)
print(pMS4A1_bar)
dev.off()

plot_title = "L1_CXCL13_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 3)
print(pCXCL13_bar)
dev.off()

plot_title = "L1_CXCR5_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 3)
print(pCXCR5_bar)
dev.off()

plot_title = "L1_CCL19_bar.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8,
    height = 3)
print(pCCL19_bar)
dev.off()










