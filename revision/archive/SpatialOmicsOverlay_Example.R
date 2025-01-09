library(dplyr)
library(SpatialExperiment)
library(ggplot2)
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Final_Qcd/breast_qcd.rds")

# names(colData(geo))

CD <- as.data.frame(colData(geo))
table(CD$roi, CD$cell_fraction)
#                          Macro Malignant Other PanCK- T_cells
# A_Islet_2                    0         2     2      0       1
# A_Islet_4                    1         2     2      0       2
# A_Islet_TME_1                0         2     1      1       1
# A_Islet_TME_2                0         2     1      1       1
# A_Islet_TME_3                1         2     1      1       1
# A_Islet_TME_4                1         2     1      1       1
# A_TME_2                      1         0     2      0       2
# A_TME_4                      1         0     1      0       1
# A_TME_7                      2         0     2      0       2
# B_Islet_1                    1         1     1      0       0
# B_Islet_2                    0         1     1      0       0
# B_Islet_3                    0         1     1      0       1
# B_Islet_4                    0         1     1      0       0
# B_Islet_5                    0         1     0      0       0
# B_Islet_6                    0         1     1      0       0
# B_Islet_7                    1         1     1      0       0
# B_TME_2                      1         0     1      0       0
# B_TME_5                      1         0     1      0       0
# B_TME_6                      1         0     0      0       0
# B_TME_7                      0         0     1      0       0
# C_Islet_1                    0         1     1      0       0
# C_Islet_2                    0         1     1      0       0
# C_Islet_3                    0         1     1      0       0
# C_Islet_4                    0         1     1      0       0
# C_Islet_5                    0         1     1      0       1
# C_Islet_6                    0         1     1      0       0
# C_Islet_7                    0         1     1      0       0
# C_Islet_8                    0         1     1      0       0
# C_TME_1                      0         0     0      0       1
# C_TME_2                      0         0     0      1       0
# C_TME_8                      0         1     1      0       1
# D_Invasive_front_1_TME       1         0     1      0       1
# D_Invasive_front_1_tumor     1         1     1      0       1
# D_Islet_1                    1         1     1      0       0
# D_Islet_2                    1         1     1      0       1
# D_Islet_3                    1         1     1      0       0
# D_TME_1                      0         0     1      0       1
# D_TME_2                      0         0     1      0       0
# D_TME_3                      1         0     1      0       1

CD_B1_3 <- CD %>%
  dplyr::filter(section_id == "B1_3") #%>%
  #dplyr::filter(cell_fraction == "Malignant")

CD_B3_1 <- CD %>%
  dplyr::filter(section_id == "B3_1") %>%
  dplyr::filter(cell_fraction == "Malignant")

write.csv(CD_B3_1, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/CD_B3_1.csv")

ggplot(data = CD_B3_1, aes(x = roi_coordinate_x, y = roi_coordinate_y)) +
  geom_point()

scatterplot(CD_B3_1$scan_offset_x, CD_B3_1$scan_offset_y)


geo_B3_1 <- CD_B3_1 %>%
  dplyr::select(slide_name, sample_id2, pathology, roi_category, aoi, cell_fraction)

# write.csv(geo_B3_1, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Geo_B3_1.csv")

# # Example -----------------------------------------------------------------
library(SpatialOmicsOverlay)
library(GeomxTools)

tifFile <- downloadMouseBrainImage()
tifFile

muBrainLW <- system.file("extdata", "muBrain_subset_LabWorksheet.txt",
                                package = "SpatialOmicsOverlay")
muBrainLW <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/muBrain_ED_LabWorksheet.txt"

## Try to add in ROILabel metadata from .rds
path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/archive"
# unzip(file.path(path, "muBrain_GxT.zip"), exdir = path)

sumexp <- readRDS(file.path(path, "muBrain_GxT.RDS"))
sumexp
DF <- as.data.frame(pData(sumexp))

# test <- read.table(muBrainLW, sep = ",")
# test
# test[13,]
# test[14,]
# test[15,]

muBrain <- readSpatialOverlay(ometiff = tifFile, annots = muBrainLW, 
                              slideName = "D5761 (3)", image = FALSE,  
                              saveFile = FALSE, outline = FALSE)

muBrain
head(sampNames(muBrain))
slideName(muBrain)
head(meta(overlay(muBrain)))
head(coords(muBrain))


# -------------------------------------------------------------------------
# length(which(sampNames(muBrain) != sub(".dcc", "", rownames(DF))))
# all(sort(sampNames(muBrain)) == sub(".dcc", "", rownames(DF)))
# 
# newDFROILabel <- DF[paste0(sampNames(muBrain), ".dcc"), "ROILabel"]
# meta(overlay(muBrain)) <- data.frame(newDFROILabel)
# -------------------------------------------------------------------------
plotSpatialOverlay(overlay = muBrain, hiRes = FALSE, legend = FALSE)

muBrainAnnots <- readLabWorksheet(lw = muBrainLW, slideName = "D5761 (3)")
muBrainGeomxSet <- readRDS(unzip(system.file("extdata", "muBrain_GxT.zip",
                                             package = "SpatialOmicsOverlay"))) # This is for gene expression data
muBrain <- addPlottingFactor(overlay = muBrain, annots = muBrainAnnots, 
                             plottingFactor = "segment")
muBrain <- addPlottingFactor(overlay = muBrain, annots = muBrainGeomxSet, 
                             plottingFactor = "Calm1")
muBrain <- addPlottingFactor(overlay = muBrain, annots = muBrainAnnots, 
                             plottingFactor = "ROILabel")
muBrain
head(plotFactors(muBrain))

# -------------------------------------------------------------------------
plotSpatialOverlay(overlay = muBrain, hiRes = FALSE, colorBy = "Calm1", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Calm1 Expression in Mouse Brain")


checkValidRes(ometiff = tifFile)

res <- 8
muBrain <- addImageOmeTiff(overlay = muBrain, ometiff = tifFile, res = res)
muBrain

showImage(muBrain)
# # -------------------------------------------------------------------------
plotSpatialOverlay(overlay = muBrain, colorBy = "segment", corner = "topcenter", 
                   scaleBarWidth = 0.5, textDistance = 130, scaleBarColor = "cyan",
                   alpha = 0.3,
                   fluorLegend = TRUE)

Coord <- as.data.frame(meta(overlay(muBrain)))
muBrainAnnots$x <- Coord$X;
muBrainAnnots$y <- Coord$Y;
ggplot(data = muBrainAnnots, aes(x = x, y = y)) +
  geom_point() + 
  scale_y_reverse()














