## Method 1, replicate Visium's norm method -----------------------------
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(ggrepel)
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Final_Qcd/breast_qcd.rds")

# names(colData(geo))

CD <- as.data.frame(colData(geo))
table(CD$roi, CD$cell_fraction)

CD_B3_1 <- CD %>%
  dplyr::filter(section_id == "B3_1") 
length(table(CD_B3_1$roi)) # 11 ROI


CD_B3_1_Malig <- CD_B3_1 %>%
  dplyr::filter(cell_fraction == "Malignant") 
length(table(CD_B3_1_Malig$roi)) # 9 ROI

CD_B3_1_Other <- CD_B3_1 %>%
  dplyr::filter(cell_fraction == "Other") 
length(table(CD_B3_1_Other$roi)) # 9 ROI

CD_B3_1_Tcells <- CD_B3_1 %>%
  dplyr::filter(cell_fraction == "T_cells") 
length(table(CD_B3_1_Tcells$roi)) # 3 ROI

CD_B3_1_Panckm <- CD_B3_1 %>%
  dplyr::filter(cell_fraction == "PanCK-") 
length(table(CD_B3_1_Panckm$roi)) # 1 ROI


p <- ggplot(data = CD_B3_1, aes(x = roi_coordinate_y, y = -roi_coordinate_x,
                                      label = roi)) +
  # label = sample_id2)) +
  geom_point(size = 2) + 
  # scale_y_reverse() +
  ggrepel::geom_text_repel(size = 8, color = "blue", point.padding = 0.15, #lwd = 2,
                  min.segment.length = .1, box.padding = .2, max.overlaps = 100) + 
  theme_bw()

p
  
# -------------------------------------------------------------------------
# clus 14
clus14_id <- CD_B3_1_Malig$sample_id2[CD_B3_1_Malig$roi %in% c("C_Islet_6", "C_Islet_8", "C_TME_8")]

# clus 159
clus159_id <- CD_B3_1_Malig$sample_id2[!(CD_B3_1_Malig$roi %in% c("C_Islet_6", "C_Islet_8", "C_TME_8"))]


# -------------------------------------------------------------------------
# subset to only B3_1
geo_ <- geo[, geo$section_id == "B3_1"]
geo_$malig_sub <- ifelse(geo_$sample_id2 %in% clus14_id, "clus14",
                    ifelse(geo_$sample_id2 %in% clus159_id, "clus159", geo_$cell_fraction))

library(scater)
sce <- logNormCounts(geo_)
assays(geo_, withDimnames=FALSE)$quantile <- preprocessCore::normalize.quantiles(assay(geo_, "log1p")) 
sce <- geo_
# TODO: try Quantile + Seurat DE

library(scran)
dec <- modelGeneVar(sce, assay.type = "logcounts")
# dec <- modelGeneVar(sce, assay.type = "quantile")
top <- getTopHVGs(dec, n = 2000)
set.seed(100)
sce <- runPCA(sce, ncomponents = 7, subset_row = top, exprs_values = "logcounts")

set.seed(500)
# sce <- runUMAP(sce, assay.type = "logcounts") # 22 aois error free
sce <- runUMAP(sce, assay.type = "logcounts", n_neighbors = 5) 
set.seed(100)
# sce <- runTSNE(sce, assay.type = "logcounts")
sce <- runTSNE(sce, assay.type = "logcounts", n_neighbors = 5)

sce$cell_fraction <- as.factor(sce$cell_fraction)
plotDimRed(sce, type = "UMAP", annotate = "cell_fraction", text_by = "cell_fraction") |
  plotDimRed(sce, type = "TSNE", annotate = "cell_fraction", text_by = "cell_fraction")

sce$malig_sub <- as.factor(sce$malig_sub)
plotDimRed(sce, type = "UMAP", annotate = "malig_sub", text_by = "malig_sub") |
  plotDimRed(sce, type = "TSNE", annotate = "malig_sub", text_by = "malig_sub")

# -------------------------------------------------------------------------
# subset to only B3_1 malig
geo_malig <- geo[, geo$section_id == "B3_1" & geo$cell_fraction == "Malignant"]
geo_malig$malig_sub <- ifelse(geo_malig$sample_id2 %in% clus14_id, "clus14",
                              ifelse(geo_malig$sample_id2 %in% clus159_id, "clus159", NA))

library(scater)
sce <- logNormCounts(geo_malig)

library(scran)
dec <- modelGeneVar(sce, assay.type = "logcounts")
top <- getTopHVGs(dec, n = 2000)
set.seed(100)
sce <- runPCA(sce, ncomponents = 7, subset_row = top, exprs_values = "logcounts")

set.seed(100)
sce <- runUMAP(sce, assay.type = "logcounts", n_neighbors = 5) # 9 aois uwot error, change n_neighbors to 5
set.seed(100)
sce <- runTSNE(sce, assay.type = "logcounts", n_neighbors = 5)

sce$malig_sub <- as.factor(sce$malig_sub)
plotDimRed(sce, type = "UMAP", annotate = "malig_sub", text_by = "malig_sub") |
  plotDimRed(sce, type = "TSNE", annotate = "malig_sub", text_by = "malig_sub")


# -------------------------------------------------------------------------
library(Seurat)
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Manuscript_Figure/Fig5_Visium_clustering_biology/B3/00_Seurat_DE_heatmap_ClusPlot_by_Cluster_helper.R")

DE_result <- seurat_cluster_DE(sce, clusters = c("clus159", "clus14"), cluster_col="malig_sub", n_markers = 10)








