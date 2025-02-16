# Chromium -----------------------------------------------------------------
chrom_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"

chrom_breast_qcd <- readRDS(file.path(chrom_path, "chrom_breast.rds"))
chrom_lung_qcd <- readRDS(file.path(chrom_path, "chrom_lung.rds"))
chrom_dlbcl_qcd <- readRDS(file.path(chrom_path, "chrom_dlbcl.rds"))



chrom_B3 <- chrom_breast_qcd[, chrom_breast_qcd$patient == "B3"]

chrom_B3 <- SCTransform(chrom_B3, assay = "SoupX")
chrom_B3 <- RunPCA(object = chrom_B3, assay = "SCT", npcs = 50)
chrom_B3 <- FindNeighbors(object = chrom_B3, assay = "SCT", reduction = "pca", dims = 1:50)
chrom_B3 <- FindClusters(object = chrom_B3, resolution = 0.1)
chrom_B3 = RunUMAP(object = chrom_B3, assay = "SCT", reduction = "pca", dims = 1:50)

DimPlot(chrom_B3, group.by = "seurat_clusters") |
  FeaturePlot(chrom_B3, markers_159) | 
  FeaturePlot(chrom_B3, markers_14)

markers_159 <- c("NPPC", "GSTP1", "S100P", "GPRC5A")
markers_14 <- c("AZGP1", "MUCL1", "FASN", "SCD")
FeaturePlot(chrom_B3, markers_159)
FeaturePlot(chrom_B3, markers_14)

DimPlot(chrom_B3, group.by = "Harmonised_Level4")


# -------------------------------------------------------------------------
chrom_B3_tu <- chrom_B3[, chrom_B3$Level2 == "Tu_B3"]

chrom_B3_tu <- SCTransform(chrom_B3_tu, assay = "SoupX")
chrom_B3_tu <- RunPCA(object = chrom_B3_tu, assay = "SCT", npcs = 50)
chrom_B3_tu <- FindNeighbors(object = chrom_B3_tu, assay = "SCT", reduction = "pca", dims = 1:50)
chrom_B3_tu <- FindClusters(object = chrom_B3_tu, resolution = 0.4)
chrom_B3_tu = RunUMAP(object = chrom_B3_tu, assay = "SCT", reduction = "pca", dims = 1:50)

(DimPlot(chrom_B3_tu, group.by = "Harmonised_Level4") | 
    DimPlot(chrom_B3_tu, group.by = "seurat_clusters")) 

(FeaturePlot(chrom_B3_tu, markers_159) | 
    FeaturePlot(chrom_B3_tu, markers_14))

