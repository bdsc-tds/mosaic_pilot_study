disease = "breast"
# disease = "lung"
# disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/Visium/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/color_palette.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

## Get chromium by indication ------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/" # breast / # lung
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

vis_rawmat <- "counts"
chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
chrom_B3 <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast/breast_B3_hetero_owkin_annot.rds")
chrom_B3 <- chrom_B3[, chrom_B3$sample_id %in% c("B3", "B3_rep")]
chrom_B3$Harmonised_Level4 <- chrom_B3$Harmonised_Level4_DB
chrom_B3$Harmonised_Level4_DB <- NULL
           
chrom_other <- chrom[, chrom$patient != "B3"]


# -------------------------------------------------------------------------
chrom <- merge(chrom_B3, chrom_other)
tumor_classes <- c("Tu_B3", "Tu_B3_necrosis", "Tu_B3_NPPC", "Tu_B3_PLA2G2A")

chrom <- chrom[, !(chrom$Harmonised_Level4 %in% c("Fibroblast", "Tu_B1_MUCL1", "Tu_B1_MUCL1_necrosis", 
                                                  "Tu_B1_MUCL1_transcription", "Tu_B4_RHOB"))]
chrom$Harmonised_Level4 <- droplevels(as.factor(chrom$Harmonised_Level4))
table(chrom$Harmonised_Level4)

## Get Spot gene expr by sample -----------------------------------------
sce <- readRDS(file.path(save_path_qcd, paste0("B3_2", "_qcd.rds")))

save_bs_path_sample <- paste0(save_bs_path, foldername, "/", "B3_2", "/")

## Decon -------------------------------------------------------------------
# Ref: prep count matrix and cell anno
cell_types <- data.frame(
  cellID = colnames(chrom),
  cellType = as.factor(chrom@meta.data[["Harmonised_Level4"]]), 
  sampleInfo = "sample1")
rownames(cell_types) <- colnames(chrom)
# chrom_mat <- chrom@assays$SoupX@counts # decont chromium 
chrom_mat <- chrom@assays$RNA@counts # raw chromium: use raw chrom for decont

# sce: prep coords and count matrix
coords <- data.frame(x = spatialCoords(sce)[, 1], y = spatialCoords(sce)[, 2]); rownames(coords) <- colnames(sce)
counts_sce <- assay(sce, vis_rawmat)

CARD_obj = createCARDObject(
  sc_count = chrom_mat,
  sc_meta = cell_types,
  spatial_count = counts_sce,
  spatial_location = coords,
  ct.varname = "cellType",
  ct.select = unique(cell_types$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

results = data.frame(CARD_obj@Proportion_CARD)
results = results[, sort(colnames(results))]
write.csv(results, paste0(save_bs_path_sample, "B3_2", "_spot_Level4_decon_B3_newanno.csv"))
results <- read.csv(paste0(save_bs_path_sample, "B3_2", "_spot_Level4_decon_B3_newanno.csv"), row.names = 1)

# Visualize spatial pie chart -----------------------------------------
# colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
#            "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
#            "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D") # to cell types in alphabetic order

level4_decon_cols = level4_cellcolors[names(level4_cellcolors) %in% sort(colnames(results))]

# colors = level4_cellcolors
p1 <- CARD.visualize.pie(
  proportion = results,
  spatial_location = CARD_obj@spatial_location, 
  colors = level4_decon_cols, 
  radius = 9) + theme(legend.position = "right") # + # coord_flip() + 
# scale_y_reverse() +
# scale_x_reverse()
print(p1)

image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])

# # if need to coord_flip an image
# image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
# image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) -2

pdf(file.path(save_bs_path_sample, paste0("B3_2", "_spot_level4_decon_bycelltype_B3newanno.pdf")),
    height = plot_height, width = plot_width)
print(p1)
dev.off()

# Decon label ---------------------------------------------------------
# 1st dominant cell type
results$max <- colnames(results)[max.col(results, ties.method="first")]

result_max <- NULL
for(s in 1:nrow(results)){
  results_s <- results[s, ]
  if(results_s[results_s$max] < 0.2 ){
    results_s_max = "Mix"
  }else{
    results_s_max = results_s$max
  }
  result_max <- c(result_max, results_s_max)
}
results2 <- results
results2$max <- result_max

if(nrow(results2) < ncol(sce)){
  sce$Level4_decon_max <- NULL
  sce$barcode <- colnames(sce)
  results2$barcode <- rownames(results2)
  
  CD <- as.data.frame(colData(sce))
  CD <- CD %>% 
    left_join(results2 %>% select(max, barcode)) %>%
    rename(Level4_decon_max = max) %>%
    mutate(Level4_decon_max = factor(Level4_decon_max, levels = c(sort(unique(results2$max))[sort(unique(results2$max)) != "Mix"], "Mix")))
  rownames(CD) <- CD$barcode; CD$barcode <- NULL
  
  colData(sce) <- as(CD, "DFrame")
}else{
  sce[["Level4_decon_max"]] <- factor(results2$max, levels = c(sort(unique(results2$max))[sort(unique(results2$max)) != "Mix"], "Mix"))
}

level4_decon_cols_max <- plyr::mapvalues(levels(sce[["Level4_decon_max"]]), from = br_lu_level4_cellnames, to = level4_cellcolors)

source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/08_Baye_New_Pipeline_Decon/helpers.R")

sce_plot <- sce[, !is.na(sce$Level4_decon_max)]
p2 <- plotSpots_deconmax(sce_plot, 
                         annotate = "Level4_decon_max",
                         pt.size = 5,
                         palette = level4_decon_cols_max) # + 
# coord_flip() + 
# scale_y_reverse() +
# scale_x_reverse()

pdf(file.path(save_bs_path_sample, paste0("B3_2", "_spot_level4_decon_bycelltype_majorityvote_B3newanno.pdf")),
    height = plot_height, width = plot_width)
print(p2)
dev.off()

# -------------------------------------------------------------------------
all(rownames(results) == colnames(sce))

sce$Tu_B3 <- results$Tu_B3
sce$Tu_B3_necrosis <- results$Tu_B3_necrosis
sce$Tu_B3_NPPC <- results$Tu_B3_NPPC
sce$Tu_B3_PLA2G2A <- results$Tu_B3_PLA2G2A
sce$Fibroblast_B3 <- results$Fibroblast_B3

library(ggspavis)
(plotSpots(sce, annotate = "Tu_B3") | 
  plotSpots(sce, annotate = "Tu_B3_necrosis")) /
  (plotSpots(sce, annotate = "Tu_B3_NPPC") |
  plotSpots(sce, annotate = "Tu_B3_PLA2G2A"))

plotSpots(sce, annotate = "Fibroblast_B3")

# UMAP Chrom ---------------------------------------------------------------
chrom_B3 <- SCTransform(chrom_B3, assay = "SoupX")
chrom_B3 <- RunPCA(object = chrom_B3, assay = "SCT", npcs = 50)
chrom_B3 <- FindNeighbors(object = chrom_B3, assay = "SCT", reduction = "pca", dims = 1:50)
chrom_B3 <- FindClusters(object = chrom_B3, resolution = 0.4)
chrom_B3 = RunUMAP(object = chrom_B3, assay = "SCT", reduction = "pca", dims = 1:50)

chrom_B3$Level4_manuscript <- ifelse(!chrom_B3$Harmonised_Level4 %in% c("Tu_B3", "Tu_B3_necrosis", "Tu_B3_NPPC", "Tu_B3_PLA2G2A"),
                                     "Healthy cells", ifelse(chrom_B3$Harmonised_Level4 == "Tu_B3", "Tu_B3_generic", chrom_B3$Harmonised_Level4))
palette <- c("Healthy cells" = "#008080", 
             "Tu_B3_generic" = "#D4FF90",  
             "Tu_B3_necrosis" = "#C9DC84",
             "Tu_B3_NPPC" =  "#DBE153", 
             "Tu_B3_PLA2G2A" = "#9FCB00")
DimPlot(chrom_B3, group.by = "Level4_manuscript", cols = palette) 


# Show only tumor ball ----------------------------------------------------
chrom_B3_tu <- chrom_B3[, chrom_B3$Level4_manuscript != "Healthy cells"]
chrom_B3_tu <- SCTransform(chrom_B3_tu, assay = "SoupX")
chrom_B3_tu <- RunPCA(object = chrom_B3_tu, assay = "SCT", npcs = 50)
chrom_B3_tu <- FindNeighbors(object = chrom_B3_tu, assay = "SCT", reduction = "pca", dims = 1:50)
chrom_B3_tu <- FindClusters(object = chrom_B3_tu, resolution = 0.4)
chrom_B3_tu = RunUMAP(object = chrom_B3_tu, assay = "SCT", reduction = "pca", dims = 1:50)

Idents(chrom_B3_tu) <- chrom_B3_tu$Level4_manuscript

library(SeuratObject)
Tu_B3_NPPC <- WhichCells(object = chrom_B3_tu, idents = "Tu_B3_NPPC")
Tu_B3_PLA2G2A <- WhichCells(object = chrom_B3_tu, idents = "Tu_B3_PLA2G2A")

cells <- list(Tu_B3_NPPC = Tu_B3_NPPC, Tu_B3_PLA2G2A = Tu_B3_PLA2G2A)

library(scCustomize)
Cell_Highlight_Plot(seurat_object = chrom_B3_tu, cells_highlight = cells)

chrom_B3_tu_2_clus <- chrom_B3_tu[, chrom_B3_tu$Level4_manuscript %in% c("Tu_B3_NPPC", "Tu_B3_PLA2G2A")]
chrom_B3_tu_2_clus_markers <- FindAllMarkers(chrom_B3_tu_2_clus, assay = "SCT", slot='data',
                                 group.by = "Level4_manuscript",
                        only.pos = TRUE,
                        logfc.threshold = 1) %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)

top_markers <- chrom_B3_tu_2_clus_markers %>% group_by(cluster) # %>% top_n(20)

Hm <- Seurat::DoHeatmap(chrom_B3_tu_2_clus, 
                        assay = "SCT",
                        features = top_markers$gene, 
                        slot = 'scale.data',
                        group.by = "Level4_manuscript", 
                        angle=0, 
                        size=4, 
                        label = TRUE, 
                        raster=FALSE) + 
  guides(col = FALSE)


write.csv(chrom_B3_tu_2_clus_markers, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2/chrom_B3_tu_2_clus_markers.csv")




