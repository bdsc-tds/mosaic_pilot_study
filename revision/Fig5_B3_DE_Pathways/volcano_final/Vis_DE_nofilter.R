library(SingleCellExperiment)
library(Seurat)

# noSpotClean + logNorm + BayesSpace
sample_name = "B3_2"
disease = "breast"
sce <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered.rds") # noSpotClean + logNorm + BayesSpace

# Merge DE ----------------------------------------------------------------
sce$spatial.cluster_merge2 <- ifelse(sce$spatial.cluster %in% c(1, 5, 9), "1_5_9", sce$spatial.cluster)

save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/DiffDE_nofilter"


## Convert SCE to seurat object and use cluster_col as identifier
seurat <- Seurat::CreateSeuratObject(counts=assays(sce)[["logcounts"]],
                                     assay='Spatial',
                                     meta.data=as.data.frame(colData(x)))
LayerData(seurat, layer = "data", assay = "Spatial") <- LayerData(seurat, 
                                                                  layer = "counts", 
                                                                  assay = "Spatial") # Seurat DE works on data slot, which is logged.
# Here just migrating logcount as counts to data. 
# Other DE, limma voom works on raw counts

# Set cell identity by the cluster ID
seurat <- Seurat::SetIdent(seurat, value = "spatial.cluster_merge2")

## Subset to specified clusters
seurat <- subset(seurat, idents = clusters)

all_markers <- Seurat::FindMarkers(seurat, assay='Spatial', slot='data',
                                   # group.by=cluster_col,
                                   ident.1 = "1_5_9",
                                   ident.2 = "14", 
                                   logfc.threshold=-Inf, 
                                   min.pct = -Inf, 
                                   only.pos=FALSE)
dim(all_markers) # [1] 17883     5
range(all_markers$p_val) # [1] 7.254495e-119  1.000000e+00
range(all_markers$avg_log2FC) # -7.463204  3.441427
all_markers$gene <- rownames(all_markers)

table(all_markers$p_val_adj < 0.05 & abs(all_markers$avg_log2FC) > 1)
# FALSE  TRUE 
# 17743   140 
table(all_markers$p_val_adj < 0.05 & all_markers$avg_log2FC > 1)
# FALSE  TRUE 
# 17833    90 
table(all_markers$p_val_adj < 0.05 & all_markers$avg_log2FC < -1)
# FALSE  TRUE 
# 17833    50 

all_markers$cluster <- ifelse(all_markers$avg_log2FC > 0, "1_5_9", "14")

# contrast <- "1&5&9_14_lfc0"
contrast <- "1&5&9_14_lfc0_pval1_allgenes"
write.csv(all_markers, file.path(save_path_DE, paste0(contrast, ".csv")))
head(DE_result)
