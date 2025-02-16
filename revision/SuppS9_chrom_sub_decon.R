disease = "breast"
## Get chromium by indication ------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/" # breast / # lung
vis_rawmat <- "counts"
chrom <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast/breast_B3_hetero_owkin_annot.rds")
chrom_other_pt_tumor <- c("Tu_B1_MUCL1", "Tu_B1_MUCL1_necrosis", "Tu_B1_MUCL1_transcription", "Tu_B4_RHOB")
fibro_remove = "Fibroblast"
to_remove = c(chrom_other_pt_tumor, fibro_remove)

chrom <- chrom[, !(chrom$Harmonised_Level4_DB %in% to_remove)]
chrom$Harmonised_Level4_DB <- droplevels(as.factor(chrom$Harmonised_Level4_DB))
table(chrom$Harmonised_Level4_DB)

## Get Spot gene expr by sample -----------------------------------------
sce <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_qcd/breast_qcd/B3_2_qcd.rds")

save_bs_path_sample = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/"

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
write.csv(results, paste0(save_bs_path_sample, "B3_2_spot_Level4_decon_B3_newanno.csv"))

results <- read.csv(paste0(save_bs_path_sample, "B3_2_spot_Level4_decon_B3_newanno.csv"), row.names = 1)

sce$Tu_B3_NPPC = results$Tu_B3_NPPC
sce$Tu_B3_PLA2G2A = results$Tu_B3_PLA2G2A

df_save <- data.frame(
  Tu_B3_NPPC = sce$Tu_B3_NPPC,
  Tu_B3_PLA2G2A = sce$Tu_B3_PLA2G2A,
  x = sce$array_col,
  y = sce$array_row
)

write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS9e.csv")



