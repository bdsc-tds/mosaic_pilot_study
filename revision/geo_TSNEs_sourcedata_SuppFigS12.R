
disease = "dlbcl"
disease = "breast"
disease = "lung"

geo <- readRDS(paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", disease, "_spe.rds"))
geo_ruv <- readRDS(paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", disease, "_spe_ruv.rds"))

df_save_norm <- data.frame(
  reducedDim(geo, "TSNE")
  )
names(df_save_norm) <- paste0("Norm_", names(df_save_norm))

df_save_ruv <- data.frame(
  reducedDim(geo_ruv, "TSNE")
)
names(df_save_ruv) <- paste0("Batch_", names(df_save_ruv))

df_meta <- data.frame(
  slide = geo$slide_name, 
  sample = geo$section_id,
  aoi = geo$cell_fraction
)

df_save <- cbind(df_save_norm, df_save_ruv, df_meta)
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS12m-r.csv")
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS12a-f.csv")
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS12g-l.csv")

# -------------------------------------------------------------------------
library(ggspavis)
plotDimRed(geo, "TSNE", annotate = "cell_fraction")
plotDimRed(geo_ruv, "TSNE", annotate = "cell_fraction")
