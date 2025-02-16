chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
seu <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/final_owkin_annot.rds") # BrLu
seub <- readRDS(file.path(chrompath, "chrom_breast.rds"))
seul <- readRDS(file.path(chrompath, "chrom_lung.rds"))

df_save <- data.frame(CT1 = seu$annot_l1,
                      CT2 = seu$annot_l2,
                      CT3 = seu$annot_l3,
                      CT4 = c(seub$Harmonised_Level4, seul$Harmonised_Level4),
                      Barcode = colnames(seu),
                      seu@reductions$umap@cell.embeddings,
                      indi = seu$tissue
)
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS3a.csv")  

# --------------------------
disease = "dlbcl"; seu <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
df_save <- data.frame(CT1 = seu$Level1,
                      CT2 = seu$Level2,
                      CT3 = seu$Level3,
                      CT4 = seu$Harmonised_Level4,
                      Barcode = colnames(seu),
                      seu@reductions$umap@cell.embeddings,
                      indi = disease
)
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS3b.csv")  
