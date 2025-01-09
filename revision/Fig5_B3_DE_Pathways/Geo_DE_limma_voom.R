# 3 samples in area A

library(edgeR)
library(limma)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(ggrepel)
library(ggspavis)
library(patchwork)

# Prep dataset for DE between malig subtypes only --------------------------
geo_all <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds")
geo_malig <- geo_all[, geo_all$section_id == "B3_1" &
                       geo_all$cell_fraction == "Malignant" # & 
                     # geo_all$sample_id2 != "DSP-1001660013520-C-B12"
]
clus14_id <- c("DSP-1001660013520-C-B02", "DSP-1001660013520-C-B06", "DSP-1001660013520-C-B12")
clus159_id <- c("DSP-1001660013520-C-A03", "DSP-1001660013520-C-A05", "DSP-1001660013520-C-A07",
                "DSP-1001660013520-C-A09", "DSP-1001660013520-C-A12", "DSP-1001660013520-C-B04")
geo_malig$malig_sub <- ifelse(geo_malig$sample_id2 %in% clus14_id, "clus14",
                              ifelse(geo_malig$sample_id2 %in% clus159_id, "clus159", NA))

# DE ----------------------------------------------------------------------
sce <- geo_malig
# Remove negative control probe
sce <- sce[rownames(sce) != "NegProbe-WTX", ]

dge <- SE2DGEList(sce)
design <- model.matrix(~0 + malig_sub, data = colData(sce))
colnames(design) <- gsub(" ","_",colnames(design))
colnames(design)
#  "malig_subclus14"  "malig_subclus159"
table(sce$malig_sub)
# clus14 clus159 
# 3      6 

contr.matrix <- makeContrasts(
  BvT = malig_subclus159 - malig_subclus14,
  levels = colnames(design))

keep <- filterByExpr(dge, design)
table(keep)
rownames(dge)[!keep]
dge_all <- dge[keep, ]

# Voom DE
v <- voom(dge_all, design, plot = TRUE) 
fit <- lmFit(v)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit<- limma::decideTests(efit, p.value = 0.05, lfc = 1)
summary_efit <- summary(results_efit)

summary_efit
# BvT
# Down      12
# NotSig 18619
# Up        26

de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, 
                                  p.value = 0.05, lfc = 1)
results <- de_genes_toptable_BvT
save_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom"
write.csv(results, file.path(save_path, "Geo_3vs6_DE.csv"))

save_path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/SameDE"
saveRDS(v, file.path(save_path, "geo_voom.rds"))
saveRDS(contr.matrix, file.path(save_path, "geo_contr_mat.rds"))
