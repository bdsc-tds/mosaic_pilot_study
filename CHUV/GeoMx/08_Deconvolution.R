source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")
library(stringr)
library(tidyr)
library(dplyr)

# chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
sigmatpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5_immune/"
deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/For_level1_5_immune_decon_results"

###########################################################################
#                               Deconvolution                             #
###########################################################################
# disease = "breast"
# disease = "lung"
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/00_GeoMx_Paths.R")
sample_name <- ifelse(disease == "dlbcl", "patient", "section_id")
spe_ruv <- readRDS(file.path(datapath, paste0(disease, "_spe_ruv.rds")))


# Signature matrix --------------------------------------------------------
if(disease == "breast"){
  # chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
  sigmat_name = "sigmat_breast_add_lung_healthy.rds"
  CF_order <- c("Malignant", "Other", "T cells", "Macrophage")                                            # Fig 2, 3 
}else if(disease == "lung"){
  # chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds")) 
  sigmat_name = "sigmat_lung_add_breast_healthy.rds"
  CF_order <- c("Malignant", "Other", "T cells", "Macrophage")                                            # Fig 2, 3 
}else{
  # chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds")); chrom$patient <- chrom$sample_id
  sigmat_name = "sigmat_dlbcl.rds"
  # CT_order <- c("Epithelia", "Stroma", "B cells", "NK", "Myeloid else",  "Macrophage", "T cells", "Tumor") # Fig 6
  CF_order <- c("B cells", "Other", "T cells", "Macrophage")                                               # Fig 2, 3 
}
CT_order <- c("Epithelia", "Stroma", "Tumor", "Macrophage", "T cells", "B cells", "NK", "Myeloid else")  # Fig 2, 3 
new_ref_matrix <- readRDS(file.path(sigmatpath, sigmat_name))


# Prep decon --------------------------------------------------------------
# NegProbe name
negProbeName <- rownames(spe_ruv)[which(grepl("Neg", rownames(spe_ruv)))]

## Normed Batched
# force the normed, batch corrected count matrix to be not a df, but just matrix
assay(spe_ruv, "logcounts") <- as.matrix(assay(spe_ruv, "logcounts"))
assay(spe_ruv, "quantile_batched_linear_scale") <- as.matrix(exp(assay(spe_ruv, "logcounts")) - 1)
spd_normed_batched <- prepareSpatialDecon(spe_ruv, assay2use = "quantile_batched_linear_scale", negProbeName = negProbeName)


# Deconvolution on normalized and batch corrected object ------------------
library(SpatialDecon)
res_norm_batch <- spatialdecon(norm = spd_normed_batched$normCount,
                               bg = spd_normed_batched$backGround,
                               X = new_ref_matrix,
                               align_genes = TRUE)

subset_prop <- res_norm_batch$prop_of_all


# Gather decon df ---------------------------------------------------------
df <- data.frame(t(subset_prop))
df$sample <- rownames(df)

gathered_df <- gather(df, key = "column", value = "value", -sample)

CD <- as.data.frame(colData(spe_ruv)[, c("sample_id2", sample_name, "cell_fraction")])
CD_ <- CD %>% dplyr::rename(sample = sample_id2)

gathered_df_ <- gathered_df %>%
  left_join(CD_, by = "sample") %>%
  filter(cell_fraction != "PanCK-") %>%
  dplyr::rename(CellType = column,
         Fraction = value) %>%
  mutate(cell_fraction = ifelse(cell_fraction == "Macro", "Macrophage", cell_fraction),
         CellType = str_replace(CellType, "\\.", " ")) 


gathered_df_$cell_fraction <- factor(gathered_df_$cell_fraction, levels = CF_order)
gathered_df_$CellType <- factor(gathered_df_$CellType, levels = CT_order)

write.csv(gathered_df_, 
          file.path(deconresultpath, paste0(disease, "_batched_decon_long.csv")), 
          row.names = FALSE)


# Plot ------------------------------------------------------------------
p <- ggplot(gathered_df_, aes(x=CellType, y=Fraction)) +
  geom_boxplot() +
  # geom_jitter(aes(color=section_id), size=1, alpha=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust=1),
        panel.spacing=unit(1.5,"lines"),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 13.5, face = "bold"),
        strip.background=element_rect(fill="#DEDEDE")) +
  facet_wrap(~cell_fraction, ncol = 5) +
  guides(col = guide_legend(override.aes = list(size = 2))) +
  labs(color = "Section", x = "", y = "Cell type fraction") +
  ggtitle("")

p


