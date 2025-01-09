
###########################################################################
## No batch correction decon 
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



###########################################################################
## Batch corrected decon
deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/For_level1_5_immune_decon_results"
disease <- "lung"
gathered_df_ <- read.csv(file.path(deconresultpath, paste0(disease, "_batched_decon_long.csv")), row.names = FALSE)

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



