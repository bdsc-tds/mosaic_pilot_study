source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")
deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/"

# -------------------------------------------------------------------------
gathered_df_breast <- read.csv(file.path(deconresultpath, "breast_batched_decon_long.csv")) 
gathered_df_lung <- read.csv(file.path(deconresultpath, "lung_batched_decon_long.csv"))
gathered_df_dlbcl <- read.csv(file.path(deconresultpath, "dlbcl_batched_decon_long.csv"))
gathered_df_dlbcl$cell_fraction <- ifelse(gathered_df_dlbcl$cell_fraction == "B cells", "Malignant", gathered_df_dlbcl$cell_fraction)

CF_order <- c("Malignant", "T cells", "Macrophage", "Other")
# CT_order <- c("Tumor", "T cells", "NK", "Macrophage", "Myeloid else", "B cells", "Stroma", "Epithelia")

breast_CT_order <- c("B_cell", "B_plasma_IGHA1", "B_plasma_IGHG1", "B_plasma_IGHG3", "B_plasma_IGHM", "B_plasma_IGKC", "B_plasma_IGLC1", "DC_1", "DC_2", "DC_activated", "DC_pc", "Endothelia_vascular", "Fibroblast", "Fibroblast_B3", "Granulocyte", "Macrophage",               
                     "Mast_cell", "Monocyte", "Muscle_smooth", "NK", "Pericyte", "T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing", "Tu_B1_MUCL1", "Tu_B1_MUCL1_necrosis", "Tu_B1_MUCL1_transcription", "Tu_B2", "Tu_B3_CYP4F8", "Tu_B4_RHOB")

lung_CT_order <- c("B_cell", "B_plasma_IGHA1", "B_plasma_IGHG1", "B_plasma_IGHG3", "B_plasma_IGHM", "B_plasma_IGKC", "B_plasma_IGLC1", "DC_1", "DC_2", "DC_activated", "DC_pc", "Endothelia_vascular", "Epi_lung", "Fibroblast", "Granulocyte", "Macrophage", "Mast_cell", 
                   "Monocyte", "Muscle_smooth", "NK", "Pericyte", "T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing", "Tu_L1_SFTPB", "Tu_L2_FXYD2", "Tu_L3_G0S2", "Tu_L3_G0S2_immune_signature", "Tu_L4_KRT17_immune_signature", "Tu_L4_KRT17_mucous",
                   "Tu_L4_KRT17_necrosis", "Tu_L4_KRT17_neutrophil_signature")

dlbcl_CT_order <- c("B_plasma", "DC_1", "DC_2", "DC_pc", "Endothelia", "Epi_bronchus", "Epi_mucous", "Epi_Mucous_surface_gastric", "Epi_parietal_gastric", "Fibro_Muscle", "Mono_Macro", "NK", "Pericyte", "T_CD4", "T_CD4_reg", "T_CD8", "T_dividing", "Tu_D1_LMO2", 
                    "Tu_D1_RGS13", "Tu_D1_SMIM14", "Tu_D2_mito", "Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD", "Tu_D4_BCL7A", "Tu_D4_PNN", "Tu_D5_CCL22", "Tu_D6_BCL2")


# -------------------------------------------------------------------------
plotHeatmap <- function(gathered_df, CT_order){ 
  ## Per cell_fraction sum to 1
  stats <- gathered_df %>%
    group_by(cell_fraction, CellType) %>%
    summarize(mean_value = mean(Fraction), .groups = "drop") %>%
    mutate(per = round(100 * mean_value, 0))
  stats$CellType <- factor(stats$CellType, levels = rev(CT_order))
  stats$cell_fraction <- factor(stats$cell_fraction, levels = CF_order)

  p <- ggplot(stats, aes(x = cell_fraction, y = CellType, fill = per)) + 
    geom_tile() + 
    geom_text(aes(x = cell_fraction, y = CellType, label = paste0(per, "%"))) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colors = c("white", "#21751F"), transform = "sqrt") + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(angle = 45, size = 11, hjust = 1, vjust = 1)
    ) + 
    labs(fill = paste0("Avg. cell type deconvolution fraction per AOI label"), title = "", y = "") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
  
  
  
  return(p)
}

# -------------------------------------------------------------------------
p_breast_patho_decon <- plotHeatmap(gathered_df_breast %>% filter(cell_fraction != "PanCK-"), breast_CT_order) 
p_lung_patho_decon <- plotHeatmap(gathered_df_lung %>% filter(cell_fraction != "PanCK-"), lung_CT_order)
p_dlbcl_patho_decon <- plotHeatmap(gathered_df_dlbcl %>% filter(cell_fraction != "PanCK-"), dlbcl_CT_order)

p_final <- p_breast_patho_decon | p_lung_patho_decon | p_dlbcl_patho_decon + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 8),        
        legend.title = element_text(size = 10))
p_final

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_Vis_Heatmap/SuppFig"
pdf(file = file.path(fig_path, "Geo_AOI_Decon_Heatmap_level4.pdf"),
    width = 22,
    height = 13)
print(p_final)
dev.off()



