library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(dplyr)
library(tidyr)

vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Visium_Decon/Level4/C2L/" # C2L - level 4  
vis_patho_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium/Annotations/" # C2L - level 4  

breast_samples <- c("B1_2.csv", "B1_4.csv", "B2_2.csv", "B3_2.csv", "B4_2.csv")
lung_samples <- c("L1_2.csv", "L1_4.csv", "L2_2.csv", "L3_2.csv", "L4_2.csv")
dlbcl_samples <- c("DLBCL_1.csv", "DLBCL_2.csv", "DLBCL_3.csv", "DLBCL_4.csv", "DLBCL_5.csv", "DLBCL_6.csv")

getDeconPatho_DF <- function(samples){
  decon_long <- NULL
  for(i in 1:length(samples)){
    decon_i <- read.csv(file.path(vis_decon_path, samples[i]))
    patho_i <- read.csv(file.path(vis_patho_path, samples[i]))
    
    if(samples[i] %in% breast_samples & samples[i] != "B3_2.csv"){
      if("Fibroblast_B3" %in% names(decon_i)){
        decon_i <- decon_i %>%
          mutate(Fibroblast = Fibroblast + Fibroblast_B3) %>%
          select(-c(Fibroblast_B3))
      }
    }
    
    decon_i_long <- decon_i %>%
      pivot_longer(cols = -X, names_to = "CellType", values_to = "Fraction") %>%
      dplyr::rename(Barcode = X) %>%
      inner_join(patho_i %>% mutate(Region = ifelse(Region == "Vessels", "Vessel", Region)), by = "Barcode") %>%
      mutate(Sample = samples[i]) %>%
      filter(!(Region %in% c("Artefact_Fold_exclude", "Exclude")))
    
    decon_long <- rbind(decon_long, decon_i_long)
  }
  
  return(decon_long)
}

breast_decon_patho = getDeconPatho_DF(samples = breast_samples); 
breast_CT_order <- c("B_cell", "B_plasma_IGHA1", "B_plasma_IGHG1", "B_plasma_IGHG3", "B_plasma_IGHM", "B_plasma_IGKC", "B_plasma_IGLC1", "DC_1", "DC_2", "DC_activated", "DC_pc", "Endothelia_vascular", "Fibroblast", "Fibroblast_B3", "Granulocyte", "Macrophage",               
                     "Mast_cell", "Monocyte", "Muscle_smooth", "NK", "Pericyte", "T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing", "Tu_B1_MUCL1", "Tu_B1_MUCL1_necrosis", "Tu_B1_MUCL1_transcription", "Tu_B2", "Tu_B3_CYP4F8", "Tu_B4_RHOB")

lung_decon_patho = getDeconPatho_DF(samples = lung_samples)
lung_CT_order <- c("B_cell", "B_plasma_IGHA1", "B_plasma_IGHG1", "B_plasma_IGHG3", "B_plasma_IGHM", "B_plasma_IGKC", "B_plasma_IGLC1", "DC_1", "DC_2", "DC_activated", "DC_pc", "Endothelia_vascular", "Epi_lung", "Fibroblast", "Granulocyte", "Macrophage", "Mast_cell", 
                   "Monocyte", "Muscle_smooth", "NK", "Pericyte", "T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing", "Tu_L1_SFTPB", "Tu_L2_FXYD2", "Tu_L3_G0S2", "Tu_L3_G0S2_immune_signature", "Tu_L4_KRT17_immune_signature", "Tu_L4_KRT17_mucous",
                   "Tu_L4_KRT17_necrosis", "Tu_L4_KRT17_neutrophil_signature")

dlbcl_decon_patho = getDeconPatho_DF(samples = dlbcl_samples)
dlbcl_CT_order <- c("B_plasma", "DC_1", "DC_2", "DC_pc", "Endothelia", "Epi_bronchus", "Epi_mucous", "Epi_Mucous_surface_gastric", "Epi_parietal_gastric", "Fibro_Muscle", "Mono_Macro", "NK", "Pericyte", "T_CD4", "T_CD4_reg", "T_CD8", "T_dividing", "Tu_D1_LMO2", 
                    "Tu_D1_RGS13", "Tu_D1_SMIM14", "Tu_D2_mito", "Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD", "Tu_D4_BCL7A", "Tu_D4_PNN", "Tu_D5_CCL22", "Tu_D6_BCL2")

# -------------------------------------------------------------------------
plotHeatmap <- function(decon_patho_long_all, CT_order){ 
  ## Per region per cell type
  stats <- decon_patho_long_all %>%
    group_by(Region, CellType) %>%
    summarize(mean_value = mean(Fraction), .groups = "drop") %>%
    mutate(per = round(100 * mean_value, 0))
  
  stats$CellType <- factor(stats$CellType, levels = rev(CT_order))
  
  p <- ggplot(stats, aes(x = Region, y = CellType, fill = per)) + 
    geom_tile() + 
    geom_text(aes(x = Region, y = CellType, label = paste0(per, "%"), size = )) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colors = c("white", "dodgerblue"), transform = "sqrt") + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(angle = 45, size = 11, hjust = 0, vjust = 1)
    ) + 
    labs(fill = paste0("Avg. cell type deconvolution fraction per pathology label"), title = "", y = "") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
  
  return(p)
}


# -------------------------------------------------------------------------
p_breast_patho_decon <- plotHeatmap(breast_decon_patho, breast_CT_order)
p_lung_patho_decon <- plotHeatmap(lung_decon_patho, lung_CT_order)
p_dlbcl_patho_decon <- plotHeatmap(dlbcl_decon_patho, dlbcl_CT_order)

p_final <- p_breast_patho_decon | p_lung_patho_decon | p_dlbcl_patho_decon + 
  plot_layout(ncol = 2, widths = c(1.3, 1.1, 2)) & # so weird that I need to specify ncol = 2 to have 3 width correct
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 8),        
        legend.title = element_text(size = 10))
p_final

df_FigS8a <- rbind(breast_decon_patho, lung_decon_patho, dlbcl_decon_patho)
write.csv(df_FigS8a, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS8a_patho_decon_vis.csv")

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_Vis_Heatmap/SuppFig"
pdf(file.path(fig_path, "Vis_Patho_Decon_Heatmap_level4.pdf"), width = 22, height = 13)
print(p_final)
dev.off()


