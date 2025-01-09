library(dplyr)
library(tidyverse)
library(patchwork)

# Merge in patho annotation ----------------------------------------------
patho_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium/Annotations/"
files_to_read <- list.files(patho_path)

patho_all <- NULL
for(i in 1:length(files_to_read)){
  patho_i <- read.csv(file.path(patho_path, files_to_read[i]))
  patho_i$Section <- gsub(".csv", "", files_to_read[i])
  print(dim(patho_i))
  
  patho_all <- rbind(patho_all, patho_i)
}

table(patho_all$Section)
# B1_2    B1_4    B2_2    B3_2    B4_2 DLBCL_1 DLBCL_2 DLBCL_3 DLBCL_4 DLBCL_5 DLBCL_6    L1_2    L1_4    L2_2    L3_2    L4_2 
# 1897    1999    2569    2220    1850    4710    4121    4951    1460    1778    1743     874     944     842    2443    2955 


# Use CARD level 4 for now ----------------------------------------------
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/CARD"
vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Visium_Decon/Level4/C2L"

files_to_read <- list.files(vis_decon_path)

decon_all <- NULL
for(i in 1:length(files_to_read)){
  decon_i <- read.csv(file.path(vis_decon_path, files_to_read[i])) %>%
    dplyr::rename(Barcode = X)
  
  # if(files_to_read[i] == "B1_2_spot_Level4_decon.csv"){ # CARD
  if(files_to_read[i] == "B1_2.csv"){                     # C2L
    sce_qcd <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_qcd/breast_qcd/B1_2_qcd.rds")
    decon_i <- decon_i[decon_i$Barcode %in% sce_qcd$Barcode, ]
  }
  
  print(dim(decon_i))

  decon_patho_long_i <- pivot_longer(decon_i, names_to = "CellType", values_to = "Fraction", -c("Barcode"))
  # decon_patho_long_i$Section <- gsub("_spot_Level4_decon.csv", "", files_to_read[i]) # CARD
  decon_patho_long_i$Section <- gsub(".csv", "", files_to_read[i]) # C2L
  
  decon_all <- rbind(decon_all, decon_patho_long_i)
}



# Merged -------------------------------------------------------------------
decon_patho_long_all <- decon_all %>%
  left_join(patho_all, by = c("Barcode", "Section")) %>%
  mutate(Region = case_when(Region == "Vessel" ~ "Vessels",
                            Region == "Most_likely_Tumor" ~ "Most_likely_tumor",
                            TRUE ~ Region)) %>%
  filter(!is.na(Region))

# Combine Categories to Fig 2 ---------------------------------------------
## Breast and lung 
tumor_related <- c("Breast_Carcinoma_in_situ", "Tumor_pure", "Most_likely_tumor", "Most_likely_Tumor")
lymphocytes_related <- c("Lymphocytes")
macrophage_related <- c("Immune_Cell_mix")
stroma_related <- c("Large_Vessel", "Intratumoral_Stroma", "Intratumoral_Vessel", "Normal_Lung_Epithelium")
breast_lung_patho <- c(tumor_related, lymphocytes_related, macrophage_related, stroma_related)

## DLBCL
d_tumor_related <- c("Tumor")
d_lymphocytes_related <- c("Small lymphocytes")
d_macrophage_related <- c("Epithelium")
d_stroma_related <- c("Stroma", "Vessels")
dlbcl_patho <- c(d_tumor_related, d_lymphocytes_related, d_macrophage_related, d_stroma_related)


# -------------------------------------------------------------------------
decon_patho_long_all_breast <- decon_patho_long_all %>%
  mutate(CellType = case_when(CellType %in% c("B_cell", "B_plasma", "B_plasma_IGHA1", "B_plasma_IGHG1", "B_plasma_IGHG3", "B_plasma_IGHM", "B_plasma_IGKC", "B_plasma_IGLC1") ~ "B cells",
                              CellType %in% c("Endothelia_vascular", "Fibroblast", "Fibroblast_B3", "Muscle_smooth", "Pericyte") ~ "Stroma",
                              CellType %in% c("Tu_B1_MUCL1", "Tu_B1_MUCL1_necrosis", "Tu_B1_MUCL1_transcription", "Tu_B2", "Tu_B3_CYP4F8", "Tu_B4_RHOB") ~ "Tumor",
                              CellType %in% c("DC_1", "DC_2", "DC_activated", "DC_pc", "Granulocyte", "Mast_cell") ~ "Myeloid else",
                              CellType == "NK" ~ "NK",
                              CellType %in% c("Macrophage", "Monocyte") ~ "Macrophage",
                              CellType %in% c("T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing") ~ "T cells")) %>%
  filter(Region %in% breast_lung_patho) %>%
  mutate(Region = case_when(Region %in% tumor_related ~ "Tumor-enriched",
                                Region %in% lymphocytes_related ~ "Lymphocytes-enriched",
                                Region %in% macrophage_related ~ "Immune Cell Mix", 
                                Region %in% stroma_related ~ "Stroma-enriched"))






CT_order <- c("Tumor", "Stroma", "Epithelia", "Macrophage", "T cells", "B cells", "NK", "Myeloid else") 
PA_order <- c("Tumor-enriched", "Stroma-enriched", "Lymphocytes-enriched", "Immune Cell Mix") 

CT_order <- c("Tumor", "Stroma", "Macrophage", "T cells", "B cells", "NK", "Myeloid else") 
PA_order <- c("Tumor-enriched", "Stroma-enriched", "Lymphocytes-enriched", "Immune Cell Mix") 


decon_patho_long_all_dlbcl <- decon_patho_long_all %>%
  mutate(CellType = case_when(CellType == "B_plasma" ~ "B cells",
                              CellType %in% c("Epi_bronchus", "Epi_lung", "Epi_mucous", "Epi_Mucous_surface_gastric", "Epi_parietal_gastric") ~ "Stroma",
                              CellType %in% c( "Tu_D1_LMO2", "Tu_D1_RGS13", "Tu_D1_SMIM14", "Tu_D2_mito", "Tu_D3_BCL2A1", "Tu_D3_dividing", 
                                               "Tu_D3_FAM3C", "Tu_D3_IGHD", "Tu_D4_BCL7A", "Tu_D4_PNN", "Tu_D5_CCL22", "Tu_D6_BCL2") ~ "Tumor",
                              CellType %in% c("DC_1", "DC_2", "DC_pc") ~ "Myeloid else",
                              CellType == "NK" ~ "NK",
                              CellType == "Mono_Macro" ~ "Macrophage",
                              CellType %in% c("T_CD4", "T_CD4_reg", "T_CD8", "T_dividing") ~ "T cells")) %>%
  filter(Region %in% breast_lung_patho) %>%
  mutate(Region = case_when(Region %in% tumor_related ~ "Tumor-enriched",
                            Region %in% lymphocytes_related ~ "Lymphocytes-enriched",
                            Region %in% macrophage_related ~ "Epithelium", 
                            Region %in% stroma_related ~ "Stroma-enriched"))


CT_order <- c("Tumor", "Stroma", "Epithelia", "Macrophage", "T cells", "B cells", "NK", "Myeloid else") 
PA_order <- c("Tumor-enriched", "Stroma-enriched", "Lymphocytes-enriched", "Epithelium") 

decon_patho_long_all_test$CellType <- factor(decon_patho_long_all_test$CellType, levels = CT_order)
decon_patho_long_all_test$Region <- factor(decon_patho_long_all_test$Region, levels = PA_order)


# -------------------------------------------------------------------------
plotHeatmap <- function(decon_patho_long_all, disease = "breast"){ 
  if(disease == "breast"){
    sample_list = c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2")
  }else if(disease == "lung"){
    sample_list = c("L1_2", "L1_4", "L2_2", "L3_2", "L4_2")
  }else{ # dlbcl
    sample_list = c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
  }
  
  ## Per region per cell type
  stats <- decon_patho_long_all %>%
    filter(Section %in% sample_list) %>%
    group_by(Region, CellType) %>%
    summarize(mean_value = mean(Fraction), .groups = "drop") %>%
    mutate(per = round(100 * mean_value, 0))
  
  ## Per cell type sum to 1
  stats_ <- stats %>%
    group_by(CellType) %>%
    mutate(countT= sum(mean_value)) %>%
    mutate(per=round(100*mean_value/countT,0))
  
  stats_$CellType <- factor(stats_$CellType, levels = rev(CT_order))
  stats_$Region <- factor(stats_$Region, levels = PA_order)
  
  # Create the heatmap Region per decon
  p <- ggplot(stats_, aes(x = CellType, y = Region, fill = per)) + 
    geom_tile() + 
    geom_text(aes(x = CellType, y = Region, label = paste0(per, "%"), size = )) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colors = c("white", "dodgerblue")) + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(angle = 45, size = 11, hjust = 0, vjust = 1)
    ) + 
    labs(fill = paste0("Average deconvolution fraction per-region weighted by per-cell-type"), title = "", y = "") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
  
  return(p)
}

# -------------------------------------------------------------------------
p_breast_patho_decon <- plotHeatmap(decon_patho_long_all_breast, disease = "breast")
p_lung_patho_decon <- plotHeatmap(decon_patho_long_all_lung, disease = "lung")
p_dlbcl_patho_decon <- plotHeatmap(decon_patho_long_all_dlbcl, disease = "dlbcl")

p_final <- ((p_breast_patho_decon + theme(axis.title.x = element_blank())) | 
  (p_lung_patho_decon + theme(axis.title.x = element_blank())) | 
  (p_dlbcl_patho_decon + theme(axis.title.x = element_blank()))) + 
  plot_layout(# guides = "collect", 
    widths = c(1.7, 1.5, 0.9)) & 
  theme(legend.position = "bottom" , 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.margin = margin(1, 20, 1, 1))



figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Patho_Decon_agreement"

pdf(file.path(figpath, "Vis_Patho_Decon_Heatmap_final_fraction_percelltype_final_C2L.pdf"), width = 20, height = 7)
print(p_final)
dev.off()


