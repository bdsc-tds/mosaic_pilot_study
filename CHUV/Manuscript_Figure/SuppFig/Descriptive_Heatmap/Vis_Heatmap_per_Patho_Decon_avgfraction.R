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


# -------------------------------------------------------------------------
plotHeatmap <- function(decon_patho_long_all, disease = "breast"){ 
  if(disease == "breast"){
    sample_list = c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2")
  }else if(disease == "lung"){
    sample_list = c("L1_2", "L1_4", "L2_2", "L3_2", "L4_2")
  }else{ # dlbcl
    sample_list = c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
  }
  
  ## Per region sum to 1
  stats <- decon_patho_long_all %>%
    filter(Section %in% sample_list) %>%
    group_by(Region, CellType) %>%
    summarize(mean_value = mean(Fraction), .groups = "drop") %>%
    mutate(per = round(100 * mean_value, 0))
  
  ## Per cell type sum to 1?
  stats_ <- stats %>%
    group_by(CellType) %>%
    mutate(countT= sum(mean_value)) %>%
    mutate(per=round(100*mean_value/countT,0))
  
  stats_$CellType <- factor(stats_$CellType, levels = rev(levels(factor(stats_$CellType))))
  
  p <- ggplot(stats_, aes(x = Region, y = CellType, fill = per)) + 
    geom_tile() + 
    geom_text(aes(x = Region, y = CellType, label = paste0(per, "%"))) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left") +
    scale_fill_viridis_c(option = "mako") + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(angle = 45, size = 11, hjust = 0, vjust = 1)
    ) + 
    labs(fill = paste0("Average deconvolution fraction per-region per-cell-type"), title = "", y = "") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
  
  return(p)
}

# -------------------------------------------------------------------------
p_breast_patho_decon <- plotHeatmap(decon_patho_long_all, disease = "breast")
p_lung_patho_decon <- plotHeatmap(decon_patho_long_all, disease = "lung")
p_dlbcl_patho_decon <- plotHeatmap(decon_patho_long_all, disease = "dlbcl")

p_final <- ((p_breast_patho_decon + theme(axis.title.x = element_blank())) | 
  (p_lung_patho_decon + theme(axis.title.x = element_blank())) | 
  (p_dlbcl_patho_decon + theme(axis.title.x = element_blank()))) + 
  plot_layout(# guides = "collect", 
    widths = c(1.7, 1.5, 0.9)) & 
  theme(legend.position = "bottom" , 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.margin = margin(1, 20, 1, 1))



figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Gallery_Patho_Decon_Pt_Heatmap"

pdf(file.path(figpath, "Vis_Patho_Decon_Heatmap_final_fraction_percelltype_final_C2L.pdf"), width = 20, height = 14)
print(p_final)
dev.off()
