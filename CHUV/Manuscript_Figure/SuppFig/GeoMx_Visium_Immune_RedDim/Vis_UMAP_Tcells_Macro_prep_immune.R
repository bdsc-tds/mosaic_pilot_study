library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(tidyr)

# -------------------------------------------------------------------------
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2_final_immune" # CARD
vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L/level1_5_immune_long" # C2L - level 4 combined to level 1.5 immunne 


# After spotclean ---------------------------------------------------------
get_Decon_DF_to_Merge <- function(vis_decon_path, vis_sample_list, filename){
  decon_gather <- read.csv(file.path(vis_decon_path, filename))
  
  vis_decon <- decon_gather %>%
    mutate(Section_idx = case_when(Section == vis_sample_list[1] ~ "1",
                                   Section == vis_sample_list[2] ~ "2",
                                   Section == vis_sample_list[3] ~ "3",
                                   Section == vis_sample_list[4] ~ "4",
                                   Section == vis_sample_list[5] ~ "5"
                                   ,Section == vis_sample_list[6] ~ "6"
    ),
    Barcode = paste0(Barcode, "_", Section_idx)) %>%
    filter(CellType %in% c("T cells", "Macrophage")) %>%
    select(Barcode, CellType, Fraction) %>%
    pivot_wider(names_from = CellType, values_from = Fraction)
  
  return(vis_decon)
}


saveToSeuFracImmune <- function(vis_decon_path, disease){
  datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
  savepath.spotclean = paste0(datapath.spotclean, "/Results")
  save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
  
  seu <- readRDS(paste0(savepath.spotclean, save_rds_name))
  
  
  # -------------------------------------------------------------------------
  if(disease == "breast"){
    vis_sample_list <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2")
  }else if(disease == "lung"){
    vis_sample_list <- c("L1_2", "L1_4", "L2_2", "L3_2", "L4_2")
  }else{ # dlbcl
    vis_sample_list <- c("D1", "D2", "D3", "D4", "D5", "D6")
  }
  

  decon_df <- get_Decon_DF_to_Merge(vis_decon_path, 
                                    vis_sample_list,
                                    filename = paste0("vis_", disease, "_decon_immune_long.csv"))
  
  
  # Merge in to seu -------------------------------
  seu <- seu[, colnames(seu) %in% decon_df$Barcode]
  seuCD <- data.frame(Barcode = colnames(seu))
  
  seuCD_ <- seuCD %>%
    left_join(decon_df)
  
  seu$Macrophage_frac <- seuCD_$Macrophage
  seu$Tcells_frac <- seuCD_$`T cells`
  
  
  saveRDS(seu, paste0(savepath.spotclean, save_rds_name))
}

saveToSeuFracImmune(vis_decon_path, disease = "breast")
saveToSeuFracImmune(vis_decon_path, disease = "lung")
saveToSeuFracImmune(vis_decon_path, disease = "dlbcl")




