## sample example 
# # subspot joint decont all
# test <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_joint_decont/Lung/L1_2/L1_2_baye_clustered_all_enhanced_expr.rds")
# 
# # spot joint decont all
# test <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_joint_decont/Lung/L1_2/L1_2_baye_clustered.rds")

# -------------------------------------------------------------------------
disease = "lung"
# disease = "breast"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
samples_with_anno <- c("B1_4", "B2_2", "B3_2", "B4_2", "L1_2", "L1_4", "L2_2", "L3_2", "L4_2")

for(i in 1:nsamples){
  print(i)
  if(save_names[i] %in% samples_with_anno){
    ## Get Spot SCE by sample ---------------------------------------------
    filepath <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_joint_decont/", foldername, "/", save_names[i])
    sce <- readRDS(file.path(filepath, paste0(save_names[i], "_baye_clustered.rds")))
    sce$barcode <- colnames(sce)

    # Map in Pathology annotation -----------------------------------------
    if(disease == "breast"){
    vis_anno <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/", patho_filenames[i]))
    vis_anno_sample <- vis_anno %>%
      select(Barcode, CS) %>%
      rename(barcode = Barcode,
             Region = CS) # 1999   17
    
    }else if(disease == "lung"){
      if(save_names[i] == "L3_2"){
        vis_anno <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/", "L3_2.csv"))
        vis_anno_sample <- vis_anno %>%
          rename(barcode = Barcode,
                 Region = KvL)  %>%
          mutate(barcode = paste0("L3_2", "_", barcode)) %>%
          select(barcode, Region) 
      }else{
        vis_anno <- read_excel("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/Visium_combined_annotations_CS_ED.xlsx")
        vis_anno_sample <- vis_anno %>%
          filter(Section == save_names[i]) %>%
          rename(barcode = Barcode) %>%
          mutate(barcode = paste0(Section, "_", barcode)) %>%
          select(barcode, Region) 
      }
    }


    # Merge in spot -------------------------------------------------------
    CD_spot <- data.frame(colData(sce))

    CD_spot <- CD_spot %>%
      left_join(vis_anno_sample, by = "barcode")
    
    colData(sce) <- as(CD_spot, "DFrame")
    colnames(sce) <- CD_spot$barcode

    # print(length(which(is.na(sce$Region))))
    # table(sce$Region)
    if(save_names[i] == "B2_2"){
      sce <- sce[, sce$Region != "Fold_Exclude"] # B2_2, since new patho annotation update ...
    }else if(save_names[i] %in% c("B3_2", "B4_2")){
      sce <- sce[, sce$Region != "Exclude"] # B3_2, B4_2 since new patho annotation update ...
    }else if(save_names[i] == "L3_2"){
      sce <- sce[, sce$Region != "Artefact_Fold_exclude"] # L3_2 since new patho annotation ...
    }
    
    saveRDS(sce, file.path(filepath, paste0(save_names[i], "_baye_clustered.rds")))
  }
}