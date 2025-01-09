# Breast Lung ----------------------------------------
br_lu_level4_cellnames <- c(
  "B_cell",                           "B_plasma_IGHA1",                   "B_plasma_IGHG1",                   "B_plasma_IGHG3",
  "B_plasma_IGHM",                    "B_plasma_IGKC",                    "B_plasma_IGLC1",                   "DC_1",
  "DC_2",                             "DC_activated",                     "DC_pc",                            "Endothelia_lymphatic",
  "Endothelia_vascular",              "Epi_lung",                         "Fibroblast",                       "Fibroblast_B3",
  
  "Granulocyte",                      "Macrophage",                       "Mast_cell",                        "Monocyte",
  "Muscle_smooth",                    "NK",                               "Pericyte",                         "T_CD4",
  "T_CD8_exhausted",                  "T_CTL",                            "T_CXCL13",                         "T_reg",
  "TNK_dividing",
  
  "Tu_B1_MUCL1",                      "Tu_B1_MUCL1_necrosis",             "Tu_B1_MUCL1_transcription",        
  #"Tu_B3_CYP4F8", 
  "Tu_B3",                            # "Tu_B3_necrosis",                   "Tu_B3_NPPC",                       "Tu_B3_PLA2G2A",
  "Tu_B4_RHOB",                       "Tu_B2",
  
  "Tu_L1_SFTPB",                      "Tu_L2_FXYD2",                      "Tu_L3_G0S2",                       "Tu_L3_G0S2_immune_signature",
  "Tu_L4_KRT17_immune_signature",     "Tu_L4_KRT17_mucous",               "Tu_L4_KRT17_necrosis",             "Tu_L4_KRT17_neutrophil_signature",
  "Mix"
)

br_lu_level4_cellcolors <- c(
  "#EEEE00", "#FF00FF", "#CD00CD", "#DA70D6", "#FF83FA", "#8B4789", "#B452CD", "#C71585", "#FFAEB9", "#CD8C95", "#FF6347", "#FFDAB9","#FF7256", "#BC8F8F", "#388E8E","#008080",
  
  "#FFC125", "#FF9900", "#FFD700", "#FF8000", "#668B8B", "#79CDCD", "#A0522D", "#4169E1", "#00BFFF", "#009ACD", "#3D59AB", "#87CEFF", "#97FFFF",
  
  "#A2CD5A", "#00EE76", "#ADFF2F", # Breast
  "#EEEE00", 
  #"#D4FF90", "#C9DC84", "#DBE153", "#9FCB00", 
  "#FFD700", "#66FF66", 
  
  "#EEB4B4", "#FFA54F", "#CDAF95", "#CD853F", "#8B4513", "#CD661D", "#F4A460", "#8B2500", # Lung
  "#5A5A5A" # Mix
)

names(br_lu_level4_cellcolors) <- br_lu_level4_cellnames


# -------------------------------------------------------------------------
level1_5_cellnames <- c("B", "Myeloid", "Vessel", "Fibro_muscle", "Granulocyte", "Epithelia", "T_NK", "Tumor")
br_level1_5_cellcolors <- c("#EEEE00", "#9A32CD", "#FF7256", "#388E8E", "#FFC125", "#BC8F8F", "#4169E1", "#00C78C")
lu_level1_5_cellcolors <- c("#EEEE00", "#9A32CD", "#FF7256", "#388E8E", "#FFC125", "#BC8F8F", "#4169E1", "#FFA54F")
names(br_level1_5_cellcolors) <- level1_5_cellnames
names(lu_level1_5_cellcolors) <- level1_5_cellnames


# DLBCL -----------------------------------------------
dlbcl_level4_cellnames <- c(
  "B_plasma",     "DC_1",           "DC_2",        "DC_pc",
  "Endothelia",   "Epi_bronchus",   "Epi_mucous",  "Epi_Mucous_surface_gastric", "Epi_parietal_gastric", "Fibro_Muscle",
  "Mono_Macro",   "NK",             "Pericyte",    "T_CD4",                      "T_CD4_reg",            "T_CD8",         "T_dividing",
  "Tu_D1_LMO2",   "Tu_D1_RGS13",  "Tu_D1_SMIM14",
  "Tu_D2_mito",
  "Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD",
  "Tu_D4_BCL7A",  "Tu_D4_PNN",
  "Tu_D5_CCL22",
  "Tu_D6_BCL2",
  "Mix"
)

dlbcl_level4_cellcolors <- c(
  "#FF00FF", "#C71585", "#FFAEB9", "#FF6347",                                  # B, DCs
  "#FF7256", "#BC8F8F", "#A52A2A", "#FFC1C1", "#CD5C5C", "#388E8E",            # Endo, Epi, Muscles
  "#9A32CD", "#79CDCD", "#A0522D", "#4169E1", "#6495ED", "#009ACD", "#97FFFF", # Mostly immunes
  "#FF8C00", "#E3A869", "#CD6600",                                            # Tu_D1
  "#EEEE00",                                                                  # Tu_D2
  "#CDC673", "#006400", "#FFD700", "#00CD00",                                 # Tu_D3
  "#A2CD5A", "#6B8E23",                                                       # Tu_D4
  "#00EE76",                                                                  # Tu_D5
  "#ADFF2F",                                                                  # Tu_D6
  "#5A5A5A"                                                                   # Mix
)

names(dlbcl_level4_cellcolors) <- dlbcl_level4_cellnames


# Put in script to decide palette after specify a disease ----------------
if(disease %in% c("breast", "lung")){
  level4_cellcolors <- br_lu_level4_cellcolors

  if(disease == "breast"){
    level1_5_cellcolors <- br_level1_5_cellcolors
  }else{
    level1_5_cellcolors <- lu_level1_5_cellcolors
  }

}else{
  level4_cellcolors <- dlbcl_level4_cellcolors
}


#####################################
#          Cell fraction            #
#####################################

# Breast and Lung
brlu_cellfraction_color_palette <- data.frame(cell_fraction =       c("Macro", "Malignant", "Other", "T cells"),
                                              cell_fraction_color = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1"))

cell_fraction_order = c("Macro", "Malignant", "Other", "T cells", "PanCK-") # note with PanCK-
cell_fraction_color = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1", "#09D0EF")
names(cell_fraction_color) <- cell_fraction_order

# DLBCL
dlbcl_cellfraction_color_palette <- data.frame(cell_fraction =       c("B cells", "Macro", "Other", "T cells"),
                                               cell_fraction_color = c("#FFD700", "#9A32CD", "#388E8E", "#4169E1"))



#####################################
#              Pathology            #
#####################################

# Breast pathology ------------------
breast_patho_anno <- c(
  "Breast_Carcinoma_in_situ",
  "Tumor_pure",
  "Most_likely_tumor",
  "Tumor_TIL",
  # "Tumor_Immune_mix",
  "Lymphocytes",
  "Immune_Cell_mix",
  "Tumor_Stroma_mix",
  "Large_Vessel",
  "Intratumoral_Stroma",
  "Intratumoral_Vessel",
  'Normal_Lung_Epithelium',
  "Acellular mucin",
  "Necrosis_Debris",
  "Artefact_Fold_exclude"
)

breast_patho_color <- c(
  "#00E5A1",
  "#03C78C",
  "#62e1bb",
  "#2da3f2",
  # (empty)
  "#4169E1",
  "#b2ff00",
  "#00c3c3",
  "#FF7256",
  "#388E8E",
  '#ba5440',
  "#BC8F8F",
  "#FFC1C1",
  "#666666",
  "#cccccc"
)

names(breast_patho_color) <- breast_patho_anno

# Lung pathology -----------------
lung_patho_anno <- c(
  "Tumor_pure",
  "Most_likely_Tumor",
  "Tumor_TIL",
  "Tumor_Immune_mix",
  "Lymphocytes",
  "Immune_Cell_mix",
  "Tumor_Stroma_mix",
  "Large_Vessel",
  "Intratumoral_Stroma",
  "Intratumoral_Vessel",
  "Normal_Lung_Epithelium",
  "Acellular mucin",
  "Necrosis_Debris",
  "Artefact_Fold_exclude"
)

lung_patho_color <- c(
  "#FFA54F",
  "#c47f3d",
  "#2da3f2",
  "#00e04e",
  "#4169E1",
  "#b2ff00",
  "#00c3c3",
  "#FF7256",
  "#388E8E",
  "#ba5440",
  '#BC8F8F',
  "#FFC1C1",
  "#666666",
  "#cccccc"
)

names(lung_patho_color) <- lung_patho_anno

# DLBCL pathology ----------------
dlbcl_patho_anno <- c(
  "Tumor",
  "Small lymphocytes",
  "Stroma",
  "Vessels",
  "Epithelium",
  "Hemorrhage",
  "Empty",
  "Necrosis"
)

dlbcl_patho_color <- c(
  "#FFD700",
  "#4169E1",
  "#388E8E",
  "#FF7256",
  "#BC8F8F",
  "#ea9999",
  "#cccccc",
  "#666666"
)

names(dlbcl_patho_color) <- dlbcl_patho_anno


# Level 1.5 color ---------------------------------------------------------
# Breast 
breast_1_5ct_names <- c(      "B", "Fibro_muscle", "Granulocyte", "Myeloid",    "T_NK",   "Tumor", "Vessel")
breast_1_5ct_color <- c("#EEEE00",      "#388E8E",     "#FFC125", "#9A32CD", "#4169E1", "#00C78C", "#FF7256")
names(breast_1_5ct_color) <- breast_1_5ct_names

# Lung 
lung_1_5ct_names <- c(      "B", "Epithelia", "Fibro_muscle", "Granulocyte", "Myeloid",    "T_NK",   "Tumor", "Vessel")
lung_1_5ct_color <- c("#EEEE00",   "#BC8F8F",      "#388E8E",     "#FFC125", "#9A32CD", "#4169E1", "#FFA54F", "#FF7256")
names(lung_1_5ct_color) <- lung_1_5ct_names

# DLBCL 
dlbcl_1_5ct_names <- c(      "B", "Epithelia", "Fibro_Muscle", "Myeloid", "T_NK",    "Tumor",   "Vessel")
dlbcl_1_5ct_color <- c("#EEEE00", "#BC8F8F",   "#388E8E",      "#9A32CD", "#4169E1", "#FFD700", "#FF7256")
names(dlbcl_1_5ct_color) <- dlbcl_1_5ct_names



# Patient palette ----------------------------------------------------------
breast_pt_palette <- c("#D9636CFF", "#F588AFFF", "#E6A0C4FF", "#F26386FF")


lung_pt_palette <- c("#FFA54F", "#FF8C00", "#CD853F", "#8B4513")


dlbcl_pt_palette <- c("#A4D984FF","#EEEE00","#FFD700","#A2CD5A","#00EE76","#ADFF2F")







