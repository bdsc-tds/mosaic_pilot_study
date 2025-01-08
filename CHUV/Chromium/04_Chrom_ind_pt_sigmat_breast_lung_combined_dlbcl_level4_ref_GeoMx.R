library(SpatialDecon)
library(dplyr)
library(Seurat)
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
sigmatpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level4/"

healthy_cells_brlu <- c(
  "B_cell", "B_plasma_IGHA1", "B_plasma_IGHG1", "B_plasma_IGHG3", "B_plasma_IGHM", "B_plasma_IGKC", 
  "B_plasma_IGLC1", "DC_1", "DC_2", "DC_activated", "DC_pc", "Endothelia_vascular", "Fibroblast", 
  "Fibroblast_B3", "Granulocyte", "Macrophage", "Mast_cell", "Monocyte", "Muscle_smooth", "NK", 
  "Pericyte", "T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing"
  )

healthy_cells_dlbcl <- c(
  "B_plasma", "DC_1", "DC_2", "DC_pc", "Endothelia", "Epi_bronchus", "Epi_mucous", 
  "Epi_Mucous_surface_gastric", "Epi_parietal_gastric", "Fibro_Muscle", "Mono_Macro",
  "NK", "Pericyte", "T_CD4", "T_CD4_reg", "T_CD8", "T_dividing"
  )

# helper --------------------------------------------------------------------
getSigmat <- function(chrom, disease, pt_name = "L1"){
  if(pt_name == "B2"){ # Only B2, take tumor as all other pts tumor combined.
    chrom$Harmonised_Level4 <- ifelse(grepl("Tu_B1|Tu_B3|Tu_B4", chrom$Harmonised_Level4), "Tu_B2", chrom$Harmonised_Level4)
  }
  
  if(disease %in% c("breast", "lung")){
    if(pt_name == "B3"){
      fibro_remove = "Fibroblast"
    }else{ # all other samples
      fibro_remove = "Fibroblast_B3"
    }
  }
  
  if(disease == "breast"){
    chrom_pt <- chrom[, chrom$Harmonised_Level4 %in% c(healthy_cells_brlu, unique(chrom$Harmonised_Level4)[grepl(paste0("Tu_", pt_name), unique(chrom$Harmonised_Level4))])]
    chrom_pt <- chrom_pt[, !(chrom_pt$Harmonised_Level4 %in% fibro_remove)]
  }else if(disease == "lung"){
    chrom_pt <- chrom[, chrom$Harmonised_Level4 %in% c(healthy_cells_brlu, "Epi_lung", unique(chrom$Harmonised_Level4)[grepl(paste0("Tu_", pt_name), unique(chrom$Harmonised_Level4))])]
    chrom_pt <- chrom_pt[, !(chrom_pt$Harmonised_Level4 %in% fibro_remove)]
  }else{
    chrom_pt <- chrom[, chrom$Harmonised_Level4 %in% c(healthy_cells_dlbcl, unique(chrom$Harmonised_Level4)[grepl(paste0("Tu_", pt_name), unique(chrom$Harmonised_Level4))])]
  }
  
  print(table(chrom_pt$Harmonised_Level4))
  
  annot_df <- data.frame(
    CellID = colnames(chrom_pt),
    LabeledCellType = chrom_pt$Harmonised_Level4
  )
  
  #TODO: add DE gene selection (For Sabrina)
  
  custom_mtx <- create_profile_matrix(mtx = chrom_pt@assays$RNA@counts,            # cell x gene count matrix
                                      cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
                                      cellTypeCol = "LabeledCellType",  # column containing cell type
                                      cellNameCol = "CellID",           # column containing cell ID/name
                                      matrixName = "custom_mini_colon", # name of final profile matrix
                                      outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                      normalize = FALSE,                # Should data be normalized? 
                                      minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                      minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                      scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                      discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.
  
  saveRDS(custom_mtx, file.path(sigmatpath, paste0("sigmat_", pt_name, ".rds")))
  p <- heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"), labRow = NA, margins = c(10, 5), cexCol = 0.7)
  print(p)
  return(custom_mtx)
}


# ## Breast only add lung immune stroma -------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds")) # table(chrom$patient, chrom$Harmonised_Level4)
#    B_cell B_plasma_IGHA1 B_plasma_IGHG1 B_plasma_IGHG3 B_plasma_IGHM B_plasma_IGKC B_plasma_IGLC1 DC_1 DC_2 DC_activated DC_pc Endothelia_vascular Fibroblast
# B1     51             82             51             54            98           356             22    3   32            8     6                  88        536
# B2      0              0              0              1             0             4              0    0    2            0     0                   3          2
# B3     18              2              2              0             0             4              1    1   20            3     1                 361          2
# B4     10              2              0              0             0             1              0    3   12            1     0                 130        192
# L1     24              7              0              1             1             5              0    4   26            3     3                  35         12
# L2     84             92             48             10            32           129             16    1   22            0     2                  31         44
# L3    354            276            274             52            28           489             94   90  488          164   215                 176       1871
# L4   1553            673           1188            323           193          1386            266  138  564          117   197                 412        518
# 
#    Fibroblast_B3 Granulocyte Macrophage Mast_cell Monocyte Muscle_smooth   NK Pericyte T_CD4 T_CD8_exhausted T_CTL T_CXCL13 T_reg TNK_dividing Tu_B1_MUCL1
# B1             0           2        103        98       95           282   95       52   161               6   411       23   270            8         254
# B2             0           0          5         0       17            13    0        3     0               0     0        0     0            0           0
# B3           780           1         17        30      242             9    3      123    78               1    90        0    41            2           0
# B4             0           0        112         0       75           382    3      105    13               2    47        2    28            3           0
# L1             0           0         20        11       34            39   15       21    61               0   100        0    26            0           0
# L2             0           2         65         5       17            36   58        8    48              38    78       35   131           15           0
# L3             0         135       3879       107      733           818  199       52   818            1151  2107      317  2098          195           0
# L4             0          47        647       316      559           578  335      120  1126              37   582      348  1276          204           0
# 
#    Tu_B1_MUCL1_necrosis Tu_B1_MUCL1_transcription Tu_B3_CYP4F8 Tu_B4_RHOB
# B1                  298                       392            0          0
# B2                    0                         0            0          0
# B3                    0                         0         2897          0
# B4                    0                         0            0        850
# L1                    0                         0            0          0
# L2                    0                         0            0          0
# L3                    0                         0            0          0
# L4                    0                         0            0          0

sigmatB1 <- getSigmat(chrom, disease = "breast", pt_name = "B1")
sigmatB2 <- getSigmat(chrom, disease = "breast", pt_name = "B2")
sigmatB3 <- getSigmat(chrom, disease = "breast", pt_name = "B3")
sigmatB4 <- getSigmat(chrom, disease = "breast", pt_name = "B4")


## Lung only add breast immune stroma -------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds")) # table(chrom$patient, chrom$Harmonised_Level4)
#    B_cell B_plasma_IGHA1 B_plasma_IGHG1 B_plasma_IGHG3 B_plasma_IGHM B_plasma_IGKC B_plasma_IGLC1 DC_1 DC_2 DC_activated DC_pc Endothelia_vascular Epi_lung Fibroblast
# B1     51             82             51             54            98           356             22    3   32            8     6                  88        0        536
# B2      0              0              0              1             0             4              0    0    2            0     0                   3        0          2
# B3     18              2              2              0             0             4              1    1   20            3     1                 361        0          2
# B4     10              2              0              0             0             1              0    3   12            1     0                 130        0        192
# L1     24              7              0              1             1             5              0    4   26            3     3                  35       48         12
# L2     84             92             48             10            32           129             16    1   22            0     2                  31       72         44
# L3    354            276            274             52            28           489             94   90  488          164   215                 176       46       1871
# L4   1553            673           1188            323           193          1386            266  138  564          117   197                 412      482        518
# 
#    Fibroblast_B3 Granulocyte Macrophage Mast_cell Monocyte Muscle_smooth   NK Pericyte T_CD4 T_CD8_exhausted T_CTL T_CXCL13 T_reg TNK_dividing Tu_L1_SFTPB Tu_L2_FXYD2
# B1             0           2        103        98       95           282   95       52   161               6   411       23   270            8           0           0
# B2             0           0          5         0       17            13    0        3     0               0     0        0     0            0           0           0
# B3           780           1         17        30      242             9    3      123    78               1    90        0    41            2           0           0
# B4             0           0        112         0       75           382    3      105    13               2    47        2    28            3           0           0
# L1             0           0         20        11       34            39   15       21    61               0   100        0    26            0         306           0
# L2             0           2         65         5       17            36   58        8    48              38    78       35   131           15           0         637
# L3             0         135       3879       107      733           818  199       52   818            1151  2107      317  2098          195           0           0
# L4             0          47        647       316      559           578  335      120  1126              37   582      348  1276          204           0           0
# 
#    Tu_L3_G0S2 Tu_L3_G0S2_immune_signature Tu_L4_KRT17_immune_signature Tu_L4_KRT17_mucous Tu_L4_KRT17_necrosis Tu_L4_KRT17_neutrophil_signature
# B1          0                           0                            0                  0                    0                                0
# B2          0                           0                            0                  0                    0                                0
# B3          0                           0                            0                  0                    0                                0
# B4          0                           0                            0                  0                    0                                0
# L1          0                           0                            0                  0                    0                                0
# L2          0                           0                            0                  0                    0                                0
# L3        311                         267                            0                  0                    0                                0
# L4          0                           0                          404                 81                  600                              322
sigmatL1 <- getSigmat(chrom, disease = "lung", pt_name = "L1")
sigmatL2 <- getSigmat(chrom, disease = "lung", pt_name = "L2")
sigmatL3 <- getSigmat(chrom, disease = "lung", pt_name = "L3")
sigmatL4 <- getSigmat(chrom, disease = "lung", pt_name = "L4")


# DLBCL  ------------------------------------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds")); chrom$patient <- chrom$sample_id # table(chrom$patient, chrom$Harmonised_Level4)
#         B_plasma DC_1 DC_2 DC_pc Endothelia Epi_bronchus Epi_mucous Epi_Mucous_surface_gastric Epi_parietal_gastric Fibro_Muscle Mono_Macro   NK Pericyte T_CD4 T_CD4_reg T_CD8 T_dividing Tu_D1_LMO2
# DLBCL_1       13   58   81    29         67            1          3                          1                    0         1112        694   13       23   207         0   312          3       3659
# DLBCL_2        0   17   63    21         12            1          1                          3                    0          191        647   10       23    12         0   472          0          0
# DLBCL_3       27  200  237   611         96            2          9                          1                    0         1050        828  319       35   855         0  1223        241          0
# DLBCL_4       32    2   33     0         76           19        941                       1232                  349          295        269    7        3     6         9    57          0          0
# DLBCL_5        2   12   26     0         94           97         15                         15                    0          287        600    0        4     3         0   311          0          0
# DLBCL_6        3    2    5     1          0            1        102                        123                   29           19         37    0        1     1         0     9          0          0
# 
#         Tu_D1_RGS13 Tu_D1_SMIM14 Tu_D2_mito Tu_D3_BCL2A1 Tu_D3_dividing Tu_D3_FAM3C Tu_D3_IGHD Tu_D4_BCL7A Tu_D4_PNN Tu_D5_CCL22 Tu_D6_BCL2
# DLBCL_1         767         1893          0            0              0           0          0           0         0           0          0
# DLBCL_2           0            0       5624            0              0           0          0           0         0           0          0
# DLBCL_3           0            0          0          325           1814        5795        707           0         0           0          0
# DLBCL_4           0            0          0            0              0           0          0        1470       535           0          0
# DLBCL_5           0            0          0            0              0           0          0           0         0        1019          0
# DLBCL_6           0            0          0            0              0           0          0           0         0           0       1152
sigmatD1 <- getSigmat(chrom, disease = "dlbcl", pt_name = "D1")
sigmatD2 <- getSigmat(chrom, disease = "dlbcl", pt_name = "D2")
sigmatD3 <- getSigmat(chrom, disease = "dlbcl", pt_name = "D3")
sigmatD4 <- getSigmat(chrom, disease = "dlbcl", pt_name = "D4")
sigmatD5 <- getSigmat(chrom, disease = "dlbcl", pt_name = "D5")
sigmatD6 <- getSigmat(chrom, disease = "dlbcl", pt_name = "D6")


# # Sanity checks -----------------------------------------------------------
# list.files(sigmatpath)
# sigmat_dlbcl_combined <- readRDS(file.path(sigmatpath, "sigmat_dlbcl.rds"))
# heatmap(sweep(sigmat_dlbcl_combined, 1, apply(sigmat_dlbcl_combined, 1, max), "/"), labRow = NA, margins = c(10, 5), cexCol = 0.7) # Tumor not so clear as breast/lung, whether combine or pt-specific



