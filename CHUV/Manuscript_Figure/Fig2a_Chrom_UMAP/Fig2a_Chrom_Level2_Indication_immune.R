library(Seurat)
library(dplyr)

# Breast and lung -------------------------------------------------------
chrom_brlu_palette <- data.frame(
  rbind(
    c("Epithelia", "#BC8F8F"), 
    c("Stroma", "#388E8E"), 
    c("B cells", "#EEEE00"), 
    c("T cells", "#4169E1"),
    c("NK", "#79CDCD"),
    c("Macrophage", "#9A32CD"),
    c("Myeloid else", "#C71585"),
    c("Tu_B1", "#00C78C"), 
    c("Tu_B3", "#66CDAA"), 
    c("Tu_B4", "#3CB371"), 
    c("Tu_L1", "#EEB4B4"), 
    c("Tu_L2", "#FFA54F"), 
    c("Tu_L3", "#CD853F"), 
    c("Tu_L4", "#8B4513")
  )
)
color <- chrom_brlu_palette$X2
names(color) <- chrom_brlu_palette$X1

# DLBCL ----------------------------------------------------------------
chrom_dlbcl_palette <- data.frame(
  rbind(
    c("Epithelia", "#BC8F8F"), 
    c("Stroma", "#388E8E"), 
    c("B cells", "#EEEE00"), 
    c("T cells", "#4169E1"),
    c("NK", "#79CDCD"),
    c("Macrophage", "#9A32CD"),
    c("Myeloid else", "#C71585"),
    c("Tu_D1", "#FF8C00"),
    c("Tu_D2", "#EEEE00"),
    c("Tu_D3", "#FFD700"),
    c("Tu_D4", "#A2CD5A"),
    c("Tu_D5", "#00EE76"),
    c("Tu_D6", "#ADFF2F")
  )
)

color <- chrom_dlbcl_palette$X2
names(color) <- chrom_dlbcl_palette$X1

#########################################################################
# -------------------------------------------------------------------------
chrom_prev_brlu <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/final_owkin_annot.rds")
# breast / lung
chrom_prev_brlu$level1_5_immune <- case_when(chrom_prev_brlu$annot_l1 == "B cells" ~ "B cells",   
                                             chrom_prev_brlu$annot_l1 %in% c("Fibro_muscle", "Endothelia") ~ "Stroma",
                                             chrom_prev_brlu$annot_l1 == "Lung_epithelia" ~ "Epithelia",
                                             chrom_prev_brlu$annot_l1 == "tumour_L1" ~ "Tu_L1", 
                                             chrom_prev_brlu$annot_l1 == "tumour_L2" ~ "Tu_L2", 
                                             chrom_prev_brlu$annot_l1 == "tumour_L3" ~ "Tu_L3",
                                             chrom_prev_brlu$annot_l1 == "tumour_L4" ~ "Tu_L4", 
                                             chrom_prev_brlu$annot_l1 == "tumour_B1" ~ "Tu_B1", 
                                             chrom_prev_brlu$annot_l1 == "tumour_B3" ~ "Tu_B3", 
                                             chrom_prev_brlu$annot_l1 == "tumour_B4" ~ "Tu_B4", 
                                             chrom_prev_brlu$annot_l2 == "NK" ~ "NK",
                                             chrom_prev_brlu$annot_l2 %in% c("DC_1", "DC_2", "DC_activated", "DC_plasmacytoid", "Granulocytes", "Mast cells") ~ "Myeloid else",
                                             chrom_prev_brlu$annot_l2 %in% c("Macrophages", "Monocytes") ~ "Macrophage",
                                             chrom_prev_brlu$annot_l2 %in% c("T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_regs", "Proliferating_TNK") ~ "T cells") 

chrom_prev_brlu$level1_5_immune <- factor(chrom_prev_brlu$level1_5_immune, levels = chrom_brlu_palette$X1)
table(chrom_prev_brlu$level1_5_immune)
# Epithelia       Stroma      B cells      T cells           NK   Macrophage Myeloid else        Tu_B1        Tu_B3        Tu_B4        Tu_L1        Tu_L2        Tu_L3        Tu_L4 
#       648         7834         8357        11977          708         6617         2880          944         2900          850          306          637          578         1407 

Idents(chrom_prev_brlu) <- chrom_prev_brlu$level1_5_immune
p <- DimPlot(chrom_prev_brlu, cols = color)

red_df <- data.frame(chrom_prev_brlu@reductions$umap@cell.embeddings)
red_df$celltype <- chrom_prev_brlu$level1_5_immune
red_df$indication <- chrom_prev_brlu$tissue
write.csv(red_df, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig2a_breast_lung.csv")

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2"
pdf(file.path(figpath, "/Chrom_UMAP_BrLu_level1_5_immune.pdf"), width = 9, height = 6)
print(p)
dev.off()


# -------------------------------------------------------------------------
color <- c("#00c78cff", "#ffa54fff")
names(color) <- c("Breast", "Lung")
chrom_prev_brlu$indication <- 
  case_when(chrom_prev_brlu$sample_id %in% c("B3", "B2", "B1", "B4", "B3_rep", "B4_rep") ~ "Breast",
            chrom_prev_brlu$sample_id %in% c("L3", "L2", "L1", "L4") ~ "Lung")
Idents(chrom_prev_brlu) <- factor(chrom_prev_brlu$indication, levels = c("Breast", "Lung"))
p <- DimPlot(chrom_prev_brlu, cols = color)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2"
pdf(file.path(figpath, "/Chrom_UMAP_BrLu_Indication.pdf"), width = 9, height = 6)
print(p)
dev.off()


# -------------------------------------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom_dlbcl <- readRDS(file.path(chrompath, "chrom_dlbcl.rds")) # 16066 39713
chrom_dlbcl$level2_immune <- ifelse(as.character(chrom_dlbcl$level1_5_immune) == "Tumor", as.character(chrom_dlbcl$Level2), as.character(chrom_dlbcl$level1_5_immune))
table(chrom_dlbcl$level2_immune)
# B cells    Epithelia   Macrophage Myeloid else           NK       Stroma      T cells        Tu_D1        Tu_D2        Tu_D3        Tu_D4        Tu_D5        Tu_D6 
#      77         2945         3075         1398          349         3388         3721         6319         5624         8641         2005         1019         1152 

Idents(chrom_dlbcl) <- factor(chrom_dlbcl$level1_5_immune, levels = chrom_dlbcl_palette$X1)
p <- DimPlot(chrom_dlbcl, cols = color)

red_df <- data.frame(chrom_dlbcl@reductions$umap@cell.embeddings)
red_df$celltype <- chrom_dlbcl$level1_5_immune
write.csv(red_df, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig2a_dlbcl.csv")

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2"
pdf(file.path(figpath, "/Chrom_UMAP_DLBCL_level1_5_immune.pdf"), width = 9, height = 6)
print(p)
dev.off()





