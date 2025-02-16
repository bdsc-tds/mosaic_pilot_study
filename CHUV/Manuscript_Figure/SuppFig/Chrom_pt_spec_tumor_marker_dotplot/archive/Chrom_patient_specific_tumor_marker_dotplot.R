library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
# library(SCpubr)

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrompath_br_lung <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/final_owkin_annot.rds"
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Chrom_Pt_Spec_Tumor_DotPlot"

## My order
# tu_breast <- c(
#   "NAMPT", "LTF", "STAC2",  "ANKRD30A",                    # Tu_B1
#   "FKBP5", "MUCL1", "FASN", "FABP7", "ALDH3B2", "ALOX15B", # Tu_B1 Tu_B3, more Tu_B1
#   "ALCAM",                                                 # Tu_B1 Tu_B3, both
#   "ACSM1", "APOD",                                         # Tu_B1 Tu_B3, more Tu_B3
#   "CYP4F8", "CLU",                                         # Tu_B3
#   "AZGP1", "MLPH",                                         # all patients, more B1 B3
#   "MGP",                                                   # Tu_B3 Tu_B4, more Tu_B3
#   "GATA3", "EVL", "IER2", "NTN4", "ESR1", "SHROOM1",       # Tu_B4
#   "KRT19", "RHOB"                                          # all patients, more B4
# )
# 
# tu_lung <- c(
#   "SFTPB", "AKR1C2", "CPS1", "PCSK2", "AKR1C3", "HOPX", "TCIM",                         # Tu_L1
#   "F3", "WFDC2", "MUC1", "MSLN", "FTH1", "ITGA3",                                       # Tu_L1 & other 
#   "B2M", "CD74", "DCBLD2", "LAMB3", "CAV2", "TGFBI", "MDK", "COL17A1", "G0S2", "FAM3C", # Tu_L2
#   "MYH14", "GPRC5A", "ANXA4", "CD24", "MYO1D", "DDX52", "PKHD1", "FXYD2",               # Tu_L3
#   "FXYD3", "DSP", "GSTP1", "AQP3", "FURIN", "SERPINB1",                                 # Tu_L4 & other 
#   "KRT17", "S100A2", "S100A9", "TRIM29"                                                 # Tu_L4
# )

tu_breast <- c(
  "NAMPT", "FKBP5", "LTF", "STAC2", "ALOX15B", "ANKRD30A", "ALCAM", "FABP7",
  "CYP4F8", "AZGP1", "CLU", "ALDH3B2", "FASN", "MGP", "ACSM1", "MUCL1", "APOD",
  "GATA3", "EVL", "IER2", "NTN4", "ESR1", "SHROOM1", "KRT19", "MLPH", "RHOB"
)

tu_lung <- c(
  "WFDC2", "F3", "TCIM", "HOPX", "AKR1C3", "PCSK2", "CPS1", "AKR1C2", "SFTPB", "B2M",
  "FTH1", "CD74", "DCBLD2", "LAMB3", "CAV2", "TGFBI", "MDK", "COL17A1", "G0S2", "FAM3C",
  "ITGA3", "GPRC5A", "PKHD1", "MYH14", "MSLN", "FXYD2", "ANXA4", "CD24", "MYO1D",
  "DDX52", "MUC1", "SERPINB1", "FURIN", "FXYD3", "AQP3", "DSP", "GSTP1", "TRIM29",
  "S100A9", "S100A2", "KRT17"
)

tu_br_ind_markers_ <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/dbuszta/breast_tum_de_genes.csv")
tu_br_ind_markers <- tu_br_ind_markers_ %>%
  group_by(tumour) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n = 50)

tu_lu_ind_markers_ <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/dbuszta/lung_tum_de_genes.csv")
tu_lu_ind_markers <- tu_lu_ind_markers_ %>%
  group_by(tumour) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n = 50)

# helper  -----------------------------------------------------------------
plot_save_pt_spec_tumor_dotplot <- function(seu, genelist #, 
                                            #savefig_width = 3.5, savefig_height = 8, 
                                            #save_title = "Breast_Chrom_DotPlot.pdf"
                                            ){
  
  seu@assays$RNA@counts <- seu@assays$RNA@counts[rownames(seu@assays$SCT@counts), ] # DLBCL bug
  seu@assays$RNA@data <- seu@assays$RNA@counts[rownames(seu@assays$SCT@data), ]     # DLBCL bug
  seu@assays$RNA@meta.features <- seu@assays$RNA@meta.features[rownames(seu@assays$SCT@data), ]  # DLBCL bug
  
  DefaultAssay(seu) <- "RNA" # SCT smooth out signal
  seu <- NormalizeData(seu, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000) 
  # seu <- FindVariableFeatures(seu)
  # seu <- ScaleData(seu)    # scale or not does not matter, but for dotplot, it takes log transformed data from the "data" slot of any assay
  
  p <- SCpubr::do_DotPlot(sample = seu, features = rev(genelist), 
                          assay = "RNA",
                          axis.text.x.angle = 0,
                          flip = TRUE, legend.position = "bottom"
  ) + 
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.title = element_text(size = 10),
          axis.text.x = element_text(size = 15), 
          legend.box = "vertical") +
    scale_fill_viridis_c(option = "rocket", direction = -1) + 
    scale_x_discrete(position = "top") + 
    guides(fill = guide_colorbar(title = "Avg. Expression", title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.5, "cm"))
    )
  return(p)
  
  # pdf(file.path(figpath, save_title), width = savefig_width, height = savefig_height)
  # print(p)
  # dev.off()
}


# Patient specific tumor marker -----------------------------------------
seu <- readRDS(file.path(chrompath, "chrom_breast.rds"))
# seu <- seu[rownames(seu) %in% tu_breast, ]
seu <- seu[, seu$Harmonised_Level4 %in% c("Tu_B1_MUCL1", "Tu_B1_MUCL1_necrosis", 
                                          "Tu_B1_MUCL1_transcription",
                                          "Tu_B3_CYP4F8", "Tu_B4_RHOB")]
# Idents(seu) <- factor(seu$patient, levels = paste0("B", 1:4))
Idents(seu) <- factor(seu$Level2)
  
p_breast <- plot_save_pt_spec_tumor_dotplot(seu, genelist = tu_breast# ,
                                #savefig_width = 3.5, savefig_height = 8,
                                # save_title = "breast_Chrom_DotPlot_RNA_lognorm_full.pdf"
                                )

df_breast <- as.data.frame(t(seu@assays$RNA@counts[tu_breast, ]))
df_breast$tu_group <- seu$Level2

write.csv(df_breast, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS4b_dot.csv")
# -------------------------------------------------------------------------
seu <- readRDS(file.path(chrompath, "chrom_lung.rds"))

# seu <- NormalizeData(seu, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000) 
# test <- rowSums(seu@assays$RNA$data)
# plot(density(test))

seu <- seu[, seu$Harmonised_Level4 %in% c("Tu_L1_SFTPB", "Tu_L2_FXYD2", "Tu_L3_G0S2", "Tu_L3_G0S2_immune_signature", 
                                          "Tu_L4_KRT17_immune_signature", "Tu_L4_KRT17_mucous",
                                          "Tu_L4_KRT17_necrosis", "Tu_L4_KRT17_neutrophil_signature")]

# seu <- seu[rownames(seu) %in% tu_lung, ]
# Idents(seu) <- factor(seu$patient, levels = paste0("L", 1:4))
Idents(seu) <- factor(seu$Level2)

# plot_save_pt_spec_tumor_dotplot(seu, genelist = tu_lung,
#                                 # savefig_width = 3.5, savefig_height = 10.5,
#                                 # save_title = "lung_Chrom_DotPlot_RNA_lognorm_full.pdf"
#                                 )

p_lung <- plot_save_pt_spec_tumor_dotplot(seu, genelist = tu_lung)

df_lung <- as.data.frame(t(seu@assays$RNA@counts[tu_lung, ]))
df_lung$tu_group <- seu$Level2

write.csv(df_lung, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS4c_dot.csv")

# Try annotating breast and lung together --------------------------------
genes_selected <- c(tu_breast, tu_lung)

seu <- readRDS(chrompath_br_lung)
# seu <- seu[rownames(seu) %in% genes_selected, ]
Idents(seu) <- factor(seu$patient, levels = c(paste0("B", 1:4), paste0("L", 1:4)))

plot_save_pt_spec_tumor_dotplot(seu, genelist = genes_selected,
                                savefig_width = 3.8, savefig_height = 17,
                                save_title = "breast_lung_Chrom_DotPlot_RNA_lognorm_full.pdf")


# -------------------------------------------------------------------------
tu_dlbcl <- c("NME2", "FCRL1", "LMO2", "IGKC", "MT-CO3", "TMSB4X", "GRHPR",
              "MPEG1", "CD1C", "CD52", "LTB", "SMIM14", "TNFRSF13C", "CD79A",
              "CD22", "FCRL5",                                    # D1 clus 2, 5, 13, 19:
              "MT-CYB", "NIBAN3", "SPIB", "MT-ATP6", "HIST1H1C", "IGHM",
              "MT-ND4", "TCL1A", "MT-ND4L",                     # D2 clus 0:
              "CD83", "SWAP70", "SEL1L3", "FCRL3", "ACTB", "IL4I1", "FAM3C",
              "ACTG2", "CD74", "MS4A1", "CNN2", "LRMP", "LCP1",  # D3 clus 1, 6, 26:
              "GGA2", "HELLS", "RASGRP2", "NKX6-3", "DTX1", "BCL7A", "MAT2A",
              "ARGLU1", "ALOX5", "AKNA", "ARHGEF1", "PNN", "TMC8",         # D4 clus 9&22:
              "TSPAN33", "SPATC1", "NFKB2", "BCL2L1", "CCL17",
              "CD40", "DUSP2", "ITPKB", "CCL22",                           # D5 clus 17:
              "POU2F2", "CCDC88A", "LENG8", "KLHL6", "TNFRSF13B", "BCL2",
              "NFATC1", "MT-ND6", "FCRL2"                        # D6 clus 12:
)


seu <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
seu <- seu[, grepl("Tu_D", seu$Level2)]
Idents(seu) <- factor(seu$Level2)

p <- plot_save_pt_spec_tumor_dotplot(seu, genelist = tu_dlbcl
                                #      ,
                                # savefig_width = 5, savefig_height = 17,
                                # save_title = "dlbcl_Chrom_DotPlot_RNA_lognorm_full.pdf"
                                )

df_dlbcl <- as.data.frame(t(seu@assays$RNA@counts[tu_dlbcl, ]))
df_dlbcl$tu_group <- seu$Level2

write.csv(df_dlbcl, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS4d_dot.csv")






