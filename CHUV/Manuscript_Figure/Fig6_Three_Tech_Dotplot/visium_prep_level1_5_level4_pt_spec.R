library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggridges)
library(patchwork)

figpath_ridge <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Fig6_Ridge"

DLBCLnChromium_Marker_Gene_List <- c(
  "MS4A1", "TNFRSF13C", "CD79B", "CD37", "PSMB8", "CD19", "TYMS",
  "TUBB", "TOP2A", "POLD4", "CD47", "CD52", "BLK", "CD38", "MAP2K1",
  "CD40", "BCL2L1", "TNFRSF8", "SMO", "RARA", "TYK2", "TNFRSF10B"
)

vispath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/dlbcl/spotclean/Results/"
vis <- readRDS(file.path(vispath, "Dlbcl-merge-SCTpostSpotClean.rds"))


# -------------------------------------------------------------------------
# "Tu_D1" = c("Tu_D1_LMO2", "Tu_D1_RGS13", "Tu_D1_SMIM14");
# "Tu_D2" = "Tu_D2_mito";
# "Tu_D3" = c("Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD");
# "Tu_D4" = c("Tu_D4_BCL7A", "Tu_D4_PNN"); 
# "Tu_D5" = "Tu_D5_CCL22";
# "Tu_D6" = "Tu_D6_BCL2"

decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Visium_Decon/Level4/C2L"

max_tumerged_all <- NULL
for(i in 1:6){
  # C2L level 4, healthy tumor together -------------------------------------
  
  decon <- read.csv(paste0(decon_path, "/DLBCL_", i, ".csv"))
  
  test_merged <- decon %>%
    rowwise() %>%
    mutate(Epithelia = sum(c_across(starts_with("Epi_"))), 
           Stroma = sum(across(c(Endothelia, Fibro_Muscle, Pericyte))),
           Plasma = B_plasma, 
           Myeloid = sum(c_across(c(DC_1, DC_2, DC_pc, Mono_Macro))),
           `T/NK` = sum(c_across(c(T_CD4, T_CD4_reg, T_CD8, T_dividing, NK))),
           Tu_D = sum(c_across(starts_with("Tu_D")))) %>%
    mutate(X = paste0(X, "_", i)) %>%
    select(X, Epithelia, Stroma, Plasma, Myeloid, `T/NK`, Tu_D)
  
  decon_result_tumerged <- test_merged %>%
    pivot_longer(cols = c(Epithelia, Stroma, Plasma, Myeloid, `T/NK`, Tu_D),
                 names_to = "CellType",
                 values_to = "Fraction") %>%
    dplyr::rename(Barcode = X) %>%
    group_by(Barcode) %>%
    arrange(desc(Fraction), .by_group = TRUE) %>%
    dplyr::slice(1) %>%               # majority voted cell type per spot
    ungroup()
  decon_result_tumerged$CellType <- ifelse(decon_result_tumerged$CellType == "Tu_D", paste0("Tu_D", i),
                                           decon_result_tumerged$CellType)
  
  
  # ## Healthy ------------------------------------------------------------
  # healthy <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/DLBCL/DLBCL_", i, 
  #                            "/DLBCL_", i, "_spot_level1_5_decon.csv"))
  # 
  # healthy_sub <- healthy %>%
  #   mutate(Stroma = Fibro_Muscle + Vessel) %>%
  #   select(X, B, Epithelia, Myeloid, Stroma, T_NK) %>%
  #   mutate(X = paste0(X, "_", i)) 
  # 
  # ## Tumor ------------------------------------------------------------
  # tumor <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/DLBCL/DLBCL_", i, 
  #                          "/DLBCL_", i, "_spot_Level4_decon.csv"))
  # 
  # # merged
  # tumor_sub_ <- tumor %>%
  #   select(X, all_of(get(paste0("Tu_D", i)))) %>%
  #   column_to_rownames("X") %>%
  #   mutate(Tu_D = rowSums(.)) %>%
  #   select(Tu_D) %>%
  #   rownames_to_column("X") %>%
  #   mutate(X = paste0(X, "_", i))
  # 
  # print(all(rownames(healthy_sub) == rownames(tumor_sub_)))
  # 
  # 
  # decon_result_tumerged <- cbind(healthy_sub %>% dplyr::rename(Barcode = X), tumor_sub_ %>% select(-X)) %>%
  #   pivot_longer(cols = c(B, Epithelia, Myeloid, Stroma, T_NK, Tu_D),
  #                names_to = "CellType",
  #                values_to = "Fraction") %>%
  #   group_by(Barcode) %>%
  #   arrange(desc(Fraction), .by_group = TRUE) %>%
  #   slice(1) %>%
  #   ungroup()
  # decon_result_tumerged$CellType <- ifelse(decon_result_tumerged$CellType == "Tu_D", paste0("Tu_D", i), 
  #                                          decon_result_tumerged$CellType)
  # 
  max_tumerged_all <- rbind(max_tumerged_all, decon_result_tumerged)
}

max_tumerged_all <- max_tumerged_all[max_tumerged_all$Barcode %in% colnames(vis), ]

table(max_tumerged_all$CellType)
# Epithelia   Myeloid    Stroma      T/NK     Tu_D1     Tu_D2     Tu_D3     Tu_D4     Tu_D5     Tu_D6 
# 1851       775       314       897      4552      3427      4095      1093       484      1092 

# max_tumerged_all$CellType <- ifelse(max_tumerged_all$CellType == "B", "B cells",
#                                     ifelse(max_tumerged_all$CellType == "T_NK", "T/NK", max_tumerged_all$CellType))

palette_ridge <- c("#BC8F8F", "#388E8E", "#EEEE00", "#9A32CD", "#4169E1", "#FF8C00", "#EEEE00", "#FFD700", "#A2CD5A", "#00EE76", "#ADFF2F")
ridge_order <- c("Epithelia", "Stroma", "Plasma", "Myeloid", "T/NK", "Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6")
names(palette_ridge) <- ridge_order

max_tumerged_all$CellType <- factor(max_tumerged_all$CellType, levels = rev(ridge_order))

# A ridge plot with all spots, even max fraction of a spot is less than 50% ---- 
p <- ggplot(max_tumerged_all,
            aes(x = Fraction, y = CellType, fill = CellType)) +
  geom_density_ridges(scale = 4, alpha = 0.5) + 
  scale_y_discrete(expand = c(0, 0)) +     
  scale_x_continuous(expand = c(0, 0)) +   
  coord_cartesian(clip = "off") + 
  xlim(c(0, 1)) +
  theme_ridges() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = palette_ridge)

write.csv(p$data, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS15a.csv")

sample_size_df <- data.frame(table(max_tumerged_all$CellType)) %>%
  dplyr::rename(n = Freq,
                ylabel = Var1) %>%
  mutate(n = paste0("n = ", n))

p <- p + geom_text(data=sample_size_df,
                   aes(label = n, x= 1, y=ylabel, vjust = -1, hjust=-0.1),
                   position = position_stack(),
                   inherit.aes = FALSE) 

plot_title = "Vis_ridge_decon_max_level1_5_level4_pt_spec_C2L.pdf"
pdf(file = file.path(figpath_ridge, plot_title),
    width = 12.5,
    height = 8)
print(p)
dev.off()


# Now a ridge only if the max fraction is > 0.5 ---------------------------
max_tumerged_all_50 <- max_tumerged_all %>% filter(Fraction > 0.5)

p <- ggplot(max_tumerged_all_50,
            aes(x = Fraction, y = CellType, fill = CellType)) +
  geom_density_ridges(scale = 4, alpha = 0.5) + 
  scale_y_discrete(expand = c(0, 0)) +     
  scale_x_continuous(expand = c(0, 0)) +   
  coord_cartesian(clip = "off") + 
  xlim(c(0, 1)) +
  theme_ridges() + 
  theme(legend.position = "none",
    plot.margin = margin(1, 1, 0, 0, "cm")) + 
  scale_fill_manual(values = palette_ridge)

write.csv(p$data, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS15b.csv")

sample_size_df <- data.frame(table(max_tumerged_all_50$CellType)) %>%
  dplyr::rename(n = Freq,
                ylabel = Var1) %>%
  mutate(n = paste0("n = ", n))

p <- p + geom_text(data=sample_size_df,
                   aes(label = n, x= 1, y=ylabel, vjust = -1, hjust=-0.1),
                   position = position_stack(),
                   inherit.aes = FALSE) 

plot_title = "Vis_ridge_decon_max_level1_5_level4_pt_spec_50_C2L.pdf"
pdf(file = file.path(figpath_ridge, plot_title),
    width = 12.5,
    height = 8)
print(p)
dev.off()


# -------------------------------------------------------------------------
vis_small <- vis[, colnames(vis) %in% max_tumerged_all_50$Barcode]

CD <- data.frame(Barcode = colnames(vis_small))

decon_max_label <- data.frame(
  Barcode = max_tumerged_all_50$Barcode,
  decon_max = max_tumerged_all_50$CellType)

CD <- CD %>%
  left_join(decon_max_label)

vis_small$new_annot <- CD$decon_max

Idents(vis_small) <- factor(vis_small$new_annot, levels = ridge_order)

saveRDS(vis_small, file.path(vispath, "Dlbcl-merge-SCTpostSpotClean_small_fig6e_level1_5_level4_pt_spec_C2L.rds"))

p <- DotPlot(
  vis_small,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5))























