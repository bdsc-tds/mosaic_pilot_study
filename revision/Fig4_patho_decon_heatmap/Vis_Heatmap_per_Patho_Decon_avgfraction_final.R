library(dplyr)
library(ggplot2)
library(patchwork)

# -------------------------------------------------------------------------
CT_order <- c("Epithelia", "Stroma", "Tumor", "Macrophage", "T cells", "B cells", "NK", "Myeloid else") # Fig 2, 3 

# -------------------------------------------------------------------------
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2_final_immune" # CARD
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/cell2location_vis_long" # C2L
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/RCTD_vis_long" # RCTD - level 1.5 immunne 
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD/level1_5_immune_long" # RCTD - level 4 combined to level 1.5 immunne 
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD_DEgenes/level1_5_immune_long" # RCTD - level 4 DEgenes combined to level 1.5 immunne 
vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L/level1_5_immune_long" # C2L - level 4 combined to level 1.5 immunne 
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L_DEgenes/level1_5_immune_long" # C2L - level 4 DEgenes combined to level 1.5 immunne 
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/AnchorIntegration/level1_5_immune_long" # AnchorIntegration - level 4 DEgenes combined to level 1.5 immunne 

decon_patho_gather_breast <- read.csv(file.path(vis_decon_path, "vis_breast_decon_immune_long.csv"))
decon_patho_gather_lung <- read.csv(file.path(vis_decon_path, "vis_lung_decon_immune_long.csv"))

# table(c(decon_patho_gather_breast$Region, decon_patho_gather_lung$Region))

tumor_related <- c("Breast_Carcinoma_in_situ",
                   "Tumor_pure",
                   "Most_likely_tumor",
                   "Most_likely_Tumor")

lymphocytes_related <- c("Lymphocytes")

macrophage_related <- c("Immune_Cell_mix")

stroma_related <- c("Large_Vessel",
                    "Intratumoral_Stroma",
                    "Intratumoral_Vessel",
                    "Normal_Lung_Epithelium")

breast_lung_patho <- c(tumor_related, lymphocytes_related, macrophage_related, stroma_related)
facet_name <- c("Tumor-enriched", "Stroma-enriched", "Lymphocytes-enriched", "Immune Cell Mix")

decon_patho_gather_breast_subset <- decon_patho_gather_breast %>%
  filter(Region %in% breast_lung_patho) %>%
  mutate(Region = case_when(Region %in% tumor_related ~ "Tumor-enriched",
                            Region %in% lymphocytes_related ~ "Lymphocytes-enriched",
                            Region %in% macrophage_related ~ "Immune Cell Mix",
                            Region %in% stroma_related ~ "Stroma-enriched")) %>%
  mutate(Region = factor(Region, levels = facet_name)) 

decon_patho_gather_lung_subset <- decon_patho_gather_lung %>%
  filter(Region %in% breast_lung_patho) %>%
  mutate(Region = case_when(Region %in% tumor_related ~ "Tumor-enriched",
                            Region %in% lymphocytes_related ~ "Lymphocytes-enriched",
                            Region %in% macrophage_related ~ "Immune Cell Mix",
                            Region %in% stroma_related ~ "Stroma-enriched")) %>%
  mutate(Region = factor(Region, levels = facet_name))


# DLBCL -------------------------------------------------------------------
decon_patho_gather_dlbcl <- read.csv(file.path(vis_decon_path, "vis_dlbcl_decon_immune_long.csv"))

d_tumor_related <- c("Tumor")

d_lymphocytes_related <- c("Small lymphocytes")

d_macrophage_related <- c("Epithelium")

d_stroma_related <- c("Stroma",
                      "Vessels")

dlbcl_patho <- c(d_tumor_related, d_lymphocytes_related, d_macrophage_related, d_stroma_related)
facet_name <- c("Tumor-enriched", "Stroma-enriched", "Lymphocytes-enriched", "Epithelium")

decon_patho_gather_dlbcl_subset <- decon_patho_gather_dlbcl %>%
  filter(Region %in% dlbcl_patho) %>%
  mutate(Region = case_when(Region %in% d_tumor_related ~ "Tumor-enriched",
                            Region %in% d_lymphocytes_related ~ "Lymphocytes-enriched",
                            Region %in% d_macrophage_related ~ "Epithelium",
                            Region %in% d_stroma_related ~ "Stroma-enriched")) %>%
  mutate(Region = factor(Region, levels = facet_name))



# -------------------------------------------------------------------------
plotHeatmap <- function(decon_patho_long_all){ 
  ## Per region per cell type
  stats <- decon_patho_long_all %>%
    group_by(Region, CellType) %>%
    summarize(mean_value = mean(Fraction), .groups = "drop") %>%
    mutate(per = round(100 * mean_value, 0))
  stats$CellType <- factor(stats$CellType, levels = rev(CT_order))
  stats$Region <- factor(stats$Region, levels = PA_order)
  
  # stats_color_rank <- stats %>%
  #   group_by(CellType) %>%
  #   arrange(CellType, desc(mean_value)) %>%
  #   mutate(index = 5 - row_number())
  # 
  # stats_color_rank$CellType <- factor(stats_color_rank$CellType, levels = rev(CT_order))
  # stats_color_rank$Region <- factor(stats_color_rank$Region, levels = PA_order)
  
  
  # Create the heatmap Region per decon
  p <- ggplot(stats, aes(x = Region, y = CellType, fill = per)) + 
  # p <- ggplot(stats_color_rank, aes(x = Region, y = CellType, fill = index)) + 
    geom_tile() + 
    geom_text(aes(x = Region, y = CellType, label = paste0(per, "%"), size = )) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colors = c("white", "dodgerblue"), transform = "sqrt") + 
    # scale_fill_gradientn(colors = c("white", "dodgerblue")) + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(angle = 45, size = 11, hjust = 0, vjust = 1)
    ) + 
    labs(fill = paste0("Rank of cell type abundance by pathology class"), title = "", y = "") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
  
  return(p)
}

# -------------------------------------------------------------------------
CT_order <- c("Tumor", "T cells", "NK",  "Macrophage", "Myeloid else", "B cells", "Stroma", "Epithelia") 
PA_order <- c("Tumor-enriched", "Lymphocytes-enriched", "Immune Cell Mix", "Stroma-enriched")
p_breast_patho_decon <- plotHeatmap(decon_patho_gather_breast_subset) 
p_lung_patho_decon <- plotHeatmap(decon_patho_gather_lung_subset)

PA_order <- c("Tumor-enriched", "Lymphocytes-enriched", "Epithelium", "Stroma-enriched") 
p_dlbcl_patho_decon <- plotHeatmap(decon_patho_gather_dlbcl_subset)

vis_heatmap <- p_breast_patho_decon | p_lung_patho_decon | p_dlbcl_patho_decon

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_Vis_Heatmap"
pdf(file = file.path(fig_path, "vis_heatmap_sqrt.pdf"),
# pdf(file = file.path(fig_path, "vis_heatmap_rank.pdf"),
    width = 10,
    height = 7)
print(vis_heatmap)
dev.off()

