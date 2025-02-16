source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2"
deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/For_level1_5_immune_decon_results/"

deconpath = deconresultpath

# -------------------------------------------------------------------------
gathered_df_breast <- read.csv(file.path(deconpath, "breast_batched_decon_long.csv"))
gathered_df_lung <- read.csv(file.path(deconpath, "lung_batched_decon_long.csv"))
gathered_df_dlbcl <- read.csv(file.path(deconpath, "dlbcl_batched_decon_long.csv"))
gathered_df_dlbcl$cell_fraction <- ifelse(gathered_df_dlbcl$cell_fraction == "B cells", "Malignant", gathered_df_dlbcl$cell_fraction)

CF_order <- c("Malignant", "T cells", "Macrophage", "Other")
CT_order <- c("Tumor", "T cells", "NK", "Macrophage", "Myeloid else", "B cells", "Stroma", "Epithelia")

# -------------------------------------------------------------------------
plotHeatmap <- function(gathered_df){ 
  ## Per cell_fraction sum to 1
  stats <- gathered_df %>%
    group_by(cell_fraction, CellType) %>%
    summarize(mean_value = mean(Fraction), .groups = "drop") %>%
    mutate(per = round(100 * mean_value, 0))
  stats$CellType <- factor(stats$CellType, levels = rev(CT_order))
  stats$cell_fraction <- factor(stats$cell_fraction, levels = CF_order)
  
  # stats_color_rank <- stats %>%
  #   group_by(CellType) %>%
  #   arrange(CellType, desc(mean_value)) %>%
  #   mutate(index = 5 - row_number())
  # 
  # stats_color_rank$CellType <- factor(stats_color_rank$CellType, levels = rev(CT_order))
  # stats_color_rank$cell_fraction <- factor(stats_color_rank$cell_fraction, levels = CF_order)

  p <- ggplot(stats, aes(x = cell_fraction, y = CellType, fill = per)) + 
  # p <- ggplot(stats_color_rank, aes(x = cell_fraction, y = CellType, fill = index)) + 
    geom_tile() + 
    geom_text(aes(x = cell_fraction, y = CellType, label = paste0(per, "%"))) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left") +
    scale_fill_gradientn(colors = c("white", "#21751F"), transform = "sqrt") + 
    # scale_fill_gradientn(colors = c("white", "#21751F")) +  # for the scale
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
p_breast_patho_decon <- plotHeatmap(gathered_df_breast) 
p_lung_patho_decon <- plotHeatmap(gathered_df_lung)
p_dlbcl_patho_decon <- plotHeatmap(gathered_df_dlbcl)

colnames(gathered_df_dlbcl)[colnames(gathered_df_dlbcl) == "patient"] <- "section_id"

df_Fig4f <- rbind(gathered_df_breast, gathered_df_lung, gathered_df_dlbcl)
write.csv(df_Fig4f, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig4f_patho_decon_geo.csv")

geo_heatmap <- p_breast_patho_decon | p_lung_patho_decon | p_dlbcl_patho_decon


fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_Vis_Heatmap"
# pdf(file = file.path(fig_path, "geomx_heatmap_orig.pdf"),
pdf(file = file.path(fig_path, "geomx_heatmap_sqrt.pdf"),
# pdf(file = file.path(fig_path, "geomx_heatmap_rank.pdf"),
    width = 10,
    height = 7)
print(geo_heatmap)
dev.off()



