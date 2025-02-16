source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2"
deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/For_level1_5_immune_decon_results"

deconpath = deconresultpath; plot_title = "Geo_breast_lung_decon_sbs_box_outline_color_immune_ylim01.pdf"

# -------------------------------------------------------------------------
CT_order <- c("Epithelia", "Stroma", "Tumor", "Macrophage", "T cells", "B cells", "NK", "Myeloid else") # Fig 2, 3 

# -------------------------------------------------------------------------
gathered_df_breast <- read.csv(file.path(deconpath, "breast_batched_decon_long.csv"))
gathered_df_lung <- read.csv(file.path(deconpath, "lung_batched_decon_long.csv"))

# table(c(gathered_df_breast$cell_fraction, gathered_df_breast$cell_fraction))

facet_name <- c("Malignant", "Other", "T cells", "Macrophage")


# Combine breast and lung into one df -------------------------------------
gathered_df_breast$Indication <- "Breast"
gathered_df_lung$Indication <- "Lung"

gathered_df_breast_lung <- rbind(gathered_df_breast, gathered_df_lung) %>%
  mutate(cell_fraction = factor(cell_fraction, levels = facet_name))

gathered_df_breast_lung <- rbind(gathered_df_breast_lung,
                                 data.frame(sample = rep(NA, 4),
                                            CellType = rep("Epithelia", 4),
                                            Fraction = rep(0, 4),
                                            section_id = rep(NA, 4),
                                            cell_fraction = facet_name,
                                            Indication = rep("Breast", 4)))

gathered_df_breast_lung$CellType <- factor(gathered_df_breast_lung$CellType, levels = CT_order)

## Combined ------------------------------------------------------------
plot_each_patho_facet <- function(plt_df, facet_name_i = "Tumor_related"){
  plt_df <- plt_df %>% filter(cell_fraction == facet_name_i)
  
  p <- ggplot(plt_df, aes(x=CellType, y=Fraction, color=Indication)) +
    geom_boxplot(outlier.size = 1, outlier.alpha = 0.9) +
    theme_classic() +
    ylim(c(0, 1)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_line(size = 0.2),
          axis.text.y = element_text(size = 8),
          panel.grid = element_blank(),
          legend.position = "none") # + 
  # facet_wrap(~cell_fraction, ncol = 4) 
  p
}
# c("#00c78c2e", "#ffa54f31", "#ffd70033")
# c("#00c78cff", "#ffa54fff", "#ffd700ff")
# palette = c("#00c78cff", "#ffa54fff")
palette = c("#00a071ff", "#df9146ff")
p1 <- plot_each_patho_facet(plt_df = gathered_df_breast_lung, facet_name_i = facet_name[1]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2 <- plot_each_patho_facet(plt_df = gathered_df_breast_lung, facet_name_i = facet_name[2]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), 
        axis.text.x = element_blank())
p3 <- plot_each_patho_facet(plt_df = gathered_df_breast_lung, facet_name_i = facet_name[3]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank())
p4 <- plot_each_patho_facet(plt_df = gathered_df_breast_lung, facet_name_i = facet_name[4]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank())

p_breast_lung <- (p1 + theme(plot.margin = unit(c(0,20,0,10), "pt"))) | 
  (p2 + theme(plot.margin = unit(c(0,20,0,0), "pt"))) | 
  (p3 + theme(plot.margin = unit(c(0,20,0,0), "pt"))) | 
  (p4 + theme(plot.margin = unit(c(0, 10 ,0,0), "pt"))) 
p_breast_lung

write.csv(gathered_df_breast_lung, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig2b_geo_breast_lung.csv")

# DLBCL -------------------------------------------------------------------
gathered_df_dlbcl <- read.csv(file.path(deconpath, "dlbcl_batched_decon_long.csv"))

gathered_df_dlbcl <- gathered_df_dlbcl %>%
  mutate(cell_fraction = ifelse(cell_fraction == "B cells", "Malignant",
                                ifelse(cell_fraction == "Macro", "Macrophage", cell_fraction))) %>%
  mutate(cell_fraction = factor(cell_fraction, levels = facet_name))

gathered_df_dlbcl$CellType <- factor(gathered_df_dlbcl$CellType, levels = CT_order)
gathered_df_dlbcl$Indication <- "DLBCL"

palette = "black"
p5 <- plot_each_patho_facet(plt_df = gathered_df_dlbcl, facet_name_i = facet_name[1]) + 
  scale_color_manual(values = palette)
p6 <- plot_each_patho_facet(plt_df = gathered_df_dlbcl, facet_name_i = facet_name[2]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p7 <- plot_each_patho_facet(plt_df = gathered_df_dlbcl, facet_name_i = facet_name[3]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p8 <- plot_each_patho_facet(plt_df = gathered_df_dlbcl, facet_name_i = facet_name[4]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_dlbcl <- (p5 + theme(plot.margin = unit(c(10,20,0,10), "pt"))) | 
  (p6 + theme(plot.margin = unit(c(10,20,0,0), "pt"))) | 
  (p7 + theme(plot.margin = unit(c(10,20,0,0), "pt"))) | 
  (p8 + theme(plot.margin = unit(c(10, 10 ,0,0), "pt"))) 
p_dlbcl

write.csv(gathered_df_dlbcl, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig2b_geo_dlbcl.csv")

p_geo_decon <- p_breast_lung / 
  p_dlbcl 


pdf(file = file.path(figpath, plot_title),
    width = 7,
    height = 4.4)
print(p_geo_decon)
dev.off()
