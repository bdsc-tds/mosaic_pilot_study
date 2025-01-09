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
facet_name <- c("Tumor_related", "Stroma_related", "Lymphocytes_related", "Macrophage_related")

decon_patho_gather_breast_subset <- decon_patho_gather_breast %>%
  filter(Region %in% breast_lung_patho) %>%
  mutate(Region = case_when(Region %in% tumor_related ~ "Tumor_related",
                            Region %in% lymphocytes_related ~ "Lymphocytes_related",
                            Region %in% macrophage_related ~ "Macrophage_related",
                            Region %in% stroma_related ~ "Stroma_related")) %>%
  mutate(Region = factor(Region, levels = facet_name)) 

decon_patho_gather_lung_subset <- decon_patho_gather_lung %>%
  filter(Region %in% breast_lung_patho) %>%
  mutate(Region = case_when(Region %in% tumor_related ~ "Tumor_related",
                            Region %in% lymphocytes_related ~ "Lymphocytes_related",
                            Region %in% macrophage_related ~ "Macrophage_related",
                            Region %in% stroma_related ~ "Stroma_related")) %>%
  mutate(Region = factor(Region, levels = facet_name))


# Combine breast and lung into one df -------------------------------------
decon_patho_gather_breast_subset$Indication <- "Breast"
decon_patho_gather_lung_subset$Indication <- "Lung"

decon_patho_gather_breast_lung_subset <- rbind(decon_patho_gather_breast_subset, 
                                               decon_patho_gather_lung_subset) 

decon_patho_gather_breast_lung_subset <- rbind(decon_patho_gather_breast_lung_subset,
                                               data.frame(Barcode = rep(NA, 4),
                                                          CellType = rep("Epithelia", 4),
                                                          Fraction = rep(0, 4),
                                                          Section = rep(NA, 4),
                                                          Region = facet_name,
                                                          Indication = rep("Breast", 4)))

decon_patho_gather_breast_lung_subset$CellType <- factor(decon_patho_gather_breast_lung_subset$CellType, levels = CT_order)
## Combined ------------------------------------------------------------
plot_each_patho_facet <- function(plt_df, facet_name_i = "Tumor_related"){
  plt_df <- plt_df %>% filter(Region == facet_name_i)
  
  p <- ggplot(plt_df, aes(x=CellType, y=Fraction, color=Indication)) +
    geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.1) +
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
  # facet_wrap(~Region, ncol = 4) 
  p
}
# c("#00c78c2e", "#ffa54f31", "#ffd70033")
# c("#00c78cff", "#ffa54fff", "#ffd700ff")
# palette = c("#00c78cff", "#ffa54fff")
palette = c("#00a071ff", "#df9146ff")
p1 <- plot_each_patho_facet(plt_df = decon_patho_gather_breast_lung_subset, facet_name_i = facet_name[1]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2 <- plot_each_patho_facet(plt_df = decon_patho_gather_breast_lung_subset, facet_name_i = facet_name[2]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), 
        axis.text.x = element_blank())
p3 <- plot_each_patho_facet(plt_df = decon_patho_gather_breast_lung_subset, facet_name_i = facet_name[3]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank())
p4 <- plot_each_patho_facet(plt_df = decon_patho_gather_breast_lung_subset, facet_name_i = facet_name[4]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank())

p_breast_lung <- (p1 + theme(plot.margin = unit(c(0,20,0,10), "pt"))) | 
  (p2 + theme(plot.margin = unit(c(0,20,0,0), "pt"))) | 
  (p3 + theme(plot.margin = unit(c(0,20,0,0), "pt"))) | 
  (p4 + theme(plot.margin = unit(c(0, 10 ,0,0), "pt"))) 
p_breast_lung


# DLBCL -------------------------------------------------------------------
decon_patho_gather_dlbcl <- read.csv(file.path(vis_decon_path, "vis_dlbcl_decon_immune_long.csv"))

d_tumor_related <- c("Tumor")

d_lymphocytes_related <- c("Small lymphocytes")

# d_macrophage_related <- c("Hemorrhage")

# d_stroma_related <- c("Stroma",
#                       "Vessels",
#                       "Epithelium")

d_macrophage_related <- c("Epithelium")

d_stroma_related <- c("Stroma",
                      "Vessels")

dlbcl_patho <- c(d_tumor_related, d_lymphocytes_related, d_macrophage_related, d_stroma_related)
facet_name <- c("Tumor_related", "Stroma_related", "Lymphocytes_related", "Macrophage_related")

decon_patho_gather_dlbcl_subset <- decon_patho_gather_dlbcl %>%
  filter(Region %in% dlbcl_patho) %>%
  mutate(Region = case_when(Region %in% d_tumor_related ~ "Tumor_related",
                            Region %in% d_lymphocytes_related ~ "Lymphocytes_related",
                            Region %in% d_macrophage_related ~ "Macrophage_related",
                            Region %in% d_stroma_related ~ "Stroma_related")) %>%
  mutate(Region = factor(Region, levels = facet_name))



decon_patho_gather_dlbcl_subset$CellType <- factor(decon_patho_gather_dlbcl_subset$CellType, levels = CT_order)
decon_patho_gather_dlbcl_subset$Indication <- "DLBCL"

palette = "black"
p5 <- plot_each_patho_facet(plt_df = decon_patho_gather_dlbcl_subset, facet_name_i = facet_name[1]) + 
  scale_color_manual(values = palette)
p6 <- plot_each_patho_facet(plt_df = decon_patho_gather_dlbcl_subset, facet_name_i = facet_name[2]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p7 <- plot_each_patho_facet(plt_df = decon_patho_gather_dlbcl_subset, facet_name_i = facet_name[3]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
p8 <- plot_each_patho_facet(plt_df = decon_patho_gather_dlbcl_subset, facet_name_i = facet_name[4]) + 
  scale_color_manual(values = palette) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_dlbcl <- (p5 + theme(plot.margin = unit(c(10,20,0,10), "pt"))) | 
  (p6 + theme(plot.margin = unit(c(10,20,0,0), "pt"))) | 
  (p7 + theme(plot.margin = unit(c(10,20,0,0), "pt"))) | 
  (p8 + theme(plot.margin = unit(c(10, 10 ,0,0), "pt"))) 
p_dlbcl

p_vis_decon <- p_breast_lung / 
  p_dlbcl 


# plot_title = "Vis_breast_lung_decon_sbs_box_outline_color_immune.pdf"
plot_title = "Vis_breast_lung_decon_sbs_box_outline_color_immune_c2l_final.pdf"
pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2/", plot_title),
    width = 7,
    height = 4.4)
print(p_vis_decon)
dev.off()
