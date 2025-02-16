library(dplyr)
library(ggplot2)

breast_pt_palette <- c("#D9636CFF", "#F588AFFF", "#E6A0C4FF", "#F26386FF")
lung_pt_palette <- c("#FFA54F", "#FF8C00", "#CD853F", "#8B4513")
dlbcl_pt_palette <- c("#A4D984FF","#EEEE00","#FFD700","#A2CD5A","#00EE76","#ADFF2F")


# -------------------------------------------------------------------------
chrom_sample_order <- c("B1", "B2", "B3", "B3_rep", "B4", "B4_rep",
                        "L1", "L2", "L3", "L4",
                        "D1", "D2", "D3", "D4", "D5", "D6")
chrom_df <- data.frame(
  sample = chrom_sample_order,
  n = c(3937, 50, 3477, 1252, 991, 982,
        802, 1756, 17804, 15592,
        8936, 7097, 14375, 5335, 2485, 1485)
)
chrom_df$indication <- c(rep("breast", 6), rep("lung", 4), rep("dlbcl", 6))
chrom_df$indication <- factor(chrom_df$indication, levels = c("breast", "lung", "dlbcl"))
chrom_palette <- c("#D9636CFF", "#F588AFFF", "#E6A0C4FF", "#E6A0C4FF", "#F26386FF", "#F26386FF",
                   lung_pt_palette, 
                   dlbcl_pt_palette)
names(chrom_palette) <- chrom_df$sample
chrom_df$sample <- factor(chrom_df$sample, levels = chrom_sample_order)


# -------------------------------------------------------------------------
vis_sample_order <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2",
                      "L1_2", "L1_4", "L2_2", "L3_2", "L4_2",
                      "D1", "D2", "D3", "D4", "D5", "D6")
vis_df <- data.frame(
  sample = vis_sample_order,
  n = c(1730, 1988, 2560, 2156, 1652,
        870, 944, 822, 2378, 2750,
        4704, 4111, 4951, 1436, 1729, 1649)
  )
vis_df$indication <- c(rep("breast", 5), rep("lung", 5), rep("dlbcl", 6))
vis_df$indication <- factor(vis_df$indication, levels = c("breast", "lung", "dlbcl"))
vis_palette <- c("#D9636CFF", "#D9636CFF", "#F588AFFF", "#E6A0C4FF", "#F26386FF",
                 "#FFA54F", "#FFA54F", "#FF8C00", "#CD853F", "#8B4513", 
                 dlbcl_pt_palette)
names(vis_palette) <- vis_df$sample
vis_df$sample <- factor(vis_df$sample, levels = vis_sample_order)

# -------------------------------------------------------------------------
geo_sample_order <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1",
  "L1_1", "L2_1", "L3_1", "L3_3", "L4_3",
  "D1", "D2", "D3", "D4", "D5", "D6")
geo_df <- data.frame(
  sample = geo_sample_order,
  n = c(24, 20, 22, 22, 24,
        25, 24, 23, 20, 26,
        24, 23, 23, 21, 26, 19)
)
geo_df$indication <- c(rep("breast", 5), rep("lung", 5), rep("dlbcl", 6))
geo_df$indication <- factor(geo_df$indication, levels = c("breast", "lung", "dlbcl"))
geo_palette <- c("#D9636CFF", "#D9636CFF", "#F588AFFF", "#E6A0C4FF", "#F26386FF",
                 "#FFA54F", "#FF8C00", "#CD853F", "#CD853F", "#8B4513", 
                 dlbcl_pt_palette)
names(geo_palette) <- geo_df$sample
geo_df$sample <- factor(geo_df$sample, levels = geo_sample_order)


write.csv(chrom_df, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS2a.csv")     
write.csv(vis_df, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS2b.csv")  
write.csv(geo_df, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS2d.csv")  

# -------------------------------------------------------------------------
plot_bar_each_indication <- function(plt_df, lung_dodge, dlbcl_dodge, facet_gap){
  bar_width <- 0.8
  p <- ggplot(data = plt_df, aes(x = sample, y = n, fill = sample)) +
    geom_bar(stat = "identity", width = bar_width) +
    geom_text(data = plt_df %>% filter(indication == "breast"), 
              aes(x = as.numeric(sample), y = n, label = as.character(n)),
              vjust = -0.5,
              color = "black",
              size = 3.5) +
    geom_text(data = plt_df %>% filter(indication == "lung"), 
              aes(x = as.numeric(sample) - lung_dodge, y = n, label = as.character(n)),
              vjust = -0.5,
              color = "black",
              size = 3.5) +
    geom_text(data = plt_df %>% filter(indication == "dlbcl"), 
              aes(x = as.numeric(sample) - dlbcl_dodge, y = n, label = as.character(n)),
              vjust = -0.5,
              color = "black",
              size = 3.5) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5, hjust=0.5),
          panel.border = element_blank(), 
          legend.position =  "none",
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(facet_gap, "cm"),
          plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 14)) + 
    labs(y = "Counts", x = "Sample") + 
    facet_grid(.~indication, scales = "free_x", space="free_x")
  
  return(p)
}


p_chrom <- plot_bar_each_indication(chrom_df, lung_dodge = 6, dlbcl_dodge = 10, facet_gap = 4) + scale_fill_manual(values = chrom_palette)
p_vis <- plot_bar_each_indication(vis_df, lung_dodge = 5, dlbcl_dodge = 10, facet_gap = 4) + scale_fill_manual(values = vis_palette)
p_geo <- plot_bar_each_indication(geo_df, lung_dodge = 5, dlbcl_dodge = 10, facet_gap = 4) + scale_fill_manual(values = geo_palette)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Descriptive/sample_size"

pdf(file.path(figpath, "chrom_sample_size.pdf"), width = 13, height = 4.5)
print(p_chrom)
dev.off()

pdf(file.path(figpath, "vis_sample_size.pdf"), width = 13, height = 4.5)
print(p_vis)
dev.off()

pdf(file.path(figpath, "geo_sample_size.pdf"), width = 13, height = 4.1)
print(p_geo)
dev.off()


