library(dplyr)
library(ggplot2)

contrast <- "1&5&9_14" # "ConsensusA_B"
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2/archive"
DE_result <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))

results <- as.data.frame(DE_result)

# results$avg_log2FC <- ifelse(results$cluster == "4_14_Consensus", -1 * results$avg_log2FC, results$avg_log2FC) # upper tissue location one
results$avg_log2FC <- ifelse(results$cluster == "14", -1 * results$avg_log2FC, results$avg_log2FC) # upper tissue location one

results = results %>% arrange(p_val_adj)
head(results)

## Categorize results based on P-value & FDR for plotting
# results$Color[results$cluster == "1_5_7_9_Consensus" & results$p_val_adj < 0.05] <- "Clus159_05"
# results$Color[results$cluster == "1_5_7_9_Consensus" & results$p_val_adj < 0.01] <- "Clus159_01"
# results$Color[results$cluster == "4_14_Consensus" & results$p_val_adj < 0.05] <- "Clus14_05"
# results$Color[results$cluster == "4_14_Consensus" & results$p_val_adj < 0.01] <- "Clus14_01"

results$Color[results$cluster == "1_5_9" & results$p_val_adj < 0.05] <- "Clus159_05"
results$Color[results$cluster == "1_5_9" & results$p_val_adj < 0.01] <- "Clus159_01"
results$Color[results$cluster == "14" & results$p_val_adj < 0.05] <- "Clus14_05"
results$Color[results$cluster == "14" & results$p_val_adj < 0.01] <- "Clus14_01"

results$Color <- factor(results$Color,
                        levels = c("Clus159_05", "Clus159_01",
                                   "Clus14_05", "Clus14_01"))

top_n_genes = 10
# results$gene <- rownames(results)
# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$p_val)) * sign(results$avg_log2FC)
top_g = c(results[, 'gene'][
  order(results[, "invert_P"], decreasing = TRUE)[1:top_n_genes]],
  results[, 'gene'][
    order(results[, 'invert_P'], decreasing = FALSE)[1:top_n_genes]])
top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix

blue_vis_geo <- c("SULT1C3", "PEG10", "PLA2G2A", "AZGP1", "MUCL1", "NR2F1", "HIPK2", "ACADM", "OLFML3")
pink_vis_geo <- c("NPPC", "AKR1C2", "GPRC5A", "GSTP1", "PGC", "LBP", "THRSP", "S100P", "ST6GALNAC4", "BST2", "CRIP1", 
                  "ACSL3", "CAPS", "FAM234B", "AP1S3", "SH3BGRL", "PKIB", "ATP2A3", "SGK1", "MBOAT2", "FNIP2")

results$DEvisgeo <- ifelse(results$gene %in% c(blue_vis_geo, pink_vis_geo), TRUE, FALSE)

# Graph results
p <- ggplot(results,
            aes(x = avg_log2FC, y = -log10(p_val),
                color = Color, label = gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point(size = 7) +
  geom_point(data = results[results$DEvisgeo, ],     # Highlighted points
             aes(x = avg_log2FC, y = -log10(p_val)),
             size = 4.5, color = "black", fill = NA, shape = 21, stroke = 2.5) +
  labs(x = expression("log"[2]*"FoldChange"),
       y = expression("Significance, -log"[10]*"p-value"),
       color = "Significance") +
  scale_color_manual(values = c("Clus159_05" = "#fc98fc",
                                "Clus159_01" = "#ff00ff",
                                "Clus14_05" = "#9393ff",
                                "Clus14_01" = "#0000ff"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  ggrepel::geom_text_repel(data = subset(results, (gene %in% top_g & p_val_adj < 0.001)), # | (DEvisgeo & -log10(p_val) >= 21)),
                           size = 10, point.padding = 0.15, color = "black",
                           min.segment.length = .1, box.padding = .2, lwd = 2,
                           max.overlaps = 50) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20))

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Volcano_final"
# figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5/B3"
# figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Poster_Figures/"
  
plot_title = "Vis_B3.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 20,
    height = 8.5)
print(p)
dev.off()


