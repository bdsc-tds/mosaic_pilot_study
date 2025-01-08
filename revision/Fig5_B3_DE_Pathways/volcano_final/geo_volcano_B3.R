# de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05,
#                                   lfc = 1)
# results <- de_genes_toptable_BvT

library(dplyr)
library(ggplot2)

# contrast <- "Geo_3vs6_DE_lfc0"
contrast <- "Geo_3vs6_DE_lfc0_pval1"
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/DiffDE_nofilter"
results <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))

results = results %>% arrange(adj.P.Val)
head(results)

## Categorize results based on P-value & FDR for plotting
results$Color[results$logFC > 0 & results$adj.P.Val < 0.05] <- "Clus159_05"     # color by FDR 
results$Color[results$logFC > 0 & results$adj.P.Val < 0.01] <- "Clus159_01"
results$Color[results$logFC < 0 & results$adj.P.Val < 0.05] <- "Clus14_05"
results$Color[results$logFC < 0 & results$adj.P.Val < 0.01] <- "Clus14_01"
results$Color[abs(results$logFC) < 1 | results$adj.P.Val >= 0.05] <- paste0(expression("log"[2]*"FoldChange"), " < 1 or NS")
results$Color[results$P.Value < 0.05 & results$adj.P.Val >= 0.05 & abs(results$logFC) > 1] <- "Pval < 0.05 & Pvaladj >= 0.05"

results$Color <- factor(results$Color,
                        levels = c("Clus159_05", "Clus159_01",
                                   "Clus14_05", "Clus14_01", 
                                   "Pval < 0.05 & Pvaladj >= 0.05",
                                   paste0(expression("log"[2]*"FoldChange"), " < 1 or NS")))

top_n_genes = 10
results$invert_P <- (-log10(results$adj.P.Val)) * sign(results$logFC)           # label by FDR
top_g = c(results[, 'TargetName'][
  order(results[, "invert_P"], decreasing = TRUE)[1:top_n_genes]],
  results[, 'TargetName'][
    order(results[, 'invert_P'], decreasing = FALSE)[1:top_n_genes]])
top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix

blue_vis_geo <- c("SULT1C3", "PEG10", "PLA2G2A", "AZGP1", "MUCL1", "NR2F1", "HIPK2", "ACADM", "OLFML3")
pink_vis_geo <- c("NPPC", "AKR1C2", "GPRC5A", "GSTP1", "PGC", "LBP", "THRSP", "S100P", "ST6GALNAC4", "BST2", "CRIP1", 
                  "ACSL3", "CAPS", "FAM234B", "AP1S3", "SH3BGRL", "PKIB", "ATP2A3", "SGK1", "MBOAT2", "FNIP2")

results$DEvisgeo <- ifelse(results$TargetName %in% c(blue_vis_geo, pink_vis_geo), TRUE, FALSE)

results <- results %>%
  dplyr::rename(p_val_adj = adj.P.Val,
                p_val = P.Value, 
                avg_log2FC = logFC,
                gene = TargetName)

plot_volcano <- function(results, legendpos = "none"){
  p <- ggplot(results,
              aes(x = avg_log2FC, y = -log10(p_val),                                     # y axis by p value 
                  color = Color, label = gene)) +
    geom_vline(xintercept = c(1, -1), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point(data = results[abs(results$avg_log2FC) <= 1, ],
               aes(size = 8)) +
    geom_point(data = results[abs(results$avg_log2FC) > 1, ],
               aes(size = 8)) +
    geom_point(data = results[results$DEvisgeo, ],     
               aes(x = avg_log2FC, y = -log10(p_val)),                                   # point location by p value
               size = 3.5, color = "black", fill = NA, shape = 21, stroke = 2.5) +
    labs(x = expression("log"[2]*"FoldChange"),
         y = expression("Significance, -log"[10]*"p-value"),
         color = "Significance") +
    scale_color_manual(values = c("Clus159_01" = "#ff00ff",
                                  "Clus14_05" = "#9393ff",
                                  "Clus14_01" = "#0000ff",
                                  "Pval < 0.05 & Pvaladj >= 0.05" = "darkgrey",
                                  "\"log\"[2] * \"FoldChange\" < 1 or NS" = "lightgrey"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    ggrepel::geom_text_repel(data = subset(results, (gene %in% top_g & p_val_adj < 0.001 & abs(avg_log2FC) > 1) | gene == "ADRA2A"),  # label by FDR
                             size = 10, point.padding = 0.15, color = "black",
                             min.segment.length = .1, box.padding = .2, lwd = 2,
                             max.overlaps = 50) +
    theme_bw() +
    theme(legend.position = legendpos, 
          panel.grid = element_blank(),
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 20)) + 
    ylim(c(0, 25)) + 
    xlim(c(-2.1, 4.9))
  
  return(p)
}


# range(-log10(results$p_val))
# [1]  3.56972 24.42781

# range(results$avg_log2FC)
# -1.463069  4.754574

p <- plot_volcano(results, legendpos = "none")
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Volcano_final"
plot_title = "Geo_B3.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 20,
    height = 8.5)
print(p)
dev.off()

