DE_result <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2/archive/chrom_B3_tu_2_clus_markers.csv")

results <- as.data.frame(DE_result)

results$avg_log2FC <- ifelse(results$cluster == "Tu_B3_PLA2G2A", -1 * results$avg_log2FC, results$avg_log2FC) # upper tissue location one
results = results %>% arrange(p_val_adj)
head(results)

# Categorize results based on P-value & FDR for plotting
results$Color[results$cluster == "Tu_B3_NPPC" & results$p_val_adj < 0.05] <- "ClusNPPC_05"
results$Color[results$cluster == "Tu_B3_NPPC" & results$p_val_adj < 0.01] <- "ClusNPPC_01"
results$Color[results$cluster == "Tu_B3_PLA2G2A" & results$p_val_adj < 0.05] <- "ClusPLA2G2A_05"
results$Color[results$cluster == "Tu_B3_PLA2G2A" & results$p_val_adj < 0.01] <- "ClusPLA2G2A_01"

results$Color <- factor(results$Color,
                        levels = c("ClusNPPC_05", "ClusNPPC_01",
                                   "ClusPLA2G2A_05", "ClusPLA2G2A_01"))

top_n_genes = 12
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

# Graph results
p <- ggplot(results,
            aes(x = avg_log2FC, y = -log10(p_val),
                color = Color, label = gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point(size = 5) +
  labs(x = expression("log"[2]*"FoldChange"),
       y = expression("Significance, -log"[10]*"p-value"),
       color = "Significance") +
  scale_color_manual(values = c("ClusNPPC_05" = "#fc98fc",
                                "ClusNPPC_01" = "#ff00ff",
                                "ClusPLA2G2A_05" = "#9393ff",
                                "ClusPLA2G2A_01" = "#0000ff"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  ggrepel::geom_text_repel(data = subset(results, gene %in% top_g & p_val_adj < 0.001),
                           size = 8, point.padding = 0.15, color = "black",
                           min.segment.length = .1, box.padding = .2, lwd = 2,
                           max.overlaps = 50) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

write.csv(p$data, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS9f.csv")

# plot_title = "Chrom_B3_DE_Volcano_NPPC_PLA2G2A.pdf"
# pdf(file = file.path("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5/B3", plot_title),
#     width = 10,
#     height = 10)
# print(p)
# dev.off()


