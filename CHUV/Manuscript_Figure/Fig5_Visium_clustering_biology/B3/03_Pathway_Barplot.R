# install.packages("enrichR")
# devtools::install_github("wjawaid/enrichR")
library(enrichR)
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5/B3"

# contrast <- "ConsensusA_B"
# contrast <- "1&5&9_14"
# contrast <- "chrom_B3_tu_2_clus_markers"
contrast <- "Geo_3vs6_DE"
# save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2"
# save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2/archive"
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/DiffDE/"
DE_result <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))
results <- as.data.frame(DE_result)

# results <- results %>% filter(cluster == "1_5_7_9_Consensus")
# results <- results %>% filter(cluster == "1_5_9") %>%
# results <- results %>% filter(cluster == "Tu_B3_NPPC") %>%
#   arrange(p_val_adj)
# gene_list <- results$gene # [1:100]

results <- results %>% filter(logFC > 0) %>%
  arrange(adj.P.Val)
gene_list <- results$gene # 26

databases <- c("KEGG_2021_Human", "MSigDB_Hallmark_2020", "MSigDB_Computational", "MSigDB_Oncogenic_Signatures")
enrichment_results <- enrichr(gene_list, databases)
enrichment_results <- enrichment_results$MSigDB_Hallmark_2020 
head(enrichment_results)

enrichment_results_plt_df <- data.frame(enrichment_results) %>%
  arrange(Odds.Ratio) %>%
  tail(8)

head(enrichment_results_plt_df)

enrichment_results_plt_df$Term <- factor(enrichment_results_plt_df$Term,
                                         levels = enrichment_results_plt_df$Term)

pathway_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/enrichR"
# write.csv(enrichment_results_plt_df, file.path(pathway_path, "pink_vis_pathway.csv"))
# write.csv(enrichment_results_plt_df, file.path(pathway_path, "pink_chrom_pathway.csv"))
write.csv(enrichment_results_plt_df, file.path(pathway_path, "pink_geo_pathway.csv")) # no need as nothing

# Plot barplot of adjusted p-values
p <- ggplot(data = enrichment_results_plt_df, aes(x=Term, y=Odds.Ratio, fill=-log10(P.value))) +
  geom_bar(stat="identity") + 
  ggtitle("Enrichment Analysis Results") + 
  labs(fill = expression("-log"[10]*"p-value")) + 
  xlab("Enriched Terms") + 
  ylab("Odds Ratio") + coord_flip() + 
  scale_fill_gradient(low =  "#fc98fc",
                      high = "#ff00ff" 
  ) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        plot.title = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "bottom"
  ) # + ylim(0,15.5)

plot_title = "Vis_B3_2_DE_Pathway_159_consensus_final.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 8, 
    height = 5.5)
print(p)
dev.off()


# -------------------------------------------------------------------------
DE_result <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))
results <- as.data.frame(DE_result)

# results <- results %>% filter(cluster == "4_14_Consensus")
# results <- results %>% filter(cluster == "14") %>%
# results <- results %>% filter(cluster == "Tu_B3_PLA2G2A") %>%
#   arrange(p_val_adj)
# gene_list <- results$gene # [1:100]
results <- results %>% filter(logFC > 0) %>%
  arrange(adj.P.Val)
gene_list <- results$gene 

databases <- c("KEGG_2021_Human", "MSigDB_Hallmark_2020", "MSigDB_Computational", "MSigDB_Oncogenic_Signatures")
enrichment_results <- enrichr(gene_list, databases)
enrichment_results <- enrichment_results$MSigDB_Hallmark_2020 
head(enrichment_results)

enrichment_results_plt_df <- data.frame(enrichment_results) %>%
  arrange(Odds.Ratio) %>%
  tail(8)

head(enrichment_results_plt_df)

enrichment_results_plt_df$Term <- factor(enrichment_results_plt_df$Term,
                                         levels = enrichment_results_plt_df$Term)

pathway_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/enrichR"
# write.csv(enrichment_results_plt_df, file.path(pathway_path, "blue_vis_pathway.csv"))
# write.csv(enrichment_results_plt_df, file.path(pathway_path, "blue_chrom_pathway.csv"))
write.csv(enrichment_results_plt_df, file.path(pathway_path, "blue_geo_pathway.csv"))

# Plot barplot of adjusted p-values
p <- ggplot(data = enrichment_results_plt_df, aes(x=Term, y=Odds.Ratio, fill=-log10(P.value))) +
  geom_bar(stat="identity") + 
  ggtitle("Enrichment Analysis Results") + 
  labs(fill = expression("-log"[10]*"p-value")) + 
  xlab("Enriched Terms") + 
  ylab("Odds Ratio") + coord_flip() + 
  scale_fill_gradient(low = "#9393ff",  
                      high = "#0000ff" 
  ) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        plot.title = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "bottom"
  ) # + ylim(0,15.5)


plot_title = "Vis_B3_2_DE_Pathway_14_consensus_final.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 7,
    height = 5.5)
print(p)
dev.off()

