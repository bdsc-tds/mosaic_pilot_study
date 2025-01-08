# install.packages("enrichR")
library(enrichR)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5/B3"

DE_result <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2/chrom_B3_tu_2_clus_markers.csv")
results <- as.data.frame(DE_result)

results <- results %>% filter(cluster == "Tu_B3_NPPC")
gene_list <- results$gene # [1:100]

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
  ) # + ylim(0,30)

plot_title = "Vis_B3_2_DE_Pathway_159.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 8, 
    height = 5.5)
print(p)
dev.off()


# -------------------------------------------------------------------------
DE_result <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2/chrom_B3_tu_2_clus_markers.csv")
results <- as.data.frame(DE_result)

results <- results %>% filter(cluster == "Tu_B3_PLA2G2A")
gene_list <- results$gene # [1:100]

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
  ) + ylim(0,30)


plot_title = "Vis_B3_2_DE_Pathway_14.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 7,
    height = 5.5)
print(p)
dev.off()

