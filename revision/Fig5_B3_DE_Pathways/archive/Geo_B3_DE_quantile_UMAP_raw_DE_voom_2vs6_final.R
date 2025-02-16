# 2 samples in area A

library(edgeR)
library(limma)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(ggrepel)
library(ggspavis)
library(patchwork)
library(standR)
library(enrichR)

###########################################################################
# UMAP --------------------------------------------------------------------
geo_all <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds")
# subset to only B3_1 malig
geo <- geo_all[, geo_all$section_id == "B3_1"]
clus14_id <- c("DSP-1001660013520-C-B02", "DSP-1001660013520-C-B06") # , "DSP-1001660013520-C-B12"
clus159_id <- c("DSP-1001660013520-C-A03", "DSP-1001660013520-C-A05", "DSP-1001660013520-C-A07",
                "DSP-1001660013520-C-A09", "DSP-1001660013520-C-A12", "DSP-1001660013520-C-B04")
geo$malig_sub <- ifelse(geo$sample_id2 %in% clus14_id, "Malignant Area A",
                        ifelse(geo$sample_id2 %in% clus159_id, "Malignant Area B", "Healthy"))
geo$malig_sub <- factor(geo$malig_sub, levels = c("Malignant Area A", "Malignant Area B", "Healthy"))

## Quantile Norm
assays(geo, withDimnames=FALSE)$quantile <- preprocessCore::normalize.quantiles(assay(geo, "log1p"))
assays(geo)

assay = 8
which.assay = "quantile"

set.seed(100)
geo <- scater::runPCA(geo, assay.type = which.assay, ncomponents = 11)
set.seed(500)
geo <- scater::runUMAP(geo, dimred = "PCA")

p <- plotDimRed(geo, type = "UMAP", annotate = "malig_sub", text_by = "malig_sub", pt.size = 2)
p

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3"
plot_title = "Geo_B3_DE_UMAP.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 13,
    height = 4)
print(p)
dev.off()

###########################################################################
# Prep dataset for DE between malig subtypes only --------------------------
geo_all <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds")
geo_malig <- geo_all[, geo_all$section_id == "B3_1" &
                       geo_all$cell_fraction == "Malignant" & 
                       geo_all$sample_id2 != "DSP-1001660013520-C-B12"]
clus14_id <- c("DSP-1001660013520-C-B02", "DSP-1001660013520-C-B06") # , "DSP-1001660013520-C-B12"
clus159_id <- c("DSP-1001660013520-C-A03", "DSP-1001660013520-C-A05", "DSP-1001660013520-C-A07",
                "DSP-1001660013520-C-A09", "DSP-1001660013520-C-A12", "DSP-1001660013520-C-B04")
geo_malig$malig_sub <- ifelse(geo_malig$sample_id2 %in% clus14_id, "clus14",
                              ifelse(geo_malig$sample_id2 %in% clus159_id, "clus159", NA))

# # library(scater)
# # sce <- logNormCounts(geo_malig)
# # BvT
# # Down      45
# # NotSig 18545
# # Up        70
# 
# assays(geo_malig, withDimnames=FALSE)$quantile <- preprocessCore::normalize.quantiles(assay(geo_malig, "log1p"))
# # BvT
# # Down      45
# # NotSig 18545
# # Up        70


# Voom works on raw counts, so quantile is only for UMAP on full sample B3
# Here on directly work on malig AOIs only

# DE ----------------------------------------------------------------------
sce <- geo_malig
# counts(sce) <- assay(sce, "quantile")
dge <- SE2DGEList(sce)
design <- model.matrix(~0 + malig_sub, data = colData(sce))
colnames(design) <- gsub(" ","_",colnames(design))
colnames(design)
#  "malig_subclus14"  "malig_subclus159"
table(sce$malig_sub)
# clus14 clus159 
# 2       6 

contr.matrix <- makeContrasts(
  BvT = malig_subclus159 - malig_subclus14,
  levels = colnames(design))

keep <- filterByExpr(dge, design)
table(keep)
rownames(dge)[!keep]
dge_all <- dge[keep, ]

# Voom DE
v <- voom(dge_all, design, plot = TRUE) 
fit <- lmFit(v)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit<- decideTests(efit, p.value = 0.05, lfc = 1
                           # , adjust.method = "bonferroni"
                           )
summary_efit <- summary(results_efit)

summary_efit
## BvT.  # 3 vs 6 # BH, lfc = 0
## Down      45  
## NotSig 18545
## Up        70

# BvT.  # 2 vs 6 # BH, lfc = 0
# Down      79
# NotSig 18476
# Up       110

# BvT.      # 2 vs 6 # bonferroni, lfc = 1
# Down       8
# NotSig 18632
# Up        25

# BvT.      # 2 vs 6 # BH, lfc = 1
# Down      12
# NotSig 18626
# Up        27

# Volcano ----------------------------------------------------------------
de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05,
                                  lfc = 1,
                                  # adjust.method = "bonferroni"
                                  adjust.method = "BH"
                                  )
results <- de_genes_toptable_BvT

results = results %>% arrange(adj.P.Val)
head(results)

## Categorize results based on P-value & FDR for plotting
results$Color[results$logFC > 0 & results$adj.P.Val < 0.05] <- "Clus159_05"
results$Color[results$logFC > 0 & results$adj.P.Val < 0.01] <- "Clus159_01"
results$Color[results$logFC < 0 & results$adj.P.Val < 0.05] <- "Clus14_05"
results$Color[results$logFC < 0 & results$adj.P.Val < 0.01] <- "Clus14_01"

results$Color <- factor(results$Color, levels = c("Clus159_05", "Clus159_01", "Clus14_05", "Clus14_01"))

top_n_genes = 12
# results$TargetName <- rownames(results)
# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$P.Value)) * sign(results$logFC)
top_g = c(results[, 'TargetName'][
  order(results[, "invert_P"], decreasing = TRUE)[1:top_n_genes]],
  results[, 'TargetName'][
    order(results[, 'invert_P'], decreasing = FALSE)[1:top_n_genes]])
top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix


# Graph results
p <- ggplot(results,
            aes(x = logFC, y = -log10(P.Value),
                color = Color, label = TargetName)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point(size = 5) +
  labs(x = expression("log"[2]*"FoldChange"),
       y = expression("Significance, -log"[10]*"p-value"),
       color = "Significance") +
  scale_color_manual(values = c("Clus159_05" = "#fc98fc",
                                "Clus159_01" = "#ff00ff",
                                "Clus14_05" = "#9393ff",
                                "Clus14_01" = "#0000ff"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  ggrepel::geom_text_repel(data = subset(results, TargetName %in% top_g ), # & adj.P.Val < 0.01
                           size = 8, point.padding = 0.15, color = "black",
                           min.segment.length = .1, box.padding = .2, lwd = 2,
                           max.overlaps = 50) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

save_path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Volcano_Pathway"
# write.csv(results, file.path(save_path, "voom_B3malig_DE_2vs6.csv"))
# write.csv(results, file.path(save_path, "voom_B3malig_DE_2vs6_lfc1_BF.csv"))
write.csv(results, file.path(save_path, "voom_B3malig_DE_2vs6_lfc1_BH.csv"))

# pdf(file = file.path(save_path, "voom_B3malig_volcano_2vs6.pdf"),
# pdf(file = file.path(save_path, "voom_B3malig_volcano_2vs6_lfc1_BF.pdf"),
pdf(file = file.path(save_path, "voom_B3malig_volcano_2vs6_lfc1_BH.pdf"),
    width = 11.5,
    height = 12)
print(p)
dev.off()


###########################################################################
# Pathway -----------------------------------------------------------------
save_path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Volcano_Pathway"
DE_results <- read.csv(file.path(save_path, 
                                 # "voom_B3malig_DE_2vs6.csv"
                                 # "voom_B3malig_DE_2vs6_lfc1_BF.csv"
                                 "voom_B3malig_DE_2vs6_lfc1_BH.csv"
                                 ),
                       row.names = 1)

results <- DE_results %>% filter(Color %in% c("Clus159_05", "Clus159_01"))
gene_list <- results$TargetName # [1:100]

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
  ) + ylim(0, 34.2)

# plot_title = "Geo_B3_DE_Pathway_159_pink_2vs6.pdf"
# plot_title = "Geo_B3_DE_Pathway_159_pink_2vs6_lfc1_BF.pdf"
plot_title = "Geo_B3_DE_Pathway_159_pink_2vs6_lfc1_BH.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 8, 
    height = 5.5)
print(p)
dev.off()


# -------------------------------------------------------------------------
results <- DE_results %>% filter(Color %in% c("Clus14_01", "Clus14_05"))
gene_list <- results$TargetName # [1:100]

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
  ) + ylim(0, 34.2)


# plot_title = "Geo_B3_DE_Pathway_14_blue_2vs6.pdf"
# plot_title = "Geo_B3_DE_Pathway_14_blue_2vs6_lfc1_BF.pdf"
plot_title = "Geo_B3_DE_Pathway_14_blue_2vs6_lfc1_BH.pdf"
pdf(file = file.path(save_path, plot_title),
    width = 7,
    height = 5.5)
print(p)
dev.off()



