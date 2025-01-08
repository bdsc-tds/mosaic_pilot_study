## Method 2, GeoMx's quantile no RUV4 -----------------------------------
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(ggrepel)
library(ggspavis)
library(patchwork)
library(standR)

# Use quantile norm on one sample
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds")
# Use logCPM + edgeR limma voom combo
# geo_norm <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe.rds")

clus14_id <- c("DSP-1001660013520-C-B02", "DSP-1001660013520-C-B06", "DSP-1001660013520-C-B12")
clus159_id <- c("DSP-1001660013520-C-A03", "DSP-1001660013520-C-A05", "DSP-1001660013520-C-A07",
                "DSP-1001660013520-C-A09", "DSP-1001660013520-C-A12", "DSP-1001660013520-C-B04")

geo <- geo[, geo$section_id == "B3_1"] # quantile
# geo <- geo_norm[, geo_norm$section_id == "B3_1"] # logCPM
# geo <- geo[, geo$section_id == "B3_1" & geo$cell_fraction == "Malignant"] # Malig only
geo$malig_sub <- ifelse(geo$sample_id2 %in% clus14_id, "Malignant Area A",
                         ifelse(geo$sample_id2 %in% clus159_id, "Malignant Area B", "Healthy"))
geo$malig_sub <- factor(geo$malig_sub, levels = c("Malignant Area A", "Malignant Area B", "Healthy"))


# Sample spec quantile norm, no need batch correct, PCA, UMAP, TSNE -------
## Quantile Norm
assays(geo, withDimnames=FALSE)$quantile <- preprocessCore::normalize.quantiles(assay(geo, "log1p"))
assays(geo)

assay = 8
which.assay = "quantile"


# Sample spec logCPM norm, no need batch correct, PCA, UMAP, TSNE -------
geo <- geomxNorm(geo, method = "CPM")
assay = 2
which.assay = "logcounts"

set.seed(100)
geo <- scater::runPCA(geo, assay.type = which.assay, ncomponents = 11)
# geo <- scater::runPCA(geo, assay.type = which.assay, ncomponents = 4)
set.seed(500)
geo <- scater::runUMAP(geo, dimred = "PCA")
# geo <- scater::runUMAP(geo, dimred = "PCA", n_neighbors = 4)
# set.seed(100)
# geo <- scater::runTSNE(geo, dimred = "PCA")
# 
# geo$malig_sub <- as.factor(geo$malig_sub)
# plotDimRed(geo, type = "UMAP", annotate = "malig_sub", text_by = "malig_sub") |
#   plotDimRed(geo, type = "TSNE", annotate = "malig_sub", text_by = "malig_sub")

p <- plotDimRed(geo, type = "UMAP", annotate = "malig_sub", text_by = "malig_sub", pt.size = 2) 
# p <- plotDimRed(geo, type = "UMAP", annotate = "cell_fraction", text_by = "cell_fraction", pt.size = 2) 
# p <- plotDimRed(geo, type = "UMAP", annotate = "roi", text_by = "roi", pt.size = 2) 
# + scale_color_manual(values = c("Malignant Area A" = "#0000FF", "Malignant Area B" = "#FF00FF", "Healthy" = "grey")) 
p

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3"
plot_title = "Geo_B3_DE_UMAP.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 13,
    height = 4)
print(p)
dev.off()


# # Geo_malig ---------------------------------------------------------------
# geo_malig <- geo[, geo$section_id == "B3_1" & geo$cell_fraction == "Malignant"]
# 
# # Sample spec quantile norm, no need batch correct, PCA, UMAP, TSNE -------
# ## Quantile Norm
# assays(geo_malig, withDimnames=FALSE)$quantile <- preprocessCore::normalize.quantiles(assay(geo_malig, "log1p")) 
# assays(geo_malig)
# 
# assay = 8
# which.assay = "quantile"
# 
# set.seed(100)
# geo_malig <- scater::runPCA(geo_malig, assay.type = which.assay, ncomponents = 4)
# set.seed(600)
# geo_malig <- scater::runUMAP(geo_malig, dimred = "PCA", n_neighbors = 4)
# 
# plotDimRed(geo_malig, type = "UMAP", annotate = "malig_sub", text_by = "malig_sub", pt.size = 2) 


# edgeR DE ---------------------------------------------------------------
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds")
geo <- geo[, geo$section_id == "B3_1"]

assays(geo, withDimnames=FALSE)$quantile <- preprocessCore::normalize.quantiles(assay(geo, "log1p"))
assays(geo)

clus14_id <- c("DSP-1001660013520-C-B02", "DSP-1001660013520-C-B06", "DSP-1001660013520-C-B12")
clus159_id <- c("DSP-1001660013520-C-A03", "DSP-1001660013520-C-A05", "DSP-1001660013520-C-A07",
                "DSP-1001660013520-C-A09", "DSP-1001660013520-C-A12", "DSP-1001660013520-C-B04")

geo$malig_sub <- ifelse(geo$sample_id2 %in% clus14_id, "clus14",
                        ifelse(geo$sample_id2 %in% clus159_id, "clus159", geo$cell_fraction))
geo$malig_sub <- factor(geo$malig_sub, levels = c("clus14", "clus159", "Other", "T cells"))
sce <- geo
sce <- sce[, sce$cell_fraction != "PanCK-"]
library(edgeR)
library(limma)

counts(sce) <- assay(sce, "meannormcounts")
dge <- SE2DGEList(sce)
design <- model.matrix(~0 + malig_sub, data = colData(sce))
colnames(design) <- gsub(" ","_",colnames(design))
colnames(design)
# "malig_subclus14"  "malig_subclus159" "malig_subOther"   "malig_subT_cells"
# "malig_subMalignant_Area_A" "malig_subMalignant_Area_B" "malig_subHealthy"
table(sce$malig_sub)
# clus14 clus159   Other T cells 
# 3       6      10       3 
# Malignant Area A Malignant Area B          Healthy 
#                3                6               13 

contr.matrix <- makeContrasts(
  BvT = malig_subclus159 - malig_subclus14,
  levels = colnames(design))

keep <- filterByExpr(dge, design)
table(keep)
rownames(dge)[!keep]
dge_all <- dge[keep, ]
# dge_all <- dge

# Voom DE
v <- voom(dge_all, design, plot = TRUE) 
fit <- lmFit(v)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit<- decideTests(efit, p.value = 0.05, lfc = 1)
summary_efit <- summary(results_efit)

summary_efit
#          BvT
# Down       4
# NotSig 18664
# Up         9

library(ggrepel)
library(tidyverse)

# de_results_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf)

de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05, lfc = 1)

de_results_BvT %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                     ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
  ggplot(aes(AveExpr, logFC, col = DE)) + 
  geom_point(shape = 1, size = 1) + 
  geom_text_repel(data = de_genes_toptable_BvT %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                                       ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("Clus 159 vs. Clus14") +
  scale_color_manual(values = c("blue","gray","deeppink1")) +
  theme(text = element_text(size=15))


# Previous visualization --------------------------------------------------
results <- de_genes_toptable_BvT

results = results %>% arrange(adj.P.Val)
head(results)

## Categorize results based on P-value & FDR for plotting
results$Color[results$logFC > 0 & results$adj.P.Val < 0.05] <- "Clus159_05"
results$Color[results$logFC > 0 & results$adj.P.Val < 0.01] <- "Clus159_01"
results$Color[results$logFC < 0 & results$adj.P.Val < 0.05] <- "Clus14_05"
results$Color[results$logFC < 0 & results$adj.P.Val < 0.01] <- "Clus14_01"

results$Color <- factor(results$Color, levels = c("Clus159_05", "Clus159_01", "Clus14_05", "Clus14_01"))

top_n_genes = 40
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
write.csv(results, file.path(save_path, "quantile_DE.csv"))

pdf(file = file.path(save_path, "quantile_volcano.pdf"),
    width = 11.5,
    height = 12)
print(p)
dev.off()

# -------------------------------------------------------------------------
DE_results <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Volcano_Pathway/voom_B3malig_DE_3vs6_lfc1_BH.csv",
                        row.names = 1)

## EnrichR ----------------------------------------------------------------
# results <- DE_results %>% filter(Color %in% c("Clus159_05", "Clus159_01")) %>%
#   arrange(adj.P.Val) # %>%
#   # top_n(9)
# gene_list <- results$gene # [1:100]
# 
# databases <- c("KEGG_2021_Human", "MSigDB_Hallmark_2020", "MSigDB_Computational", "MSigDB_Oncogenic_Signatures")
# enrichment_results <- enrichr(gene_list, databases)
# enrichment_results <- enrichment_results$MSigDB_Hallmark_2020 
# head(enrichment_results)
# 
# enrichment_results_plt_df <- data.frame(enrichment_results) %>%
#   arrange(Odds.Ratio) %>%
#   tail(8)
# 
# head(enrichment_results_plt_df)
# 
# enrichment_results_plt_df$Term <- factor(enrichment_results_plt_df$Term,
#                                          levels = enrichment_results_plt_df$Term)


# Define the gene set collections you want to use
# KEGG, Gene Ontology (GO), Reactome, etc.
gs.annots <- buildIdx(entrezIDs = rownames(v$E), 
                      species = "human", # Adjust for species (e.g., mouse = "mouse")
                      msigdb.gsets = c("c2", "c5"), # Molecular Signature Database collections
                      go.gsets = TRUE, 
                      kegg.gsets = TRUE)

# Perform EGSEA analysis
egsea.results <- egsea(fit = fit, 
                       contrasts = colnames(design)[2], # Specify the contrast of interest
                       gs.annots = gs.annots,
                       symbolsMap = v$genes, # Optional: Map to gene symbols
                       sort.by = "p.adj", 
                       report = TRUE, # Generates an HTML report
                       num.threads = 4) # Adjust threads for parallel processing



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
results <- DE_results %>% filter(Color %in% c("Clus14_01", "Clus14_05")) %>%
  arrange(adj.P.Val) %>%
  top_n(9)
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
  ) + ylim(0,15.5)


plot_title = "Vis_B3_2_DE_Pathway_14_consensus_final.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 7,
    height = 5.5)
print(p)
dev.off()









