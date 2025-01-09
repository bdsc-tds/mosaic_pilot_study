library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(ggrepel)
library(ggspavis)
library(patchwork)
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(msigdbr)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/DiffDE_GSEA"

###########################################################################
contrast <- "Geo_3vs6_DE"
DE_result <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))
results <- as.data.frame(DE_result)
# Geo GSEA --------------------------------------------------------------------

ranked_genes <- results$logFC
names(ranked_genes) <- results$TargetName    # Use gene symbols
ranked_genes <- ranked_genes %>% sort(., decreasing = TRUE)

gsea_result <- GSEA(geneList = ranked_genes,
                    TERM2GENE = hallmark_sets[, c("gs_id", "gene_symbol")],
                    TERM2NAME = hallmark_sets[, c("gs_id", "gs_name")],
                    pvalueCutoff = 1,
                    minGSSize = 1,
                    verbose = TRUE)
gsea_result <- gsea_result %>% as.data.frame()

write.csv(gsea_result, file.path(save_path, "Geo_3vs6_GSEA.csv"))

table(gsea_result$enrichmentScore > 0)
# FALSE  TRUE 
# 10     8 


# Visium GSEA -------------------------------------------------------------
contrast <- "1&5&9_14"
DE_result <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))
results <- as.data.frame(DE_result)
results$avg_log2FC <- ifelse(results$cluster == "14", results$avg_log2FC * -1, results$avg_log2FC)

ranked_genes <- results$avg_log2FC
names(ranked_genes) <- results$gene    # Use gene symbols
ranked_genes <- ranked_genes %>% sort(., decreasing = TRUE)

gsea_result <- GSEA(geneList = ranked_genes,
                    TERM2GENE = hallmark_sets[, c("gs_id", "gene_symbol")],
                    TERM2NAME = hallmark_sets[, c("gs_id", "gs_name")],
                    pvalueCutoff = 1,
                    minGSSize = 1,
                    verbose = TRUE)
gsea_result <- gsea_result %>% as.data.frame()

write.csv(gsea_result, file.path(save_path, "Vis_GSEA.csv"))

table(gsea_result$enrichmentScore > 0)
# FALSE  TRUE 
# 17    21 

# Chromium GSEA -----------------------------------------------------------
contrast <- "chrom_B3_tu_2_clus_markers"
DE_result <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))
results <- as.data.frame(DE_result)
results$avg_log2FC <- ifelse(results$cluster == "Tu_B3_PLA2G2A", results$avg_log2FC * -1, results$avg_log2FC)

ranked_genes <- results$avg_log2FC
names(ranked_genes) <- results$gene    # Use gene symbols
ranked_genes <- ranked_genes %>% sort(., decreasing = TRUE)

gsea_result <- GSEA(geneList = ranked_genes,
                    TERM2GENE = hallmark_sets[, c("gs_id", "gene_symbol")],
                    TERM2NAME = hallmark_sets[, c("gs_id", "gs_name")],
                    pvalueCutoff = 1,
                    minGSSize = 1,
                    verbose = TRUE)
gsea_result <- gsea_result %>% as.data.frame()

write.csv(gsea_result, file.path(save_path, "Chrom_GSEA.csv"))

table(gsea_result$enrichmentScore > 0)
# FALSE  TRUE 
# 8     9 


