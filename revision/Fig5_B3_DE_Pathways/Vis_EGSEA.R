library(edgeR)
library(limma)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(ggrepel)
library(ggspavis)
library(patchwork)

# Prep dataset for DE between malig subtypes only --------------------------
vis_all <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered.rds")
vis_malig <- vis_all[, vis_all$spatial.cluster %in% c("1", "5", "9", "14")]
vis_malig$malig_sub <- ifelse(vis_malig$spatial.cluster == "14", "clus14",
                              ifelse(vis_malig$spatial.cluster %in% c("1", "5", "9"), "clus159", NA))

vis_DE <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/DiffDE/1&5&9_14.csv")
vis_DE <- vis_DE %>%
  mutate(avg_log2FC = ifelse(cluster == "14", -1 * avg_log2FC, avg_log2FC))


# Map to entrezID ---------------------------------------------------------
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# DE genes
gene_mapping_DE <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                      filters = "hgnc_symbol", 
                      values = vis_DE$gene,    
                      mart = mart) # 140
gene_mapping_DE <- gene_mapping_DE %>%
  filter(!is.na(entrezgene_id)) # 140

# all vis QCed genes
gene_mapping_vis <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                         filters = "hgnc_symbol", 
                         values = rownames(vis_malig),  # 17883
                         mart = mart) # 17551
gene_mapping_vis <- gene_mapping_vis %>%
  filter(!is.na(entrezgene_id)) # 17516


gs.annots = buildIdx(entrezIDs = vis_DE$entrezgene_id, species = "human",
                     msigdb.gsets = "h")

vis_DE <- vis_DE %>%
  dplyr::rename(hgnc_symbol = gene) %>%
  inner_join(gene_mapping_DE, by = "hgnc_symbol")
deGenes = vis_DE$entrezgene_id
logFC = vis_DE$avg_log2FC
names(logFC) = deGenes

egsea = egsea.ora(geneIDs = deGenes, 
                universe = gene_mapping_vis$entrezgene_id, # all genes mapping
                logFC = logFC, gs.annots = gs.annots,
                symbolsMap = gene_mapping_DE, # de genes only
                display.top = 5, 
                num.threads = 4, report = FALSE)

t = topSets(egsea, contrast = 1, gs.label = "h", 
            number = 40, names.only = FALSE)
head(t)
t$Description = rownames(t)
table(t$direction)
# Down   Up 
# 14   17 
# Paths ----------------------------------------------------------------------
save_path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/EGSEA"
write.csv(t, file.path(save_path, "vis_EGSEA.csv"))
