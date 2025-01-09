library(edgeR)
library(limma)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Prep dataset for DE between malig subtypes only --------------------------
chrom_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"

chrom_breast_qcd <- readRDS(file.path(chrom_path, "chrom_breast.rds"))
chrom_B3_tu <- chrom_breast_qcd[, chrom_breast_qcd$patient == "B3" & chrom_breast_qcd$Level2 == "Tu_B3"] # 17170


chrom_DE <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/DiffDE/chrom_B3_tu_2_clus_markers.csv")
chrom_DE <- chrom_DE %>%
  mutate(avg_log2FC = ifelse(cluster == "Tu_B3_PLA2G2A", -1 * avg_log2FC, avg_log2FC))


# Map to entrezID ---------------------------------------------------------
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# DE genes
gene_mapping_DE <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                         filters = "hgnc_symbol", 
                         values = chrom_DE$gene,    
                         mart = mart) # 40
gene_mapping_DE <- gene_mapping_DE %>%
  filter(!is.na(entrezgene_id)) # 40

# all chrom QCed genes
gene_mapping_chrom <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                          filters = "hgnc_symbol", 
                          values = rownames(chrom_B3_tu),  # 17170
                          mart = mart) # 16873
gene_mapping_chrom <- gene_mapping_chrom %>%
  filter(!is.na(entrezgene_id)) # 16837


gs.annots = buildIdx(entrezIDs = chrom_DE$entrezgene_id, species = "human",
                     msigdb.gsets = "h")

chrom_DE <- chrom_DE %>%
  dplyr::rename(hgnc_symbol = gene) %>%
  inner_join(gene_mapping_DE, by = "hgnc_symbol")
deGenes = chrom_DE$entrezgene_id
logFC = chrom_DE$avg_log2FC
names(logFC) = deGenes

egsea = egsea.ora(geneIDs = deGenes, 
                  universe = gene_mapping_chrom$entrezgene_id, # all genes mapping
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
# 1    4 
# Paths ----------------------------------------------------------------------
save_path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/EGSEA"
write.csv(t, file.path(save_path, "chrom_EGSEA.csv"))




