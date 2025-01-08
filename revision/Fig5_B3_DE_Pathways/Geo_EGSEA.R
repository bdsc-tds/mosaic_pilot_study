# 3 samples in area A

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
# BiocManager::install("EGSEA")

###########################################################################
# Prep dataset for DE between malig subtypes only --------------------------
geo_all <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds")
geo_malig <- geo_all[, geo_all$section_id == "B3_1" &
                       geo_all$cell_fraction == "Malignant" # & 
                       # geo_all$sample_id2 != "DSP-1001660013520-C-B12"
                     ]
clus14_id <- c("DSP-1001660013520-C-B02", "DSP-1001660013520-C-B06", "DSP-1001660013520-C-B12")
clus159_id <- c("DSP-1001660013520-C-A03", "DSP-1001660013520-C-A05", "DSP-1001660013520-C-A07",
                "DSP-1001660013520-C-A09", "DSP-1001660013520-C-A12", "DSP-1001660013520-C-B04")
geo_malig$malig_sub <- ifelse(geo_malig$sample_id2 %in% clus14_id, "clus14",
                              ifelse(geo_malig$sample_id2 %in% clus159_id, "clus159", NA))

# DE ----------------------------------------------------------------------
sce <- geo_malig
# Remove negative control probe
sce <- sce[rownames(sce) != "NegProbe-WTX", ]

dge <- SE2DGEList(sce)
design <- model.matrix(~0 + malig_sub, data = colData(sce))
colnames(design) <- gsub(" ","_",colnames(design))
colnames(design)
#  "malig_subclus14"  "malig_subclus159"
table(sce$malig_sub)
# clus14 clus159 
# 3      6 

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

results_efit<- limma::decideTests(efit, p.value = 0.05, lfc = 1)
summary_efit <- summary(results_efit)

summary_efit
# BvT
# Down      12
# NotSig 18619
# Up        26


# EGSEA -------------------------------------------------------------------
## Map gene symbols to Entrez IDs (not needed)
# library(biomaRt) 
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# # Define your gene symbols
# gene_symbols <- rownames(v$E)
# 
# # Map gene symbols to Entrez IDs
# gene_mapping <- getBM(attributes = c("hgnc_symbol", "entrezgene_id", "external_gene_name"),
#                       filters = "hgnc_symbol",  # Input type: Gene symbols
#                       values = gene_symbols,    # Your input genes
#                       mart = mart)
# 
# # View the results
# head(gene_mapping)

# Define the gene set collections you want to use
# KEGG, Gene Ontology (GO), Reactome, etc.
# gs.annots <- buildIdx(entrezIDs = gene_mapping$entrezgene_id, # 1-1 mapping
rownames(v) <- v$genes$GeneID

gs.annots <- buildIdx(entrezIDs = rownames(v), # includes pseudo-genes, like ABCC6P2, a non-coding part of ABCC6
                      species = "human", # Adjust for species (e.g., mouse = "mouse")
                      msigdb.gsets = "h" # Molecular Signature Database collections
                      # kegg.exclude = "Metabolism"
                      )

names(gs.annots)
# [1] "h"   "kegg"
summary(gs.annots$kegg)
summary(gs.annots$h)

egsea.sort()
# Perform EGSEA analysis
egsea.results <- egsea(voom.results = v,
                       contrasts = contr.matrix, 
                       gs.annots = gs.annots,
                       logFC.cutoff = 1,
                       # symbolsMap = v$genes$TargetName, # Optional: Map to gene symbols
                       sort.by = "p.adj", 
                       report = FALSE, 
                       num.threads = 4) 


t = topSets(egsea.results, contrast = 1, gs.label = "h", # sort.by = "ora",
            number = 50, names.only = FALSE)
head(t, 10)
#                                    Rank      p.value      p.adj vote.rank avg.rank med.rank   min.pvalue min.rank avg.logfc avg.logfc.dir direction significance
# HALLMARK_ANDROGEN_RESPONSE            1 0.0009919477 0.02371646         5 11.83333      9.5 8.269991e-05        1  1.372742     -1.430874      Down     94.50649
# HALLMARK_CHOLESTEROL_HOMEOSTASIS      2 0.1338648260 0.23080142         5 11.58333      6.0 1.190476e-02        2  1.313288     -1.313288      Down     34.80458
# HALLMARK_IL6_JAK_STAT3_SIGNALING      3 0.1338648260 0.23080142         5 11.75000      8.5 1.190476e-02        1  1.450419     -1.450419      Down     38.54324
# HALLMARK_FATTY_ACID_METABOLISM        4 0.1338648260 0.23080142        15 20.33333     19.0 1.190476e-02        4  1.137642     -1.137642      Down     30.01585
# HALLMARK_MTORC1_SIGNALING             5 0.0014448255 0.02371646        50 21.08333     12.0 1.204819e-04        2  1.451701     -1.619116      Down    100.00000
# HALLMARK_ESTROGEN_RESPONSE_EARLY      6 0.1338648260 0.23080142        25 18.58333     19.0 1.190476e-02        5  1.097768     -1.097768      Down     28.92874
# HALLMARK_TGF_BETA_SIGNALING           7 0.0004284493 0.02142247        10 12.00000      7.0 3.571112e-05        1  1.079633     -1.079633      Down     76.15606
# HALLMARK_INTERFERON_ALPHA_RESPONSE    8 0.1338648260 0.23080142        35 23.75000     27.0 1.190476e-02        4  1.475737      1.475737        Up     39.23351
# HALLMARK_SPERMATOGENESIS              9 0.2367987534 0.38047815        20 20.16667     19.5 2.226779e-02        2  1.229596     -1.229596      Down     21.09379
# HALLMARK_UV_RESPONSE_UP              10 0.4431625818 0.55395323        25 27.75000     27.5 4.761905e-02       10  1.229596     -1.229596      Down     12.50490

table(t$direction)
# Down   Up 
# 42    8 

# t_top8_blue <- t %>%
#   mutate(pathway = rownames(t)) %>%
#   filter(direction == "Down") %>%
#   slice_head(n = 8)
# 
# t_top8_pink <- t %>%
#   mutate(pathway = rownames(t)) %>%
#   filter(direction == "Up") %>%
#   slice_head(n = 8)
t$Description = rownames(t)
save_path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/EGSEA"
write.csv(t, file.path(save_path, "geo_EGSEA.csv"))


