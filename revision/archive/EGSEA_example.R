# EGSEA on a count matrix
# EGSEA on a list of genes
library(EGSEA)
library(EGSEAdata)
data(il13.data)
voom.results = il13.data$voom
contrast = il13.data$contra
# find Differentially Expressed genes
library(limma)
##
## Attaching package: ’limma’
## The following object is masked from ’package:BiocGenerics’:
##
## plotMA
vfit = lmFit(voom.results, voom.results$design)
vfit = contrasts.fit(vfit, contrast)
vfit = eBayes(vfit)
# select DE genes (Entrez IDs and logFC) at p-value <= 0.05
# and |logFC| >= 1
top.Table = topTable(vfit, coef = 1, number = Inf, p.value = 0.05, lfc = 1)
# > head(top.Table)
# FeatureID Symbols     logFC  AveExpr         t      P.Value    adj.P.Val        B
# 30835     30835   CD209  3.882492 4.818066  25.93935 1.598970e-10 2.018738e-06 14.52844
# 3554       3554   IL1R1  2.303161 5.728493  22.98881 5.267500e-10 2.339320e-06 13.62705
# 2208       2208   FCER2  5.818129 3.356316  24.97210 2.328014e-10 2.018738e-06 13.25669
# 2675       2675   GFRA2  3.277526 3.319754  22.93290 5.395423e-10 2.339320e-06 13.13380
# 55022     55022    PID1 -3.454681 4.012016 -21.56122 9.900617e-10 3.434128e-06 12.71054
# 7850       7850   IL1R2  3.336171 4.158437  20.19777 1.880207e-09 5.434739e-06 12.26003
deGenes = as.character(top.Table$FeatureID)
# length(deGenes)
# # [1] 453
# head(deGenes)
# # [1] "30835" "3554"  "2208"  "2675"  "55022" "7850" 
logFC = top.Table$logFC
names(logFC) = deGenes
# > head(logFC)
# 30835      3554      2208      2675     55022      7850 
# 3.882492  2.303161  5.818129  3.277526 -3.454681  3.336171 
# > table(logFC > 0)
# 
# FALSE  TRUE 
# 158   295 

# build the gene set collection index
gs.annots = buildIdx(entrezIDs = deGenes, species = "human",
                     msigdb.gsets = "none", kegg.exclude = c("Metabolism"))
# $kegg
# An object of class "GSCollectionIndex"
# Number of gene sets: 164
# Annotation columns: ID, GeneSet, NumGenes, Type
# Total number of indexing genes: 453
# Species: Homo sapiens
# Collection name: KEGG Pathways
# Collection uniqe label: kegg
# Database version: NA
# Database update date: 07 March 2017

## [1] "Building KEGG pathways annotation object ... "
# dim(voom.results$genes)
# [1] 17343     2

# > head(voom.results$genes)
# FeatureID      Symbols
# 1          1         A1BG
# 4       1000         CDH2
# 5      10000         AKT3
# 14 100009605       TRNAF1
# 17 100009613     ANO1-AS2
# 18 100009676 LOC100009676

# > head(top.Table[, c(1, 2)])
# FeatureID Symbols
# 30835     30835   CD209
# 3554       3554   IL1R1
# 2208       2208   FCER2
# 2675       2675   GFRA2
# 55022     55022    PID1
# 7850       7850   IL1R2
# > dim(top.Table)
# [1] 453   8
gsa = egsea.ora(geneIDs = deGenes, 
                universe = as.character(voom.results$genes[, 1]), # all genes mapping
                logFC = logFC, title = "X24IL13-X24", gs.annots = gs.annots,
                symbolsMap = top.Table[, c(1, 2)], #de genes only
                display.top = 5, 
                num.threads = 4, report = FALSE)

