library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SCpubr)
library(patchwork)

library(tidyverse)
library(presto)
library(ComplexHeatmap)
library(circlize)

library(scran)
library(RColorBrewer)


disease = "dlbcl"
# spe3 <- readRDS(file.path(data_path, paste0("dlbcl", "_spe_ruv.rds")))
# disease = "breast"
# spe2 <- readRDS(file.path(data_path, paste0("breast", "_spe_ruv.rds")))


# -------------------------------------------------------------------------
# disease = "lung"

source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/GeoMx/00_GeoMx_Paths.R")
data_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
spe <- readRDS(file.path(data_path, paste0(disease, "_spe_ruv.rds"))) # expression on normalized and batch corrected counts 

if(class(assay(spe, "logcounts"))[1] != "dgCMatrix"){assay(spe, "logcounts") <- as(as.matrix(assay(spe, "logcounts")), "dgCMatrix")}

set.seed(123)
nn.clusters <- scran::clusterCells(spe, use.dimred="PCA")
table(nn.clusters)
spe$cluster <- nn.clusters


# -------------------------------------------------------------------------
seu <- as.Seurat(spe)
seu <- ScaleData(seu)

Idents(seu) <- seu$cell_fraction

all.markers <- FindAllMarkers(object = seu, only.pos = TRUE)
all.markers %>%
  group_by(cluster) %>%
  arrange(avg_log2FC) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 400) %>%
  ungroup() -> top2000


# cate_var <- c(cluster, pathology, roi_type, patient_id, section_id, cell_fraction)
# cont_var <- c(area, nuclei, genes_detected)

#########################################################################
#                             Complex Heatmap                           #
#########################################################################
exp_mat_test <- as(logcounts(spe), "matrix") 
exp_mat <- exp_mat_test[rownames(exp_mat_test) %in% top2000$gene, ]

dim(exp_mat)

df_save <- data.frame(
  t(exp_mat),
  # cluster = spe$cluster,
  area = spe$area, 
  nuclei = spe$nuclei, 
  genes_detected = spe$genes_detected, 
  pathology = spe$pathology, 
  roi_type = spe$roi_type, 
  patient_id = spe$patient_id, 
  section_id = spe$section_id, 
  cell_fraction = spe$cell_fraction)

write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS5c.csv")

# Brewer palettes ----------------------------------------------------------
# Note: does not help to order the color, whatever...
paired_palette <- brewer.pal(12, "Paired")
set3_palette <- brewer.pal(12, "Set3")
set2_palette <- brewer.pal(8, "Set2")
set1_palette <- brewer.pal(9, "Set1")
dark2_palette <- brewer.pal(8, "Dark2")
accent_palette <- brewer.pal(8, "Accent")

# cont_col_fun1 <- colorRamp2(range(spe$area), c("white", "darkorchid"))
# cont_col_fun2 <- colorRamp2(range(spe$nuclei), c("white", "deeppink"))
# cont_col_fun3 <- colorRamp2(range(spe$genes_detected), c("white", "darkblue"))
#########################################################################
# Add col annotations
#########################################################################
if(disease %in% c("breast", "lung")){
  column_ha <- HeatmapAnnotation(
    cluster = spe$cluster,
    area = spe$area, # anno_simple(spe$area, col = cont_col_fun1),
    nuclei = spe$nuclei, # anno_simple(, col = cont_col_fun2),
    genes_detected = spe$genes_detected, # anno_simple(spe$genes_detected, col = cont_col_fun3),
    pathology = spe$pathology, 
    roi_type = spe$roi_type, 
    patient_id = spe$patient_id, 
    section_id = spe$section_id, 
    cell_fraction = spe$cell_fraction,
    col = list(cluster = setNames(paired_palette[1:length(unique(spe$cluster))], unique(spe$cluster)),
               pathology = setNames(set3_palette[1:length(unique(spe$pathology))], unique(spe$pathology)), 
               roi_type = setNames(set2_palette[1:length(unique(spe$roi_type))], unique(spe$roi_type)), 
               patient_id = setNames(set1_palette[1:length(unique(spe$patient_id))], unique(spe$patient_id)), 
               section_id = setNames(dark2_palette[1:length(unique(spe$section_id))], unique(spe$section_id)), 
               cell_fraction = setNames(accent_palette[1:length(unique(spe$cell_fraction))], unique(spe$cell_fraction))
    ),
    na_col = "grey",
    annotation_name_side = "left", 
    show_annotation_name = TRUE
  )
}else{ # dlbcl 
  column_ha <- HeatmapAnnotation(
    cluster = spe$cluster,
    area = spe$area, 
    genes_detected = spe$genes_detected, 
    patient = spe$patient, 
    cell_fraction = spe$cell_fraction,
    col = list(cluster = setNames(paired_palette[1:length(unique(spe$cluster))], unique(spe$cluster)),
               patient = setNames(set1_palette[1:length(unique(spe$patient))], unique(spe$patient)), 
               cell_fraction = setNames(accent_palette[1:length(unique(spe$cell_fraction))], unique(spe$cell_fraction))
    ),
    na_col = "grey",
    annotation_name_side = "left", 
    show_annotation_name = TRUE
  )
}


ht <- Heatmap(exp_mat,
              name = "expression",
              row_title = "# DE Gene = 2000", 
              column_names_rot = 45,
              row_names_rot = 0,
              row_title_rot = 90, 
              column_names_gp = gpar(fontsize = 10),
              row_title_side = "left",
              row_dend_side = "left",
              cluster_rows = TRUE,
              show_row_dend = TRUE, 
              cluster_columns = TRUE,
              show_column_dend = TRUE, 
              show_row_names = FALSE, 
              show_column_names = FALSE, 
              row_title_gp = gpar(fontsize = 10),
              top_annotation = column_ha
) 
 
p <- draw(ht, merge_legend = TRUE)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/GeoMx_Expr_Heatmap"
pdf(file.path(figpath, paste0(disease, "_geo_Complex_Heatmap.pdf")), width = 17, height = 12)
print(p)
dev.off()
