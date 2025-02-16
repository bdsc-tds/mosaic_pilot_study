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

disease = "dlbcl"

# After spotclean -------------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
savepath = savepath.spotclean

# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(savepath, save_rds_name))

Idents(seu) <- seu$seurat_clusters

# Key tumor and annotation markers DLBCL
tu_dlbcl <- c(# "CD19", "CD79A", "BCL6", "IRF4", "MYC", "BCL2", "MKI67", "TP53", # Tumor
  "IGKC", "IGHG3", "JCHAIN", "IGHG1",                      # Plasma
  "MUCL3", "LCN2", "TFF1", "MUC5AC",                               # Epithelium
  "THBS1", "POSTN", "COL3A1", "COL1A2", "COL1A1",                  # Stroma
  "ACKR1", "MEOX2", "CCL14", "TSPAN7", "PLA1A", "VWF", "GIMAP4",    # Vessels
  "CD4", "CD28", "CD68", "CD163", "CD14", "SIRPA"            # Immune cells
)

p <- SCpubr::do_DotPlot(sample = seu, features = tu_dlbcl, 
                        axis.text.x.angle = 0,
                        flip = TRUE# , legend.position = "right"
) + 
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8), 
        axis.line = element_blank()) + 
  scale_fill_viridis_c(option = "mako") + 
  scale_x_discrete(position = "top") +
  guides(fill = guide_colorbar(title = "Avg. Expression", title.position = "top", title.hjust = 0.5, 
                               barwidth = unit(5, "cm"), barheight = unit(0.5, "cm")))
p


#########################################################################
#                             Complex Heatmap                           #
#########################################################################
# Pull the data from a ggplot object --------------------------------------
df<- p$data
head(df)

df_clus_grp <- data.frame(
  annot = seu$annot,
  id = seu$seurat_clusters
  ) %>%
  unique()

Plasma = c("IGKC", "IGHG3", "JCHAIN", "IGHG1")
Epithelium = c("MUCL3", "LCN2", "TFF1", "MUC5AC")
Stroma = c("THBS1", "POSTN", "COL3A1", "COL1A2", "COL1A1") 
Vessels = c("ACKR1", "MEOX2", "CCL14", "TSPAN7", "PLA1A", "VWF", "GIMAP4") 
Immune_cells = c("CD4", "CD28", "CD68", "CD163", "CD14", "SIRPA")

df_save <- df %>%
  left_join(df_clus_grp) %>%
  mutate(gene_grp = case_when(features.plot %in% Plasma ~ "Plasma",
                              features.plot %in% Epithelium ~ "Epithelium",
                              features.plot %in% Stroma ~ "Stroma",
                              features.plot %in% Vessels ~ "Vessels",
                              features.plot %in% Immune_cells ~ "Immune cells"
                              ))

write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS13d.csv")


### the matrix for the scaled expression ---------------------------------- 
# Note: Seurat plots scaled expression already wow, 
# so now scanpy standard_scale = "var" is the same as Seurat
exp_mat<-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

head(exp_mat)


## the matrix for the percentage of cells express a gene ------------------
percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

head(percent_mat)


# the range is from 0 - 100 ----------------------------------------------
range(percent_mat)

## these two matrix have the same dimension
dim(exp_mat)
dim(percent_mat)

# Palette ---------------------------------------------------------------
library(viridis)
library(Polychrome)

Polychrome::swatch(viridis(20, option = "mako"))


# -------------------------------------------------------------------------
## get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

  #        10%        50%        90%        99% 
  # 0.09311773 0.68158384 2.10158590 4.95134357 

## any value that is greater than 2 will be mapped to yellow -------------
col_fun = circlize::colorRamp2(c(0, 2, 4.5), viridis(20, option = "mako", direction = -1)[c(1, 10, 20)])
# > viridis(20)[c(1,10, 20)]
# [1] "#440154FF" "#238A8DFF" "#FDE725FF"

cell_fun = function(j, i, x, y, w, h, fill){
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}


#########################################################################
# Add col annotations
#########################################################################
colnames(exp_mat) 

cluster_anno<-  c("Tumor", "Stroma", "Necrosis", "Tumor", "Tumor", "Tumor", 
                  "Epithelium", "Vessels/Immune", "Tumor", "Plasma", "Plasma", "Tumor",
                  "Necrosis", "Plasma", "Necrosis", "Tumor", "Necrosis",
                  "Tumor", "Tumor")

# unique(cluster_anno) # "Tumor", "Stroma", "Necrosis", "Epithelium", "Vessels/Immune", "Plasma"
column_ha <- HeatmapAnnotation(
  cluster_anno = cluster_anno,
  col = list(cluster_anno = setNames(c("#FFD700", "#388E8E", "#666666", "#BC8F8F", "#FF7256", "#CD00CD"), unique(cluster_anno))), 
  na_col = "grey",
  show_annotation_name = FALSE
)


rownames(exp_mat)
gene_grp = c(rep("Plasma", 4), rep("Epithelium", 4), rep("Stroma", 5), rep("Vessels", 7), rep("Immune", 6))
# row_ha <- rowAnnotation(gene_grp = gene_grp, 
#                         col = list(gene_grp = setNames(brewer.pal(6, "Dark2"), unique(gene_grp))),
#                         na_col = "grey")
# 
# row_ha <- rowAnnotation(gene_grp = anno_block(gp = gpar(fill = gene_grp, 
#                                                         col = list(gene_grp = setNames(brewer.pal(8, "Dark2"), unique(gene_grp))),
#                                                         na_col = "grey")),
#                         width = unit(2, "mm"))

row_ha_test <- rowAnnotation(block = anno_block(gp = gpar(fill = c("Epithelium" = "#BC8F8F",  
                                                                   "Immune" = "orange", 
                                                                   "Plasma" = "#CD00CD",
                                                                   "Stroma" = "#388E8E", 
                                                                   "Vessels" = "#FB9A99"), col = NA)), 
                             width = unit(2, "mm"))

ht <- Heatmap(exp_mat,
              heatmap_legend_param=list(title="expression"),
              column_title = "Integrated DLBCL Visium Seurat Clusters", 
              column_title_gp = gpar(fontface = "bold"),
              column_names_rot = 45,
              row_names_rot = 0,
              row_title_rot = 0, 
              col=col_fun,
              rect_gp = gpar(type = "none", col = "white"),
              cell_fun = cell_fun,
              column_names_gp = gpar(fontsize = 10),
              row_names_gp = gpar(fontsize = 10),
              row_names_side = "right",
              row_title_side = "left",
              column_names_side = "bottom",
              row_dend_side = "right",
              cluster_rows = FALSE,
              show_row_dend = FALSE, 
              cluster_columns = TRUE,
              show_column_dend = TRUE, 
              row_title_gp = gpar(type = "none"),
              row_split = gene_grp,
              border = "black",
              top_annotation = column_ha
) 
draw(ht)

draw(row_ha_test + ht)

#########################################################################
# Add legend
#########################################################################
lgd_list1 = 
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "% expressed",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0  * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                             gp = gpar(fill = "black")))
  )

p <- draw(row_ha_test + ht, 
     annotation_legend_list = lgd_list1 ,
     merge_legend = TRUE
)



# -------------------------------------------------------------------------
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "/DLBCL_Vis_Marker_Complex_DotPlot.pdf"), width = 10, height = 7)
print(p)
dev.off()