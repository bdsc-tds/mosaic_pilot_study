library(BayesSpace)
library(patchwork)
library(ggplot2)
library(patchwork)

# disease = "DLBCL"
# sample = "DLBCL_6"
# disease = "Lung"
# sample = "L1_2"
disease = "Breast"
sample = "B1_2"
# resolution = "spots"
# resolution = "subspots"
read_path = paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/", disease, "/", sample)

# if(resolution == "spot"){ # spot
sce <- readRDS(file.path(read_path, paste0(sample, "_baye_clustered.rds"))) # spot
assay(sce, "log1p") <- log1p(counts(sce))
# raw_assay_name = "counts"

# }else{ # subspot
sce_enhanced <- readRDS(file.path(read_path, paste0(sample, "_baye_clustered_all_enhanced_expr.rds"))) # subspot
# raw_assay_name = "decont_subspot"
# }

dim(rowData(sce_enhanced))

# -------------------------------------------------------------------------
markers <- c("ITGB6")

plot_expression <- function(sce_obj = sce, marker = "CD4"){
  featurePlot(sce_obj, marker, assay.type = "log1p", color=NA) +
    scale_x_reverse() + 
    viridis::scale_fill_viridis(option="A") +
    labs(title=marker, fill="Log\nexpression") +
    theme(plot.title = element_blank())
  # theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
}

# With legend
feature.plots1 <- plot_expression(sce, markers)
enhanced.feature.plots1 <- plot_expression(sce_enhanced, markers)

(feature.plots1 | enhanced.feature.plots1) + 
  plot_annotation(title = paste0(sample, ": spot vs. subspot ITGB6 expression"))


