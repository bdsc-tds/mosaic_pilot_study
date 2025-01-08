library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Gallery_Patho_Decon_Pt_Heatmap"

# disease = "breast"
# disease = "lung"
# disease = "dlbcl"

disease_list <- c("breast", "lung", "dlbcl")

for(disease in disease_list){
  # After spotclean -----------------------------------------------------
  datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean/Results")
  save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
  
  assign(paste0("seu_", disease), readRDS(paste0(datapath.spotclean, save_rds_name)))
}


# -------------------------------------------------------------------------
plotHeatmap <- function(seu, text = "region"){ # scale to per region percentage, or just counts
  
  df <- as.data.frame(table(seu@meta.data[["Region"]], seu@meta.data[["Level4_decon_max"]]))
  df2 <- df %>%
    group_by(Var1) %>%
    mutate(countT= sum(Freq)) %>%
    mutate(per=round(100*Freq/countT,0))
  
  data_frame <- data.frame(patho = df2$Var1, decon = df2$Var2, per = df2$per)
  # data_frame <- data.frame(patho = df2$Var1, decon = df2$Var2, per = df2$Freq)
  data_frame$patho <- factor(data_frame$patho, levels = rev(levels(factor(data_frame$patho))))
  
  # Create the heatmap patho per decon
  p <- ggplot(data_frame, aes(x = decon, y = patho, fill = per)) + 
    geom_tile() + 
    geom_text(aes(x = decon, y = patho, label = paste0(per, "%"))) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left") +
    scale_fill_viridis_c(option = "mako") + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5)
    ) + 
    labs(fill = paste0("Per-", text, " Percentage"), title = "", y = "Region") + 
    # labs(fill = paste0("Counts"), title = "", y = "Region") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
  
  return(p)
}


# -------------------------------------------------------------------------
p_breast_patho_decon <- plotHeatmap(seu_breast, text = "region")
p_lung_patho_decon <- plotHeatmap(seu_lung, text = "region")
p_dlbcl_patho_decon <- plotHeatmap(seu_dlbcl, text = "region")

p_final <- (p_breast_patho_decon + theme(axis.title.x = element_blank())) / 
  (p_lung_patho_decon + theme(axis.title.x = element_blank())) / 
  (p_dlbcl_patho_decon + theme(axis.title.x = element_blank())) + 
  plot_layout(# guides = "collect", 
              heights = c(1.2, 1.1, 0.7)) & 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 8),        
        legend.title = element_text(size = 10))



pdf(file.path(figpath, "Vis_Patho_Decon_Heatmap_region.pdf"), width = 12, height = 17)
# pdf(file.path(figpath, "Vis_Patho_Decon_Heatmap_counts.pdf"), width = 10, height = 17)
print(p_final)
dev.off()


# ##########################################################################
# plotHeatmap <- function(seu, text = "cell-type"){ # scale to per decon cell type percentage
#   
#   df <- as.data.frame(table(seu@meta.data[["Region"]], seu@meta.data[["Level4_decon_max"]]))
#   df2 <- df %>%
#     group_by(Var2) %>%
#     mutate(countT= sum(Freq)) %>%
#     mutate(per=round(100*Freq/countT,0))
#   
#   data_frame <- data.frame(patho = df2$Var1, decon = df2$Var2, per = df2$per)
#   data_frame$patho <- factor(data_frame$patho, levels = rev(levels(factor(data_frame$patho))))
#   
#   # Create the heatmap patho per decon
#   p <- ggplot(data_frame, aes(x = decon, y = patho, fill = per)) + 
#     geom_tile() + 
#     scale_x_discrete(position = "bottom") +
#     scale_y_discrete(position = "left") +
#     scale_fill_viridis_c(option = "mako") + 
#     theme_minimal() +
#     theme(legend.position = "bottom",
#           axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5)
#     ) + 
#     labs(fill = paste0("Per-", text, " Percentage"), title = "", y = "Region") + 
#     guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
#                                  barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
#   
#   return(p)
# }
# 
# 
# # -------------------------------------------------------------------------
# p_breast_patho_decon <- plotHeatmap(seu_breast, text = "cell-type")
# p_lung_patho_decon <- plotHeatmap(seu_lung, text = "cell-type")
# p_dlbcl_patho_decon <- plotHeatmap(seu_dlbcl, text = "cell-type")
# 
# p_final <- (p_breast_patho_decon + theme(axis.title.x = element_blank())) / 
#   (p_lung_patho_decon + theme(axis.title.x = element_blank())) / 
#   (p_dlbcl_patho_decon + theme(axis.title.x = element_blank())) + 
#   plot_layout(guides = "collect", 
#     heights = c(1.2, 1.1, 0.7)) & 
#   theme(legend.position = "bottom", 
#         legend.text = element_text(size = 8),        
#         legend.title = element_text(size = 10))
# 
# 
# 
# pdf(file.path(figpath, "Vis_Patho_Decon_Heatmap_cell_type.pdf"), width = 10, height = 17)
# print(p_final)
# dev.off()


