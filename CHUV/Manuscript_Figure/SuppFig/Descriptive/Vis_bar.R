library(ggridges)
library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)

# disease = "breast"
# disease = "lung"
disease = "dlbcl"
# After spotclean -----------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
seu <- readRDS(paste0(datapath.spotclean, save_rds_name))

plt_df <- as.data.frame(table(seu$Region)) %>%
  dplyr::rename(Region = Var1,
                count = Freq) %>%
  mutate(per = count/sum(count)) 
plt_df_dlbcl <- plt_df; plt_df_dlbcl$indication = "DLBCL"
plt_df_breast <- plt_df; plt_df_breast$indication = "Breast"
plt_df_lung <- plt_df; plt_df_lung$indication = "Lung"

write.csv(rbind(plt_df_breast, 
                plt_df_lung,
                plt_df_dlbcl),
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS2c.csv")  


source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")
if(disease == "breast"){
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  patho_color <- breast_patho_color
}else if(disease == "lung"){
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  patho_color <- lung_patho_color
}else{
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  patho_color <- dlbcl_patho_color
}


# -------------------------------------------------------------------------
plotBar_seu <- function(seu, decreasing = TRUE, ylabel = "Percentage", annotate = "Region"
                        # ,title = "Pathology Region Frequency in Annotated Visium - Breast"
){
  CD <- seu@meta.data
  cnt <- plyr::count(CD[[annotate]])
  CD[["annotate"]] <- factor(CD[[annotate]],
                             levels = cnt$x[order(cnt$freq, decreasing = decreasing)])

  p <- ggplot(data = CD, aes(x = annotate, fill = annotate, color = NULL)) +
    geom_bar(aes(y = after_stat(count)/sum(after_stat(count)))) +
    theme_ridges() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    ylab(ylabel) + xlab("") +
    # ggtitle(title) +
    scale_fill_manual(values = patho_color[names(table(CD[[annotate]]))],
                      na.value = "#d3d3d3")

  p
}

p <- plotBar_seu(seu, decreasing = TRUE, ylabel = "Percentage", annotate = "Region")

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Descriptive/Vis_bar"
pdf(file = file.path(fig_path, "Vis_breast_bar.pdf"),
    width = 10,
    height = 10)
print(p)
dev.off()

pdf(file = file.path(fig_path, "Vis_lung_bar.pdf"),
    width = 10,
    height = 6)
print(p)
dev.off()

pdf(file = file.path(fig_path, "Vis_dlbcl_bar.pdf"),
    width = 5,
    height = 12)
print(p)
dev.off()

