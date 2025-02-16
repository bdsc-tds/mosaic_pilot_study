library(ggridges)
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
# disease = "breast"; seu <- readRDS(file.path(chrompath, "chrom_breast.rds"))
# disease = "lung"; seu <- readRDS(file.path(chrompath, "chrom_lung.rds"))
disease = "dlbcl"; seu <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

if(disease == "breast"){
  source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/Owkin_Pilot/Code/color_palette.R")
  source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  CT_color <- breast_1_5ct_color
  CT_color <- breast_level4_cellcolors
}else if(disease == "lung"){
  source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/Owkin_Pilot/Code/color_palette.R")
  source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  CT_color <- lung_1_5ct_color
  CT_color <- lung_level4_cellcolors
}else{
  source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/Owkin_Pilot/Code/color_palette.R")
  source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  CT_color <- dlbcl_1_5ct_color
  CT_color <- dlbcl_level4_cellcolors
}

plotBar_seu <- function(seu, decreasing = TRUE, ylabel = "Percentage", annotate = "Harmonised_Level4"
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
    scale_fill_manual(values = CT_color[names(table(CD[[annotate]]))],
                      na.value = "#d3d3d3")
  
  p
}

df_save <- data.frame(CT4 = seu$Harmonised_Level4,
                 CT1_5 = seu$Level1,
                 Barcode = colnames(seu)
                 )
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS4_dlbcl.csv")     
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS4_breast.csv")  
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS4_lung.csv")  

p <- plotBar_seu(seu, decreasing = TRUE, ylabel = "Percentage", annotate = "Harmonised_Level4")
p <- plotBar_seu(seu, decreasing = TRUE, ylabel = "Percentage", annotate = "level1_5")

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Descriptive/Chrom_bar/"
pdf(file = file.path(fig_path, "Chrom_breast_bar_level4.pdf"),
    width = 15,
    height = 5)
print(p)
dev.off()

pdf(file = file.path(fig_path, "Chrom_breast_bar_level1_5.pdf"),
    width = 4,
    height = 5)
print(p)
dev.off()


# -------------------------------------------------------------------------
pdf(file = file.path(fig_path, "Chrom_lung_bar_level4.pdf"),
    width = 15,
    height = 5)
print(p)
dev.off()

pdf(file = file.path(fig_path, "Chrom_lung_bar_level1_5.pdf"),
    width = 4,
    height = 4)
print(p)
dev.off()


# -------------------------------------------------------------------------
pdf(file = file.path(fig_path, "Chrom_dlbcl_bar_level4.pdf"),
    width = 15,
    height = 5)
print(p)
dev.off()

pdf(file = file.path(fig_path, "Chrom_dlbcl_bar_level1_5.pdf"),
    width = 4,
    height = 6.5)
print(p)
dev.off()
