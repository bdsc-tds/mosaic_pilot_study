library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(patchwork)

save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/DiffDE_GSEA"
resultfilegeo <- "Geo_3vs6_GSEA.csv"
resultfilevis <- "Vis_GSEA.csv"
resultfilechrom <- "Chrom_GSEA.csv"

# Venn diagram ------------------------------------------------------------
gsea_resultgeo <- read.csv(file.path(save_path_DE, resultfilegeo))
gsea_resultvis <- read.csv(file.path(save_path_DE, resultfilevis))
gsea_resultchrom <- read.csv(file.path(save_path_DE, resultfilechrom))

get_top8_path_df <- function(gsea_resultdf){
  gsea_resultdf_ <- gsea_resultdf %>%
    mutate(group = ifelse(enrichmentScore > 0, "pink", "blue")) %>%
    group_by(group) %>%
    arrange(group, desc(abs(enrichmentScore))) %>%
    slice_head(n = 8)
  return(gsea_resultdf_)
}

gsea_resultgeo_ <- get_top8_path_df(gsea_resultgeo)
gsea_resultvis_ <- get_top8_path_df(gsea_resultvis)
gsea_resultchrom_ <- get_top8_path_df(gsea_resultchrom)

# install.packages("VennDiagram")
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
dir = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/results/DiffDE"
setwd(dir)
plot_venn <- function(set1, set2, set3, name){
  VennDiagram::venn.diagram(
    x = list(set1, set2, set3),
    category.names = c("Geo_Paths" , "Vis_Paths " , "Chrom_Paths"),
    filename = name,
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol
  )
}

plot_venn(set1 = gsea_resultgeo$Description[gsea_resultgeo$enrichmentScore < 0], 
          set2 = gsea_resultvis$Description[gsea_resultvis$enrichmentScore < 0],
          set3 = gsea_resultchrom$Description[gsea_resultchrom$enrichmentScore < 0], 
          name = "blue_paths.png")
plot_venn(set1 = gsea_resultgeo$Description[gsea_resultgeo$enrichmentScore > 0], 
          set2 = gsea_resultvis$Description[gsea_resultvis$enrichmentScore > 0],
          set3 = gsea_resultchrom$Description[gsea_resultchrom$enrichmentScore > 0], 
          name = "pink_paths.png")

plot_venn(set1 = gsea_resultgeo_$Description[gsea_resultgeo_$group == "blue"], 
          set2 = gsea_resultvis_$Description[gsea_resultvis_$group == "blue"],
          set3 = gsea_resultchrom_$Description[gsea_resultchrom_$group == "blue"], 
          name = "blue_paths_top8.png")
plot_venn(set1 = gsea_resultgeo_$Description[gsea_resultgeo_$group == "pink"], 
          set2 = gsea_resultvis_$Description[gsea_resultvis_$group == "pink"],
          set3 = gsea_resultchrom_$Description[gsea_resultchrom_$group == "pink"], 
          name = "pink_paths_top8.png")


# enrichment barplots -----------------------------------------------------
library(stringr)
plot_top8_pathways <- function(plt_df, collow, colhigh, dir, name){
  plt_df_ <- plt_df %>%
    mutate(Description = str_to_title(gsub("HALLMARK", "", gsub("_", " ", Description)))) %>%
    arrange(desc(-log10(pvalue)))
  plt_df_$Description = factor(plt_df_$Description, levels = rev(plt_df_$Description))
  
  p <- ggplot(data = plt_df_, aes(x=Description, y=-log10(pvalue), fill=-log10(pvalue))) +
    geom_bar(stat="identity") + 
    ggtitle("Enrichment Analysis Results") + 
    labs(fill = expression("-log"[10]*"p-value")) + 
    xlab("Enriched Terms") + 
    ylab(expression("-log"[10]*"p-value")) + coord_flip() + 
    scale_fill_gradient(low = collow, high = colhigh) + 
    theme_minimal() + 
    theme(panel.grid = element_blank(), 
          plot.title = element_blank(), 
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.position = "bottom"
    ) # + ylim(0,15.5)

  ggsave(filename = file.path(dir, name), plot = p, bg = "white", width = 30, height = 18, units = "cm")
  return(p)
}

plot_top8_pathways(gsea_resultgeo_ %>% filter(group == "pink"), collow = "#fc98fc", colhigh = "#ff00ff", dir, "geo_pink_paths.png")
plot_top8_pathways(gsea_resultgeo_ %>% filter(group == "blue"), collow = "#9393ff", colhigh = "#0000ff", dir, "geo_blue_paths.png")

plot_top8_pathways(gsea_resultvis_ %>% filter(group == "pink"), collow = "#fc98fc", colhigh = "#ff00ff", dir, "vis_pink_paths.png")
plot_top8_pathways(gsea_resultvis_ %>% filter(group == "blue"), collow = "#9393ff", colhigh = "#0000ff", dir, "vis_blue_paths.png")

plot_top8_pathways(gsea_resultchrom_ %>% filter(group == "pink"), collow = "#fc98fc", colhigh = "#ff00ff", dir, "chrom_pink_paths.png")
plot_top8_pathways(gsea_resultchrom_ %>% filter(group == "blue"), collow = "#9393ff", colhigh = "#0000ff", dir, "chrom_blue_paths.png")

