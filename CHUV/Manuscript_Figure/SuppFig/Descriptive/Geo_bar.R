library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")


# helper ------------------------------------------------------------------
plot_bar_each_indication <- function(plt_df){
  bar_width <- 0.8
  p <- ggplot(data = plt_df, aes(x = cell_fraction, y = count, fill = cell_fraction)) +
    geom_bar(stat = "identity", width = bar_width) +
    geom_text(aes(x = as.numeric(cell_fraction), y = count, label = paste0(count, " (", round(per, 2)* 100, "%)")),
              vjust = -0.5,
              color = "black",
              size = 3.5) +
    scale_fill_manual(values = cell_fraction_color) + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5, hjust=0.5),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(), 
          legend.position =  "none",
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 14)) + 
    labs(y = "Counts", x = "Cell fraction") + 
    ggtitle(str_to_title(disease))
  
  return(p)
}


# -------------------------------------------------------------------------
disease_list <- c("breast", "lung", "dlbcl")
# disease = "breast"

for(disease in disease_list){
  ## Geo
  datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
  spe_ruv <- readRDS(file.path(datapath, paste0(disease, "_spe_ruv.rds")))

  plt_df <- as.data.frame(table(spe_ruv$cell_fraction)) %>%
    dplyr::rename(cell_fraction = Var1,
                  count = Freq) %>%
    mutate(cell_fraction = as.character(cell_fraction),
           cell_fraction = ifelse(cell_fraction == "Macro", "Macrophage", cell_fraction),
           per = count/sum(count)) 
  
  if(disease %in% c("breast", "lung")){
    plt_df$cell_fraction <- factor(plt_df$cell_fraction, levels = c("Malignant", "Other", "T cells", "Macrophage", "PanCK-"))
    cell_fraction_order = c("Macrophage", "Malignant", "Other", "T cells", "PanCK-") # note with PanCK-
    cell_fraction_color = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1", "#09D0EF")
    names(cell_fraction_color) <- cell_fraction_order
  }else{
    plt_df$cell_fraction <- factor(plt_df$cell_fraction, levels = c("B cells", "Other", "T cells", "Macrophage"))
    cell_fraction_order = c("B cells", "Macrophage", "Other", "T cells") 
    cell_fraction_color = c("#FFD700", "#9A32CD", "#388E8E", "#4169E1")
    names(cell_fraction_color) <- cell_fraction_order
  }

  assign(paste0("plt_df_", disease), plt_df)
  
  # p <- plot_bar_each_indication(plt_df)
  # assign(paste0("p_", disease), p)
}

plt_df_breast$indication = "Breast" 
plt_df_lung$indication = "Lung"
plt_df_dlbcl$indication = "DLBCL"

write.csv(rbind(plt_df_breast, 
                plt_df_lung,
                plt_df_dlbcl),
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS2e.csv")     


p_breast | p_lung | p_dlbcl


fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Descriptive/Geo_bar/"
pdf(file = file.path(fig_path, "Geo_breast_bar.pdf"),
    width = 4,
    height = 6.2)
print(p_breast)
dev.off()

pdf(file = file.path(fig_path, "Geo_lung_bar.pdf"),
    width = 4,
    height = 6)
print(p_lung)
dev.off()

pdf(file = file.path(fig_path, "Geo_dlbcl_bar.pdf"),
    width = 3.2,
    height = 7.5)
print(p_dlbcl)
dev.off()


# Pie chart ---------------------------------------------------------------
# Load ggplot2
library(ggplot2)
library(dplyr)

# Create Data
data <- data.frame(
  group=LETTERS[1:5],
  value=c(13,7,9,21,2)
)

# Compute the position of labels
data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  
  geom_text(aes(y = ypos, label = group), color = "white", size=6) +
  scale_fill_brewer(palette="Set1")