library(stringr)
library(patchwork)

scaleFUN <- function(x) sprintf("%.0f", x) # have no decimal on x-axis

## GeoMx
read_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/archive/GeoMx_test_norm/"

plt_df_per_indication <- function(disease){
  spe <- readRDS(file.path(read_path, paste0(disease, ".rds")))
  
  if(disease == "dlbcl"){spe$section_id <- spe$patient}
  plt_df <- data.frame(section = spe$section_id, 
                       # loglibsize = spe$log_lib_size,
                       # ngenedetected = spe$gene_detection_rate,
                       ngenedetected = spe$genes_detected, 
                       indication = ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease)),
                       cell_fraction = spe$cell_fraction)
  return(plt_df)
}

plt_df_geo_breast <- plt_df_per_indication(disease = "breast")
plt_df_geo_lung <- plt_df_per_indication(disease = "lung")
plt_df_geo_dlbcl <- plt_df_per_indication(disease = "dlbcl")

write.csv(rbind(plt_df_geo_breast, plt_df_geo_lung, plt_df_geo_dlbcl), 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig1c_geo.csv")
plot_violin_per_indication_geo <- function(plt_df = plt_df_geo_dlbcl, 
                                           section_order = c("D1", "D2", "D3", "D4", "D5", "D6"), 
                                           cell_fraction_order = c("B cells", "Macro", "Other", "T cells"),
                                           color_palette_order = c("#FFD700", "#9A32CD", "#388E8E", "#4169E1"),
                                           indication_title = "DLBCL",
                                           legen.pos = c(0.75, 0.835),
                                           legen.nrow = 2,
                                           legen.text.size = 5,
                                           legen.title.size = 5.5,
                                           legen.key.size = 0.5, 
                                           strip_color = "#00c78cff",
                                           strip_fill = "#00c78c2e"){
  
  plt_df$section <- factor(plt_df$section, levels = section_order)
  plt_df$cell_fraction <- factor(plt_df$cell_fraction, cell_fraction_order)
  
  plt_df <- plt_df %>%
    filter(cell_fraction != "PanCK-") 
  
  p <- ggplot(plt_df, aes(x = section, y = ngenedetected)) + 
    geom_violin(trim=FALSE, alpha = 0.5) + 
    geom_jitter(aes(color = cell_fraction), size=1.5, width = 0.2) + 
    scale_colour_manual(values = color_palette_order) +
    theme_bw() +
    xlab("") + ylab("")
  
  p <- p + 
    theme(legend.position = legen.pos, 
          legend.key = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = legen.text.size),
          legend.title = element_text(size = legen.title.size, face = "bold"),
          legend.key.size = unit(legen.key.size, "lines"),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(face="bold", size=13.5), 
          strip.background=element_rect(color = strip_color, 
                                        fill = strip_fill)) +
    guides(color = guide_legend(nrow = legen.nrow, byrow = TRUE)) +
    labs(color = "Cell fraction") +
    facet_wrap(~ indication)
  
  p <- p + scale_y_continuous(labels=scaleFUN)
  p
}


geo_plot <- 
  (plot_violin_per_indication_geo(plt_df = plt_df_geo_breast, 
                                 section_order = c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1"), 
                                 cell_fraction_order = c("Macro", "Malignant", "Other", "T cells"),
                                 color_palette_order = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1"),
                                 indication_title = "Breast",
                                 legen.pos = c(0.75, 0.775),
                                 legen.nrow = 2,
                                 legen.text.size = 5,
                                 legen.title.size = 5.5,
                                 legen.key.size = 0.5, 
                                 strip_color = "#00c78cff",
                                 strip_fill = "#00c78c2e")) / # + ylim(0, 1))
  plot_spacer() /
  (plot_violin_per_indication_geo(plt_df = plt_df_geo_lung, 
                                 section_order = c("L1_1", "L2_1", "L3_1", "L3_3", "L4_3"), 
                                 cell_fraction_order = c("Macro", "Malignant", "Other", "T cells"),
                                 color_palette_order = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1"),
                                 indication_title = "Lung",
                                 legen.pos = c(0.75, 0.775),
                                 legen.nrow = 2,
                                 legen.text.size = 5,
                                 legen.title.size = 5.5,
                                 legen.key.size = 0.5, 
                                 strip_color = "#ffa54fff",
                                 strip_fill = "#ffa54f31")) / #  + ylim(0, 1))
  plot_spacer() /
  (plot_violin_per_indication_geo(plt_df = plt_df_geo_dlbcl, 
                                 section_order = c("D1", "D2", "D3", "D4", "D5", "D6"), 
                                 cell_fraction_order = c("B cells", "Macro", "Other", "T cells"),
                                 color_palette_order = c("#FFD700", "#9A32CD", "#388E8E", "#4169E1"),
                                 indication_title = "DLBCL",
                                 legen.pos = c(0.75, 0.885), # c(0.25, 0.15),
                                 legen.nrow = 2, # 4,
                                 legen.text.size = 5, # 3,
                                 legen.title.size = 5.5, # 3.5,
                                 legen.key.size = 0.4, 
                                 strip_color = "#ffd700ff",
                                 strip_fill = "#ffd70033") + theme(legend.margin = margin(0, 0, 0, 0))) + # + ylim(0, 1)
  plot_layout(heights = c(4, -0.95, 4, -0.95, 4))
geo_plot

pdf(file = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig1/Fig1c/geo_ngenedetected_violin_bold_colorfacet_notGDR.pdf",
      # "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1c/geo_ngenedetected_violin_bold.pdf",
      # "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1c/geo_ngenedetected_violin.pdf",
    width = 3.3,
    height = 7)
print(geo_plot)
dev.off()



## Visium ------------------------------------------------------------------------------------------------------------------------
plt_df_per_indication <- function(disease){
  foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
  source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/Visium/01_params.R")
  save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
  plt_df <- NULL
  for(i in 1:nsamples){
    sce <- readRDS(file.path(save_bs_path, paste0(save_names[i], "_qcd.rds")))
    plt_df_i <- data.frame(section = save_names[i],
                           # loglibsize = log1p(sce$sum),
                           ngenedetected = colSums(counts(sce) > 0)
                           )
    plt_df <- rbind(plt_df, plt_df_i)                     
  }
  plt_df$indication = foldername
  
  return(plt_df)
}

plt_df_vis_breast <- plt_df_per_indication(disease = "breast")
plt_df_vis_lung <- plt_df_per_indication(disease = "lung")
plt_df_vis_dlbcl <- plt_df_per_indication(disease = "dlbcl")
plt_df_vis_dlbcl$section <- paste0("D", substr(plt_df_vis_dlbcl$section, 7, 7))

write.csv(rbind(plt_df_vis_breast, plt_df_vis_lung, plt_df_vis_dlbcl), 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig1c_vis.csv")

plot_violin_per_indication <- function(plt_df, 
                                       section_order, 
                                       indication_title = "Breast", 
                                       strip_color = "#00c78cff",
                                       strip_fill = "#00c78c2e"){
  
  n_jitters <- 150 # Each sample contributes 100 jitters
  # Subsample a portion of the data for jittering
  if(any(table(plt_df$section) < n_jitters)){
    plt_df_smallclass <- plt_df %>% filter(section == names(which(table(plt_df$section) < n_jitters))) # small sample all in
    plt_df_otherclass <- plt_df %>% filter(section != names(which(table(plt_df$section) < n_jitters)))
    
    sampled_data <- plt_df_otherclass %>%  
      group_by(section) %>%
      sample_n(size = n_jitters, replace = FALSE)
    
    sampled_data <- rbind(sampled_data, plt_df_smallclass) %>% arrange(section)
    
  }else{
    sampled_data <- plt_df %>%  
      group_by(section) %>%
      sample_n(size = n_jitters, replace = FALSE) %>% arrange(section)
  }
  
  sampled_data$section <- factor(sampled_data$section, levels = section_order)
  plt_df$section <- factor(plt_df$section, levels = section_order)
  
  # Basic violin plot
  p <- ggplot(plt_df, aes(x = section, y = ngenedetected)) + 
    geom_jitter(data = sampled_data, size=0.05, alpha=0.8, width = 0.2) + 
    geom_violin(trim=FALSE, alpha = 0.5) + 
    theme_bw() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 10),
          panel.grid = element_blank(),
          strip.text = element_text(face="bold", size=13.5), 
          strip.background=element_rect(color = strip_color, 
                                        fill = strip_fill)) + 
    xlab("") + ylab("") + 
    facet_wrap(~ indication)

  p
}

vis_plot <- (plot_violin_per_indication(plt_df_vis_breast, c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2"), "Breast", 
                                        strip_color = "#00c78cff",
                                        strip_fill = "#00c78c2e") + ylim(0, 11500)) /
  plot_spacer() /
  (plot_violin_per_indication(plt_df_vis_lung, c("L1_2", "L1_4", "L2_2", "L3_2", "L4_2"), "Lung", 
                              strip_color = "#ffa54fff",
                              strip_fill = "#ffa54f31") + ylim(0, 11500)) /
  plot_spacer() /
  (plot_violin_per_indication(plt_df_vis_dlbcl, c("D1", "D2", "D3", "D4", "D5", "D6"), "DLBCL", 
                              strip_color = "#ffd700ff",
                              strip_fill = "#ffd70033") + ylim(0, 11500)) + 
  plot_layout(guides = "collect", heights = c(4, -0.95, 4, -0.95, 4))
vis_plot

pdf(file = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig1/Fig1c/vis_ngenedetected_violin_bold_colorfacet.pdf",
    # "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1c/vis_ngenedetected_violin_bold.pdf",
    # "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1c/vis_ngenedetected_violin.pdf",
    width = 3.3,
    height = 7)
print(vis_plot)
dev.off()


## Chromium ------------------------------------------------------------------------------------------------------------------------
chrom_breast <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_breast.rds")
chrom_lung <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_lung.rds")
chrom_dlbcl <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/DLBCL/dlbcl_final_owkin_annot_pub.rds")
chrom_dlbcl$sample_id <- paste0("D", substr(chrom_dlbcl$sample_id, 7, 7))

plot_violin_per_indication <- function(obj = chrom_breast, 
                                       section_order = c("B1", "B2", "B3", "B3_rep", "B4", "B4_rep"), 
                                       indication_title = "Breast", 
                                       strip_color = "#00c78cff",
                                       strip_fill = "#00c78c2e"){
  
  plt_df <- data.frame(section = obj$sample_id,
                       # loglibsize = log1p(obj$nCount_RNA),
                       ngenedetected = obj$nFeature_RNA,
                       indication = rep(indication_title, ncol(obj)))
  
  n_jitters <- 150 # Each sample contributes 100 jitters
  # Subsample a portion of the data for jittering
  if(any(table(plt_df$section) < n_jitters)){
    plt_df_smallclass <- plt_df %>% filter(section == names(which(table(plt_df$section) < n_jitters))) # small sample all in
    plt_df_otherclass <- plt_df %>% filter(section != names(which(table(plt_df$section) < n_jitters)))
    
    sampled_data <- plt_df_otherclass %>%  
      group_by(section) %>%
      sample_n(size = n_jitters, replace = FALSE)
    
    sampled_data <- rbind(sampled_data, plt_df_smallclass) %>% arrange(section)
    
  }else{
    sampled_data <- plt_df %>%  
      group_by(section) %>%
      sample_n(size = n_jitters, replace = FALSE) %>% arrange(section)
  }
  
  sampled_data$section <- factor(sampled_data$section, levels = section_order)
  plt_df$section <- factor(plt_df$section, levels = section_order)
  
  # Basic violin plot
  p <- ggplot(plt_df, aes(x = section, y = ngenedetected)) + 
    geom_jitter(data = sampled_data, size=0.0001, alpha=0.8, width = 0.2) + 
    # geom_jitter(size=0.001, alpha=0.8, width = 0.2) + 
    geom_violin(trim=FALSE, alpha = 0.5) + 
    theme_bw() +
    theme(axis.text.y = element_text(size = 6), 
          axis.text.x = element_text(size = 10),
          panel.grid = element_blank(),
          strip.text = element_text(face="bold", size=13.5), 
          strip.background=element_rect(color = strip_color, 
                                        fill = strip_fill)) + 
    xlab("") + ylab("") + 
    facet_wrap(~ indication)

  p
}

write.csv(rbind(plt_df_breast, plt_df_lung, plt_df_dlbcl), 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig1c_chrom.csv")

chrom_plot <- (plot_violin_per_indication(chrom_breast, c("B1", "B2", "B3", "B3_rep", "B4", "B4_rep"), "Breast", 
                                          strip_color = "#00c78cff",
                                          strip_fill = "#00c78c2e") + ylim(0, 11500)) /
  plot_spacer() /
  (plot_violin_per_indication(chrom_lung, c("L1", "L2", "L3", "L4"), "Lung", 
                              strip_color = "#ffa54fff",
                              strip_fill = "#ffa54f31") + ylim(0, 11500)) /
  plot_spacer() /
  (plot_violin_per_indication(chrom_dlbcl, c("D1", "D2", "D3", "D4", "D5", "D6"), "DLBCL", 
                              strip_color = "#ffd700ff",
                              strip_fill = "#ffd70033") + ylim(0, 11500)) + 
  plot_layout(guides = "collect", heights = c(4, -0.95, 4, -0.95, 4))
chrom_plot

pdf(file = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig1/Fig1c/chrom_ngenedetected_violin_bold_colorfacet.pdf", 
    # "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1c/chrom_ngenedetected_violin.pdf",
    width = 3.3,
    height = 7)
print(chrom_plot)
dev.off()


