library(SpatialExperiment)
# Recreate Fig 1b ---------------------------------------------------------
dv200_blockage <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200_BlockAgeMonth.csv")
dv200_blockage$Histology <- factor(dv200_blockage$Histology, 
                                   levels=c("Breast ductal carcinoma", "Breast lobular carcinoma", "Breast mucinous carcinoma", 
                                            "Lung adeno-carcinoma", "Lung squamous carcinoma", "DLBCL"))
dv200_blockage$Sampling_Type <- factor(dv200_blockage$Sampling_Type, 
                                   levels=c("Resection", "Biopsy"))

write.csv(dv200_blockage, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig1b.csv")
library(ggplot2)
ggplot(dv200_blockage, aes(x=Block_age_month, y=DV200, shape=Sampling_Type, color=Histology, label=IDCHUV)) + 
  geom_point(size=6) +
  ggtitle("Block characteristics") +
  theme(legend.position="bottom") +
  ggrepel::geom_text_repel(size=8, color="black") +
  # scale_fill_manual(values=c("DLBCL"="#FFD700",
  #                   "Lung squamous carcinoma"="#FFA44F",
  #                   "Lung adeno-carcinoma"="#B37337",
  #                   "Breast lobular carcinoma"="#00A674",
  #                   "Breast mucinous carcinoma"="#004D36",
  #                   "Breast ductal carcinoma"="#00FFB3")) +
  theme_minimal()


# DV200 vs ngene_detected -------------------------------------------------
library(stringr)
library(patchwork)

###########################################################################
## GeoMx
read_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/archive/GeoMx_test_norm/"

plt_df_per_indication <- function(disease){
  spe <- readRDS(file.path(read_path, paste0(disease, ".rds")))
  
  if(disease == "dlbcl"){spe$section_id <- spe$patient}
  plt_df <- data.frame(section = spe$section_id, 
                       loglibsize = spe$log_lib_size,
                       # ngenedetected = spe$gene_detection_rate,
                       ngenedetected = spe$genes_detected, 
                       indication = ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease)),
                       cell_fraction = spe$cell_fraction)
  return(plt_df)
}
df1 <- plt_df_per_indication("breast")
df2 <- plt_df_per_indication("lung")
df3 <- plt_df_per_indication("dlbcl")

df <- rbind(df1, df2, df3)

library(dplyr)
df_ <- df %>%
  mutate(IDCHUV=substr(section, 1,2)) %>%
  left_join(dv200_blockage)

library(patchwork)
df_$IDCHUV <- factor(df_$IDCHUV, levels=c("B1", "B2", "B3", "B4", 
                                          "L1", "L2", "L3", "L4", 
                                          "D1", "D2", "D3", "D4", "D5", "D6"))
df_$cell_fraction <- ifelse(df_$cell_fraction == "Macro", "Macrophage", df_$cell_fraction)
df_$cell_fraction <- as.factor(df_$cell_fraction)

names(df_)[names(df_) == "cell_fraction"] <- "AOI label"

minand1 <- function(x){
  sort(x)[5] + 1000
}


p <- ggplot(df_, aes(x=DV200, y=ngenedetected)) + 
  geom_jitter(aes(shape=Histology, color=`AOI label`), size=2, width = 0.2, alpha = 0.7) + 
  geom_boxplot(aes(fill=IDCHUV), alpha = 0.5, width=0.5, show.legend = FALSE, outlier.shape = NA) + 
  scale_color_manual(values = c("Macrophage" = "#9A32CD", "Malignant" = "#BC8F8F", "Other" = "#388E8E",
                                "PanCK-" = "#09D0EF", "T cells" = "#4169E1", "B cells" = "#FFD700")) +
  ylab("Number of detected genes per AOI") + 
  ggrepel::geom_text_repel(aes(label=IDCHUV), 
                           stat = "summary", 
                           fun = "minand1", 
                           size = 5, 
                           box.padding = 0.5) +  
  theme_minimal()

write.csv(df_, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFig1de.csv")
plot_width = 15
plot_height = 8
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200"
pdf(file = file.path(fig_path, "Geo_Ngenedetected_DV200.pdf"), 
    width = plot_width,
    height = plot_height)
print(p)
dev.off()


p <- ggplot(df_, aes(x=Block_age_month, y=ngenedetected)) + 
  geom_jitter(aes(shape=Histology, color=`AOI label`), size=2, width = 0.2, alpha = 0.7) + 
  geom_boxplot(aes(fill=IDCHUV), alpha = 0.5, width=0.5, show.legend = FALSE, outlier.shape = NA) + 
  scale_color_manual(values = c("Macrophage" = "#9A32CD", "Malignant" = "#BC8F8F", "Other" = "#388E8E",
                                "PanCK-" = "#09D0EF", "T cells" = "#4169E1", "B cells" = "#FFD700")) +
  ylab("Number of detected genes per AOI") + 
  xlab("Block age month") + 
  ggrepel::geom_text_repel(aes(label=IDCHUV), 
                           stat = "summary", 
                           fun = "minand1", 
                           size = 5, 
                           box.padding = 0.5) +  
  theme_minimal()


plot_width = 15
plot_height = 8
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200"
pdf(file = file.path(fig_path, "Geo_Ngenedetected_Block_age.pdf"), 
    width = plot_width,
    height = plot_height)
print(p)
dev.off()


## Loglibsize -------------------
minand1 <- function(x){
  sort(x)[5] + 1
}

p <- ggplot(df_, aes(x=DV200, y=loglibsize)) + 
  geom_jitter(aes(shape=Histology, color=`AOI label`), size=2, width = 0.2, alpha = 0.7) + 
  geom_boxplot(aes(fill=IDCHUV), alpha = 0.5, width=0.5, show.legend = FALSE, outlier.shape = NA) + 
  scale_color_manual(values = c("Macrophage" = "#9A32CD", "Malignant" = "#BC8F8F", "Other" = "#388E8E",
                                "PanCK-" = "#09D0EF", "T cells" = "#4169E1", "B cells" = "#FFD700")) +
  ggrepel::geom_text_repel(aes(label=IDCHUV), 
                           stat = "summary", 
                           fun = "minand1", 
                           size = 5, 
                           box.padding = 0.5) +  
  theme_minimal()


plot_width = 15
plot_height = 8
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200"
pdf(file = file.path(fig_path, "Geo_Loglibsize_DV200.pdf"), 
    width = plot_width,
    height = plot_height)
print(p)
dev.off()

p <- ggplot(df_, aes(x=Block_age_month, y=loglibsize)) + 
  geom_jitter(aes(shape=Histology, color=`AOI label`), size=2, width = 0.2, alpha = 0.7) + 
  geom_boxplot(aes(fill=IDCHUV), alpha = 0.5, width=0.5, show.legend = FALSE, outlier.shape = NA) + 
  scale_color_manual(values = c("Macrophage" = "#9A32CD", "Malignant" = "#BC8F8F", "Other" = "#388E8E",
                                "PanCK-" = "#09D0EF", "T cells" = "#4169E1", "B cells" = "#FFD700")) +
  ggrepel::geom_text_repel(aes(label=IDCHUV), 
                           stat = "summary", 
                           fun = "minand1", 
                           size = 5, 
                           box.padding = 0.5) +  
  theme_minimal()


plot_width = 15
plot_height = 8
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200"
pdf(file = file.path(fig_path, "Geo_Loglibsize_Block_age.pdf"), 
    width = plot_width,
    height = plot_height)
print(p)
dev.off()



########################################################################
## Visium
plt_df_per_indication <- function(disease){
  foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
  source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/Visium/01_params.R")
  save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
  plt_df <- NULL
  for(i in 1:nsamples){
    sce <- readRDS(file.path(save_bs_path, paste0(save_names[i], "_qcd.rds")))
    plt_df_i <- data.frame(section = save_names[i],
                           patho = sce$Region, 
                           loglibsize = log1p(sce$sum),
                           ngenedetected = colSums(counts(sce) > 0)
    )
    plt_df <- rbind(plt_df, plt_df_i)                     
  }
  plt_df$indication = foldername
  
  return(plt_df)
}

disease = "breast"
df1 <- plt_df_per_indication(disease = "breast")
disease = "lung"
df2 <- plt_df_per_indication(disease = "lung")
disease = "dlbcl"
df3 <- plt_df_per_indication(disease = "dlbcl")
df3$section <- paste0(substr(df3$section,1,1), substr(df3$section,7,7))

df <- rbind(df1, df2, df3)

df_ <- df %>%
  mutate(IDCHUV=substr(section, 1,2)) %>%
  left_join(dv200_blockage)

library(patchwork)
df_$IDCHUV <- factor(df_$IDCHUV, levels=c("B1", "B2", "B3", "B4", 
                                          "L1", "L2", "L3", "L4", 
                                          "D1", "D2", "D3", "D4", "D5", "D6"))
df_$patho <- ifelse(df_$patho == "Most_likely_Tumor", "Most_likely_tumor", df_$patho)
df_$patho <- as.factor(df_$patho)

names(df_)[names(df_) == "patho"] <- "Pathology label"

# minand1 <- function(x){
#   sort(x)[5] + 1000
# }

source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/color_palette.R")
plt_df <- df_
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

write.csv(df_, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFig1bc.csv")
p <- ggplot(df_, aes(x=DV200, y=ngenedetected)) + 
  geom_jitter(data = sampled_data, aes(shape=Histology), size=1, width = 0.2, alpha = 0.7) + # , color=`Pathology label`
  geom_boxplot(aes(fill=IDCHUV), alpha = 0.5, width=0.5, show.legend = FALSE, outlier.shape = NA) + 
  scale_color_manual(values = c(breast_patho_color, lung_patho_color, dlbcl_patho_color)) +
  ylab("Number of detected genes per spot") + 
  ggrepel::geom_text_repel(aes(label=IDCHUV), 
                           stat = "summary", 
                           fun = "median", 
                           size = 5, 
                           box.padding = 0.5) +  
  theme_minimal()


plot_width = 15
plot_height = 8
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200"
pdf(file = file.path(fig_path, "Vis_Ngenedetected_DV200.pdf"), 
    width = plot_width,
    height = plot_height)
print(p)
dev.off()


p <- ggplot(df_, aes(x=Block_age_month, y=ngenedetected)) + 
  geom_jitter(data = sampled_data, aes(shape=Histology), size=1, width = 0.2, alpha = 0.7) + # , color=`Pathology label`
  geom_boxplot(aes(fill=IDCHUV), alpha = 0.5, width=0.5, show.legend = FALSE, outlier.shape = NA) + 
  scale_color_manual(values = c(breast_patho_color, lung_patho_color, dlbcl_patho_color)) +
  ylab("Number of detected genes per spot") + 
  xlab("Block age month") + 
  ggrepel::geom_text_repel(aes(label=IDCHUV), 
                           stat = "summary", 
                           fun = "median", 
                           size = 5, 
                           box.padding = 0.5) +  
  theme_minimal()


plot_width = 15
plot_height = 8
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200"
pdf(file = file.path(fig_path, "Vis_Ngenedetected_Block_age.pdf"), 
    width = plot_width,
    height = plot_height)
print(p)
dev.off()


## Loglibsize -------------------
p <- ggplot(df_, aes(x=DV200, y=loglibsize)) + 
  geom_jitter(data = sampled_data, aes(shape=Histology, color=`Pathology label`), size=1, width = 0.2, alpha = 0.7) + 
  geom_boxplot(aes(fill=IDCHUV), alpha = 0.5, width=0.5, show.legend = FALSE, outlier.shape = NA) + 
  scale_color_manual(values = c(breast_patho_color, lung_patho_color, dlbcl_patho_color)) +
  ggrepel::geom_text_repel(aes(label=IDCHUV), 
                           stat = "summary", 
                           fun = "median", 
                           size = 5, 
                           box.padding = 0.5) +  
  theme_minimal()


plot_width = 15
plot_height = 8
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200"
pdf(file = file.path(fig_path, "Vis_Loglibsize_DV200.pdf"), 
    width = plot_width,
    height = plot_height)
print(p)
dev.off()


p <- ggplot(df_, aes(x=Block_age_month, y=loglibsize)) + 
  geom_jitter(data = sampled_data, aes(shape=Histology), size=2, width = 0.2, alpha = 0.7) +  # , color=`Pathology label`
  geom_boxplot(aes(fill=IDCHUV), alpha = 0.5, width=0.5, show.legend = FALSE, outlier.shape = NA) + 
  scale_color_manual(values = c(breast_patho_color, lung_patho_color, dlbcl_patho_color)) +
  ggrepel::geom_text_repel(aes(label=IDCHUV), 
                           stat = "summary", 
                           fun = "median", 
                           size = 5, 
                           box.padding = 0.5) +  
  theme_minimal()


plot_width = 15
plot_height = 8
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/DV200"
pdf(file = file.path(fig_path, "Vis_Loglibsize_Block_age.pdf"), 
    width = plot_width,
    height = plot_height)
print(p)
dev.off()


# (ggplot(df_, aes(x=DV200, y=ngenedetected, shape=Histology, color=IDCHUV)) + 
#     geom_point(size=6) + 
#     theme_minimal()) |
#   (ggplot(df_, aes(x=Block_age_month, y=ngenedetected, shape=Histology, color=IDCHUV)) + 
#      geom_point(size=6) + 
#      theme_minimal())


