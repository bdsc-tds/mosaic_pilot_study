library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggridges)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig6/"
figpath_ridge <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Fig6_Ridge"

# Dot plot ----------------------------------------------------------------
datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched"
geo_spe <- readRDS(file.path(datapath, "dlbcl_spe_ruv.rds"))
geo <- as.Seurat(geo_spe)

geo <- readRDS("dlbcl_seu_ruv.rds")
# head(Idents(geo)) # Levels: T cells Macrophage Other Tu_D1 Tu_D2 Tu_D3 Tu_D4 Tu_D5 Tu_D6

geo <- geo[rownames(geo) %in% DLBCLnChromium_Marker_Gene_List, ]

DLBCLnChromium_Marker_Gene_List <- c(
  "MS4A1", "TNFRSF13C", "CD79B", "CD37", "PSMB8", "CD19", "TYMS",
  "TUBB", "TOP2A", "POLD4", "CD47", "CD52", "BLK", "CD38", "MAP2K1",
  "CD40", "BCL2L1", "TNFRSF8", "SMO", "RARA", "TYK2", "TNFRSF10B"
)

p <- DotPlot(
  geo,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none")

plot_title = "Geo_DLBLCnChromium_dotplot_by_AOI.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 2,
    height = 4)
print(p)
dev.off()

##########################################################################
# Introducing GeoMx
##########################################################################
geo_decon_result <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results/dlbcl_batched_decon_long.csv") # not patient specific
# geo_decon_result <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results_pt_specific/dlbcl_batched_decon_long.csv") # patient specific


geo_decon_result_ <- geo_decon_result %>%
  # mutate(Cell_Type = paste0(substr(CellType, 1, 2), "_", patient),
  #        CellType = ifelse(CellType == "Tumor", Cell_Type, CellType)) %>%
  spread(CellType, Fraction) %>%
  mutate(Other = B + Epithelia + Fibro_muscle + Vessel,
         `T cells` = T_NK,
         Macrophage = Myeloid) %>%
  select(sample, cell_fraction, Other, `T cells`, Macrophage, Tumor, patient) 

df <- geo_decon_result_ %>% 
  select(Other, `T cells`, Macrophage, Tumor)

decon_max_val <- apply(df, 1, max)
decon_max_val_indices <- max.col(df, ties.method = "first")
decon_max_val_names <- colnames(df)[decon_max_val_indices]

geo_decon_result_$decon_max_val <- decon_max_val
geo_decon_result_$decon_max_val_names <- decon_max_val_names


# If set threshold and find consensus -------------------------------------
geo_decon_result_$decon_max_val_names_threshold <- ifelse(geo_decon_result_$decon_max_val >= 0.5, geo_decon_result_$decon_max_val_names, "Mix")
#> table(geo_decon_result_$decon_max_val_names_threshold)
# Macrophage        Mix      Other    T cells      Tumor 
#        13         79         12         13         19 

table(geo_decon_result_$decon_max_val_names_threshold, geo_decon_result_$patient)
table(geo_decon_result_$decon_max_val_names_threshold, geo_decon_result_$cell_fraction)

geo_decon_result_ %>%
  filter(decon_max_val_names_threshold == "Tumor") %>%
  group_by(patient) %>%
  summarise(n = n())

geo_decon_result_ %>%
  filter(decon_max_val_names_threshold == "Tumor" & cell_fraction == "B cells") %>%
  group_by(patient) %>%
  summarise(n = n())

geo_decon_result_$decon_max_val_names_threshold2 <- ifelse(geo_decon_result_$decon_max_val > 0.7, geo_decon_result_$decon_max_val_names, "Mix")
#> table(geo_decon_result_$decon_max_val_names_threshold2)
#> Macrophage        Mix      Other    T cells 
# 6        120          9          1 


# If do not set threshold and just find consensus -------------------------
table(geo_decon_result_$decon_max_val_names, geo_decon_result_$patient)
table(geo_decon_result_$decon_max_val_names, geo_decon_result_$cell_fraction)

geo_decon_result_ %>%
  filter(decon_max_val_names == "Tumor") %>%
  group_by(patient) %>%
  summarise(n = n())

geo_decon_result_ %>%
  filter(decon_max_val_names == "Tumor" & cell_fraction == "B cells") %>%
  group_by(patient) %>%
  summarise(n = n())

#   patient     n
# <chr>   <int>
# 1 D1          8
# 2 D2          9
# 3 D3          4
# 4 D4          4
# 5 D5          9
# 6 D6          8

# -------------------------------------------------------------------------
geo_decon_result_1 <- geo_decon_result_ %>%
  gather(key = CellType, value = Fraction, -c(sample, cell_fraction, patient, decon_max_val, decon_max_val_names, decon_max_val_names_threshold, decon_max_val_names_threshold2)) %>%
  mutate(cell_fraction = ifelse(cell_fraction == "B cells", "Malignant",
                                ifelse(cell_fraction == "Macro", "Macrophage", cell_fraction)),
         cell_fraction_pt = ifelse(cell_fraction == "Malignant", paste0("Tu_", patient), cell_fraction),
         CellType_pt = ifelse(CellType == "Tumor", paste0(substr(CellType, 1, 2), "_", patient), CellType),
         decon_max_val_names_pt = ifelse(decon_max_val_names == "Tumor", paste0(substr(decon_max_val_names, 1, 2), "_", patient), decon_max_val_names),
         decon_max_val_names_consensus = case_when(decon_max_val_names == "Tumor" & cell_fraction == "Malignant" ~ "Tumor",
                                                   decon_max_val_names == "T cells" & cell_fraction == "T cells" ~ "T cells",
                                                   decon_max_val_names == "Other" & cell_fraction == "Other" ~ "Other",
                                                   decon_max_val_names == "Macrophage" & cell_fraction == "Macrophage" ~ "Macrophage",
                                                   .default = "Mix"),
         decon_max_val_names_consensus_pt = ifelse(decon_max_val_names_consensus == "Tumor", paste0(substr(decon_max_val_names_consensus, 1, 2), "_", patient), decon_max_val_names_consensus),
         decon_max_val_names_threshold_pt = ifelse(decon_max_val_names_threshold == "Tumor", paste0(substr(decon_max_val_names_threshold, 1, 2), "_", patient), decon_max_val_names_threshold)
         )

palette_ridge <- c("#388E8E", "#9A32CD", "#4169E1", "#FF8C00", "#EEEE00", "#FFD700", "#A2CD5A", "#00EE76", "#ADFF2F")
ridge_order <- c("Other", "Macrophage", "T cells", "Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6")
names(palette_ridge) <- ridge_order

# Color by AOI -------------------------------------------------------------
test_ <- geo_decon_result_1 %>%  # Keep expected decon proportion in each consensus region only
  filter((cell_fraction_pt == "Tu_D1" & CellType == "Tumor") |
           (cell_fraction_pt == "Tu_D2" & CellType == "Tumor") |
           (cell_fraction_pt == "Tu_D3" & CellType == "Tumor") |
           (cell_fraction_pt == "Tu_D4" & CellType == "Tumor") |
           (cell_fraction_pt == "Tu_D5" & CellType == "Tumor") |
           (cell_fraction_pt == "Tu_D6" & CellType == "Tumor") |
           (cell_fraction_pt == "Other" & CellType == "Other") |
           (cell_fraction_pt == "Macrophage" & CellType == "Macrophage") |
           (cell_fraction_pt == "T cells" & CellType == "T cells"))

table(test_$cell_fraction_pt)
# Macrophage      Other    T cells      Tu_D1      Tu_D2      Tu_D3      Tu_D4      Tu_D5      Tu_D6 
#         19          9         52          8         12          8          7         13          8 

# # GeoMx by AOI ------------------------------------------------------
test_$cell_fraction_pt <- factor(test_$cell_fraction_pt, levels = rev(ridge_order))
p <- ggplot(test_, aes(x = Fraction, y = cell_fraction_pt, fill = cell_fraction_pt)) +
  geom_density_ridges(scale = 4, # stat = "binline", binwidth=0.01, draw_baseline = F
                      alpha = 0.5) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  theme(legend.position = "none") + 
  xlim(c(0, 1)) + 
  scale_fill_manual(values = palette_ridge) + 
  ylab("AOI label")

sample_size_df <- data.frame(table(test_$cell_fraction_pt)) %>%
  dplyr::rename(n = Freq,
         ylabel = Var1) %>%
  mutate(n = paste0("n = ", n))

p <- p + geom_text(data=sample_size_df,
              aes(label = n, x= 1, y=ylabel, vjust = -2, hjust=-0.1),
              position = position_stack(),
              inherit.aes = FALSE)

plot_title = "Geo_ridge_AOI.pdf"
pdf(file = file.path(figpath_ridge, plot_title),
    width = 9,
    height = 6)
print(p)
dev.off()

############################################################################
# Color by consensus (to improve purity) -----------------------------------
test <- geo_decon_result_1 %>%  # Keep expected decon proportion in each consensus region only
  filter((decon_max_val_names_consensus_pt == "Tu_D1" & CellType == "Tumor") |
           (decon_max_val_names_consensus_pt == "Tu_D2" & CellType == "Tumor") |
           (decon_max_val_names_consensus_pt == "Tu_D3" & CellType == "Tumor") |
           (decon_max_val_names_consensus_pt == "Tu_D4" & CellType == "Tumor") |
           (decon_max_val_names_consensus_pt == "Tu_D5" & CellType == "Tumor") |
           (decon_max_val_names_consensus_pt == "Tu_D6" & CellType == "Tumor") |
           (decon_max_val_names_consensus_pt == "Other" & CellType == "Other") |
           (decon_max_val_names_consensus_pt == "Macrophage" & CellType == "Macrophage") |
           (decon_max_val_names_consensus_pt == "T cells" & CellType == "T cells"))

table(test$decon_max_val_names_consensus_pt)
# Macrophage      Other    T cells      Tu_D1      Tu_D2      Tu_D3      Tu_D4      Tu_D5      Tu_D6 
# 16          9         30          8          9          4          4          9          8 

# GeoMx by consensus ------------------------------------------------------
test$decon_max_val_names_consensus_pt <- factor(test$decon_max_val_names_consensus_pt, levels = rev(ridge_order))
p <- ggplot(test, aes(x = Fraction, y = decon_max_val_names_consensus_pt, fill = decon_max_val_names_consensus_pt)) +
  geom_density_ridges(scale = 4, # , stat = "binline", binwidth=0.01, draw_baseline = F
                      alpha=0.5) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  theme(legend.position = "none") + 
  xlim(c(0, 1)) + 
  scale_fill_manual(values = palette_ridge) + 
  ylab("Consensus label by AOI & Deconvolution majority vote")

sample_size_df <- data.frame(table(test$decon_max_val_names_consensus_pt)) %>%
  dplyr::rename(n = Freq,
         ylabel = Var1) %>%
  mutate(n = paste0("n = ", n))

p <- p + geom_text(data=sample_size_df,
              aes(label = n, x= 1, y=ylabel, vjust = -2, hjust=-0.1),
              position = position_stack(),
              inherit.aes = FALSE)

plot_title = "Geo_ridge_consensus.pdf"
pdf(file = file.path(figpath_ridge, plot_title),
           width = 9,
           height = 6)
print(p)
dev.off()


geo$consensus <- geo$sample_id2 %in% test$sample
saveRDS(geo, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/dlbcl_seu_ruv.rds")

#################################################################
DLBCLnChromium_Marker_Gene_List <- c(
  "MS4A1", "TNFRSF13C", "CD79B", "CD37", "PSMB8", "CD19", "TYMS",
  "TUBB", "TOP2A", "POLD4", "CD47", "CD52", "BLK", "CD38", "MAP2K1",
  "CD40", "BCL2L1", "TNFRSF8", "SMO", "RARA", "TYK2", "TNFRSF10B"
)

geo <- geo[rownames(geo) %in% DLBCLnChromium_Marker_Gene_List, ]

geo_small <- geo[, geo$sample_id2 %in% test$sample]

p <- DotPlot(
  geo_small,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5))

plot_title = "Geo_DLBLCnChromium_dotplot.pdf"
pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig6_drug_target_dlbcl_geo/", plot_title),
    width = 5,
    height = 20)
print(p)
dev.off()