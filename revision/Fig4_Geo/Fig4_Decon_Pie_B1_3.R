library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(tidyr)
library(patchwork)

## Level 4 GeoMx
geo_level4breast <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/breast_batched_decon_long.csv")

# Color
disease = "breast"
source("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/env/ydong/mosaic_pilot_study/CHUV/color_palette.R")
# colors = level4_cellcolors[names(level4_cellcolors) %in% sort(geo_level4breast_B1_Bcell_ROI$CellType)]

# Plot function
plot_AOI_pie <- function(roi_df, cf_type){
  p <- ggplot(roi_df[roi_df$cell_fraction == cf_type, ], 
              aes(fill=CellType, y=Fraction, x="")) + 
    geom_bar(position="fill", stat="identity") +
    # scale_fill_manual(values = colors) +
    coord_polar("y", start=0) + 
    theme_void() + 
    theme(legend.position="none")
  return(p)
}

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_pie_charts"


# Sanity ------------------------------------------------------------------
## QCed
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds") # 24

## Raw
geo_raw <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/GeoMx/breast_raw/B1_3_OPHI.rds") # 24
table(geo_raw$segment)
# CD3+ CD68+   CK+ Other 
# 7     3     6     8 
table(geo_raw$Cell_fraction)
# Macro Malignant     Other    PanCK-   T_cells 
# 3         6         4         4         7
table(geo_raw$roi)

# Stacked bar plot for all 24 AOIs
geo_level4breast_B1_3 <- geo_level4breast %>%
  filter(section_id == "B1_3")
length(unique(geo_level4breast_B1_3$sample)) # 24

geo_level4breast_B1_3_merged <- geo_level4breast_B1_3 %>%
  pivot_wider(values_from = Fraction, names_from = CellType) %>%
  rowwise() %>% 
  mutate(
    B_cells = B_plasma_IGHA1 + B_plasma_IGHG1 + B_plasma_IGHG3 + B_plasma_IGHM + B_plasma_IGKC + B_plasma_IGLC1,
    Granulocyte = Granulocyte + Mast_cell, 
    Myeloid = DC_1 + DC_2 + DC_activated + DC_pc + Macrophage + Monocyte,
    Stroma = Endothelia_vascular + Fibroblast_B3 + Muscle_smooth + Pericyte, 
    T_NK = T_CD4 + T_CD8_exhausted + T_CTL + T_CXCL13 + T_reg + TNK_dividing + NK,
    Tu_B1 = Tu_B1_MUCL1 + Tu_B1_MUCL1_necrosis + Tu_B1_MUCL1_transcription) %>%
  select(sample, section_id, cell_fraction, B_cells, Granulocyte, Myeloid, Stroma, T_NK, Tu_B1) %>%
  pivot_longer(cols = c(B_cells, Granulocyte, Myeloid, Stroma, T_NK, Tu_B1), 
               values_to = "Fraction", names_to = "CellType")


level1_5_cellnames <- c(    "B_cells", "Granulocyte", "Myeloid", "Stroma",  "T_NK",    "Tu_B1")
br_level1_5_cellcolors <- c("#FFDAB9", "#FFC125",     "#B452CD", "#388E8E", "#00BFFF", "#00EE76")
names(br_level1_5_cellcolors) <- level1_5_cellnames

geo_level4breast_B1_3_merged$sample_cf <- paste0(geo_level4breast_B1_3_merged$sample, "_", geo_level4breast_B1_3_merged$cell_fraction)

p <- ggplot(geo_level4breast_B1_3_merged, aes(fill=CellType, y=Fraction, x=sample_cf)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = br_level1_5_cellcolors) +
  coord_flip() + 
  theme_minimal()

write.csv(geo_level4breast_B1_3_merged, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig4bg_decon_vis.csv")

######################################################################
# B1_3 (roi 1) -------------------------------------------------------
geo_level4breast_B1_roi1 <- geo_level4breast %>%  # A_islet_2 (done)
  filter(sample %in% c(
                       "DSP-1001660018473-B-A03", # Malig
                       "DSP-1001660018473-B-A04")) %>% # Other
  arrange(cell_fraction) 

p_B1_3_roi1_pies <- plot_AOI_pie(geo_level4breast_B1_roi1, "Malignant") |
  plot_AOI_pie(geo_level4breast_B1_roi1, "Other") 

pdf(file = file.path(fig_path, "B1_3_roi1.pdf"),
    width = 4,
    height = 3)
print(p_B1_3_roi1_pies)
dev.off()


# B1_3 (roi2) ---------------------------------------------------------------
geo_level4breast_B1_roi2 <- geo_level4breast %>%
  filter(sample %in% c(
    "DSP-1001660018473-B-B03", # Malig
    "DSP-1001660018473-B-B04", # PanCK-
    "DSP-1001660018473-B-A02", # T cells
    "DSP-1001660018473-B-B05")) %>% # Macrophage
  arrange(cell_fraction) 

p_B1_3_roi2_pies <- plot_AOI_pie(geo_level4breast_B1_roi2, "Malignant") |
  plot_AOI_pie(geo_level4breast_B1_roi2, "PanCK-") | 
  plot_AOI_pie(geo_level4breast_B1_roi2, "T cells") | 
  plot_AOI_pie(geo_level4breast_B1_roi2, "Macrophage")

pdf(file = file.path(fig_path, "B1_3_roi2.pdf"),
    width = 8,
    height = 3)
print(p_B1_3_roi2_pies)
dev.off()

write.csv(geo_level4breast, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig4dg_decon_geo.csv") 

# B1_3 (roi3) ---------------------------------------------------------------
geo_level4breast_B1_roi3 <- geo_level4breast %>%
  filter(sample %in% c(
    "DSP-1001660018473-B-B07", # Other
    "DSP-1001660018473-B-B06" # T cells
  )) %>%
  arrange(cell_fraction) 

p_B1_3_roi3_pies <- 
  plot_AOI_pie(geo_level4breast_B1_roi3, "Other") | 
  plot_AOI_pie(geo_level4breast_B1_roi3, "T cells") 

pdf(file = file.path(fig_path, "B1_3_roi3.pdf"),
    width = 4,
    height = 3)
print(p_B1_3_roi3_pies)
dev.off()


# B1_3 (roi4) ---------------------------------------------------------------
geo_level4breast_B1_roi4 <- geo_level4breast %>%
  filter(sample %in% c(
    "DSP-1001660018473-B-B01", # Malig
    "DSP-1001660018473-B-B02", # PanCK-
    "DSP-1001660018473-B-B09", # T cells
    "DSP-1001660018473-B-B08")) %>% # Macrophage
  arrange(cell_fraction) 

p_B1_3_roi4_pies <- plot_AOI_pie(geo_level4breast_B1_roi4, "Malignant") |
  plot_AOI_pie(geo_level4breast_B1_roi4, "PanCK-") | 
  plot_AOI_pie(geo_level4breast_B1_roi4, "T cells") | 
  plot_AOI_pie(geo_level4breast_B1_roi4, "Macrophage")

pdf(file = file.path(fig_path, "B1_3_roi4.pdf"),
    width = 8,
    height = 3)
print(p_B1_3_roi4_pies)
dev.off()


# B1_3 (roi5) ---------------------------------------------------------------
geo_level4breast_B1_roi5 <- geo_level4breast %>% # A_islet_tme_1 (done) # A09, A10
  filter(sample %in% c(
    "DSP-1001660018473-B-A09", # Malig
    "DSP-1001660018473-B-A10", # PanCK-
    "DSP-1001660018473-B-B10" # Other->T cells
    )) %>% 
  arrange(cell_fraction) 

p_B1_3_roi5_pies <- plot_AOI_pie(geo_level4breast_B1_roi5, "Malignant") |
  plot_AOI_pie(geo_level4breast_B1_roi5, "PanCK-") | 
  plot_AOI_pie(geo_level4breast_B1_roi5, "Other") 

pdf(file = file.path(fig_path, "B1_3_roi5.pdf"),
    width = 6,
    height = 3)
print(p_B1_3_roi5_pies)
dev.off()


# B1_3 below pie (roi 6) --------------------------------------------------------
geo_level4breast_B1_roi6 <- geo_level4breast %>%
  filter(sample %in% c("DSP-1001660018473-B-B11",
                       "DSP-1001660018473-B-C01",
                       "DSP-1001660018473-B-B12")) %>%
  arrange(cell_fraction)

p_B1_3_roi6_pies <- plot_AOI_pie(geo_level4breast_B1_roi6, "Other") |
  plot_AOI_pie(geo_level4breast_B1_roi6, "T cells") | 
  plot_AOI_pie(geo_level4breast_B1_roi6, "Macrophage")

pdf(file = file.path(fig_path, "B1_3_roi6.pdf"),
    width = 6,
    height = 3)
print(p_B1_3_roi6_pies)
dev.off()


# B1_3 (roi7) ---------------------------------------------------------------
geo_level4breast_B1_roi7 <- geo_level4breast %>%
  filter(sample %in% c(
    "DSP-1001660018473-B-A11", # Malig
    "DSP-1001660018473-B-A12", # PanCK-
    "DSP-1001660018473-B-A05" # Macro->T cells
  )) %>% 
  arrange(cell_fraction) 

p_B1_3_roi7_pies <- plot_AOI_pie(geo_level4breast_B1_roi7, "Malignant") |
  plot_AOI_pie(geo_level4breast_B1_roi7, "PanCK-") | 
  plot_AOI_pie(geo_level4breast_B1_roi7, "Macrophage") 

pdf(file = file.path(fig_path, "B1_3_roi7.pdf"),
    width = 6,
    height = 3)
print(p_B1_3_roi7_pies)
dev.off()


# B1_3 (roi8) ---------------------------------------------------------------
geo_level4breast_B1_roi8 <- geo_level4breast %>% # A_islet_4 (done)
  filter(sample %in% c(
    "DSP-1001660018473-B-A07", # Malig
    "DSP-1001660018473-B-A08", # Other
    "DSP-1001660018473-B-A06" # T cells
  )) %>% 
  arrange(cell_fraction) 

p_B1_3_roi8_pies <- plot_AOI_pie(geo_level4breast_B1_roi8, "Malignant") |
  plot_AOI_pie(geo_level4breast_B1_roi8, "Other") | 
  plot_AOI_pie(geo_level4breast_B1_roi8, "T cells") 

pdf(file = file.path(fig_path, "B1_3_roi8.pdf"),
    width = 6,
    height = 3)
print(p_B1_3_roi8_pies)
dev.off()

