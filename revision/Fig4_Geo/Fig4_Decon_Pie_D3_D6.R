library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(tidyr)
library(patchwork)

# Sanity ------------------------------------------------------------------
## QCed
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/dlbcl_spe_ruv.rds") # 24
geo_D3 <- geo[, geo$patient == "D3"] # 23 AOIs
geo_D6 <- geo[, geo$patient == "D6"] # 19 AOIs

CD_D3 <- as.data.frame(colData(geo_D3))
CD_D3_roi <- CD_D3 %>% filter(roi == "C_7")

CD_D6 <- as.data.frame(colData(geo_D6))
CD_D6_roi <- CD_D6 %>% filter(roi %in% c("C_1", "C_2", "C_4", "C_5"))


## Level 4 GeoMx
geo_level4dlbcl <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/dlbcl_batched_decon_long.csv")


# Color
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/color_palette.R")
colors = level4_cellcolors[names(level4_cellcolors) %in% sort(geo_level4dlbcl$CellType)]

# Plot function
plot_AOI_pie <- function(roi_df, cf_type){
  p <- ggplot(roi_df[roi_df$cell_fraction == cf_type, ], 
              aes(fill=CellType, y=Fraction, x="")) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values =  colors) +
    coord_polar("y", start=0) + 
    theme_void() + 
    theme(legend.position="none")
  return(p)
}

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig4_Geo_pie_charts"

write.csv(geo_level4dlbcl, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig4hi_decon_geo.csv")

# D3 -------------------------------------------------------
geo_level4dlbcl_D3 <- geo_level4dlbcl %>%
  filter(sample %in% CD_D3_roi$sample_id2) %>% # "DSP-1001250001495-E-D08" "DSP-1001250001495-E-D09" "DSP-1001250001495-E-D10"
  arrange(cell_fraction) # "Mostly macrophages"   "Mostly T lymphocytes" "B lymphocytes" 

p_D3_pies <- plot_AOI_pie(geo_level4dlbcl_D3, "B cells") |
  plot_AOI_pie(geo_level4dlbcl_D3, "T cells") |
  plot_AOI_pie(geo_level4dlbcl_D3, "Macrophage") 

pdf(file = file.path(fig_path, "D3_roi.pdf"),
    width = 6,
    height = 3)
print(p_D3_pies)
dev.off()


# Sanity D6 -------------------------------------------------------------------------
geo_level4dlbcl_D6_sub <- geo_level4dlbcl %>%
  filter(patient == "D6") %>%
  left_join(CD_D6 %>% mutate(sample = sample_id2), by = "sample") %>%
  filter(roi %in% c("C_1", "C_2", "C_4", "C_5"))
length(unique(geo_level4dlbcl_D6_sub$sample)) # also 10, because decon run on QCed object. 

# Now check raw D6
geo_raw <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/GeoMx/dlbcl_raw/D6_spe.rds")
geo_raw_roi <- geo_raw[, geo_raw$roi %in% c("C_1", "C_2", "C_4", "C_5")]
df <- as.data.frame(colData(geo_raw_roi))

# QCed out aois
# DSP-1001660013971-F-C07 # C_1, CD3
# DSP-1001660013971-F-D10 # C_5, CD68


# D6 (roi1) -------------------------------------------------------
geo_level4dlbcl_D6 <- geo_level4dlbcl %>%
  filter(sample %in% CD_D6_roi$sample_id2[CD_D6_roi$roi == "C_1"]) %>% # "DSP-1001660013971-F-C08" "DSP-1001660013971-F-C09"
  arrange(cell_fraction) # "B lymphocytes"             "Epithelial cells (glands)"

p_D6_roi1_pies <- plot_AOI_pie(geo_level4dlbcl_D6, "B cells") |
  plot_AOI_pie(geo_level4dlbcl_D6, "Other")|
  plot_AOI_pie(geo_level4dlbcl_D6, "T cells")

pdf(file = file.path(fig_path, "D6_roi1.pdf"),
    width = 6,
    height = 3)
print(p_D6_roi1_pies)
dev.off()

# D6 (roi2) -------------------------------------------------------
geo_level4dlbcl_D6 <- geo_level4dlbcl %>%
  filter(sample %in% CD_D6_roi$sample_id2[CD_D6_roi$roi == "C_2"]) %>% # "DSP-1001660013971-F-D06" "DSP-1001660013971-F-D07"
  arrange(cell_fraction) # "B lymphocytes"             "Epithelial cells (glands)"

p_D6_roi2_pies <- plot_AOI_pie(geo_level4dlbcl_D6, "B cells") |
  plot_AOI_pie(geo_level4dlbcl_D6, "Other")

pdf(file = file.path(fig_path, "D6_roi2.pdf"),
    width = 4,
    height = 3)
print(p_D6_roi2_pies)
dev.off()

# D6 (roi3) -------------------------------------------------------
geo_level4dlbcl_D6 <- geo_level4dlbcl %>%
  filter(sample %in% CD_D6_roi$sample_id2[CD_D6_roi$roi == "C_4"]) %>% # DSP-1001660013971-F-D08 DSP-1001660013971-F-D09 
  arrange(cell_fraction) # "B lymphocytes"             "T cells"

p_D6_roi3_pies <- plot_AOI_pie(geo_level4dlbcl_D6, "B cells") |
  plot_AOI_pie(geo_level4dlbcl_D6, "T cells") 

pdf(file = file.path(fig_path, "D6_roi3.pdf"),
    width = 4,
    height = 3)
print(p_D6_roi3_pies)
dev.off()

# D6 (roi4) -------------------------------------------------------
geo_level4dlbcl_D6 <- geo_level4dlbcl %>%
  filter(sample %in% CD_D6_roi$sample_id2[CD_D6_roi$roi == "C_5"]) %>% # "DSP-1001660013971-F-C08" "DSP-1001660013971-F-C09"
  arrange(cell_fraction) # "B lymphocytes"             "T cells"

p_D6_roi4_pies <- plot_AOI_pie(geo_level4dlbcl_D6, "B cells") |
  plot_AOI_pie(geo_level4dlbcl_D6, "T cells") |
  plot_AOI_pie(geo_level4dlbcl_D6, "Macrophage")

pdf(file = file.path(fig_path, "D6_roi4.pdf"),
    width = 6,
    height = 3)
print(p_D6_roi4_pies)
dev.off()


# B_plasma                       DC_1                       DC_2                      DC_pc                 Endothelia 
# "#FF00FF"                  "#C71585"                  "#FFAEB9"                  "#FF6347"                  "#FF7256" 
# Epi_bronchus                 Epi_mucous Epi_Mucous_surface_gastric       Epi_parietal_gastric               Fibro_Muscle 
# "#BC8F8F"                  "#A52A2A"                  "#FFC1C1"                  "#CD5C5C"                  "#388E8E" 
# Mono_Macro                         NK                   Pericyte                      T_CD4                  T_CD4_reg 
# "#9A32CD"                  "#79CDCD"                  "#A0522D"                  "#4169E1"                  "#6495ED" 
# T_CD8                 T_dividing                 Tu_D1_LMO2                Tu_D1_RGS13               Tu_D1_SMIM14 
# "#009ACD"                  "#97FFFF"                  "#FF8C00"                  "#E3A869"                  "#CD6600" 
# Tu_D2_mito               Tu_D3_BCL2A1             Tu_D3_dividing                Tu_D3_FAM3C                 Tu_D3_IGHD 
# "#EEEE00"                  "#CDC673"                  "#006400"                  "#FFD700"                  "#00CD00" 

