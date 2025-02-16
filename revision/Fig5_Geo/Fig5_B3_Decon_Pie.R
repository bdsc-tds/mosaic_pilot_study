library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(tidyr)
library(patchwork)

## Level 4 GeoMx
geo_level4breast <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/breast_batched_decon_long.csv")

geo_level4breast_B3_ROI <- geo_level4breast %>% 
  filter(section_id == "B3_1")

## Decon ordering the bars
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds")

CD <- as.data.frame(colData(geo))
table(CD$roi, CD$cell_fraction)

CD_B3_1 <- CD %>%
  dplyr::filter(section_id == "B3_1") 
length(table(CD_B3_1$roi)) # 11 ROI


CD_B3_1_order <- CD_B3_1 %>%
  mutate(AOIn = case_when(roi == "C_Islet_7" ~ "11",
                          roi == "C_Islet_3" ~ "10",
                          roi == "C_Islet_4" ~ "9",
                          roi == "C_Islet_5" ~ "8",
                          roi == "C_TME_2" ~ "7",
                          roi == "C_Islet_2" ~ "6",
                          roi == "C_Islet_1" ~ "5",
                          roi == "C_TME_1" ~ "4",
                          roi == "C_Islet_6" ~ "3",
                          roi == "C_TME_8" ~ "2",
                          roi == "C_Islet_8" ~ "1"
  )) %>%
  mutate(AOIn_cf = paste0(AOIn, "_", cell_fraction)) %>%
  select(sample_id2, AOIn, AOIn_cf) %>%
  arrange(as.numeric(AOIn)) %>%
  dplyr::rename(sample = sample_id2)


# Merge in  ---------------------------------------------------------------
head(geo_level4breast_B3_ROI)
geo_level4breast_B3_ROI_sum <- geo_level4breast_B3_ROI %>%
  pivot_wider(values_from = Fraction, names_from = CellType) %>%
  rowwise() %>% 
  mutate(
    B_cells = B_plasma_IGHA1 + B_plasma_IGHG1 + B_plasma_IGHG3 + B_plasma_IGHM + B_plasma_IGKC + B_plasma_IGLC1,
    Granulocyte = Granulocyte + Mast_cell, 
    Myeloid = DC_1 + DC_2 + DC_activated + DC_pc + Macrophage + Monocyte,
    Stroma = Endothelia_vascular + Fibroblast_B3 + Muscle_smooth + Pericyte, 
    T_NK = T_CD4 + T_CD8_exhausted + T_CTL + T_CXCL13 + T_reg + TNK_dividing + NK) %>%
  select(sample, section_id, cell_fraction, B_cells, Granulocyte, Myeloid, Stroma, T_NK, Tu_B3_CYP4F8) %>%
  dplyr::rename(Tu_B3 = Tu_B3_CYP4F8) %>%
  pivot_longer(cols = c(B_cells, Granulocyte, Myeloid, Stroma, T_NK, Tu_B3), 
               values_to = "Fraction", names_to = "CellType")


level1_5_cellnames <- c(    "B_cells", "Granulocyte", "Myeloid", "Stroma",  "T_NK",    "Tu_B3")
br_level1_5_cellcolors <- c("#FFDAB9", "#FFC125",     "#FF9900", "#388E8E", "#00BFFF", "#EEEE00")
names(br_level1_5_cellcolors) <- level1_5_cellnames

## Merge in and Order
geo_level4breast_B3_ROI_plot <- geo_level4breast_B3_ROI_sum %>%
  left_join(CD_B3_1_order) %>%
  mutate(AOIn = factor(AOIn, levels = as.character(1:11)),
         cell_fraction = factor(cell_fraction, levels = c("Malignant", "Other", "T cells", "PanCK-")),
         CellType = factor(CellType, levels = level1_5_cellnames)) %>%
  arrange(as.numeric(AOIn), cell_fraction, CellType) %>%
  filter(!is.na(AOIn))


unique(geo_level4breast_B3_ROI_plot$AOIn_cf)
geo_level4breast_B3_ROI_plot$AOIn_cf = factor(geo_level4breast_B3_ROI_plot$AOIn_cf,
                                              levels = rev(c("1_Malignant",  "1_Other",      "2_Malignant",  "2_Other",      "2_T_cells",    "3_Malignant",  "3_Other",      "4_Other",      "4_T_cells",    "5_Malignant",  "5_Other",      "6_Malignant",
                                                             "6_Other",      "7_PanCK-",     "8_Malignant",  "8_Other",      "8_T_cells",    "9_Malignant",  "9_Other",      "10_Malignant", "10_Other",     "11_Malignant", "11_Other")))


# Plot function
plot_AOI_pie <- function(roi_df, cf_type){
  p <- ggplot(roi_df[roi_df$cell_fraction == cf_type, ], 
              aes(fill=CellType, y=Fraction, x="")) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values =  br_level1_5_cellcolors) +
    coord_polar("y", start=0) + 
    theme_void() + 
    theme(legend.position="none")
  return(p)
}

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3"


write.csv(geo_level4breast_B3_ROI_plot, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig5c_decon.csv")

######################################################################
# B3 (roi 1) -------------------------------------------------------
geo_level4breast_B3_roi1 <- geo_level4breast_B3_ROI_plot %>%  
  filter(AOIn == 1) %>% 
  arrange(cell_fraction) 

p_B3_roi1_pies <- plot_AOI_pie(geo_level4breast_B3_roi1, "Malignant") |
  plot_AOI_pie(geo_level4breast_B3_roi1, "Other") 

pdf(file = file.path(fig_path, "B3_roi1.pdf"),
    width = 4,
    height = 3)
print(p_B3_roi1_pies)
dev.off()


# B3 (roi2) ---------------------------------------------------------------
geo_level4breast_B3_roi2 <- geo_level4breast_B3_ROI_plot %>%
  filter(AOIn == 2) %>% 
  arrange(cell_fraction) 

p_B3_roi2_pies <- plot_AOI_pie(geo_level4breast_B3_roi2, "Malignant") |
  plot_AOI_pie(geo_level4breast_B3_roi2, "Other") | 
  plot_AOI_pie(geo_level4breast_B3_roi2, "T cells")

pdf(file = file.path(fig_path, "B3_roi2.pdf"),
    width = 6,
    height = 3)
print(p_B3_roi2_pies)
dev.off()



# B3 (roi3) ---------------------------------------------------------------
geo_level4breast_B3_roi3 <- geo_level4breast_B3_ROI_plot %>%
  filter(AOIn == 3) %>% 
  arrange(cell_fraction) 

p_B3_roi3_pies <- plot_AOI_pie(geo_level4breast_B3_roi3, "Malignant") | 
  plot_AOI_pie(geo_level4breast_B3_roi3, "Other") 

pdf(file = file.path(fig_path, "B3_roi3.pdf"),
    width = 4,
    height = 3)
print(p_B3_roi3_pies)
dev.off()


# B3 (roi4) ---------------------------------------------------------------
geo_level4breast_B3_roi4 <- geo_level4breast_B3_ROI_plot %>%
  filter(AOIn == 4) %>% 
  arrange(cell_fraction) 

p_B3_roi4_pies <-  plot_AOI_pie(geo_level4breast_B3_roi4, "Other") |
  plot_AOI_pie(geo_level4breast_B3_roi4, "T cells") 

pdf(file = file.path(fig_path, "B3_roi4.pdf"),
    width = 4,
    height = 3)
print(p_B3_roi4_pies)
dev.off()


# B3 (roi5) ---------------------------------------------------------------
geo_level4breast_B3_roi5 <- geo_level4breast_B3_ROI_plot %>% 
  filter(AOIn == 5) %>% 
  arrange(cell_fraction) 

p_B3_roi5_pies <- plot_AOI_pie(geo_level4breast_B3_roi5, "Malignant") |
  plot_AOI_pie(geo_level4breast_B3_roi5, "Other") 

pdf(file = file.path(fig_path, "B3_roi5.pdf"),
    width = 4,
    height = 3)
print(p_B3_roi5_pies)
dev.off()


# B3 below pie (roi 6) --------------------------------------------------------
geo_level4breast_B3_roi6 <- geo_level4breast_B3_ROI_plot %>%
  filter(AOIn == 6) %>% 
  arrange(cell_fraction) 

p_B3_roi6_pies <- plot_AOI_pie(geo_level4breast_B3_roi6, "Malignant") |
  plot_AOI_pie(geo_level4breast_B3_roi6, "Other") 

pdf(file = file.path(fig_path, "B3_roi6.pdf"),
    width = 4,
    height = 3)
print(p_B3_roi6_pies)
dev.off()


# B3 (roi7) ---------------------------------------------------------------
geo_level4breast_B3_roi7 <- geo_level4breast_B3_ROI_plot %>%
  filter(AOIn == 7) %>% 
  arrange(cell_fraction) 

p_B3_roi7_pies <- plot_AOI_pie(geo_level4breast_B3_roi7, "PanCK-") 

pdf(file = file.path(fig_path, "B3_roi7.pdf"),
    width = 2,
    height = 3)
print(p_B3_roi7_pies)
dev.off()


# B3 (roi8) ---------------------------------------------------------------
geo_level4breast_B3_roi8 <- geo_level4breast_B3_ROI_plot %>% 
  filter(AOIn == 8) %>% 
  arrange(cell_fraction) 

p_B3_roi8_pies <- plot_AOI_pie(geo_level4breast_B3_roi8, "Malignant") |
  plot_AOI_pie(geo_level4breast_B3_roi8, "Other") | 
  plot_AOI_pie(geo_level4breast_B3_roi8, "T cells") 

pdf(file = file.path(fig_path, "B3_roi8.pdf"),
    width = 6,
    height = 3)
print(p_B3_roi8_pies)
dev.off()


# B3 (roi 9) --------------------------------------------------------
geo_level4breast_B3_roi9 <- geo_level4breast_B3_ROI_plot %>%
  filter(AOIn == 9) %>% 
  arrange(cell_fraction) 

p_B3_roi9_pies <- plot_AOI_pie(geo_level4breast_B3_roi9, "Malignant") |
  plot_AOI_pie(geo_level4breast_B3_roi9, "Other") 

pdf(file = file.path(fig_path, "B3_roi9.pdf"),
    width = 4,
    height = 3)
print(p_B3_roi9_pies)
dev.off()


# B3 (roi 10) --------------------------------------------------------
geo_level4breast_B3_roi10 <- geo_level4breast_B3_ROI_plot %>%
  filter(AOIn == 10) %>% 
  arrange(cell_fraction) 

p_B3_roi10_pies <- plot_AOI_pie(geo_level4breast_B3_roi10, "Malignant") |
  plot_AOI_pie(geo_level4breast_B3_roi10, "Other") 

pdf(file = file.path(fig_path, "B3_roi10.pdf"),
    width = 4,
    height = 3)
print(p_B3_roi10_pies)
dev.off()


# B3 (roi 11) --------------------------------------------------------
geo_level4breast_B3_roi11 <- geo_level4breast_B3_ROI_plot %>%
  filter(AOIn == 11) %>% 
  arrange(cell_fraction) 

p_B3_roi11_pies <- plot_AOI_pie(geo_level4breast_B3_roi11, "Malignant") |
  plot_AOI_pie(geo_level4breast_B3_roi11, "Other") 

pdf(file = file.path(fig_path, "B3_roi11.pdf"),
    width = 4,
    height = 3)
print(p_B3_roi11_pies)
dev.off()




