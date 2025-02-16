library(ggplot2)
library(dplyr)
library(stringr)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig3"


# deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results_pt_specific"
# prev_deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results/"
# 
# deconpath = deconresultpath; plot_title = "Geo_breast_lung_decon_combined_box_outline_color_pt_specific.pdf"
# deconpath = prev_deconresultpath; plot_title = "Geo_breast_lung_decon_combined_box_outline_color.pdf"


# Only breast and lung ----------------------------------------------------
# Mapped spot and AOI IDs (spot)
mapped_geo_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Registration/GeoMx_mapped_spot_LW"
geo_samples <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1", "L1_1", "L2_1", "L3_3", "L4_1")
vis_samples <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2", "L1_2", "L2_2", "L3_2", "L4_2")
disease_list <- c(rep("breast", 5), rep("lung", 4))

# mapped_all <- NULL
# for(i in 1:length(geo_samples)){
#   mapped_i <- read.csv(file.path(mapped_geo_path, paste0(geo_samples[i], "_mapped_LW.csv")), row.names = 1) %>% 
#     select(barcode, Section_ID, Sample_ID, Cell_fraction) %>%
#     dplyr::rename(cell_fraction = Cell_fraction) %>%
#     mutate(Section_ID_Vis = vis_samples[i])
#   
#   # Merge in pathology ------------------------------------------------------
#   sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease_list[i], "_qcd")
#   sce <- readRDS(file.path(sce_path, paste0(vis_samples[i], "_qcd.rds")))
#   
#   mapped_i <- mapped_i %>%
#     left_join(data.frame(barcode = colnames(sce),
#                          Region = sce$Region))
#   
#   mapped_all <- rbind(mapped_all, mapped_i)
# }
# 
# write.csv(mapped_all, 
#           file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration", "spot_mapped_cf_region.csv"), row.names = FALSE)


# -------------------------------------------------------------------------
mapped_all <- read.csv(file.path("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Registration", "spot_mapped_cf_region.csv"))

# Geo decon results (spot)
# geo_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results"
deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/For_level1_5_immune_decon_results"

geo_decon <- rbind(read.csv(file.path(deconresultpath, "breast_batched_decon_long.csv")),
                   read.csv(file.path(deconresultpath, "lung_batched_decon_long.csv"))
                   # read.csv(file.path(geo_decon_path, "dlbcl_batched_decon_long.csv")) %>% rename(section_id = patient)
)
# geo_macro <- read.csv(file.path(geo_decon_path, "geo_batched_macro_long.csv"))
# geo_macro <- geo_macro %>%
#   filter(section_id == "L3_1") # 237 -> 32

geo_decon_mapped <- rbind(
  geo_decon[geo_decon$sample %in% unique(mapped_all$Sample_ID), ] #,
  # geo_macro
  ) %>%
  dplyr::rename(sample_id = sample) # 431 -> 668 by adding 237 macros 
geo_decon_mapped$cell_fraction = ifelse(geo_decon_mapped$cell_fraction == "Macro", "Macrophage", geo_decon_mapped$cell_fraction)  


# > head(geo_decon_mapped) # 431   5
#                     sample CellType   Fraction section_id cell_fraction
# 14 DSP-1001660018473-B-D03     T_NK 0.18231696       B4_1         Other
# 16 DSP-1001660018473-B-D05     T_NK 0.17330728       B4_1     Malignant
# 17 DSP-1001660018473-B-D06     T_NK 0.19302811       B4_1         Other
# 19 DSP-1001660018473-B-D08     T_NK 0.26285257       B4_1         Other
# 21 DSP-1001660018473-B-D10     T_NK 0.14996134       B4_1         Other
# 24 DSP-1001660018473-B-E01     T_NK 0.09200995       B4_1         Other

# # Vis decon results (spot) --------------------------------------------
# # vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2c"
# # ## CARD
# vis_decon <- rbind(read.csv(file.path(vis_decon_path, "vis_breast_decon_long.csv")),
#                    read.csv(file.path(vis_decon_path, "vis_lung_decon_long.csv"))
# )
# # ## RCTD
# # vis_decon <- rbind(read.csv(file.path(vis_decon_path, "vis_breast_decon_long_RCTD.csv")),
# #                    read.csv(file.path(vis_decon_path, "vis_lung_decon_long_RCTD.csv"))
# # )

# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2_final_immune" # CARD
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/cell2location_vis_long" # C2L
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/RCTD_vis_long" # RCTD - level 1.5 immune
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD/level1_5_immune_long" # RCTD - level 4 combined to level 1.5 immunne 
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD_DEgenes/level1_5_immune_long" # RCTD - level 4 DEgenes combined to level 1.5 immunne 
vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L/level1_5_immune_long" # C2L - level 4 combined to level 1.5 immunne 
# vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L_DEgenes/level1_5_immune_long" # C2L - level 4 DEgenes combined to level 1.5 immunne 
decon_patho_gather_breast <- read.csv(file.path(vis_decon_path, "vis_breast_decon_immune_long.csv"))
decon_patho_gather_lung <- read.csv(file.path(vis_decon_path, "vis_lung_decon_immune_long.csv"))

vis_decon <- rbind(decon_patho_gather_breast,
                   decon_patho_gather_lung)

# vis_macro <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2d/vis_spot_70pct_myeloid_breast_lung_decon_long.csv")
# vis_macro <- vis_macro %>%
#   mutate(Sample_ID = "-",
#          Section_ID = "L3_1",
#          Section = "L3_2",
#          cell_fraction = "Macrophage")

vis_decon_mapped_ <- mapped_all %>%
  dplyr::rename(Barcode = barcode,
                Section = Section_ID_Vis) %>%
  filter(cell_fraction != "PanCK-") %>%
  left_join(vis_decon, by = c("Barcode", "Section"))

  
vis_decon_mapped <- rbind(vis_decon_mapped_
                          #, vis_macro
                          ) %>% # 3385 -> 3425 # increase 40 rows, 5 spots for macro
  dplyr::rename(section_id = Section_ID,
       sample_id = Sample_ID) %>%
  mutate(cell_fraction = ifelse(cell_fraction == "T_cells", "T cells", cell_fraction)) %>%
  mutate(CellType = ifelse(CellType == "Fibro_Muscle", "Fibro_muscle", CellType))

# >  head(vis_decon_mapped)   # no macrophage: 3385    8
#              Barcode Section_ID               Sample_ID cell_fraction Section Region     CellType     Fraction
# 1 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE         T_NK 0.0002814366
# 2 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE            B 0.0248407509
# 3 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE        Tumor 0.1323340298
# 4 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE       Vessel 0.0085140267
# 5 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE Fibro_Muscle 0.8316757981
# 6 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE      Myeloid 0.0004073472
  
vis_geo_spot_decon <- rbind(
  geo_decon_mapped %>% select(CellType, Fraction, cell_fraction, section_id, sample_id) %>%
    mutate(Platform = "GeoMx"),
  vis_decon_mapped %>% select(CellType, Fraction, cell_fraction, section_id, sample_id) %>%
    mutate(Platform = "Visium")
)

vis_geo_spot_decon <- vis_geo_spot_decon %>%
  filter(!is.na(CellType)) # B1_2 after patho qcd out

vis_geo_spot_decon$CellType <- factor(
  vis_geo_spot_decon$CellType,
  levels = c("Epithelia", "Stroma", "Tumor", "Macrophage", "T cells", "B cells", "NK", "Myeloid else")) # Fig 2, 3 

# -------------------------------------------------------------------------
# All mapped T cells spots are from breast, so no epithelia cell type, so remove it from T cells facet x axis by facet_wrap2(scales = "free_x")
# check <- gathered_df_breast_lung %>% filter(cell_fraction == "T cells")
# table(check$CellType)
# Epithelia       Stroma        Tumor   Macrophage      T cells      B cells           NK Myeloid else 
#   0            8            8            8            8            8            8            8 
# vis_geo_spot_decon_ <- vis_geo_spot_decon
# 
# CellType 
# Fraction = c(0, 0)
# cell_fraction = rep("")
# section_id = rep(NA,2)
# sample_id = rep(NA,2)
# Platform = c("GeoMx", "Visium")

# vis_geo_spot_decon_test <- vis_geo_spot_decon %>%
#   filter(is.na(CellType))
# 
# library(rstatix)
# stat.test <- vis_geo_spot_decon %>%
#   filter(cell_fraction == "Malignant") %>%
#   group_by(CellType) %>%
#   t_test(Fraction ~ Platform) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance("p.adj")
# stat.test
# #     CellType     .y.      group1 group2    n1    n2 statistic    df        p    p.adj p.adj.signif
# #     <chr>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>       
# #   1 B            Fraction GeoMx  Visium    28   117    -8.53  125.  4.12e-14 3.30e-13 ****       #  
# #   2 Epithelia    Fraction GeoMx  Visium    17    73    -4.30   30.2 1.63e- 4 1.30e- 3 **         # 
# #   3 Fibro_muscle Fraction GeoMx  Visium    28   117    -6.45  136.  1.79e- 9 1.43e- 8 ****       #  
# #   4 Granulocyte  Fraction GeoMx  Visium    28   117    -1.37   31.4 1.79e- 1 1   e+ 0 ns          
# #   5 Myeloid      Fraction GeoMx  Visium    28   117     0.438  59.2 6.63e- 1 1   e+ 0 ns          
# #   6 T_NK         Fraction GeoMx  Visium    28   117     4.40   33.3 1.04e- 4 8.32e- 4 ***        #  
# #   7 Tumor        Fraction GeoMx  Visium    28   117     5.26   74.1 1.35e- 6 1.08e- 5 ****       # 
# #   8 Vessel       Fraction GeoMx  Visium    28   117    -0.497  47.0 6.21e- 1 1   e+ 0 ns     
# 
# stat.test <- vis_geo_spot_decon %>%
#   filter(cell_fraction == "Other") %>%
#   filter(CellType != "Epithelia") %>% # just one sample per platform
#   group_by(CellType) %>%
#   t_test(Fraction ~ Platform) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance("p.adj")
# stat.test
# #     CellType     .y.      group1 group2    n1    n2 statistic    df          p     p.adj p.adj.signif
# #    <chr>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>      <dbl>     <dbl> <chr>       
# #   1 B            Fraction GeoMx  Visium    28   350    -2.50   43.9 0.0164     0.115     ns          
# #   2 Fibro_muscle Fraction GeoMx  Visium    28   350     1.19   71.7 0.237      1         ns          
# #   3 Granulocyte  Fraction GeoMx  Visium    28   350    -5.69   33.9 0.00000219 0.0000153 ****        
# #   4 Myeloid      Fraction GeoMx  Visium    28   350    -0.574  33.5 0.57       1         ns          
# #   5 T_NK         Fraction GeoMx  Visium    28   350    -4.46   34.8 0.0000824  0.000577  ***         
# #   6 Tumor        Fraction GeoMx  Visium    28   350    -1.44   76.4 0.155      1         ns          
# #   7 Vessel       Fraction GeoMx  Visium    28   350     3.32   60.4 0.00151    0.0106    *      # 
# 
# stat.test <- vis_geo_spot_decon %>%
#   filter(cell_fraction == "T cells") %>%
#   group_by(CellType) %>%
#   t_test(Fraction ~ Platform) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance("p.adj")
# stat.test
# # CellType     .y.      group1 group2    n1    n2 statistic    df        p   p.adj p.adj.signif
# #  <chr>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>   <dbl> <chr>       
# # 1 B            Fraction GeoMx  Visium     3     6     0.362  4.04 0.735    1       ns          
# # 2 Fibro_muscle Fraction GeoMx  Visium     3     6    -0.823  6.86 0.438    1       ns          
# # 3 Granulocyte  Fraction GeoMx  Visium     3     6    -4.66   5.82 0.00376  0.0263  *           
# # 4 Myeloid      Fraction GeoMx  Visium     3     6     0.260  3.32 0.81     1       ns          
# # 5 T_NK         Fraction GeoMx  Visium     3     6     0.261  5.56 0.804    1       ns          
# # 6 Tumor        Fraction GeoMx  Visium     3     6    -0.630  6.91 0.549    1       ns          
# # 7 Vessel       Fraction GeoMx  Visium     3     6     6.32   7.00 0.000397 0.00278 **           #
# 
# annotation_df <- data.frame(
#   color = c("E", "H"),
#   start = c("Good", "Fair"),
#   end = c("Very Good", "Good"),
#   y = c(3.6, 4.7),
#   label = c("Comp. 1", "Comp. 2")
# )

# -------------------------------------------------------------------------
library(ggh4x)

# cell_fraction_order = c("Macrophage", "Malignant", "Other", "T cells", "PanCK-") # note with PanCK-
# cell_fraction_color = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1", "#09D0EF")
# c("#bc8f8f55" "#388e8e26", "#4169e11e")

facet_colors <- c("Malignant" = "#BC8F8F", "Other" = "#388E8E", "T cells" = "#4169E1")
facet_fills <- c("Malignant" = "#bc8f8f55", "Other" = "#388e8e26", "T cells" = "#4169e11e")


p <- ggplot(vis_geo_spot_decon, aes(x=CellType, y=Fraction, fill= Platform)) +
  geom_boxplot(outlier.size = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c("#92d051ff", "#257cf2ff")) + 
  theme(axis.text.x = element_text(size = 12.5, angle = 90, vjust = 0.5, hjust=1),
        panel.spacing=unit(1.5,"lines"),
        panel.grid = element_blank(), 
        axis.title = element_blank(),
        strip.text.x = element_text(size = 13.5, face = "bold"), 
        strip.background=element_rect(color = "#BC8F8F", 
                                      fill = "#bc8f8f55")) + 
  ggh4x::facet_wrap2(~cell_fraction, 
                     scales = "free_x",
                     strip = strip_themed(background_x = elem_list_rect(color = facet_colors,
                                                                        fill = facet_fills)),
                     ncol = 4)

write.csv(vis_geo_spot_decon, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig3c.csv")

library(rstatix)
stat.test <- vis_geo_spot_decon %>%
  filter(cell_fraction == "Malignant") %>%
  group_by(CellType) %>%
  t_test(Fraction ~ Platform) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# CellType     .y.      group1 group2    n1    n2 statistic    df        p    p.adj p.adj.signif
# <fct>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>       
#   1 Epithelia    Fraction GeoMx  Visium    17    73     -1.50  39.5 1.41e- 1 1   e+ 0 ns          
#   2 Stroma       Fraction GeoMx  Visium    28   116     -3.39  74.4 1.12e- 3 8.96e- 3 **          
#   3 Tumor        Fraction GeoMx  Visium    28   116      3.54  60.7 7.69e- 4 6.15e- 3 **          
#   4 Macrophage   Fraction GeoMx  Visium    28   116     -1.07  80.7 2.88e- 1 1   e+ 0 ns          
#   5 T cells      Fraction GeoMx  Visium    28   116      4.36  30.7 1.36e- 4 1.09e- 3 **          
#   6 B cells      Fraction GeoMx  Visium    28   116     -8.04 139.  3.57e-13 2.86e-12 ****        
#   7 NK           Fraction GeoMx  Visium    28   116     -4.03  73.1 1.37e- 4 1.10e- 3 **          
#   8 Myeloid else Fraction GeoMx  Visium    28   116     -2.07  74.9 4.17e- 2 3.34e- 1 ns  

stat.test <- vis_geo_spot_decon %>%
  filter(cell_fraction == "Other") %>%
  filter(CellType != "Epithelia") %>% # just one sample per platform
  group_by(CellType) %>%
  t_test(Fraction ~ Platform) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# CellType     .y.      group1 group2    n1    n2 statistic    df        p    p.adj p.adj.signif
# <fct>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>       
#   1 Stroma       Fraction GeoMx  Visium    28   342     16.1   64.4 3.27e-24 2.29e-23 ****        
#   2 Tumor        Fraction GeoMx  Visium    28   342     -4.63  78.8 1.44e- 5 1.01e- 4 ***         
#   3 Macrophage   Fraction GeoMx  Visium    28   342     -2.60  40.0 1.31e- 2 9.17e- 2 ns          
#   4 T cells      Fraction GeoMx  Visium    28   342    -10.6   38.8 5.28e-13 3.70e-12 ****        
#   5 B cells      Fraction GeoMx  Visium    28   342     -3.33  42.2 1.83e- 3 1.28e- 2 *           
#   6 NK           Fraction GeoMx  Visium    28   342    -16.4  298.  1.72e-43 1.20e-42 ****        
#   7 Myeloid else Fraction GeoMx  Visium    28   342    -12.9   49.8 1.61e-17 1.13e-16 ****    

stat.test <- vis_geo_spot_decon %>%
  filter(cell_fraction == "T cells") %>%
  group_by(CellType) %>%
  t_test(Fraction ~ Platform) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# CellType     .y.      group1 group2    n1    n2 statistic    df      p p.adj p.adj.signif
# <fct>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>  <dbl> <dbl> <chr>       
# 1 Stroma       Fraction GeoMx  Visium     3     5     1.48   5.96 0.189  1     ns          
# 2 Tumor        Fraction GeoMx  Visium     3     5     0.221  5.99 0.833  1     ns          
# 3 Macrophage   Fraction GeoMx  Visium     3     5    -1.18   5.09 0.291  1     ns          
# 4 T cells      Fraction GeoMx  Visium     3     5    -0.411  5.91 0.695  1     ns          
# 5 B cells      Fraction GeoMx  Visium     3     5     0.454  3.84 0.674  1     ns          
# 6 NK           Fraction GeoMx  Visium     3     5    -0.888  4.25 0.422  1     ns          
# 7 Myeloid else Fraction GeoMx  Visium     3     5    -2.73   5.86 0.0349 0.244 ns    


fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig3"
pdf(file = file.path(fig_path, "Vis_Geo_Regi_sbs_AOI_facet_Cell2location_final.pdf"),
    width = 10,
    height = 4.5)
print(p)
dev.off()





