library(readxl)
library(dplyr)
library(writexl)
library(SingleCellExperiment)

# ## LabWorkSheet ----------------------------------------------------
# owkinLW <- read_excel("~/Desktop/ProjectWork/Owkin_Pilot/data/owkinLW_update.xlsx", sheet = "merged_LabWorkSheet")
# # owkinLW <- read_excel("~/Desktop/ProjectWork/Owkin_Pilot/data/owkinLW.xlsx", sheet = "merged_LabWorkSheet")
# # owkinLW$roi[owkinLW$Section_ID == "L1_1" & owkinLW$`slide name` == "#CH_L_p002_WTA" & owkinLW$roi == "A_Islet_1"] <- "A_Islet_1_"
# # write_xlsx(list(merged_LabWorkSheet = owkinLW), "~/Desktop/ProjectWork/Owkin_Pilot/data/owkinLW_update.xlsx")
# 
# 
# spotmatch <- read.csv("~/Desktop/ProjectWork/Owkin_Pilot/data/spots_aois_matching.csv")
# 
# # Owkin LW small ----------------------------------------------------------
# owkinLW <- owkinLW %>%
#   select(Section_ID, everything()) %>%
#   arrange(Section_ID, Sample_ID, roi) 
# 
# owkinLW_small <- owkinLW %>%
#   select(Sample_ID, Section_ID, segment, Cell_fraction, Cell_fraction_location, `slide name`, roi, ROI_category) %>%
#   arrange(Section_ID, roi) %>%
#   filter(!is.na(ROI_category))
# head(owkinLW_small[, c("Section_ID", "roi", "segment")])
# 
# owkinLW_small <- owkinLW_small %>%
#   mutate(segment = case_when(segment %in% c("PanCK", "CK+", "CK") ~ "PanCK+",
#                              segment %in% c("CD3", "CD3+") ~ "CD3+",
#                              segment %in% c("other", "Other") ~ "Other",
#                              segment %in% c("CD68", "CD68+") ~ "CD68+",
#                              segment == "Full ROI" ~ "Full ROI"
#   ))                                                         # dim 242 8
# 
# 
# # Spot match small ---------------------------------------------------------
# spotmatch_small <- spotmatch
# colnames(spotmatch_small)[colnames(spotmatch_small) == "Sample.ID"] <- "Section_ID_Vis"
# colnames(spotmatch_small)[colnames(spotmatch_small) == "ROI.name"] <- "roi"
# colnames(spotmatch_small)[colnames(spotmatch_small) == "AOI"] <- "segment"
# spotmatch_small$Section_ID_Vis <- gsub("-", "_", spotmatch_small$Section_ID_Vis)
# 
# spotmatch_small <- spotmatch_small %>%
#   arrange(Section_ID_Vis, roi) %>%
#   mutate(Section_ID = case_when(Section_ID_Vis == "B1_2" ~ "B1_1",
#                                 Section_ID_Vis == "B1_4" ~ "B1_3",
#                                 Section_ID_Vis == "B2_2" ~ "B2_1",
#                                 Section_ID_Vis == "B3_2" ~ "B3_1",
#                                 Section_ID_Vis == "B4_2" ~ "B4_1",
#                                 Section_ID_Vis == "L1_2" ~ "L1_1",
#                                 Section_ID_Vis == "L1_4" ~ "None",
#                                 Section_ID_Vis == "L2_2" ~ "L2_1",
#                                 Section_ID_Vis == "L3_2" ~ "L3_3",
#                                 Section_ID_Vis == "L4_2" ~ "L4_1",
#   ))
# head(owkinLW_small[, c("Section_ID", "roi", "segment")])
# 
# spotmatch_small <- spotmatch_small %>%
#   mutate(segment = case_when(segment %in% c("PanCK", "CK+") ~ "PanCK+",
#                              segment %in% c("CD3", "CD3+") ~ "CD3+",
#                              segment %in% c("other", "Other") ~ "Other",
#                              segment == "Full ROI" ~ "Full ROI"
#   )) %>%
#   filter(Total.intersection.area..px./Spot.area..px. >= 0.7) # dim 627 9
# 
# 
# 
# final <- spotmatch_small %>%
#   left_join(owkinLW_small, by = c("Section_ID", "roi", "segment")) %>%
#   arrange(Section_ID, roi)
# 
# write_xlsx(final, "~/Desktop/Owkin/data/GeoMx/spot_aoi_matched.xlsx")


# Plot matched spots ------------------------------------------------------
final_spot <- read_xlsx("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Vis_Geo_Match/spot_aoi_matched.xlsx")
final_spot <- final_spot %>%
  filter(overlap_pct > 0.7)

# disease = "breast"
# source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
# samples_for_registration <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2") # Visium
# geo_reg_names <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1")

disease = "lung"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
samples_for_registration <- c("L1_2", "L2_2", "L3_2", "L4_2") # Visium
geo_reg_names <- c("L1_1", "L2_1", "L3_3", "L4_1")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

# i = 3
for(i in 1:length(samples_for_registration)){

  # Expr -----------------------------------------------------------
  # save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace/", foldername, "/", samples_for_registration[i], "/")
  # sce <- readRDS(file.path(save_bs_path, paste0(samples_for_registration[i], "_baye_clustered.rds"))) # spot
  save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/", foldername, "/", samples_for_registration[i], "/")
  sce <- readRDS(file.path(save_bs_path, paste0(samples_for_registration[i], "_baye_clustered.rds"))) # spot
  assays(sce)
  
  spotmatch_small_subset <- final_spot %>%
    filter(Section_ID_Vis == samples_for_registration[i]) %>%
    dplyr::rename(barcode = Spot.barcode,
           sample = Section_ID_Vis) #%>%
    # mutate(barcode = paste0(sample, "_", barcode))
  
  # # A tibble: 6 × 15
  # Section_ID_Vis ROI.ID roi       Spot.array.location Spot.barcode segment Total.intersection.a…¹ Spot.area..px. Section_ID overlap_pct Sample_ID Cell_fraction
  # <chr>          <chr>  <chr>     <chr>               <chr>        <chr>                    <dbl>          <dbl> <chr>            <dbl> <chr>     <chr>        
  # 1 L1_2           ROI:0  A_Islet_1 049x037             GTCAGTAGGTG… CD3+                     118.          19691. L1_1           0.00598 DSP-1001… T_cells      
  # 2 L1_2           ROI:0  A_Islet_1 050x036             TGTCTTGCGAC… CD3+                     292.          19691. L1_1           0.0148  DSP-1001… T_cells      
  # 3 L1_2           ROI:0  A_Islet_1 050x040             CTCAGGTCACT… CD3+                      63.3         19691. L1_1           0.00321 DSP-1001… T_cells      
  # 4 L1_2           ROI:0  A_Islet_1 050x042             CTTCGTCAAGA… CD3+                     103.          19691. L1_1           0.00522 DSP-1001… T_cells      
  # 5 L1_2           ROI:0  A_Islet_1 051x037             TACGACGCGAT… CD3+                    1120.          19691. L1_1           0.0569  DSP-1001… T_cells      
  # 6 L1_2           ROI:0  A_Islet_1 052x036             GTTCCAGTGCC… CD3+                      67.0         19691. L1_1           0.00340 DSP-1001… T_cells 
  
  sce$matched_spots <- colnames(sce) %in% spotmatch_small_subset$barcode
  sce$sample <- samples_for_registration[i]
  
  spatialCoords(sce) <- NULL
  plotSpotQC(sce, annotate = "matched_spots", x_coord = "pxl_col_in_fullres", y_coord= "pxl_row_in_fullres", 
             in_tissue = NULL, type = "spots") 
  
  CD <- as.data.frame(colData(sce))
  CD$barcode <- colnames(sce)
  
  # in_tissue array_row array_col sample_id sum_counts sum_decont   sum detected subsets_Mito_sum subsets_Mito_detected subsets_Mito_percent
  # AACACTTGGCAAGGAA-1      TRUE        47        71  sample01      23760  26841.527 23760     7521              576                    11             2.424242
  # AACAGGATTCATAGTT-1      TRUE        49        43  sample01      23033  25411.331 23033     7706              321                    12             1.393653
  # AACAGGTTATTGCACC-1      TRUE        28        86  sample01       7143   8051.770  7143     4084              133                    10             1.861963
  # AACCACTGCCATAGCC-1      TRUE        29        49  sample01       5340   5056.173  5340     3219               95                    11             1.779026
  # AACCGCCAGACTACTT-1      TRUE        39        45  sample01      18990  20506.662 18990     6638              368                    11             1.937862
  # AACCTACTGTAACTCA-1      TRUE        35        17  sample01      19287  23944.442 19287     6794              308                    11             1.596931
  # total mito_drop libsize_drop  edge patho_exclude all_drop pxl_col_in_fullres pxl_row_in_fullres cluster.init spatial.cluster
  # AACACTTGGCAAGGAA-1 23760     FALSE        FALSE FALSE         FALSE    FALSE               1072               1622            5               5
  # AACAGGATTCATAGTT-1 23033     FALSE        FALSE FALSE         FALSE    FALSE                769               1659            2               2
  # AACAGGTTATTGCACC-1  7143     FALSE        FALSE FALSE         FALSE    FALSE               1236               1265            2               7
  # AACCACTGCCATAGCC-1  5340     FALSE        FALSE FALSE         FALSE    FALSE                835               1283            3               3
  # AACCGCCAGACTACTT-1 18990     FALSE        FALSE FALSE         FALSE    FALSE                791               1471            3               3
  # AACCTACTGTAACTCA-1 19287     FALSE        FALSE FALSE         FALSE    FALSE                489               1394            3               9
  # barcode              Region matched_spots
  # AACACTTGGCAAGGAA-1 AACACTTGGCAAGGAA-1    Tumor_Stroma_mix         FALSE
  # AACAGGATTCATAGTT-1 AACAGGATTCATAGTT-1 Intratumoral_Stroma          TRUE
  # AACAGGTTATTGCACC-1 AACAGGTTATTGCACC-1 Intratumoral_Stroma         FALSE
  # AACCACTGCCATAGCC-1 AACCACTGCCATAGCC-1 Intratumoral_Vessel         FALSE
  # AACCGCCAGACTACTT-1 AACCGCCAGACTACTT-1          Tumor_pure         FALSE
  # AACCTACTGTAACTCA-1 AACCTACTGTAACTCA-1    Tumor_Stroma_mix         FALSE
  
  
  CD <- CD %>%
    left_join(spotmatch_small_subset, by = c("barcode"))
  
  colData(sce) <- as(CD, "DFrame")
  
  
  # Plotting ----------------------------------------------------------------
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  n = 5
  cols = gg_color_hue(n) # "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
  
  colData(sce) <- as(CD, "DFrame")
  sce$Cell_fraction_color <- case_when(sce$Cell_fraction == "Macro" ~ "#A3A500",
                                       sce$Cell_fraction == "Malignant" ~ "#F8766D",
                                       sce$Cell_fraction == "Other" ~ "#00BF7D",
                                       sce$Cell_fraction == "PanCK-" ~ "#00B0F6",
                                       sce$Cell_fraction == "T_cells" ~ "#E76BF3",
                                       is.na(sce$Cell_fraction) ~ "#757B82")
  
  names(cols) <- c("Malignant", "Macro", "Other" ,"PanCK-", "T_cells")
  
  p_i <- plotSpots(sce, annotate = "Cell_fraction", x_coord = "pxl_col_in_fullres", y_coord= "pxl_row_in_fullres", pt.size = 1) + 
    scale_color_manual(values = c(cols[names(table(sce$Cell_fraction))], "#757B82"))
  assign(paste0("p", i), p_i)
  
  if("Region" %in% colnames(colData(sce))){t_i <- table(sce$Region, sce$Cell_fraction)}else{t_i <- table(sce$Cell_fraction)}
  assign(paste0("t", i), t_i)
  # plotSpots(sce, annotate = "ROI.ID", x_coord = "pxl_col_in_fullres", y_coord= "pxl_row_in_fullres", pt.size = 1)
  
  # Save location -----------------------------------------------------------
  regis_savepath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration/"
  saveRDS(sce, file.path(paste0(regis_savepath, "/Visium_mapped_spot_obj"),
                         paste0(samples_for_registration[i], "_mapped_spot.rds")))
  
  geomx_mapped_subspot_LW <- CD %>%
    filter(!is.na(segment))
  write.csv(geomx_mapped_subspot_LW, 
            file.path(paste0(regis_savepath, "/GeoMx_mapped_spot_LW"), 
                      paste0(geo_reg_names[i], "_mapped_LW.csv"))) # column Sample_ID to be filtered out
  
}





