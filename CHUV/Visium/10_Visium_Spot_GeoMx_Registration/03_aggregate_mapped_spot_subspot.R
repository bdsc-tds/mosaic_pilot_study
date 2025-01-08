# Only breast and lung ----------------------------------------------------

############################################################################
#                                   Spot                                   #
############################################################################

# Mapped spot and AOI IDs (spot)
mapped_geo_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration/GeoMx_mapped_spot_LW"
geo_samples <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1", "L1_1", "L2_1", "L3_3", "L4_1")
vis_samples <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2", "L1_2", "L2_2", "L3_2", "L4_2")
disease_list <- c(rep("breast", 5), rep("lung", 4))

mapped_all <- NULL
for(i in 1:length(geo_samples)){
  mapped_i <- read.csv(file.path(mapped_geo_path, paste0(geo_samples[i], "_mapped_LW.csv")), row.names = 1) %>% 
    select(barcode, Section_ID, Sample_ID, Cell_fraction) %>%
    dplyr::rename(cell_fraction = Cell_fraction) %>%
    mutate(Section_ID_Vis = vis_samples[i])
  
  # Merge in pathology ------------------------------------------------------
  sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease_list[i], "_qcd")
  sce <- readRDS(file.path(sce_path, paste0(vis_samples[i], "_qcd.rds")))
  
  mapped_i <- mapped_i %>%
    left_join(data.frame(barcode = colnames(sce),
                         Region = sce$Region))
  
  mapped_all <- rbind(mapped_all, mapped_i)
}

mapped_all <- mapped_all %>% # 10 Visium spots from B1_2 patho exclude
  filter(!is.na(Region)) 

write.csv(mapped_all, 
          file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration", "spot_mapped_cf_region.csv"), row.names = FALSE)


