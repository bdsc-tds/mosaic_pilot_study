library(tidyr)
library(dplyr)
library(SpatialExperiment)
# # Sanity (PanCK- samples included) -----------------------------------------
# breast_geo_level4_decon <- read.csv(file.path(deconresultpath, "breast_batched_decon_long.csv"))
# lung_geo_level4_decon <- read.csv(file.path(deconresultpath, "lung_batched_decon_long.csv"))
# dlbcl_geo_level4_decon <- read.csv(file.path(deconresultpath, "dlbcl_batched_decon_long.csv"))
# length(unique(breast_geo_level4_decon$sample))
# length(unique(lung_geo_level4_decon$sample))
# length(unique(dlbcl_geo_level4_decon$sample))

# -------------------------------------------------------------------------
# disease = "breast"
# disease = "lung"
disease = "dlbcl"

source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/00_GeoMx_Paths.R")

datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
spe_ruv <- readRDS(file.path(datapath, paste0(disease, "_spe_ruv.rds")))

deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/For_level1_5_immune_decon_results"

# -------------------------------------------------------------------------
get_Decon_DF_to_Merge <- function(disease){
  geo_decon <- read.csv(file.path(deconresultpath, 
                                          paste0(disease, "_batched_decon_long.csv"))) %>%
    filter(CellType %in% c("T cells", "Macrophage")) %>%
    select(sample, CellType, Fraction) %>%
    pivot_wider(names_from = CellType, values_from = Fraction)
  
  return(geo_decon)
}

geo_decon <- get_Decon_DF_to_Merge(disease = disease)

CD <- as.data.frame(colData(spe_ruv)) %>%
  select(sample_id2) %>% 
  dplyr::rename(sample = sample_id2)

CD <- CD %>%
  left_join(geo_decon)

spe_ruv$Macrophage_frac <- CD$Macrophage
spe_ruv$Tcells_frac <- CD$`T cells`

# For breast and lung, NA is for the PanCK- segments, set to 0
spe_ruv$Macrophage_frac <- ifelse(is.na(spe_ruv$Macrophage_frac), 0, spe_ruv$Macrophage_frac)
spe_ruv$Tcells_frac <- ifelse(is.na(spe_ruv$Tcells_frac), 0, spe_ruv$Tcells_frac)

which(is.na(spe_ruv$Macrophage_frac))

saveRDS(spe_ruv, file.path(datapath, paste0(disease, "_spe_ruv.rds")))






