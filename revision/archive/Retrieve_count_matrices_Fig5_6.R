chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

vispath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/dlbcl/spotclean/Results/"
vis <- readRDS(file.path(vispath, "Dlbcl-merge-SCTpostSpotClean.rds"))

geopath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched"
geo <- readRDS(file.path(geopath, "dlbcl_seu_ruv.rds"))

# In dotplot
chrom@assays$SCT@data # log1p of SCT normalized decontaminated (soupx) matrix
vis@assays$SCT@data #  log1p of SCT normalized decontaminated (spotclean) matrix
geo@assays$originalexp@data # quantile normalized on log1p, RUV4 batch corrected matrix

# Raw data
chrom@assays$RNA@counts # raw

vis@assays$Spatial@layers$counts.1 #  raw
vis@assays$Spatial@layers$counts.2 #  raw
vis@assays$Spatial@layers$counts.3 #  raw
vis@assays$Spatial@layers$counts.4 #  raw
vis@assays$Spatial@layers$counts.5 #  raw
vis@assays$Spatial@layers$counts.6 #  raw

geo@assays$originalexp@counts # raw
