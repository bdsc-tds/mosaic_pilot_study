disease_list = c("breast", "lung")
d = 1

# for(d in 1:length(disease_list)){
disease = disease_list[d]
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")
if(disease == "breast"){
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  samples_for_registration <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2") # Visium
  geo_reg_names <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1")
  patho_color <- breast_patho_color
}else if(disease == "lung"){
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  samples_for_registration <- c("L1_2", "L2_2", "L3_2", "L4_2") # Visium
  geo_reg_names <- c("L1_1", "L2_1", "L3_1", "L4_1")
  patho_color <- lung_patho_color
}

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
regis_savepath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Registration/"
plot_savepath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Breast_Lung_patho"


# Read in SCE -------------------------------------------------------------
sample = "B1_2"
sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
sce <- readRDS(file.path(sce_path, paste0(sample, "_qcd.rds")))

samples <- c(c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2"), c("L1_2", "L2_2", "L3_2", "L4_2"))
diseases <- c(rep("breast", 5), rep("lung", 4))
df_save_all <- data.frame(matrix(nrow = 0, ncol = 8)) 
names(df_save_all) <- c("barcode", "x_coord", "y_coord", "Region", "Level4_decon_max", 
                        "sample", "matched_spots", "cell_fraction")
for(i in 1:length(samples)){
  sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_qcd/", diseases[i], "_qcd")
  sce <- readRDS(file.path(sce_path, paste0(samples[i], "_qcd.rds")))
  sce_mapped <- readRDS(file.path(paste0(regis_savepath, "/Visium_mapped_spot_obj"), paste0(samples[i], "_mapped_spot.rds"))) # spot
  df_mapped <- data.frame(barcode = sce_mapped$barcode,
                          matched_spots = sce_mapped$matched_spots,
                          cell_fraction = sce_mapped$Cell_fraction)
  
  df_save <- data.frame(barcode = sce$Barcode,
                        x_coord = sce$array_col, 
                        y_coord = sce$array_row,
                        Region = sce$Region,
                        Level4_decon_max = sce$Level4_decon_max,
                        sample = samples[i]) %>%
    inner_join(df_mapped)
  
  df_save_all <- rbind(df_save_all, df_save)
}

write.csv(df_save_all, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS6bcd.csv")


samples <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
diseases <- c(rep("dlbcl", 6))
df_save_all <- data.frame(matrix(nrow = 0, ncol = 6)) 
names(df_save_all) <- c("barcode", "x_coord", "y_coord", "Region", "Level4_decon_max", "sample")
for(i in 1:length(samples)){
  sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/Visium_qcd/", diseases[i], "_qcd")
  sce <- readRDS(file.path(sce_path, paste0(samples[i], "_qcd.rds")))
  
  df_save <- data.frame(barcode = sce$Barcode,
                        x_coord = sce$array_col, 
                        y_coord = sce$array_row,
                        Region = sce$Region,
                        Level4_decon_max = sce$Level4_decon_max,
                        sample = samples[i])
  
  df_save_all <- rbind(df_save_all, df_save)
}

write.csv(df_save_all, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS7ab.csv")

y_rev = ifelse(sample == "B1_4", FALSE, TRUE)
# Patho plot --------------------------------------------------------------
point_size = ifelse(sample %in% c("L3_2", "L4_2", "B2_2", "B4_2"), 2,
                    ifelse(sample == "B3_2", 2.5, 3))
plot_patho <- plotSpots(sce, in_tissue = NULL, annotate = "Region",
                        x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", pt.size = point_size, y_reverse = y_rev) +
  scale_color_manual(values = patho_color[names(table(sce$Region))],
                     na.value = "#d3d3d3") + labs(color = "Pathology")

if(sample == "B1_4"){plot_patho <- plot_patho + coord_flip()}

pdf(file = file.path(plot_savepath, paste0(sample, "_patho.pdf")),
    width = 10,
    height = 5.5)
print(plot_patho)
dev.off()


