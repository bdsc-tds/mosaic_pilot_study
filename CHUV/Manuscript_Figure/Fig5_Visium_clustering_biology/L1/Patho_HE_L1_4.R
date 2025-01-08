library(ggspavis)

disease = "lung"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Visium/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/color_palette.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Manuscript_Figure/Fig4_Vis_patho_decon_gallery/alternative/Patho_hires_HE_helper.R")
sample = "L1_4"
if(disease == "lung"){
  patho_color <- lung_patho_color
}else if(disease == "breast"){
  patho_color <- breast_patho_color
}else{
  patho_color <- dlbcl_patho_color
}

sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
sce <- readRDS(file.path(sce_path, paste0(sample, "_qcd.rds")))

datapath_rerun <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/VisiumRerun/", sample)
coord_new <- read.csv(file.path(datapath_rerun, "tissue_positions.csv"))

CD <- data.frame(colData(sce))

CD_new <- coord_new %>% 
  mutate(Barcode = barcode) %>% 
  select(Barcode, pxl_col_in_fullres, pxl_row_in_fullres)

CD <- CD %>%
  left_join(CD_new, by = "Barcode")
rownames(CD) <- CD$Barcode

colData(sce) <- as(CD, "DFrame")


sce_ <- sce
colors <- patho_color[names(table(sce_$Region))]

scalef <- rjson::fromJSON(file = file.path(datapath_rerun, "scalefactors_json.json"))

# imgData(sce_)

sce_ <- rmvImg(sce_)
sce_ <- addImg(sce_, 
               sample_id = sample, 
               image_id = "hires",
               imageSource = file.path(datapath_rerun, "tissue_hires_image.png"), 
               scaleFactor = scalef$tissue_hires_scalef, 
               load = TRUE)
# img <- imgRaster(sce_, 
#                  sample_id = "L4_2", 
#                  image_id = "hires")
# plot(img)
# imgData(sce_)

sce_$Region <- as.factor(sce_$Region)
p2 <- plotVisium_new(
  sce_, 
  annotate = "Region", 
  palette = patho_color[names(table(sce_$Region))], 
  zoom = TRUE, 
  pt.size = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2.5
)

# Image size and save -----------------------------------------------------
image_height <- (range(spatialCoords(sce_)[, 2])[2] - range(spatialCoords(sce_)[, 2])[1])
image_width <- (range(spatialCoords(sce_)[, 1])[2] - range(spatialCoords(sce_)[, 1])[1])


plot_width <- 15
plot_height <- plot_width * (image_height/image_width) - 2

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5"
pdf(file = file.path(paste0(fig_path, "/L1/"),
                     paste0(sample, "_patho_HE.pdf")),
    width = plot_width,
    height = plot_height)
print(p2)
dev.off()

