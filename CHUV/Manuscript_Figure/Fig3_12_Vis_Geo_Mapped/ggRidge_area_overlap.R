library(ggridges)
library(dplyr)
library(readxl)
library(ggplot2)

cell_fraction_order = c("Macrophage", "Malignant", "Other", "T cells", "PanCK-") # note with PanCK-
cell_fraction_color = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1", "#09D0EF")
names(cell_fraction_color) <- cell_fraction_order

final_spot <- read_xlsx("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Vis_Geo_Match/spot_aoi_matched.xlsx")
mapped_all_spot <- read.csv(file.path("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Registration", "spot_mapped_cf_region.csv"))
mapped_all_spot <- mapped_all_spot %>%
  filter(cell_fraction != "PanCK-")

# -------------------------------------------------------------
final_spot <- final_spot %>%
  filter(Cell_fraction != "PanCK-") %>%
  mutate(Cell_fraction = case_when(Cell_fraction == "Macro" ~ "Macrophage", 
                                    Cell_fraction == "T_cells" ~ "T cells",
                                    .default = as.character(Cell_fraction)))

final_spot$Cell_fraction <- factor(final_spot$Cell_fraction, levels = rev(c("Malignant", "Other", "T cells", "Macrophage")))

# final_spot_ <- final_spot %>%
#   mutate(indication = ifelse(Section_ID_Vis %in% c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2"), "Breast", "Lung"))
# final_spot_ %>% 
#   group_by(indication, Cell_fraction) %>%
#   summarise(median_area = median(overlap_pct))

final_spot_mapped <- final_spot %>%
  filter(paste0(Section_ID_Vis, Spot.barcode) %in% paste0(mapped_all_spot$Section_ID_Vis, mapped_all_spot$barcode)) 

# -------------------------------------------------------------------------
p <- ggplot(final_spot, aes(x = overlap_pct, y = Cell_fraction, fill = Cell_fraction)) +
  geom_density_ridges(scale = 4, alpha = 0.5) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + 
  theme_ridges() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y = element_text(hjust = 0),
        axis.title = element_text(size = 11)) +
  scale_fill_manual(values = cell_fraction_color) + 
  xlab("AOI & spot area overlap percentage") + 
  ylab("AOI label") + 
  geom_vline(xintercept = 0.70, linetype = "dashed", color = "darkorange", size = 1)

p <- p + annotate("text", x = 0.75, y = 5, 
                  label = "Spots mapped: \n>70% AOI overlap", 
                  color = "darkorange", 
                  size = 4, 
                  hjust = 0)


final_spotmore70 <- final_spot_mapped %>%
  filter(overlap_pct > 0.7)
sample_size_dfmore70 <- data.frame(table(final_spotmore70$Cell_fraction)) %>%
  dplyr::rename(n = Freq,
                ylabel = Var1) %>%
  mutate(n = paste0("n = ", n))


final_spotless70 <- final_spot %>%
  filter(overlap_pct <= 0.7)
sample_size_dfless70 <- data.frame(table(final_spotless70$Cell_fraction)) %>%
  dplyr::rename(n = Freq,
                ylabel = Var1) %>%
  mutate(n = paste0("n = ", n))


p_final <- p + 
  geom_text(data=sample_size_dfmore70,
              aes(label = n, x = 1, y = ylabel, vjust = -1, hjust = -0.1),
              position = position_stack(),
              inherit.aes = FALSE) +
  geom_text(data=sample_size_dfless70,
            aes(label = n, x = 0, y = ylabel, vjust = -1, hjust = -0.1),
            position = position_stack(),
            inherit.aes = FALSE)

write.csv(final_spot, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig3b_density.csv")
write.csv(sample_size_dfless70, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig3b_stats_less.csv")
write.csv(sample_size_dfmore70, 
          "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig3b_stats_more.csv")


figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig3"
pdf(file.path(figpath, "ggRidge_Area.pdf"), width = 8, height = 5)
print(p_final)
dev.off()
