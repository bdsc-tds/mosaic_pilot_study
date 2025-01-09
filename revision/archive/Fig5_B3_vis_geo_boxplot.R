# Read decon results of Vis and Geo
decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/Visium_Decon/Level4/C2L"
vis_level4breast_B3 <- read.csv(file.path(decon_path, paste0("B3_2", ".csv")), row.names = 1)
vis_level4breast_B3 <- vis_level4breast_B3 %>%
  pivot_longer(cols = everything(),
               names_to = "CellType",
               values_to = "Fraction") %>%
  mutate(CellType = ifelse(CellType == "Tu_B3_CYP4F8", "Tu_B3", CellType))
# TODO: add column of area AB to both geo and vis

ggplot(data, aes(x = Group, y = Value, fill = Subgroup)) +
  geom_boxplot() +
  labs(title = "Grouped Box Plot", x = "Group", y = "Value") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1")

## Level 4 GeoMx
geo_level4breast <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific/breast_batched_decon_long.csv")
geo_level4breast_B3 <- geo_level4breast %>%
  filter(section_id == "B3_1") %>%
  mutate(CellType = ifelse(CellType == "Tu_B3_CYP4F8", "Tu_B3", CellType))



