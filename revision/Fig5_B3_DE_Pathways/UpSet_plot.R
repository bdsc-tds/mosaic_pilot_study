read_path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/Geo_Vis_Chrom/DiffDE"
vis_DE <- read.csv(file.path(read_path, "1&5&9_14.csv"))
geo_DE <- read.csv(file.path(read_path, "Geo_3vs6_DE.csv"))
chrom_DE <- read.csv(file.path(read_path, "chrom_B3_tu_2_clus_markers.csv"))

vis_blue_gene <- vis_DE$gene[vis_DE$cluster == "14"]
vis_pink_gene <- vis_DE$gene[vis_DE$cluster == "1_5_9"]

geo_blue_gene <- geo_DE$TargetName[geo_DE$logFC < 0]
geo_pink_gene <- geo_DE$TargetName[geo_DE$logFC > 0]

chrom_blue_gene <- chrom_DE$gene[chrom_DE$cluster == "Tu_B3_PLA2G2A"]
chrom_pink_gene <- chrom_DE$gene[chrom_DE$cluster == "Tu_B3_NPPC"]


blue_gene_list = list(Visium = vis_blue_gene, GeoMx = geo_blue_gene, Chromium = chrom_blue_gene)
pink_gene_list = list(Visium = vis_pink_gene, GeoMx = geo_pink_gene, Chromium = chrom_pink_gene)


# To label in volcano -----------------------------------------------------
intersect(vis_blue_gene, geo_blue_gene)
blue_vis_geo <- c("SULT1C3", "PEG10", "PLA2G2A", "AZGP1", "MUCL1", "NR2F1", "HIPK2", "ACADM", "OLFML3")
intersect(vis_pink_gene, geo_pink_gene)
pink_vis_geo <- c("NPPC", "AKR1C2", "GPRC5A", "GSTP1", "PGC", "LBP", "THRSP", "S100P", "ST6GALNAC4", "BST2", "CRIP1", 
                  "ACSL3", "CAPS", "FAM234B", "AP1S3", "SH3BGRL", "PKIB", "ATP2A3", "SGK1", "MBOAT2", "FNIP2")

get_m_ss_cs <- function(list_input, mode){
  m <- make_comb_mat(list_input, mode = mode)
  ss <- set_size(m)
  cs <- comb_size(m)
  
  return(list(m, ss, cs))
}

blue_dis_list <- get_m_ss_cs(blue_gene_list, mode = "distinct")
blue_int_list <- get_m_ss_cs(blue_gene_list, mode = "intersect")
pink_dis_list <- get_m_ss_cs(pink_gene_list, mode = "distinct")
pink_int_list <- get_m_ss_cs(pink_gene_list, mode = "intersect")


saveUpSetPlot <- function(list, bar_col, dir, name){
  m = list[1][[1]]; ss = list[2][[1]]; cs = list[3][[1]]
  
  ht <- UpSet(m, 
              set_order = colnames(m), 
              comb_order = order(comb_degree(m)),
              top_annotation = HeatmapAnnotation(
                "Distinct diff. genes" = anno_barplot(
                  cs, 
                  ylim = c(0, max(cs)*1.1),
                  border = FALSE, 
                  gp = gpar(fill = bar_col), 
                  height = unit(4, "cm")
                ), 
                annotation_name_side = "left", 
                annotation_name_rot = 90),
              right_annotation = HeatmapAnnotation(
                which = "row",
                "Total" = anno_barplot(
                  ss, 
                  ylim = c(0, max(ss)*1.1),
                  border = FALSE, 
                  gp = gpar(fill = bar_col), 
                  width = unit(4, "cm")
                )
              )
  )
  
  pdf(file = file.path(dir, name), width = 10, height = 4)
  ht = draw(ht)
  od = column_order(ht)
  rod = row_order(ht)
  decorate_annotation("Distinct diff. genes", {
    grid.text(cs[od], 
              x = seq_along(cs), 
              y = unit(cs[od], "native") + unit(2, "pt"), 
              default.units = "native", just = c("left", "bottom"), 
              gp = gpar(fontsize = 8, col = "#404040"), rot = 45)
  })
  decorate_annotation("Total", {
    grid.text(ss[rod], 
              x = unit(ss[rod], "native") + unit(20, "pt"), 
              y = rev(seq_along(ss)), 
              default.units = "native", just = c("right", "bottom"), 
              gp = gpar(fontsize = 8, col = "#404040"))
  })
  
  dev.off()
}

save_path = "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Manuscript_revision/Fig5_Geo_B3/UpSetPlot"
saveUpSetPlot(blue_dis_list, bar_col = "#0000ff", dir = save_path, name = "blue_dis.pdf")
saveUpSetPlot(blue_int_list, bar_col = "#9393ff", dir = save_path, name = "blue_int.pdf")
saveUpSetPlot(pink_dis_list, bar_col = "#ff00ff", dir = save_path, name = "pink_dis.pdf")
saveUpSetPlot(pink_int_list, bar_col = "#fc98fc", dir = save_path, name = "pink_int.pdf")


distinct_df <- data.frame(
  val = c(c(39, 3, 2, 4, 2, 0, 5), c(56, 5, 10, 13, 13, 0, 8)),
  group = c(rep("blue", 7), rep("pink", 7)),
  comb = rep(c("Visium", "GeoMx", "Chromium", "vis_geo", "vis_chrom", "geo_chrom", "vis_geo_chrom"), 2),
  type = rep("total", 14)
)

intersect_df <- data.frame(
  val = c(c(50, 12, 9, 9, 7, 5, 5), c(90, 26, 31, 21, 21, 8, 8)),
  group = c(rep("blue", 7), rep("pink", 7)),
  comb = rep(c("Visium", "GeoMx", "Chromium", "vis_geo", "vis_chrom", "geo_chrom", "vis_geo_chrom"), 2),
  type = rep("total", 14)
)

total_df <- data.frame(
  val = c(50, 12, 9, 90, 26, 31),
  group = c(rep("blue", 3), rep("pink", 3)),
  comb = rep(c("Visium", "GeoMx", "Chromium"), 2),
  type = rep("total", 6)
)

all_df <- rbind(distinct_df, intersect_df, total_df) 
write.csv(all_df, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/Fig5f.csv")
# %>%
#   mutate(group = factor(group, levels = c("blue", "pink")),
#          comb = factor(group, levels = c("Visium", "GeoMx", "Chromium", "vis_geo", "vis_chrom", "geo_chrom", "vis_geo_chrom")),
#          type = factor(group, levels = c("total", "distinct", "intersect")))


# side by side bar plot ---------------------------------------------------
p <- ggplot(all_df, aes(x=comb, y=val, fill=group)) +
  geom_bar(stat='identity', 
                   position='dodge2',
                   width = 0.8,
                   color = "black") + 
  geom_text(aes(label = val),  
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 4) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_fill_manual(values = c("blue" = "#0000ff", "pink" = "#ff00ff"))
p









