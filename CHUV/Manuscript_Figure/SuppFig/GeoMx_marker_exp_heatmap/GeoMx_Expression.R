geopath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched"

disease = "breast"
disease = "lung"
disease = "dlbcl"

get_plt_df <- function(disease){
  geo <- readRDS(file.path(geopath, paste0(disease, "_spe.rds")))
  
  logcounts <- assay(geo, "logcounts")
  table(geo$segment)
  if(disease %in% c("breast", "lung")){
    plt_df <- data.frame(
      cf = geo$cell_fraction,
      EPCAM = logcounts["EPCAM", ],
      KRT17 = logcounts["KRT17", ],
      PECAM1 = logcounts["PECAM1", ],
      # CD79A = logcounts["CD79A", ],
      CD68 = logcounts["CD68", ],
      CD3D = logcounts["CD3D", ]
    )
    
    plt_df_ <- plt_df %>%
      gather(CT, normexp, -cf)
    
  }else{
    plt_df <- data.frame(
      cf = geo$cell_fraction,
      Patient = geo$patient,
      PTPRC = logcounts["PTPRC", ],
      CD3D = logcounts["CD3D", ],
      CD4 = logcounts["CD4", ],
      CD8A = logcounts["CD8A", ],
      CD19 = logcounts["CD19", ],
      CD68 = logcounts["CD68", ],
      CD79A = logcounts["CD79A", ]
    )
    
    plt_df_ <- plt_df %>%
      gather(CT, normexp, -c(cf, Patient))
  }
  return(plt_df_)
}


# Breast Lung -------------------------------------------------------------
br <- get_plt_df(disease = "breast")
lu <- get_plt_df(disease = "lung")

br_lu <- rbind(br %>% mutate(Indication = "Breast"), 
               lu %>% mutate(Indication = "Lung"))
br_lu$CT <- as.factor(br_lu$CT)
br_lu <- br_lu %>% mutate(cf = ifelse(cf == "Macro", "Macrophage", cf))

write.csv(br_lu, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS5a.csv") 

bp <- ggplot(br_lu, aes(x=cf, y=normexp)) + 
  geom_violin(trim=TRUE, alpha = 0.5) + 
  geom_jitter(aes(color=Indication), size=1, alpha=0.9) +
  scale_color_manual(values = c("Breast" = "#00c78cff", "Lung" = "#ffa54fff")) + 
  scale_y_continuous(labels = function(x) round(x, 0)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(0.8, "lines")) + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  facet_wrap(~ CT, ncol = 7, scales = "free") +
  xlab("") + ylab("Normalized expression") 
bp

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/GeoMx_Expr_Heatmap"
pdf(file = file.path(fig_path, "Expr_br_lu.pdf"),
    width = 13,
    height = 4)
print(bp)
dev.off()

# DLBCL -------------------------------------------------------------------
dlbcl <- get_plt_df(disease = "dlbcl") %>% mutate(cf = ifelse(cf == "Macro", "Macrophage", cf))

write.csv(dlbcl, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppFigS5b.csv") 

# Box with jitter ---------------------------------------------------------
bp <- ggplot(dlbcl, aes(x=cf, y=normexp)) + 
  geom_violin(trim=TRUE, alpha = 0.5) + 
  geom_jitter(aes(color=Patient), size=1, alpha=0.9) +
  scale_y_continuous(labels = function(x) round(x, 0)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(0.8, "lines")) + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  facet_wrap(~ CT, ncol = 7, scales = "free") +
  xlab("") + ylab("Normalized expression") 
bp

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/GeoMx_Expr_Heatmap"
pdf(file = file.path(fig_path, "Expr_dlbcl.pdf"),
    width = 18,
    height = 4.3)
print(bp)
dev.off()

