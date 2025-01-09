seurat_cluster_DE <- function(x, 
                              cluster_col="spatial.cluster", 
                              assay_name="logcounts", 
                              n_markers=10, 
                              logfc.threshold=1, 
                              palette = NULL, 
                              clusters=NULL, 
                              plot = TRUE,
                              scale = TRUE,
                              return_all_marker = TRUE,
                              ...) {
  
  if (class(x) != "Seurat"){
    ## Convert SCE to seurat object and use cluster_col as identifier
    seurat <- Seurat::CreateSeuratObject(counts=assays(x)[[assay_name]],
                                         assay='Spatial',
                                         meta.data=as.data.frame(colData(x)))
    LayerData(seurat, layer = "data", assay = "Spatial") <- LayerData(seurat, 
                                                                      layer = "counts", 
                                                                      assay = "Spatial")
  } else {
    # generate assay called Spatial
    if (assay_name %in% Assays(x)){
      x[["Spatial"]] <- x[[assay_name]]
    } else {
      x[["Spatial"]] <- x[[DefaultAssay(x)]]
    }
    
    # Trim the Seurat data
    DefaultAssay(x) <- "Spatial"
    seurat <- Seurat::DietSeurat(x, assays = "Spatial")
    
    if (!("data" %in% Layers(seurat))){
      stop("Please make sure the seurat assay of choice has 
           normalized data in 'data' layer!")
    }
  }
  
  # Set cell identity by the cluster ID
  seurat <- Seurat::SetIdent(seurat, value = cluster_col)
  
  # Check and generate color palette
  cluNum <- seurat[[cluster_col]] %>% table() %>% length()
  if (is.null(palette)){
    color <- grDevices::colors()[grep('gr(a|e)y', 
                                      grDevices::colors(), 
                                      invert = TRUE)]
    set.seed(2023)
    palette <- sample(color, cluNum, replace = FALSE)
  } 
  
  ## Subset to specified clusters
  if (!is.null(clusters)) {
    seurat <- subset(seurat, idents = clusters)
    
    # Check palette viability
    if (length(palette) < length(clusters)){
      stop(sprintf("palette number (%s) is too short for %s clusters!", 
                   length(palette), length(clusters)))
    }
    palette <- palette[clusters]
  } else if (length(palette) < cluNum){
    stop("palette is too short for the number of clusters to plot.")
  }
  
  ## Scale data
  if (scale == TRUE & !("scale.data" %in% Layers(seurat))){
    scaledData <- LayerData(seurat, layer = "data", assay = "Spatial") %>% 
      as.matrix %>% t %>% scale %>% t %>% as("sparseMatrix")
    LayerData(seurat, layer = "scale.data", assay = "Spatial") <- scaledData
  } else if (!("scale.data" %in% Layers(seurat) )){
    LayerData(seurat, layer = "scale.data", assay = "Spatial") <- LayerData(seurat, 
                                                                            layer = "data", 
                                                                            assay = "Spatial")
  }
  
  ## Select top n significant(p<0.05) markers from each cluster (by log fold change)
  all_markers <- Seurat::FindAllMarkers(seurat, assay='Spatial', slot='data',
                                        group.by=cluster_col,
                                        logfc.threshold=logfc.threshold, 
                                        only.pos=TRUE,
                                        ...) %>% 
    dplyr::filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
  
  top_markers <- all_markers %>% group_by(cluster) %>% top_n(n_markers)
  
  if (plot){
    ## Plot expression of markers
    Hm <- Seurat::DoHeatmap(seurat, 
                            assay = "Spatial",
                            features = top_markers$gene, 
                            slot = 'scale.data',
                            group.by = cluster_col, 
                            group.colors = palette, 
                            angle=0, 
                            size=4, 
                            label = TRUE, 
                            raster=FALSE) + 
      guides(col = FALSE)
    
    ## Explicitly plot the heatmap
    print(Hm)
  }
  
  if (return_all_marker){
    return(all_markers)
  } else {
    return(top_markers)
  }
}

clusterPlot_by_cluster <- function(sce_obj, type = "subspot", sample_name = "B3_2", annotate = "spatial.cluster"){
  for (j in levels(as.factor(sce_obj[[annotate]]))){
    sce_obj[[paste0("clus", j)]] <- ifelse(sce_obj[[annotate]] == j, TRUE, FALSE)
  }
  clus <- paste0("clus", levels(as.factor(sce_obj[[annotate]])))
  feat.plots <- purrr::map(clus, function(x) clusterPlot(sce_obj, label = x, palette = c("grey", "red")) +
                             scale_y_reverse() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                             ggtitle(x))
  clusters_plot <- patchwork::wrap_plots(feat.plots, ncol=4) + plot_annotation(
    title = paste0(sample_name, " - BayesSpace ", type, " clusters")
  ) & theme(plot.title = element_text(hjust = 0.5))
  return(clusters_plot)
}