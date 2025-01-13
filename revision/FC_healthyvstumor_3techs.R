library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpmisc)

path <- "~/PilotStudy"
chrom <- readRDS(file.path(path, "chrom_dlbcl.rds"))
vis <- readRDS(file.path(path, "vis_dlbcl.rds"))
geo <- readRDS(file.path(path, "geo_dlbcl.rds"))


genes <- c("MS4A1","TNFRSF13C","CD79B","CD37","CD19","BTK","TUBB","CD47","CD52","CD38","MAP2K1","CD40","BCL2L1","TNFRSF8","RARA","TYK2","TNFRSF10B","MEF2B","NOP56","MYO1G","SEMA7A","FCRL5","FCRL3","BCL6","BCL11A","BCL2L11","BCL7A")

## Chromium
chrom2 <- CreateSeuratObject(chrom@assays$SCT@data, project = "chrom", assay = "SCT", slot="data")
chrom2@meta.data$level1_5_immune_tumor_subtypes <- chrom@meta.data$level1_5_immune_tumor_subtypes
chrom2@meta.data$sample_id <- chrom@meta.data$sample_id
chrom2@meta.data$MagvsH <- ifelse(grepl("Tu_",chrom2@meta.data$level1_5_immune_tumor_subtypes),"Tumor","Healthy")
#chrom2 <- chrom2[rownames(chrom2) %in% genes,]
Idents(chrom2) <- chrom2@meta.data$sample_id 

for(sub in unique(chrom2@meta.data$sample_id)){
  chrom3 <- subset(x = chrom2, idents = sub)  
  Idents(chrom3) <- chrom3@meta.data$MagvsH 
  print(table(Idents(chrom3)))
  MvsH <- FindMarkers(chrom3, ident.2 = "Healthy", ident.1 = "Tumor",logfc.threshold =0,min.pct =0,min.cells.feature = 0,min.cells.group = 0)
  MvsH$Ranking <- seq(1:nrow(MvsH))
  write.table(MvsH,file.path(path,paste0("chrom_",sub,"_all.txt")),sep="\t",row.names=TRUE,col.names=NA)
  MvsHsub <- MvsH[rownames(MvsH) %in% genes,]
  colnames(MvsHsub) <- paste0(colnames(MvsHsub),"_",sub)
  MvsHsub <- MvsHsub[,c(2,5,6)]
  assign(paste0("chrom_",sub),MvsHsub)
  if(sub == "DLBCL_1"){
    dat <- MvsHsub
  }else{
    dat <- merge(dat,MvsHsub,by="row.names")  
    rownames(dat) <- dat[,1]
    dat[,1] <- NULL
  }
}



## Geomx
geo@meta.data$MagvsH <- ifelse(grepl("Tu_",geo@meta.data$clusters),"Tumor","Healthy")
geo@meta.data$patient <- geo@meta.data$patient
Idents(geo) <- geo@meta.data$patient 

for(sub in unique(geo@meta.data$patient)){
  geo3 <- subset(x = geo, idents = sub)  
  Idents(geo3) <- geo3@meta.data$MagvsH 
  print(table(Idents(geo3)))
  MvsH <- FindMarkers(geo3, ident.2 = "Healthy", ident.1 = "Tumor",logfc.threshold =0,min.pct =0,min.cells.feature = 0,min.cells.group = 0)
  MvsH$Ranking <- seq(1:nrow(MvsH))
  write.table(MvsH,file.path(path,paste0("geo_",sub,"_all.txt")),sep="\t",row.names=TRUE,col.names=NA)
  MvsHsub <- MvsH[rownames(MvsH) %in% genes,]
  colnames(MvsHsub) <- paste0(colnames(MvsHsub),"_",sub)
  MvsHsub <- MvsHsub[,c(2,5,6)]
  assign(paste0("geo_",sub),MvsHsub)
  if(sub == "D2"){
    dat_geo <- MvsHsub
  }else{
    dat_geo <- merge(dat_geo,MvsHsub,by="row.names")  
    rownames(dat_geo) <- dat_geo[,1]
    dat_geo[,1] <- NULL
  }
}


# Visium
vis2 <- CreateSeuratObject(vis@assays$SCT@data, project = "vis", assay = "SCT", slot="data")
vis2@meta.data$annot <- vis@meta.data$annot
vis2@meta.data$sample_id <- vis@meta.data$sample_id
vis2$MagvsH <- ifelse(grepl("Tu_",vis2@meta.data$annot),"Tumor","Healthy") 
Idents(vis2) <- vis2@meta.data$sample_id 

for(sub in unique(vis2@meta.data$sample_id)){
  vis3 <- subset(x = vis2, idents = sub)  
  Idents(vis3) <- vis3@meta.data$MagvsH 
  print(table(Idents(vis3)))
  MvsH <- FindMarkers(vis3, ident.2 = "Healthy", ident.1 = "Tumor",logfc.threshold =0,min.pct =0,min.cells.feature = 0,min.cells.group = 0)
  MvsH$Ranking <- seq(1:nrow(MvsH))
  write.table(MvsH,file.path(path,paste0("visium_",sub,"_all.txt")),sep="\t",row.names=TRUE,col.names=NA)
  MvsHsub <- MvsH[rownames(MvsH) %in% genes,]
  colnames(MvsHsub) <- paste0(colnames(MvsHsub),"_",sub)
  MvsHsub <- MvsHsub[,c(2,5,6)]
  assign(paste0("visium_",sub),MvsHsub)
  if(sub == "DLBCL_1"){
    dat_vis <- MvsHsub
  }else{
    dat_vis <- merge(dat_vis,MvsHsub,by="row.names")  
    rownames(dat_vis) <- dat_vis[,1]
    dat_vis[,1] <- NULL
  }
}
