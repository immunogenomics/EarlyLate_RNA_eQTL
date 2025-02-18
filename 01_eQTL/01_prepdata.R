library(data.table)
meta<-read.csv("cell-annotation.n424.csv.gz",header=T)
meta$cell.type<-gsub(" ","_",meta$cell.type)

library(Seurat)
library(SeuratDisk)
library(SeuratData)
returndata<-function(celltype,filename){
  cell.in<-as.character(meta[meta$cell.type==celltype,]$barcode)
  meta2<-meta[meta$cell.type==celltype,]
  d<-LoadH5Seurat(paste0("../../data/ROSMAP/object/",filename,".h5Seurat"), assays = "counts")
  d<-NormalizeData(d,
        normalization.method = "LogNormalize",
        scale.factor = 10000)
  d<-d[["SCT"]]@data
  d<-d[,colnames(d)%in%cell.in]
  return(d)
}

df1<-returndata(celltype="Endothelial",filename="vascular.niche")
df2<-returndata(celltype="Astrocyte",filename="astrocytes")
df3<-returndata(celltype="Microglia",filename="microglia")
df4<-returndata(celltype="Inhibitory_Neurons",filename="inhibitory")
df5<-returndata(celltype="Oligodendrocytes",filename="oligodendroglia")
df6<-returndata(celltype="OPCs",filename="oligodendroglia")
df7<-returndata(celltype="Excitatory_Neurons",filename="cux2-")
df8<-returndata(celltype="Excitatory_Neurons",filename="cux2+")

genes.intersect<-Reduce(intersect,list(rownames(df1),rownames(df2),rownames(df3),rownames(df4),rownames(df5),rownames(df6),rownames(df7),rownames(df8)))

df1<-df1[genes.intersect,]
df2<-df2[genes.intersect,]
df3<-df3[genes.intersect,]
df4<-df4[genes.intersect,]
df5<-df5[genes.intersect,]
df6<-df6[genes.intersect,]
df7<-df7[genes.intersect,]
df8<-df8[genes.intersect,]

allgeneids <- genes.intersect
allcells <- c(colnames(df1),colnames(df2),colnames(df3),colnames(df4),colnames(df5),colnames(df6),colnames(df7),colnames(df8))

out <- matrix(0, nrow=length(allgeneids), ncol=length(unique(meta$individualID)) )
rownames(out) <- allgeneids
colnames(out) <- unique(meta$individualID)

for( ind in unique( meta$individualID ) ){
  cells <- subset(meta, individualID==ind )$barcode # cells from target donor
  cells <- intersect(cells, allcells)
  if( length(cells) > 10 ){ 
    d3.1 <- df1[, colnames(df1) %in% cells]
    d3.2 <- df2[, colnames(df2) %in% cells]
    d3.3 <- df3[, colnames(df3) %in% cells]
    d3.4 <- df4[, colnames(df4) %in% cells]
    d3.5 <- df5[, colnames(df5) %in% cells]
    d3.6 <- df6[, colnames(df6) %in% cells]
    d3.7 <- df7[, colnames(df7) %in% cells]
    d3.8 <- df8[, colnames(df8) %in% cells]
    
    d3 <- do.call(cbind,list(d3.1,d3.2,d3.3,d3.4,d3.5,d3.6,d3.7,d3.8))
    d3 <- as.matrix(d3)
    rowaverage <- rowSums(d3) / ncol(d3)
    out[, ind] <- rowaverage
  }
}

out <- data.frame( gene=row.names(out), out )

saveRDS(out, paste0("ROSMAP.All_Cell.norm_av.rds"))
