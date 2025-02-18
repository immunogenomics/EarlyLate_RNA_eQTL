library(coloc)
options(stringsAsFactors=F)
library(stringr)
library(dplyr)
library(data.table)

cell<-commandArgs(trailingOnly = T)[1]
chunk<-commandArgs(trailingOnly = T)[2]


n<-sample_n # number of samples

chunkfile<-read.table("/path/to/chunkfile.txt")
colnames(chunkfile)<-c("chr","gene")

data="ROSMAP"

res<-data.frame()
for(i in 1:nrow(chunkfile)){
	chr=chunkfile[i,]$chr
	gene=chunkfile[i,]$gene
	maf<-read.table(paste0("maf/ROSMAP.",cell,".chr",chr,".hg19.updated.eQTL.afreq.gz"))
	maf<-maf[,c(2,5)]
	colnames(maf)<-c("V2","af")
	maf$maf<-maf$af
	maf[maf$af>0.5,]$maf<-1-maf[maf$af>0.5,]$af
	nuc<-fread(paste0("fastqtl/",data,".",cell,".PEER_norm.chr",chr,".nominal.txt.gz"))
	nuc<-data.frame(nuc)
	nuc<-nuc[nuc$V1==gene,]
	nuc<-left_join(nuc,maf,by="V2")
	nuc<-nuc[nuc$maf>0,]
  nuc_coloc<-list(pvalues = nuc$V4, N=n, MAF= nuc$maf, type = "quant", beta=nuc$V5, snp=as.character(nuc$V2),sdY=1)
	nuc_res<- coloc::finemap.abf(dataset=nuc_coloc)
	nuc_res$gene<-gene
	nuc_res<-nuc_res[-(nrow(nuc_res)),]
	write.table(nuc_res,paste0("comparison_sc_nuc/stats/",data,".",gene,".finemap.all.txt"),sep="\t",quote=F,row=F,col=F)
  }
