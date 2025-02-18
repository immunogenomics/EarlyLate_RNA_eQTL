library(coloc)
options(stringsAsFactors=F)
library(stringr)
library(dplyr)
library(data.table)

chunk<-commandArgs(trailingOnly = T)[1]

# prep chunk file
chunkfile<-read.table(paste0("/path/to/chunkfile.txt"))
colnames(chunkfile)<-c("chr","gene")

# prep sample size
nuc_n<-sample_n_in_nucseq

bulk_n<-1092

data="ROSMAP"

res<-data.frame()
for(i in 1:nrow(chunkfile)){
	chr=chunkfile[i,]$chr
	gene=chunkfile[i,]$gene
        cell<-"All_Cell"
        tryCatch(
        {	
	# prep nuc-seq
	maf<-fread(paste0("maf/ROSMAP.",cell,".chr",chr,".hg19.updated.eQTL.afreq.gz"))
	maf<-data.frame(maf)
	maf<-maf[,c(2,5)]
	colnames(maf)<-c("ID","af")
	maf$maf<-maf$af
	maf[maf$af>0.5,]$maf<-1-maf[maf$af>0.5,]$af
	nuc<-fread(paste0("fastqtl/",data,".",cell,".PEER_norm.chr",chr,".nominal.txt.gz"))
	nuc<-data.frame(nuc)
	colnames(nuc)<-c("gene","ID","dum","p","beta")
	nuc<-nuc[nuc$gene==gene,]
	nuc<-left_join(nuc,maf,by="ID")
	nuc<-nuc[nuc$maf>0,]
	nuc<-nuc %>% group_by(ID) %>% slice(which.min(p))
  nuc_ready<-data.frame(nuc)
	nuc<-NULL

	# prep bulk
	maf<-fread(paste0("wgs/qced/NIA_JG_1898_samples_GRM_WGS_b37_",chr,".PASS.Broad_Rush.MAF0.01.afreq.gz"))
	maf<-data.frame(maf)
	maf<-maf[,c(2,5)]
	colnames(maf)<-c("ID","af")
	maf$maf<-maf$af
	maf[maf$af>0.5,]$maf<-1-maf[maf$af>0.5,]$af
	bulk<-fread(paste0("brain_bulk_eQTL/eQTLchr",chr,"_1Mb.csv.gz"))
	bulk<-data.frame(bulk)
  bulk<-bulk[bulk$geneSym==gene,]
	conv<-fread(paste0("brain_bulk_eQTL/rsID_wgsID_chr",chr,".txt.gz"))
	bulk<-merge(bulk,conv)
	bulk<-left_join(bulk,maf,by="ID")
	bulk<-bulk[bulk$maf>0,]
  bulk<-bulk %>% group_by(ID) %>% slice(which.min(p))
	bulk$alt<-str_split_fixed(bulk$ID,":",4)[,4]
	bulk[bulk$A1!=bulk$alt,]$beta <- (-1)*bulk[bulk$A1!=bulk$alt,]$beta # adjust beta for alt allele
  bulk_ready<-data.frame(bulk)
	bulk<-NULL
	
	# set shared SNPs
	snps <- intersect(nuc_ready$ID,bulk_ready$ID)
	nuc_ready <- nuc_ready[nuc_ready$ID %in% snps,]
	bulk_ready <- bulk_ready[bulk_ready$ID %in% snps,]
	
	# input for coloc
	nuc_coloc<-list(pvalues = nuc_ready$p, N=nuc_n, MAF= nuc_ready$maf, type = "quant", beta=nuc_ready$beta, snp=as.character(nuc_ready$ID),sdY=1)
	bulk_coloc <- list(varbeta = (bulk_ready$se)^2, N=bulk_n, MAF= bulk_ready$maf, type = "quant", beta=bulk_ready$beta, snp=as.character(bulk_ready$ID),sdY=1)
	
	# coloc
	coloc_res = coloc::coloc.abf(nuc_coloc, bulk_coloc)
	PP.H4.abf<-coloc_res$summary["PP.H4.abf"]
	PP.H3.abf<-coloc_res$summary["PP.H3.abf"]
	nsnps<-coloc_res$summary["nsnps"]
	tmp<-data.frame(cell=cell,gene=gene,PP.H3.abf=PP.H3.abf,PP.H4.abf=PP.H4.abf,nsnps=nsnps)
  }, error = function(cond){tmp<-data.frame(cell=cell,gene=gene,PP.H3.abf=NA,PP.H4.abf=NA,nsnps=NA)})
	res<-rbind(res,tmp)
  }

write.table(res,paste0("stats/ROSMAP.chunk_",chunk,".nuc_bulk.coloc.txt"),sep="\t",quote=F,row=F,col=T)

