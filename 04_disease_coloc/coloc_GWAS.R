library(coloc)
options(stringsAsFactors=F)
library(stringr)
library(dplyr)
library(data.table)

data<-commandArgs(trailingOnly = T)[1]
trait<-commandArgs(trailingOnly = T)[2]
chunk<-commandArgs(trailingOnly = T)[3]

subset_dataset <- function (dataset, index) 
{
    if (!length(index)) 
        return(dataset)
    if (!is.numeric(index)) 
        stop("index must be a numeric vector of indexing which snps to keep")
    d = dataset
    n = length(d$snp)
    if (n <= 1) {
        stop("not trimming length 1 dataset. check snp element exists")
    }
    if (max(index) > n) 
        stop("cannot subset to more than the number of snps in the dataset. check index")
    message("trimming dataset from ", n, " snps to ", length(index))
    for (v in names(dataset)) {
        if (is.matrix(d[[v]]) && ncol(d[[v]]) == n && nrow(d[[v]])) {
            d[[v]] = d[[v]][index, index]
            next
        }
        if (is.vector(d[[v]]) && length(d[[v]]) == n) {
            d[[v]] = d[[v]][index]
            next
        }
    }
    return(d)
}

# prep chunk file
chunkfile<-read.table(paste0("chunk/",data,"_sig_eGenes.for_GWAS_coloc_",chunk))
colnames(chunkfile)<-c("chr","gene")

# prep sample size
sampln<-read.table("exp/ROSMAP_celltype_ngenes_sample.txt")
colnames(sampln)<-c("dataset","celltype","gene","sample")
nuc_n<-sampln[sampln$celltype=="All_Cell",]$sample
bulk_n<-1092

res<-data.frame()
for(i in 1:nrow(chunkfile)){
	chr=chunkfile[i,]$chr
	gene=chunkfile[i,]$gene
  cell<-"All_Cell"
	# load GWAS coloc file
	gwas_coloc<-readRDS(paste0("GWAS/",trait,".",chr,".coloc_base.rds"))
  tryCatch({
		if(data=="ROSMAP"){
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
			snps <- intersect(nuc_ready$ID,gwas_coloc[["snp"]])
			nuc_ready <- nuc_ready[nuc_ready$ID %in% snps,]
			eqtl_coloc<-list(pvalues = nuc_ready$p, N=nuc_n, MAF= nuc_ready$maf, type = "quant", beta=nuc_ready$beta, snp=as.character(nuc_ready$ID),sdY=1)
			gwas_coloc<-subset_dataset(gwas_coloc, c(1:length(gwas_coloc[["snp"]]))[gwas_coloc[["snp"]] %in% snps])
			if(any(names(gwas_coloc)=="pvalues")){
				min_gwas_p<-min(gwas_coloc[["pvalues"]])
			}else{
				min_gwas_p<-2*pnorm(max(abs(gwas_coloc[["beta"]]/sqrt(gwas_coloc[["varbeta"]]))), lower.tail = FALSE)
			}
			}
		if(data=="ROSMAP_BulkRNA"){
			# prep bulk
			maf<-fread(paste0("wgs/NIA_JG_1898_samples_GRM_WGS_b37_",chr,".PASS.Broad_Rush.MAF0.01.afreq.gz"))
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
			snps <- intersect(bulk_ready$ID,gwas_coloc[["snp"]])
			bulk_ready <- bulk_ready[bulk_ready$ID %in% snps,]
			eqtl_coloc <- list(varbeta = (bulk_ready$se)^2, N=bulk_n, MAF= bulk_ready$maf, type = "quant", beta=bulk_ready$beta, snp=as.character(bulk_ready$ID),sdY=1)
			gwas_coloc<-subset_dataset(gwas_coloc, c(1:length(gwas_coloc[["snp"]]))[gwas_coloc[["snp"]] %in% snps])
			if(any(names(gwas_coloc)=="pvalues")){
				min_gwas_p<-min(gwas_coloc[["pvalues"]])
			}else{
				min_gwas_p<-2*pnorm(max(abs(gwas_coloc[["beta"]]/sqrt(gwas_coloc[["varbeta"]]))), lower.tail = FALSE)
			}
			}
			# coloc
			if(min_gwas_p<5e-8){
				coloc_res = coloc::coloc.abf(eqtl_coloc, gwas_coloc)
				PP.H4.abf<-coloc_res$summary["PP.H4.abf"]
				PP.H3.abf<-coloc_res$summary["PP.H3.abf"]
				nsnps<-coloc_res$summary["nsnps"]
				tmp<-data.frame(cell=cell,gene=gene,PP.H3.abf=PP.H3.abf,PP.H4.abf=PP.H4.abf,nsnps=nsnps,data=data,trait=trait)
				res<-rbind(res,tmp)
			}
			}, error = function(cond){
			tmp<-data.frame(cell=cell,gene=gene,PP.H3.abf=NA,PP.H4.abf=NA,nsnps=NA,data=data,trait=trait)
			res<-rbind(res,tmp)
			})			
  }


write.table(res,paste0("stats/",data,".chunk_",chunk,".",trait,".GWAS_coloc.txt"),sep="\t",quote=F,row=F,col=T)

