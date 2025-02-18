options(stringsAsFactors=F)
library(stringr)
library(dplyr)
library(data.table)

res<-data.frame()
for(data in c("ROSMAP_BulkRNA","ROSMAP")){
  all<-fread(paste0("stats/",data,".FDRsig_eGenes.snps.bed.gz"))
  base<-all %>% group_by(V6) %>% summarise(allPIP=sum(V5))
  colnames(base)[1]<-"gene"
  
  for(categ in c("transcript","3UTR","5UTR","promoter","braincCRE","braincCRE_wo_TSS100kb","braincCRE_wo_TSS250kb","psychEncode_cCRE")){
    tmp<-fread(paste0("stats/",data,".FDRsig_eGenes.snps_in_",categ,".withPIP.txt.gz"))
    tmp<-unique(tmp)
    tmp<-tmp %>% group_by(V3) %>% summarise(sumPIP=sum(V2))
    colnames(tmp)[1]<-"gene"
    tmp<-left_join(base,tmp,by="gene")
    tmp$sumPIP[is.na(tmp$sumPIP)]<-0
    tmp$adjPIP<-tmp$sumPIP/tmp$allPIP
    tmp$data=data
    tmp$categ=categ
    res<-rbind(res,tmp)
  }
  
  transcript<-fread(paste0("stats/",data,".snps_in_transcript.withPIP.txt.gz"))
  transcript<-unique(transcript)
  transcript_snp_gene<-unique(paste(transcript$V1,transcript$V3,sep="_"))

  for(categ in c("exon","m6A","ALLeCLIP")){
    tmp<-fread(paste0("stats/",data,".FDRsig_eGenes.snps_in_",categ,".withPIP.txt.gz"))
    tmp<-unique(tmp)
    tmp$snp_gene<-paste(tmp$V1,tmp$V3,sep="_")
    tmp<-tmp[tmp$snp_gene %in% transcript_snp_gene,]
    tmp<-tmp %>% group_by(V3) %>% summarise(sumPIP=sum(V2))
    colnames(tmp)[1]<-"gene"
    tmp<-left_join(base,tmp,by="gene")
    tmp$sumPIP[is.na(tmp$sumPIP)]<-0
    tmp$adjPIP<-tmp$sumPIP/tmp$allPIP
    tmp$data=data
    tmp$categ=categ
    res<-rbind(res,tmp)
  }
}

res<-data.frame(res)

stats<-data.frame()
for(categ in all_tested_categories){
  tmp.1<-data.frame(res[res$data=="ROSMAP_BulkRNA"&res$categ==categ,])
  rownames(tmp.1)<-tmp.1$gene
  tmp.2<-data.frame(res[res$data=="ROSMAP"&res$categ==categ,])
  rownames(tmp.2)<-tmp.2$gene
  gene_assess<-intersect(tmp.1$gene,tmp.2$gene)
  tmp.1<-tmp.1[gene_assess,]
  tmp.2<-tmp.2[gene_assess,]
  ttest<-t.test(tmp.1$adjPIP,tmp.2$adjPIP,paired=T)
  tmp<-data.frame(data="ROSMAP",categ=categ,t=ttest$statistic[[1]],diff=ttest$estimate[[1]],lower=ttest$conf.int[[1]],higher=ttest$conf.int[[2]],p=ttest$p.value)
  stats<-rbind(stats,tmp)}

write.table(stats,paste0("summary/ROSMAP.adjPIP_enrich.summary.txt"),sep="\t",quote=F,row=F)
