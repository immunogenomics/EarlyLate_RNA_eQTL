library(peer)

dataset = commandArgs(trailingOnly=TRUE)[1]
celltype = commandArgs(trailingOnly=TRUE)[2]

meta<-read.table("meta_info.n424.selected.txt",header=T,sep="\t")
colnames(meta)[19]<-"RID"

pcs <- read.table("pca.n424.evec")
pcs <- pcs[,c(-1,-13:-23)]
colnames(pcs)<-c("IID",paste0("PC",1:10))
pcs <- pcs[, 1:6 ]

cov<-merge(pcs,meta,by="IID")

# expression data
d1 <- read.csv(paste0("exp/",dataset,".",celltype,".norm_av.csv"))
d1 <- d1[,-1]

rownames(d1)<-d1$gene
d2 <- as.matrix(d1[ ,c(-1)])

cov <- cov[cov$RID %in% colnames(d2),]
row.names(cov) <- as.character(cov$RID)

# common sample
commonsample <- intersect( cov$RID, colnames(d2) )

# filtering sample and genes
d2 <- d2[, commonsample] #same order of samples
d2[is.na(d2)] <- 0 #missing values were set to 0

exprrate <- apply(d2, 1, function(x){ length( x [ x > 0 ] ) / length(x) })
d3 <- d2[ exprrate > 0.2,  ]

# inverse normal normalization
rn<-apply(d3, 1, function(x){
   qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) )
})
rn <- t(rn)

# peer
QQ <- rn

model = PEER()
PEER_setPhenoMean(model, as.matrix(t(QQ)) )
dim(PEER_getPhenoMean(model))

PEER_setAdd_mean(model, TRUE)

cov2 <- cov[commonsample, ]
cov2 <- cov2[,c("PC1","PC2","PC3","PC4","PC5","Study","msex","age_death","pmi")]
PEER_setCovariates(model, as.matrix(cov2))

K = 10
PEER_setNk(model,K)
PEER_getNk(model)
PEER_setNmax_iterations(model,10000)

# perform the inference
PEER_update(model)

# output
residuals = PEER_getResiduals(model)
colnames(residuals)  <- row.names(QQ)
row.names(residuals) <- colnames(QQ)
out <- data.frame( ID = row.names(QQ),t(residuals))
write.table(out, paste0("exp/",dataset,".",celltype,".PEER_norm.tsv"), sep = "\t", quote = F, row.names = FALSE)

# when you want to get PEER factors
factors = PEER_getX(model)
saveRDS(factors,paste0("exp/",dataset,".",celltype,".PEER_factor.rds"))

row.names(factors) <- colnames(QQ)

write.table(factors, paste0("exp/",dataset,".",celltype,".PEER_factor.tsv"), sep = "\t", quote = F)

