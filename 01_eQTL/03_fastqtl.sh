#!/bin/bash

cell="$1"
CHR="$2"

module load samtools
module load bcftools
module load python/3.8.2

orgvcf="path/to/genotype.vcf.gz"

bcftools view -S ID/${cell}.exp_bed.vcfID.order -Oz -o vcf/${cell}.chr${CHR}.hg19.updated.eQTL.vcf.gz ${orgvcf}
tabix -f -p vcf vcf/${cell}.chr${CHR}.hg19.updated.eQTL.vcf.gz

/path/to/fastQTL.static \
      --vcf  vcf/${cell}.chr${CHR}.hg19.updated.eQTL.vcf.gz \
      --bed  exp/${cell}.PEER_norm.vcfID.gene.bed.gz \
      --region  ${CHR} \
      --permute 10000 \
      --out  fastqtl/${cell}_v2.PEER_norm.chr${CHR}.perm.txt

/path/to/fastQTL.static \
      --vcf  vcf/${cell}.chr${CHR}.hg19.updated.eQTL.vcf.gz \
      --bed  exp/${cell}.PEER_norm.vcfID.gene.bed.gz \
      --region  ${CHR} \
      --out  fastqtl/${cell}_v2.PEER_norm.chr${CHR}.nominal.txt

gzip -f fastqtl/${cell}_v2.PEER_norm.chr${CHR}.nominal.txt

