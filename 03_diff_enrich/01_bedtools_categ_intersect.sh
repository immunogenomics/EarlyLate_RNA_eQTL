#!/bin/bash
data="$1"
categ="$2"

bedtools intersect -a stats/${data}.FDRsig_eGenes.snps.bed.gz -b /path/to/hg19.refGene.${categ}_per_gene.bed.gz -wa -wb | awk '{OFS="\t"}{if($6==$10)print $4,$5,$6}' | uniq | gzip -c > stats/${data}.FDRsig_eGenes.snps_in_${categ}.withPIP.txt.gz
