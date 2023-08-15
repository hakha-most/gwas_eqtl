
set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

#########
#get protein-coding genes expressed in blood

pcfile="gene_annotations/genes.protein_coding.v39.gtf"
d_pc_genes=fread(pcfile)

#expression levels in GTEx
gtexfile="colocalization/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
d_rna=fread(gtexfile,skip=2)
d_rna=d_rna %>% select("Name","Whole Blood")
colnames(d_rna)=c("Gene","tpm")
d_rna <- (d_rna) %>% separate(Gene,c("GeneSymbol","nn"),sep="\\.",remove=T) %>% select(-nn)

d_exp=d_rna[(d_rna$GeneSymbol %in% d_pc_genes$GeneSymbol) & (d_rna$tpm>0),]
d_exp=left_join(d_exp,(d_pc_genes %>% select(GeneSymbol,hgnc_id)), by="GeneSymbol")

#####
#find closest expressed gene to blood-trait GWAS SNPs

gwasfile="gwas_props/filter_indep_gwas.assoc" #GWAS results for all 44 traits
bloodfile="colocalization/blood_traits.txt"

d_gwas=fread(gwasfile)
d_gwas=d_gwas %>% distinct(`pheno`, `SNP`, .keep_all = TRUE) 

d_blood=fread(bloodfile,sep=",")
d_gwas_blood=d_gwas[(d_gwas$pheno %in% d_blood$pheno),]

d_snp=d_gwas_blood %>% select(SNP)
d_snp=d_snp[!duplicated(d_snp$SNP),]
d_snp <- (d_snp) %>% separate(SNP,c("chr","pos","A1","A2"),sep=":",remove=F) 
d_snp$chr=as.numeric(d_snp$chr)
d_snp$pos=as.numeric(d_snp$pos)

d_intervals=d_pc_genes %>% separate(chr,c("x","chr"),sep = "r",remove = T) %>% select(hgnc_id,chr,tss,tes)
d_intervals$hgnc_id=as.character(d_intervals$hgnc_id)
d_intervals$chr=as.numeric(as.character(d_intervals$chr))
d_intervals$tss=as.numeric(as.character(d_intervals$tss))
d_intervals$tes=as.numeric(as.character(d_intervals$tes))
d_intervals=d_intervals[d_intervals$hgnc_id %in% d_exp$hgnc_id,]

d_distance=d_snp %>% select(SNP,chr,pos)

d_tss_str="d_tss"
d_gene_tss="tss_gene"

d_distance[,d_tss_str]=0
d_distance[,d_gene_tss]="X"

header_d_tss=d_tss_str
header_tss_gene=d_gene_tss

#loop over GWAS SNPs
for (i in 1:nrow(d_distance)){
  print(i)
  
  chr_temp=d_distance$chr[i]
  pos=d_distance$pos[i]
  
  d_temp=d_intervals[d_intervals$chr==chr_temp,]
   
  dX=d_temp[order(abs(d_temp$tss-pos)),]
  D=(pos - dX[1,]$tss); if (dX[1,]$tss>dX[1,]$tes)(D=-D)
  d_distance[i,header_d_tss] = D; d_distance[i,header_tss_gene] = dX[1,]$hgnc_id
  
}

d_out=d_distance %>% select(SNP,d_tss,tss_gene)

outfile="colocalization/blood_trait.snps.txt"
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)

