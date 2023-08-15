library(data.table)
library(tidyverse)
library(dplyr)

##########
#code to generate data for Supplementary Fig. 1

set.seed(0)

snpfile="snp_annotations/filter_snps.txt"
gwasfile="gwas_props/gwas_hits.all_phenos.clumped_sorted.txt"
outfile="extended_gwas_analyses/ukbiobank/gwas_hits.filtered_phenos.clumped_sorted.txt"

# data from Neale lab, see Supp Data
rg_file="data/ukb_neale_lab/rg.pairwise_min_0.5.txt"
phenos_file="data/ukb_neale_lab/filtered.phenos.txt"

d_gwas=fread(gwasfile,header=F)
colnames(d_gwas)=c("pheno","SNP","Pval")
d_gwas$SNP=as.character(d_gwas$SNP)

d_rg=fread(rg_file)
rg_phenos=unique(c(d_rg$p1,d_rg$p2))

good_phenos=fread(phenos_file,header = F)
colnames(good_phenos)=c("p","p_conv") # p: trait label; p_conv: trait label modified for categorical traits

d_gwas=d_gwas[d_gwas$pheno %in% good_phenos$p_conv,]
d_gwas=d_gwas[d_gwas$pheno %in% rg_phenos,]

d_snp=fread(snpfile)
d_gwas=d_gwas[(d_gwas$SNP %in% d_snp$SNP),]

write.table(d_gwas,file=outfile,quote=F,sep="\t",row.names=F)

