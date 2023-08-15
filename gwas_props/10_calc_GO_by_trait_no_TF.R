
## Script to compute GO enrichments (excluding TF genes) of GWAS SNPs for each trait separately

# This script, as provided, is run on GWAS SNPs. 
# The same script was run using 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="gwas_props/filter_indep_gwas.assoc" 
infofile="snp_annotations/filter_snps.txt"
GO_file="genes_multiGO/genes_noTF_GO_BP.txt" #see supp data
outfile="gwas_props/gwas.GO_props_noTF.by_trait"

d_info=fread(infofile)

d_go=fread(GO_file)

d_all=left_join(d_info %>% select(SNP,gene),d_go,by="gene")

d_gwas=fread(gwasfile)
d_gwas_annots=left_join(d_gwas,d_all,by="SNP")


#####

d_test=d_gwas_annots[d_gwas_annots$pheno=="50_irnt",] %>% select(-pheno,-SNP,-Pval,-gene)
d_props_temp=apply(d_test,2,function(x) {mean(x)})
d_props_temp=data.frame(t(data.frame(d_props_temp)))
d_props_temp$pheno="50_irnt"
d_props=d_props_temp[d_props_temp$pheno=="x",]

#

phenos=sort(unique(d_gwas_annots$pheno))

for (i in 1:length(phenos)){

	pheno_temp=phenos[i]
	print(i)

	d_test=d_gwas_annots[d_gwas_annots$pheno==pheno_temp,] %>% select(-pheno,-SNP,-Pval,-gene)
	d_props_temp=apply(d_test,2,function(x) {mean(x)})
	d_props_temp=data.frame(t(data.frame(d_props_temp)))
	d_props_temp$pheno=pheno_temp
	rownames(d_props_temp)=as.character(i)
	d_props=rbind(d_props,d_props_temp)
}


write.table(d_props,file=outfile,quote=F,sep=",",row.names=F)











