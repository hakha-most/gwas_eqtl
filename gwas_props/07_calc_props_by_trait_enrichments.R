
## Script to compare basic properties for each trait and matched SNPs

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="gwas.basic_props.by_trait"
gwasfile_match="gwas_matched.basic_props.by_trait" # same as "gwasfile", but computed and concatenated 1000 times for 1000 instances of matched SNPs

d_gwas=fread(gwasfile,sep=",")
d_gwas_match=fread(gwasfile_match,sep=",")


#########

annots=unique(d_gwas$annot)
d_out=data.frame(pheno=colnames(d_gwas %>% select(-annot)))

for (annot_temp in annots){

	d_temp=d_gwas[d_gwas$annot==annot_temp,] %>% select(-annot)
	d_temp_match=d_gwas_match[d_gwas_match$annot==annot_temp,] %>% select(-annot)

	means=as.vector(unlist(lapply(d_temp_match,mean)))
	sds=as.vector(unlist(lapply(d_temp_match,sd)))

	enrich=d_temp/means
	Zs=(d_temp-means)/sds
	abs_Zs=abs(d_temp-means)/sds

	Pvals=2*pnorm(as.vector(t(abs_Zs)[,1]),mean=0,sd=1,lower.tail = F)

	d_out_temp=data.frame(data=as.vector(t(d_temp)[,1]),enrichment=as.vector(t(enrich)[,1]),Zscore=as.vector(t(Zs)[,1]),Pval=Pvals)
	colnames(d_out_temp)=paste0(colnames(d_out_temp),"_",annot_temp)

	d_out_temp$pheno=colnames(d_temp)
	d_out=left_join(d_out,d_out_temp,by="pheno")
}

#########

write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)









