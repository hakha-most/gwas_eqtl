
## Script to compare GO enrichments (excluding TF genes) for each trait and matched SNPs

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="gwas_props/gwas.GO_props_noTF.by_trait"
gwasfile_match="gwas_props/gwas_matched.GO_props_noTF.by_trait" # same as "gwasfile", but computed and concatenated 1000 times for 1000 instances of matched SNPs

d_gwas=fread(gwasfile,sep=",")
d_gwas_match=fread(gwasfile_match,sep=",")


#########

phenos=unique(d_gwas$pheno)
d_out=data.frame(pheno=as.character(),data=as.numeric(),enrich=as.numeric(),Zscore=as.numeric(),Pval=as.numeric(),annot=as.character())

annots=colnames(d_gwas %>% select(-pheno))

k=0
for (annot_temp in annots){

	k=k+1
	print(k)
	d_temp=d_gwas %>% select(annot_temp,pheno)
	d_temp_match=d_gwas_match %>% select(annot_temp,pheno)

	colnames(d_temp)=c("data","pheno")
	colnames(d_temp_match)=c("annot","pheno")
	d_means=d_temp_match[ ,list(mean=mean(annot)), by=pheno]
	d_sds=d_temp_match[ ,list(sd=sd(annot)), by=pheno]

	d1=merge(d_means,d_sds,by="pheno")
	d_sum=merge(d_temp,d1,by="pheno")

	d_sum$enrich=d_sum$data/d_sum$mean
	d_sum$Zscore=(d_sum$data-d_sum$mean)/d_sum$sd
	d_sum$abs_Z=abs(d_sum$data-d_sum$mean)/d_sum$sd

	d_sum$Pval=2*pnorm(d_sum$abs_Z,mean=0,sd=1,lower.tail = F)

	d_out_temp=d_sum %>% select(pheno,data,enrich,Zscore,Pval)

	d_out_temp$annot=annot_temp
	d_out=rbind(d_out,d_out_temp)
}

#########

write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)









