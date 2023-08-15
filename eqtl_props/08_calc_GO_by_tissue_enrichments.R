
## Script to compare GO enrichments for each tissue and matched SNPs

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

eqtlfile="eqtl_props/eqtl.GO_props.by_tissue"
eqtlfile_match="eqtl_props/eqtl_matched.GO_props.by_tissue" # same as "eqtlfile", but computed and concatenated 1000 times for 1000 instances of matched SNPs

d_eqtl=fread(eqtlfile,sep=",")
d_eqtl_match=fread(eqtlfile_match,sep=",")


#########

tissues=unique(d_eqtl$tissue)
d_out=data.frame(tissue=as.character(),data=as.numeric(),enrich=as.numeric(),Zscore=as.numeric(),Pval=as.numeric(),annot=as.character())

annots=colnames(d_eqtl %>% select(-tissue))

k=0
for (annot_temp in annots){

	k=k+1
	print(k)
	d_temp=d_eqtl %>% select(annot_temp,tissue)
	d_temp_match=d_eqtl_match %>% select(annot_temp,tissue)

	colnames(d_temp)=c("data","tissue")
	colnames(d_temp_match)=c("annot","tissue")
	d_means=d_temp_match[ ,list(mean=mean(annot)), by=tissue]
	d_sds=d_temp_match[ ,list(sd=sd(annot)), by=tissue]

	d1=merge(d_means,d_sds,by="tissue")
	d_sum=merge(d_temp,d1,by="tissue")

	d_sum$enrich=d_sum$data/d_sum$mean
	d_sum$Zscore=(d_sum$data-d_sum$mean)/d_sum$sd
	d_sum$abs_Z=abs(d_sum$data-d_sum$mean)/d_sum$sd

	d_sum$Pval=2*pnorm(d_sum$abs_Z,mean=0,sd=1,lower.tail = F)

	d_out_temp=d_sum %>% select(tissue,data,enrich,Zscore,Pval)

	d_out_temp$annot=annot_temp
	d_out=rbind(d_out,d_out_temp)
}

#########

write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)









