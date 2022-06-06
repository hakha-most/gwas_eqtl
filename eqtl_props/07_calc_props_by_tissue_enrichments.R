
## Script to compare basic properties for each tissue and matched SNPs

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

eqtlfile="eqtl.basic_props.by_tissue"
eqtlfile_match="eqtl_matched.basic_props.by_tissue" # same as "eqtlfile", but computed and concatenated 1000 times for 1000 instances of matched SNPs

d_eqtl=fread(eqtlfile,sep=",")
d_eqtl_match=fread(eqtlfile_match,sep=",")


#########

annots=unique(d_eqtl$annot)
d_out=data.frame(tissue=colnames(d_eqtl %>% select(-annot)))

for (annot_temp in annots){

	d_temp=d_eqtl[d_eqtl$annot==annot_temp,] %>% select(-annot)
	d_temp_match=d_eqtl_match[d_eqtl_match$annot==annot_temp,] %>% select(-annot)

	means=as.vector(unlist(lapply(d_temp_match,mean)))
	sds=as.vector(unlist(lapply(d_temp_match,sd)))

	enrich=d_temp/means
	Zs=(d_temp-means)/sds
	abs_Zs=abs(d_temp-means)/sds

	Pvals=2*pnorm(as.vector(t(abs_Zs)[,1]),mean=0,sd=1,lower.tail = F)

	d_out_temp=data.frame(data=as.vector(t(d_temp)[,1]),enrichment=as.vector(t(enrich)[,1]),Zscore=as.vector(t(Zs)[,1]),Pval=Pvals)
	colnames(d_out_temp)=paste0(colnames(d_out_temp),"_",annot_temp)

	d_out_temp$tissue=colnames(d_temp)
	d_out=left_join(d_out,d_out_temp,by="tissue")
}

#########

write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)









