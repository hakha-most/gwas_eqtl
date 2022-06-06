
## Script to compute GO enrichments of eQTL SNPs for each tissue separately

# This script, as provided, is run on eQTL SNPs. 
# The same script was run using 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

eqtlfile="filter_eqtls.pruned_tissues.assoc" 
infofile="snp_annots/filter_snps.txt"
GO_file="gene_annots/genes_GO_BP.txt"
outfile="eqtl.GO_props.by_tissue"

d_info=fread(infofile)
d_go=fread(GO_file)

d_all=left_join(d_info %>% select(SNP,gene),d_go,by="gene")

d_eqtl=fread(eqtlfile)
d_eqtl=d_eqtl[d_eqtl$Pval<(5e-8),] #restrict to top eQTLs using a p-value cutoff matching GWAS

d_eqtl_annots=left_join(d_eqtl,d_all,by="SNP")


#####

d_test=d_gwas_annots[d_eqtl_annots$Tissue=="Testis",] %>% select(go_annots)
d_props_temp=apply(d_test,2,function(x) {mean(x)})
d_props_temp=data.frame(t(data.frame(d_props_temp)))
d_props_temp$tissue="Testis"
d_props=d_props_temp[d_props_temp$tissue=="x",]

#

tissues=sort(unique(d_eqtl_annots$Tissue))

for (i in 1:length(tissues)){

	tissue_temp=tissues[i]
	print(i)

	d_test=d_eqtl_annots[d_eqtl_annots$Tissue==tissue_temp,] %>% select(go_annots)
	d_props_temp=apply(d_test,2,function(x) {mean(x)})
	d_props_temp=data.frame(t(data.frame(d_props_temp)))
	d_props_temp$tissue=tissue_temp
	rownames(d_props_temp)=as.character(i)
	d_props=rbind(d_props,d_props_temp)
}


write.table(d_props,file=outfile,quote=F,sep=",",row.names=F)











