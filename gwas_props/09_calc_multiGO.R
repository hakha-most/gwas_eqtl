
## Script to compute properties of GWAS SNPs wrt gene multifunctionality

# This script, as provided, is run on GWAS SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="filter_indep_gwas.assoc" 
infofile="snp_annots/filter_snps.txt"
genefile="gene_annots/genes_multiGO.txt"
outfile="gwas.multiGO_props"

d_info=fread(infofile)
d_gene=fread(genefile)

d_all=left_join(d_info,d_gene,by="gene")

d_gwas=fread(gwasfile)
d_gwas_annots=left_join(d_gwas,d_all,by="SNP")

########

d_test=d_gwas_annots

d_x1=d_test %>% group_by(GO_BP_count_decile) %>% count()
prop_BP_count_qs=t(d_x1$n/nrow(d_test))
colnames(prop_BP_count_qs)=paste0("prop_BP_count_q",as.character(1:10))

d_x1=d_test %>% group_by(GO_BP_count_400_auc_decile) %>% count()
prop_BP_count_400_auc_qs=t(d_x1$n/nrow(d_test))
colnames(prop_BP_count_400_auc_qs)=paste0("prop_BP_count_400_auc_q",as.character(1:10))

####

d_out=data.frame(prop_BP_count_qs,prop_BP_count_400_auc_qs)
write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)












