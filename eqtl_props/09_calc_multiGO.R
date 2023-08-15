
## Script to compute properties of eQTL SNPs wrt gene multifunctionality

# This script, as provided, is run on eQTL SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

eqtlfile="eqtl_props/filter_eqtls.pruned_tissues.assoc" 
infofile="snp_annotations/filter_snps.txt"
GO_file="gene_annotations/12_GO_annotations/genes_multiGO.txt"
outfile="eqtl_props/eqtl.multiGO_props"

d_info=fread(infofile)
d_gene=fread(genefile)

d_all=left_join(d_info,d_gene,by="gene")

d_eqtl=fread(eqtlfile)
d_eqtl=d_eqtl[d_eqtl$Pval<(5e-8),] #restrict to top eQTLs using a p-value cutoff matching GWAS

d_eqtl_annots=left_join(d_eqtl,d_all,by="SNP")

########

d_test=d_eqtl_annots

d_x1=d_test %>% group_by(GO_BP_count_decile) %>% count()
prop_BP_count_qs=t(d_x1$n/nrow(d_test))
colnames(prop_BP_count_qs)=paste0("prop_BP_count_q",as.character(1:10))

d_x1=d_test %>% group_by(GO_BP_count_400_auc_decile) %>% count()
prop_BP_count_400_auc_qs=t(d_x1$n/nrow(d_test))
colnames(prop_BP_count_400_auc_qs)=paste0("prop_BP_count_400_auc_q",as.character(1:10))

####

d_out=data.frame(prop_BP_count_qs,prop_BP_count_400_auc_qs)
write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)












