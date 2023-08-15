
## Script to compute promoter-enhancer enrichments of eQTL SNPs.

# This script, as provided, is run on eQTL SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

eqtlfile="eqtl_props/filter_eqtls.pruned_tissues.assoc"
infofile="snp_annotations/promoter_enhancer.counts"
outfile="eqtl_props/eqtl.promoter_enhancer.enrichments"

d_info=fread(infofile)

d_eqtl=fread(eqtlfile)
d_eqtl=d_eqtl[d_eqtl$Pval<(5e-8),] #restrict to top eQTLs using a p-value cutoff matching GWAS

d_eqtl_annots=left_join(d_eqtl,d_info,by="SNP")

p1=sapply((d_eqtl_annots %>% select(-Tissue,-SNP,-Pval,-GeneSymbol)), function(x) length(x[x>0])/length(x) ) 
p0=sapply((d_info %>% select(-SNP)), function(x) length(x[x>0])/length(x) ) 

enrichments=data.frame(t(p1/p0))

write.table(enrichments,file=outfile,quote=F,sep=",",row.names=F)
