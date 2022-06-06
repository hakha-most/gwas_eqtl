
## Script to compute promoter-enhancer enrichments of GWAS SNPs.

# This script, as provided, is run on GWAS SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="filter_indep_gwas.assoc"
infofile="snp_annots/promoter_enhancer.counts"
outfile="gwas.promoter_enhancer.enrichments"

d_info=fread(infofile)

d_gwas=fread(gwasfile)
d_gwas_annots=left_join(d_gwas,d_info,by="SNP")

p1=sapply((d_gwas_annots %>% select(-pheno,-SNP,-Pval)), function(x) length(x[x>0])/length(x) ) 
p0=sapply((d_info %>% select(-SNP)), function(x) length(x[x>0])/length(x) ) 

enrichments=data.frame(t(p1/p0))

write.table(enrichments,file=outfile,quote=F,sep=",",row.names=F)
