library(data.table)
library(tidyverse)
library(dplyr)

set.seed(0)

#get GO term count per gene
d_BP=fread("12_GO_annotations/genes_BP_min400.txt")
d_BP_x=d_BP %>% dplyr::select(-gene)

d_genes=d_BP %>% dplyr::select(gene) 
d_genes$GO_BP_count=rowSums(d_BP_x)

#get independent GO term count per gene
d_BP_indep=fread("12_GO_annotations/indep_BP_terms.txt",header = F); colnames(d_BP_indep)="GO.ID"
d_genes$GO_BP_count_400_auc=rowSums(d_BP_x %>% dplyr::select(d_BP_indep$GO.ID))
d_genes$GO_BP_count_400_auc_decile=ntile(d_genes$GO_BP_count_400_auc,10)

write.table(d_genes,file="12_GO_annotations/genes_multiGO.txt",quote=F,sep="\t",row.names=F)



