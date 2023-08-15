
set.seed(as.numeric(5))

library(data.table)
library(tidyverse)
library(dplyr)
library(glmnet)

##########

AUC_thresh=0.75

d_go=fread("12_GO_annotations/genes_BP_min400.txt")
d_go_names=fread("12_GO_annotations/GO_names_BP_min400.txt")

#results from GWAS and eQTL analyses, same as Table S6 
d_gwas=fread("12_GO_annotations/GWAS.two_tailed_enrichment.GO_props.by_trait")
d_eqtl=fread("12_GO_annotations/eQTL.two_tailed_enrichment.GO_props.by_trait")

d_all=rbind(d_gwas,d_eqtl)
d_all=d_all %>% arrange(-abs(Zscore))
d_all=d_all[d_all$GO.ID %in% d_go_names[d_go_names$Annotated<3000,]$GO.ID,] #filter GO annotations with <3K genes

d_all=left_join((d_all %>% dplyr::select(-Pval,-data,-enrich)),(d_go_names),by="GO.ID")
d_all_top_per_trait=d_all[!duplicated(d_all$pheno),] %>% dplyr::select(-pheno) #get top GO term per trait/tissue

d_top_GOs=d_all_top_per_trait %>% dplyr::select(GO.ID,Term,Annotated) #all unique top GO terms 
d_top_GOs=d_top_GOs[!duplicated(d_top_GOs$GO.ID),] %>% arrange(Annotated)

##########
#select GO terms, excluding the top GO terms above in d_top_GOs
dx=d_go %>% dplyr::select(-gene)
go_orders=data.frame(sort(colSums(dx)))
colnames(go_orders)=c("count")
go_orders$GO.ID=row.names(go_orders)
go_orders=go_orders[!(go_orders$GO.ID %in% d_top_GOs$GO.ID),]

##########
#prune top GO terms in d_top_GOs
selected_GOs=d_top_GOs$GO.ID[1:2] #pick the first 2 terms

#add to the list based on AUC
for (i in 3:nrow(d_top_GOs)){
  print(i)
  
  GO_test=d_top_GOs$GO.ID[i]
  dy=d_go %>% dplyr::select(GO_test); colnames(dy)="test"
  dx=d_go %>% dplyr::select(selected_GOs)
  
  fit <- cv.glmnet(as.matrix(dx), dy$test, family = "binomial",type.measure = "auc", keep = T)
  auc_test=max(fit$cvm)
  
  if (auc_test<AUC_thresh){selected_GOs=c(selected_GOs,GO_test)} #test whether the focal GO term can be predicted using previously selected terms
}


######
#iterate over remaining GO terms to add to the list using the same logic
for (i in 1:nrow(go_orders)){
  print(i)
  
  GO_test=go_orders$GO.ID[i]
  dy=d_go %>% dplyr::select(GO_test); colnames(dy)="test"
  dx=d_go %>% dplyr::select(selected_GOs)
  
  fit <- cv.glmnet(as.matrix(dx), dy$test, family = "binomial",type.measure = "auc", keep = T)
  auc_test=max(fit$cvm)
  
  if (auc_test<AUC_thresh){selected_GOs=c(selected_GOs,GO_test)}
  
}


######

write.table(selected_GOs,file="12_GO_annotations/indep_BP_terms.txt",quote=F,sep="\t",row.names=F,col.names = F)

