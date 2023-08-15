
set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

#========

count_file="Liu_et_al_roadmap_gene_enhancer_links/processed/all_regions.props"
pc_file="genes.protein_coding.v39.gtf"
convert_file="ensembl_versions_lookup.txt"
out_file="Liu_et_al_roadmap_gene_enhancer_links/processed/Roadmap_enh_props.txt"

####

convert_ensembl_to_hgnc <- function(gene_ids,d_pc,d_conv){
  
  dx=d_conv[ !is.na(d_conv$hgnc_id),]
  dx=dx[dx$hgnc_id!="",]
  dx=dx[dx$hgnc_id!=" ",]
  
  dx=dx[dx$ensembl_gene_id %in% gene_ids,]
  dx=dx[!duplicated(dx$ensembl_gene_id),]
  
  dx=dx[dx$hgnc_id %in% d_pc$hgnc_id,]
  dx=dx[!duplicated(dx$hgnc_id),]
  dx=dx %>% select(ensembl_gene_id,hgnc_id)
  colnames(dx)=c("GeneSymbol","hgnc_id")
  return(dx)
  
}

####

d_pc=fread(pc_file)
d_conv=fread(convert_file)

d_count=fread(count_file)
colnames(d_count)=c("GeneSymbol","Roadmap_count1","Roadmap_count2","Roadmap_count3","Roadmap_bp","Roadmap_bp_per_type")

#convert ensembl ID to HGNC ID
gene_ids=d_count$GeneSymbol
d_genes=convert_ensembl_to_hgnc(gene_ids,d_pc,d_conv)

du=left_join(d_genes,d_count) %>% select(-GeneSymbol)

#assign zero to genes without annotations in Liu et al.
dz=left_join(d_pc %>% select(hgnc_id),du)
dz[is.na(dz$Roadmap_count1),]$Roadmap_count1=0
dz[is.na(dz$Roadmap_count2),]$Roadmap_count2=0
dz[is.na(dz$Roadmap_count3),]$Roadmap_count3=0
dz[is.na(dz$Roadmap_bp),]$Roadmap_bp=0
dz[is.na(dz$Roadmap_bp_per_type),]$Roadmap_bp_per_type=0

write.table(dz,file=out_file,quote = FALSE,sep="\t",row.names = F)
