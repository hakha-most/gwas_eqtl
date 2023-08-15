library(data.table)
library(tidyverse)
library(dplyr)

###

d_pc_genes=fread("genes.protein_coding.v39.gtf")
d_conv=fread("ensembl_versions_lookup.txt")

###

d_annot=fread("ensemble_biomart_cds_length_grch38.txt")
colnames(d_annot)=c("GeneSymbol","CDS_length","Transcript_count")
d_annot=d_annot %>% arrange(GeneSymbol,-CDS_length) %>% select(-Transcript_count)

d_annot=d_annot[!duplicated(d_annot$GeneSymbol),] #get the transcript with longest CDS
d_annot_pc=d_annot[d_annot$GeneSymbol %in% d_pc_genes$GeneSymbol,]

outfile="ensemble_biomart_cds_length_grch38.processed.txt"
write.table(d_annot_pc,file=outfile,quote=F,sep="\t",row.names=F)

