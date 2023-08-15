set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

##########

d_pc_genes=fread("gene_annotations/genes.protein_coding.v39.gtf")

#######
#function to identify closest gene to the QTL
get_closest_gene <- function(chr_temp,pos_temp,d_pc_genes){
  d_temp=d_pc_genes[d_pc_genes$chr==paste0("chr",as.character(chr_temp)),] 
  d_temp$d=abs(d_temp$tss - pos_temp)
  d_temp=d_temp %>% arrange(d)
  return(d_temp$hgnc_id[1])
}

#function to compute density of TSSs per site
get_TSSD <- function(chr_temp,pos_temp,d_pc_genes){
  tss1=max(0,pos_temp-500000)
  tss2=pos_temp+500000
  d_temp=d_pc_genes[d_pc_genes$chr==paste0("chr",as.character(chr_temp)) & (d_pc_genes$tss>tss1) & (d_pc_genes$tss<tss2),]
  return(nrow(d_temp))
}

#######
load("other_qtls/methylation/cpgs2include.RData") #data from Hawe et al.

#loop over QTLs
cpgs$gene=""
cpgs$TSSD=0
for (i in 1:nrow(cpgs)){
  print(i)
  cpgs$gene[i]=get_closest_gene(cpgs$chr[i],cpgs$pos[i],d_pc_genes)
  cpgs$TSSD[i]=get_TSSD(cpgs$chr[i],cpgs$pos[i],d_pc_genes)
  
}

colnames(cpgs)[1]="cpg"

outfile="other_qtls/methylation/cpgs.txt"
write.table(cpgs,file=outfile,quote=F,sep="\t",row.names=F)


