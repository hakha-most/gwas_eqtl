library(tidyverse)
library(dplyr)
library(data.table)
library(genoppi)
library(igraph)

##########

### function to compute weighted PPI degrees 
calc_degree <- function(dg){
  
  colnames(dg)=c("ID1","ID2","weight")
  
  G=graph_from_data_frame(dg,directed = F)
  t=graph.strength(G)
  t_order=sort(t,decreasing = T)
  dt=data.frame(gene=names(t_order),degree=unname(t_order))
  return(dt)
  
}

### function to convert gene names to HGNC 
convert_name_to_hgnc <- function(gene_ids,d_pc,d_conv){
  
  dx=d_conv[d_conv$Symbol %in% gene_ids,]
  dx=dx[!duplicated(dx$Symbol),]
  
  dx=dx[dx$hgnc_id %in% d_pc$hgnc_id,]
  dx=dx[!duplicated(dx$hgnc_id),]
  
  colnames(dx)=c("gene","hgnc_id")
  return(dx)
  
}

##

d_pc_genes=fread("genes.protein_coding.v39.gtf") #protein-coding genes
d_conv=fread("ncbi_genes.hgnc_id.txt") #conversions table from gene symbol to HGNC ID

##

d_PPI=inweb_table %>% select(Gene1,Gene2,Score)
d_PPI=d_PPI[(d_PPI$Gene1 %in% d_conv$Symbol) & (d_PPI$Gene2 %in% d_conv$Symbol),]

dx=calc_degree(d_PPI %>% select(Gene1,Gene2,Score))

gene_ids=dx$gene
d_genes=convert_name_to_hgnc(gene_ids,d_pc_genes,d_conv)

dx=merge(dx,d_genes,by="gene")
dx=dx %>% arrange(-degree)
colnames(dx)[2]="PPI_degree"

dx$PPI_degree_decile=ntile(dx$PPI_degree,10)
dx$PPI_degree_quantile=ntile(dx$PPI_degree,5)

dx1=dx %>% select(gene,hgnc_id,PPI_degree,PPI_degree_decile,PPI_degree_quantile)
dx1$PPI_degree_cat=1 #indicator variable for a degree computed for the gene

#assign 0 to genes for which no degree was in the PPI data
dx2=d_pc_genes[!(d_pc_genes$hgnc_id %in% dx1$hgnc_id),] %>% select(gene,hgnc_id)
dx2$PPI_degree=0
dx2$PPI_degree_decile=0
dx2$PPI_degree_quantile=0
dx2$PPI_degree_cat=0

d_out=rbind(dx1,dx2)

outfile="PPI_degree.txt"
write.table(d_out,file=outfile,quote = FALSE,sep="\t",row.names = F)

