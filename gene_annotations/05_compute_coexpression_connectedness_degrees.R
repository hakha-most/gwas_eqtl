library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(igraph)


##saha_et_al directory is on Zenodo in dir "gene_annots"

### function to convert gene names to HGNC 
convert_name_to_hgnc <- function(gene_ids,d_pc,d_conv){
  
  dx=d_conv[d_conv$Symbol %in% gene_ids,]
  dx=dx[!duplicated(dx$Symbol),]
  
  dx=dx[dx$hgnc_id %in% d_pc$hgnc_id,]
  dx=dx[!duplicated(dx$hgnc_id),]
  
  colnames(dx)=c("gene","hgnc_id")
  return(dx)
  
}

### function to compute connectedness degrees 
calc_degree <- function(dg){
  
  colnames(dg)=c("ID1","ID2","weight")

  G=graph_from_data_frame(dg,directed = F)
  t=graph.strength(G)
  t_order=sort(t,decreasing = T)
  dt=data.frame(gene=names(t_order),degree=unname(t_order))
  return(dt)
  
}

###

indir="saha_et_al/twns/"
outdir=paste0(indir,"processed/")

d_pc_genes=fread("genes.protein_coding.v39.gtf")
d_conv=fread("ncbi_genes.hgnc_id.txt")

###

file_list=paste0(indir,"tissue_list.txt")
d_tissues=fread(file_list,header = F)
colnames(d_tissues)="tissue_id"


#####
#loop over tissues

d_genes_tissues=data.frame(gene=as.character(),rank=as.numeric(),tissue=as.character())

edge_type=1 #TE-TE interactions in Saha et al.

for (tissue in d_tissues$tissue_id){
  print(tissue)
  
  #read data from Saha et al.
  infile=paste0(indir,"raw_data/",tissue,".out.txt")
  df=fread(infile) 
  df=df[df$`Edge type`==edge_type,]
  df$Name1=as.character(df$Name1); df$Name2=as.character(df$Name2)
  
  df=df[(df$Name1 %in% d_conv$Symbol) & (df$Name2 %in% d_conv$Symbol),] #filter genes
  
  #convert gene symbol to HGNC ID
  gene_ids=sort(unique(c(df$Name1,df$Name2)))
  d_genes=convert_name_to_hgnc(gene_ids,d_pc_genes,d_conv)
  
  df=df[(df$Name1 %in% d_genes$gene) & (df$Name2 %in% d_genes$gene),]
  
  df$W1=1 #compute non-weighted degrees
  df$W2=abs(df$`Edge weight`) #compute weighted degrees to break ties

  ###
  
  d_temp1=calc_degree(df %>% select(Name1,Name2,W1)) #non-weighted
  d_temp2=calc_degree(df %>% select(Name1,Name2,W2)) #weighted
  
  d_temp=merge(d_temp1,d_temp2,by="gene")
  colnames(d_temp)=c("gene","degree1","degree2")
  
  d_temp=d_temp %>% arrange(-degree1,-degree2) #rank by non-weighted, when tied then rank by weighted
  d_temp$rank=1:nrow(d_temp)
  
  d_temp=d_temp %>% select(gene,rank)
  d_temp$tissue=tissue
  
  ###
  
  d_genes_tissues=rbind(d_genes_tissues,d_temp) #dataframe to store all tissues
  
  #save degrees per tissue
  d_temp=left_join(d_temp,d_genes,by="gene")
  d_temp=d_temp %>% select(hgnc_id,gene)

  outfile=paste0(outdir,tissue,".degree.txt")
  write.table(d_temp,file=outfile,quote=F,sep="\t",row.names=F)
  
}

#####
#combine all per-tissue ranks into one composite score

#compute last gene rank per tissue
d_last_rank=rep(0,nrow(d_tissues))
i=0
for (tissue in d_tissues$tissue_id){
  i=i+1
  d_temp=d_genes_tissues[d_genes_tissues$tissue==tissue,]
  d_last_rank[i]=tail(d_temp,1)$rank + 1
}

d_last_ranks=d_tissues; colnames(d_last_ranks)=c("tissue")
d_last_ranks$last_rank=d_last_rank

#compute prod rank
gene_list=unique(d_genes_tissues$gene)
prod_rank=rep(0,length(gene_list))
i=0
for (gene_temp in gene_list){
  i=i+1
  print(i)
  
  #gene rank across tissues
  d_temp1=d_genes_tissues[d_genes_tissues$gene==gene_temp,] %>% select(tissue,rank)
  
  #assign last tissue rank to gene if gene not ranked in that tissue 
  d_temp2=d_last_ranks[!(d_last_ranks$tissue %in% d_temp1$tissue),]; colnames(d_temp2)=c("tissue","rank") 
  
  d_temp3=rbind(d_temp1,d_temp2) #combine all tissues
  prod_rank[i]=prod(d_temp3$rank)^(1/16) #exponent added for numerical reasons
}

d_comp_rank=data.frame(gene=gene_list,Prod_Rank=prod_rank)

#convert gene symbol to HGNC ID
gene_ids=d_comp_rank$gene
d_genes=convert_name_to_hgnc(gene_ids,d_pc_genes,d_conv)

d_out=merge(d_comp_rank,d_genes,by="gene") %>% arrange(Prod_Rank) %>% select(hgnc_id,gene)
colnames(d_out)=c("hgnc_id","gene")

outfile=paste0(outdir,"genes_by_combined_rank.txt")
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)

