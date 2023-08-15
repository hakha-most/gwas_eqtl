
library(data.table)
library(tidyverse)
library(dplyr)
library(intervals)

#========

abc_dir="ABC" #directory on Zenodo in dir "gene_annots"
d_pc_genes=fread("genes.protein_coding.v39.gtf")

###########

d_abc=fread(paste0(abc_dir,"/raw_data/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"))
d_abc=d_abc[d_abc$class!="promoter",] #filter enhancers

cell_types=sort(unique(d_abc$CellType))
genes=sort(unique(d_abc$TargetGene))

counts1=rep(0,length(genes)) # active cell types count
counts2=rep(0,length(genes)) # total elements count
counts3=rep(0,length(genes)) # merged elements count
counts_bp=rep(0,length(genes)) # merged elements length
counts_bp_per_type=rep(0,length(genes)) # mean length per type

#loop over genes
for (i in 1:length(genes)){
  print(i)
  
  gene_temp=genes[i]
  d_temp=d_abc[d_abc$TargetGene==gene_temp,]
  d_temp=d_temp[d_temp$chr==d_temp$chr[1],]
  d_temp$t1=pmin(d_temp$start,d_temp$end) #interval start
  d_temp$t2=pmax(d_temp$start,d_temp$end) #interval end
  
  x=cbind(as.matrix(d_temp$t1),as.matrix(d_temp$t2))
  f=Intervals(x)
  z=as.matrix(interval_union(f)) # union of all intervals
  
  counts_bp[i]=sum(z[,2]-z[,1])
  counts1[i]=nrow(d_temp[!duplicated(d_temp$CellType),])
  counts2[i]=nrow(d_temp)
  counts3[i]=dim(z)[1]
  
  #####
  # get union of intervals per biosample
  
  z_tot=matrix(, nrow = 0, ncol = 3)
  for (k in 1:length(cell_types)){
    cell_temp=cell_types[k]
    d_temp2=d_temp[d_temp$CellType==cell_temp,]
    x_temp=cbind(as.matrix(d_temp2$t1),as.matrix(d_temp2$t2))
    f_temp=Intervals(x_temp)
    z_temp=as.matrix(interval_union(f_temp))
    z_temp=cbind(z_temp,(rep(k,nrow(z_temp))))
    z_tot=rbind(z_tot,z_temp)
  }
  
  dz=data.frame(z_tot)
  colnames(dz)=c("t1","t2","E")
  dz$d=dz$t2-dz$t1
  dk=aggregate(dz$d,by=list(dz$E), FUN=sum) #total bp per biosample
  counts_bp_per_type[i]=mean(dk$x) #mean bp across biosample
}

du=data.frame(gene=genes,ABC_count1=counts1,ABC_count2=counts2,ABC_count3=counts3,ABC_bp=counts_bp,ABC_bp_per_type=counts_bp_per_type)

outfile=paste0(abc_dir,"/processed/ABC_enh_counts.txt")
write.table(du,file=outfile,quote = FALSE,sep="\t",row.names = F)

############
#convert gene ids

### function to convert gene names to HGNC 
convert_name_to_hgnc <- function(gene_ids,d_pc,d_conv){
  
  dx=d_conv[d_conv$Symbol %in% gene_ids,]
  dx=dx[!duplicated(dx$Symbol),]
  
  dx=dx[dx$hgnc_id %in% d_pc$hgnc_id,]
  dx=dx[!duplicated(dx$hgnc_id),]
  
  colnames(dx)=c("gene","hgnc_id")
  return(dx)
  
}

####

convert_file="ncbi_genes.hgnc_id.txt"
d_conv=fread(convert_file)

d_count=fread(paste0(abc_dir,"/processed/ABC_enh_counts.txt"))

#convert gene symbol to HGNC ID
gene_ids=d_count$gene
d_genes=convert_name_to_hgnc(gene_ids,d_pc_genes,d_conv)
du=left_join(d_genes,d_count) %>% select(-gene)

#assign zero to genes without ABC annotations 
dz=left_join(d_pc_genes %>% select(hgnc_id),du)
dz[is.na(dz$ABC_count1),]$ABC_count1=0
dz[is.na(dz$ABC_count2),]$ABC_count2=0
dz[is.na(dz$ABC_count3),]$ABC_count3=0
dz[is.na(dz$ABC_bp),]$ABC_bp=0
dz[is.na(dz$ABC_bp_per_type),]$ABC_bp_per_type=0

outfile=paste0(abc_dir,"/processed/ABC_enh_counts_convert_ids.txt")
write.table(dz,file=outfile,quote = FALSE,sep="\t",row.names = F)
