args <- commandArgs(TRUE)
enh_file <- args[1]
iter <- args[2]
outfile <- args[3]

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)
library(intervals)

#========

d_ref=fread(enh_file)
colnames(d_ref)=c("chr","start","end","gene","cell_type")

cell_types=sort(unique(d_ref$cell_type))
genes=sort(unique(d_ref$gene))

##

iter=as.numeric(iter)
i_strat=(iter-1)*200+1
i_end=(iter)*200
if (i_end > length(genes)){i_end=length(genes)}

##

genes=genes[i_strat:i_end]

counts1=rep(0,length(genes)) # active cell types
counts2=rep(0,length(genes)) # total elements
counts3=rep(0,length(genes)) # merged elements
counts_bp=rep(0,length(genes)) # merged length
counts_bp_per_type=rep(0,length(genes)) # mean length per type

#loop over genes
for (i in 1:length(genes)){
  print(i)
  
  gene_temp=genes[i]
  d_temp=d_ref[d_ref$gene==gene_temp,]
  d_temp=d_temp[d_temp$chr==d_temp$chr[1],]
  d_temp$t1=pmin(d_temp$start,d_temp$end) #interval start
  d_temp$t2=pmax(d_temp$start,d_temp$end) #interval end
  
  x=cbind(as.matrix(d_temp$t1),as.matrix(d_temp$t2))
  f=Intervals(x)
  z=as.matrix(interval_union(f)) # union of all intervals
  
  counts_bp[i]=sum(z[,2]-z[,1])
  counts1[i]=nrow(d_temp[!duplicated(d_temp$cell_type),])
  counts2[i]=nrow(d_temp)
  counts3[i]=dim(z)[1]
  
  ####
  # get union of intervals per biosample
  z_tot=matrix(, nrow = 0, ncol = 3)
  for (k in 1:length(cell_types)){
    cell_temp=cell_types[k]
    d_temp2=d_temp[d_temp$cell_type==cell_temp,]
    x_temp=cbind(as.matrix(d_temp2$t1),as.matrix(d_temp2$t2))
    f_temp=Intervals(x_temp)
    z_temp=as.matrix(interval_union(f_temp))
    z_temp=cbind(z_temp,(rep(k,nrow(z_temp))))
    z_tot=rbind(z_tot,z_temp)
  }
  
  ####
  
  dz=data.frame(z_tot)
  colnames(dz)=c("t1","t2","E")
  dz$d=dz$t2-dz$t1
  dk=aggregate(dz$d,by=list(dz$E), FUN=sum) #total bp per biosample
  counts_bp_per_type[i]=mean(dk$x) #mean bp across biosample
  
  
}

du=data.frame(genes,counts1,counts2,counts3,counts_bp,counts_bp_per_type)
write.table(du,file=outfile,quote = FALSE,sep=",",row.names = F,col.names=F)
