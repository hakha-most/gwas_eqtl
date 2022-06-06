args <- commandArgs(TRUE)
infile <- args[1]
maffile <- args[2]
outfile <- args[3]

library(dplyr)
library(tidyverse)
library(data.table)

d_maf=fread(maffile,header=F); colnames(d_maf)=c("ID")
d_snp=fread(infile,header=F); colnames(d_snp)=c("pheno","ID1","P")

d_snp <- d_snp %>% 
  separate(ID1, c("chr","pos","A1","A2"), sep = c(":"),remove = FALSE) %>% 
  unite("ID2", c("chr","pos","A2","A1"), sep=":", na.rm = TRUE, remove = FALSE) 

d_maf$ID=as.character(d_maf$ID)
d_snp$ID1=as.character(d_snp$ID1)
d_snp$ID2=as.character(d_snp$ID2)

d1=d_snp[d_snp$ID1 %in% d_maf$ID,] %>% select(pheno,ID1,P)
d2=d_snp[d_snp$ID2 %in% d_maf$ID,] %>% select(pheno,ID2,P)
colnames(d1)=c("pheno","SNP","P"); colnames(d2)=c("pheno","SNP","P")
d_out=rbind(d1,d2)

write.table(d_out,file=outfile,quote=F,sep=" ",row.names=F)


