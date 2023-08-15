library(data.table)
library(tidyverse)
library(dplyr)
library(biomaRt)

d_archives=listEnsemblArchives()

d_bm=data.frame(ensembl_gene_id=as.character(),external_gene_name=as.character(),hgnc_id=as.character(),version=as.character())

#extract HGNC ID for all release versions into one dataframe
for (i in 1:nrow(d_archives)){
    
  url_temp=d_archives$url[i]
  version_temp=d_archives$name[i]
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = url_temp)
  
  d_bm_temp <- getBM(c("ensembl_gene_id","external_gene_name","hgnc_id"), mart=mart)
  d_bm_temp$version=as.character(version_temp)
  
  d_bm=rbind(d_bm,d_bm_temp)
  
}

# GRCh37 version
d_bm11=d_bm[d_bm$version=="Ensembl GRCh37" & is.na(d_bm$hgnc_id),] 
d_bm12=d_bm[d_bm$version=="Ensembl GRCh37" & !is.na(d_bm$hgnc_id),]
d_bm12$hgnc_id=paste0("HGNC:",d_bm12$hgnc_id)

# other versions
d_bm2=d_bm[d_bm$version!="Ensembl GRCh37",] 
d_bm21=d_bm2[d_bm2$version=="Ensembl 105",] # most up-to-date version at the time of study (prioritized)
d_bm22=d_bm2[d_bm2$version!="Ensembl 105",]

d_out=rbind(d_bm21,d_bm12,d_bm11,d_bm22) # in order of priority
outfile="ensembl_versions_lookup.txt"
write.table(d_out,file=outfile,quote = FALSE,sep="\t",row.names = F)