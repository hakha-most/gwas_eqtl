library(data.table)
library(tidyverse)
library(dplyr)

##########
#code to generate data for Supplementary Fig. 2

set.seed(0)

atlasdir="extended_gwas_analyses/gwas_atlas/"
d_pheno=fread(paste0(atlasdir,"atlas.phenos.txt"))
d_gwas=fread(paste0(atlasdir,"zenodo_risk_loci.txt")) %>% select(id,rsID,p)
##

infofile="snp_annotations/filter_snps.txt" #on Zenodo
d_all=fread(infofile)

d_gwas_annots=left_join(d_gwas[d_gwas$rsID %in% d_all$rsID,],d_all,by="rsID")

######################
#function to compute fraction of SNP by distance bins

calc_enrichment <- function(d_data){
  
  d_distance <- d_data %>% select(SNP,d_tss)
  
  bins1=seq(10,100,10)*1000
  bins=c(-rev(bins1),0,bins1)
  
  bins_str1=paste0("+",as.character(bins1/1000),"kb"); #bins_str1[length(bins_str1)]=">+100kb"
  bins_str2=paste0("-",as.character(bins1/1000),"kb"); #bins_str2[length(bins_str2)]="<-100kb"
  bins_str=c(rev(bins_str2),bins_str1)
  
  d_distance[,bins_str]=0
  
  for (i in 1:(length(bins)-1)){
    
    d1=bins[i]
    d2=bins[i+1]
    
    header_temp=bins_str[i]
    index <- ((d_distance[,d_tss] >= d1) & (d_distance[,d_tss] < d2))
    d_distance[index,header_temp] <- 1
    
  }
  
  d_tss_str=paste0("<-100kb")
  d_distance[,d_tss_str]=0
  index <- ((d_distance[,d_tss] < (-100000)))
  d_distance[index,d_tss_str] <- 1
  
  d_tss_str=paste0(">100kb")
  d_distance[,d_tss_str]=0
  index <- ((d_distance[,d_tss] > (100000)))
  d_distance[index,d_tss_str] <- 1
  
  cols=c("<-100kb",bins_str,">100kb")
  d_stat=(data.frame(colSums(d_distance %>% select(cols))/nrow(d_distance)))
  colnames(d_stat)=c("prop")
  rownames(d_stat)=1:nrow(d_stat)
  
  d_N=(data.frame(colSums(d_distance %>% select(cols))))
  colnames(d_N)=c("N_hits")
  rownames(d_N)=1:nrow(d_N)
  
  d_stat$prop_se=sqrt(d_stat$prop*(1-d_stat$prop)/d_N$N_hits)
  d_stat$dTSS_bin=seq(-105,105,10)
  
  
  return(d_stat)
  
}
######################

#fraction of SNP by distance bins by trait category

d_test=d_gwas_annots[d_gwas_annots$id %in% (d_pheno[d_pheno$no_rep_phenos==1,]$id),]
d_out1=calc_enrichment(d_test)
d_out1$trait_cat="no_rep_phenos"

d_test=d_gwas_annots[d_gwas_annots$id %in% (d_pheno[d_pheno$no_rep_disease==1,]$id),]
d_out2=calc_enrichment(d_test)
d_out2$trait_cat="no_rep_disease"

d_test=d_gwas_annots[d_gwas_annots$id %in% (d_pheno[d_pheno$no_rep_disease_disorder==1,]$id),]
d_out3=calc_enrichment(d_test)
d_out3$trait_cat="no_rep_disease_disorder"

d_test=d_gwas_annots[d_gwas_annots$id %in% (d_pheno[d_pheno$indep_phenos==1,]$id),]
d_out4=calc_enrichment(d_test)
d_out4$trait_cat="indep_phenos"

d_test=d_gwas_annots[d_gwas_annots$id %in% (d_pheno[d_pheno$indep_disease==1,]$id),]
d_out5=calc_enrichment(d_test)
d_out5$trait_cat="indep_disease"

d_test=d_gwas_annots[d_gwas_annots$id %in% (d_pheno[d_pheno$indep_disease_disorder==1,]$id),]
d_out6=calc_enrichment(d_test)
d_out6$trait_cat="indep_disease_disorder"

d_out=rbind(d_out1,d_out2,d_out3,d_out4,d_out5,d_out6)

#######

outfile=paste0(atlasdir,"GWAS.atlas.TSS_props.txt")
write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)

