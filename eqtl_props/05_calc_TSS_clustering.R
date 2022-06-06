
## Script to compute TSS clustering of eQTL SNPs.

# This script, as provided, is run on eQTL SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

eqtlfile="filter_eqtls.pruned_tissues.assoc"
infofile="snp_annots/filter_snps.txt"
outfile="eqtl.TSS_props"

d_info=fread(infofile)

d_eqtl=fread(eqtlfile)
d_eqtl=d_eqtl[d_eqtl$Pval<(5e-8),] #restrict to top eQTLs using a p-value cutoff matching GWAS

d_eqtl_annots=left_join(d_eqtl,d_info,by="SNP")


#########


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
  colnames(d_stat)=c("prop_hits")
  d_stat$dTSS_bin=seq(-105,105,10)
  rownames(d_stat)=1:nrow(d_stat)
  return(d_stat)
  
}


#########

d_temp=calc_enrichment(d_eqtl_annots)
d_out=data.frame(t(d_temp$prop_hits))
colnames(d_out)=d_temp$dTSS_bin

write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)











