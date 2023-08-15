
## Script to compute CI for multifunctionality properties

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="gwas_props/gwas.multiGO_props"
bootfile="gwas_props/gwas_boot.multiGO_props" # same as "gwasfile", but computed and concatenated 1000 times for 1000 bootstrapped samples
gwasfile_match="gwas_props/gwas_matched.multiGO_props" # same as "gwasfile", but computed and concatenated 1000 times for 1000 instances of matched SNPs

d_gwas=fread(gwasfile,sep=",")
d_boot=fread(bootfile,sep=",")
d_gwas_match=fread(gwasfile_match,sep=",")

#########

calc_CI <- function(d_data,d_boot){
  

  d1= apply(d_boot,2,function(x) {unname(quantile(x,probs = seq(0,1,0.001))[26])})
  d2= apply(d_boot,2,function(x) {unname(quantile(x,probs = seq(0,1,0.001))[976])})

  d1=data.frame(d1);colnames(d1)="lower"
  d2=data.frame(d2);colnames(d2)="upper"
  d_data=data.frame(t(d_data)); colnames(d_data)="mean"
  
  d_all=cbind(d_data,d1,d2)

  return(d_all)

  
}


calc_CI_mean <- function(d_data){
  

  d1= apply(d_data,2,function(x) {unname(quantile(x,probs = seq(0,1,0.001))[26])})
  d2= apply(d_data,2,function(x) {unname(quantile(x,probs = seq(0,1,0.001))[976])})
  d_mean=apply(d_data,2,function(x) {mean(x)})

  d1=data.frame(d1);colnames(d1)="lower"
  d2=data.frame(d2);colnames(d2)="upper"
  d_data=data.frame(d_mean); colnames(d_data)="mean"
  
  d_all=cbind(d_data,d1,d2)
  
  return(d_all)

  
}

#########

d_all=calc_CI(d_gwas,d_boot)
d_all_match=calc_CI_mean(d_gwas_match)

d1=d_all[rownames(d_all) %in% rownames(d_all_match),]
d1$annot=rownames(d1)
d1$cat="GWAS"
d_all_match$annot=rownames(d_all_match)
d_all_match$cat="GWAS_matched"

d_props=rbind(d1,d_all_match)

#####

write.table(d_props,file=outfile,quote=F,sep=",",row.names=F)











