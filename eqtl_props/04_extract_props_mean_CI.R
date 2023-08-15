
## Script to compute CI for basic properties

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

eqtlfile="eqtl_props/eqtl.basic_props"
bootfile="eqtl_props/eqtl_boot.basic_props" # same as "eqtlfile", but computed and concatenated 1000 times for 1000 bootstrapped samples
eqtlfile_match="eqtl_props/eqtl_matched.basic_props" # same as "eqtlfile", but computed and concatenated 1000 times for 1000 instances of matched SNPs

d_eqtl=fread(eqtlfile,sep=",")
d_boot=fread(bootfile,sep=",")
d_eqtl_match=fread(eqtlfile_match,sep=",")

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

d_all=calc_CI(d_eqtl,d_boot)
d_all_match=calc_CI_mean(d_eqtl_match)

d_all$annot=rownames(d_all)
d_all_match$annot=rownames(d_all_match)

d_all$cat="eQTL"
d_all_match$cat="eQTL_matched"

d_both=rbind(d_all,d_all_match)
rownames(d_both)=1:nrow(d_both)

d_regression=d_both[grepl("coef",d_both$annot, fixed = F),]
d_props=d_both[!grepl("coef",d_both$annot, fixed = F),]

#####

write.table(d_regression,file=outfile1,quote=F,sep=",",row.names=F)
write.table(d_props,file=outfile2,quote=F,sep=",",row.names=F)











