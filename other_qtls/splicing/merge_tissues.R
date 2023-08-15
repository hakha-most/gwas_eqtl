args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]

library(dplyr)
library(tidyverse)
library(data.table)

df=fread(infile,header=F)
colnames(df)=c("variant_id","gene_id","pval_nominal")

merge_tissues <- function(df){
  min_logP=as.numeric(sum(log(df[,"pval_nominal"])))
  return(min_logP)
}

concat_snps <- function(x){
  d_out=data.frame(names(x),unname(sapply(x, unname)))
  colnames(d_out)=c("variant_id","min_logP")
  return(d_out)
}

iterate_genes <- function(df){
  X <- by(df, df[,"variant_id"], merge_tissues)
  Y <- concat_snps(X)
  return(Y)
}

X <- by(df, df[,"gene_id"], iterate_genes)
Y <- mapply(cbind, X, "gene_id"=names(X), SIMPLIFY=F)
dg=rbindlist(Y,use.names=T)

write.table(dg,file=outfile,quote = FALSE,sep="\t",row.names = F)

