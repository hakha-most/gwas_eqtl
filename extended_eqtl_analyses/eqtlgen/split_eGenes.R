args <- commandArgs(TRUE)
outdir <- args[1]
infile <- args[2]

library(data.table)
library(tidyverse)

df=fread(infile,header=F)

batch_size=150 #eGenes per batch

k=0
for (i in seq(1, nrow(df), batch_size)) {
  k=k+1
  print(k)
  seq_size <- batch_size
  if ((i + seq_size) > nrow(df)) seq_size <- nrow(df) - i + 1
  
  d_temp=df[i:(i+seq_size-1)]
  
  outfile=paste0(outdir,"/batch_",as.character(k),".txt")
  write.table(d_temp,file=outfile,quote = FALSE,sep="\t",row.names = F, col.names = F)
  
}
