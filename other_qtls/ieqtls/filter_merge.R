args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]

library(tidyverse)
library(data.table)

set.seed(0)

d_qtl=fread(infile,header=F)
colnames(d_qtl)=c("Tissue","Cell","Gene","SNP","Pval","rank")
d_qtl <- (d_qtl) %>% separate(Gene,c("GeneSymbol","nn"),sep="\\.",remove=T) %>% select(-nn)

write.table(d_qtl,file=outfile,quote = FALSE,sep="\t",row.names = F)
