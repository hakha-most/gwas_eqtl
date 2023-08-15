args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]

library(tidyverse)
library(data.table)

set.seed(0)

d_qtl=fread(infile,header=F)
colnames(d_qtl)=c("Tissue","trait","SNP","Pval","rank")
d_qtl <- (d_qtl) %>% separate(trait,c("trait_x","nn"),sep="\\.",remove=F) %>% select(-nn)
d_qtl <- (d_qtl) %>% separate(trait_x,c("aa","bb","cc","dd","GeneSymbol"),sep=":",remove=F)
d_qtl=d_qtl %>% select(c("Tissue","GeneSymbol","SNP","Pval","rank"))

write.table(d_qtl,file=outfile,quote = FALSE,sep="\t",row.names = F)

