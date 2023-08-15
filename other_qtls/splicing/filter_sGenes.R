args <- commandArgs(TRUE)
snpfile <- args[1]
pcfile <- args[2]
infile <- args[3]
outfile1 <- args[4]
outfile2 <- args[5]

library(tidyverse)
library(data.table)

set.seed(0)

d_snp=fread(snpfile)
d_pc_genes=fread(pcfile)

d_qtl=fread(infile,header=F)
colnames(d_qtl)=c("chr","trait","SNP","Pval")
d_qtl <- (d_qtl) %>% separate(trait,c("trait_x","nn"),sep="\\.",remove=F) %>% select(-nn)
d_qtl <- (d_qtl) %>% separate(trait_x,c("aa","bb","cc","dd","GeneSymbol"),sep=":",remove=F)
d_qtl=d_qtl %>% select(c("chr","trait","GeneSymbol","SNP","Pval"))

d_qtl_filter=d_qtl[d_qtl$GeneSymbol %in% d_pc_genes$GeneSymbol,]
d_qtl_filter=d_qtl_filter[d_qtl_filter$SNP %in% d_snp$SNP,]
d_qtl_filter=d_qtl_filter %>% arrange(Pval)

d_sgenes=d_qtl_filter[!duplicated(d_qtl_filter$trait),]

write.table(d_qtl_filter,file=outfile1,quote = FALSE,sep="\t",row.names = F)
write.table(d_sgenes$trait,file=outfile2,quote = FALSE,sep="\t",row.names = F,col.names=F)
