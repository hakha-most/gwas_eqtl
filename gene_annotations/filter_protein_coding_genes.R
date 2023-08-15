args <- commandArgs(TRUE)
infile <- args[1]
outfile1 <- args[2]
outfile2 <- args[3]

library(data.table)
library(tidyverse)
library(dplyr)

d_basic=fread(infile,header = F)
colnames(d_basic)=c("chr","start","end","strand","GeneSymbol","conseq","gene","hgnc_id")

d_basic=d_basic[d_basic$hgnc_id!="",]
d_basic=d_basic[d_basic$chr!="chrX" & d_basic$chr!="chrY" & d_basic$chr!="chrM",] #filter autosomes
d_basic=d_basic[d_basic$conseq=="protein_coding",] #filter protein coding genes

d_basic$tss=d_basic$start; d_basic[d_basic$strand=="-",]$tss=d_basic[d_basic$strand=="-",]$end #add TSS coordiantes
d_basic$tes=d_basic$end; d_basic[d_basic$strand=="-",]$tes=d_basic[d_basic$strand=="-",]$start #add TES coordiantes

write.table(d_basic,file=outfile1,quote = FALSE,sep="\t",row.names = F)

####
#1Mb window around TSSs used to (i) compute gene density and (ii) filter SNPs
d_interval=d_basic %>% select(chr,tss,hgnc_id) %>% separate(chr,c("x","chr"),sep = "r",remove = T)
d_interval$chr=as.numeric(as.character(d_interval$chr))
d_interval$start=pmax((d_interval$tss-1000000),0)
d_interval$end=d_interval$tss+1000000
d_interval=d_interval %>% select(chr,start,end,hgnc_id) %>% arrange(chr,start)

write.table(d_interval,file=outfile2,quote = FALSE,sep="\t",row.names = F, col.names = F)
