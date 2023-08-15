args <- commandArgs(TRUE)
maffile <- args[1]
infile1 <- args[2]
infile2 <- args[3]
infile3 <- args[4]
outfile <- args[5]

library(tidyverse)
library(data.table)

d_var=fread(infile1)
d_pos1=fread(infile2,header=F)
d_pos2=fread(infile3,header=F)
d_maf=fread(maffile,header = F)

d_var$snp <- paste0("snp",as.character(1:nrow(d_var)))

d_pos1=d_pos1 %>% select(V1,V3,V4); colnames(d_pos1)=c("chr","pos_38","snp")
d_pos2=d_pos2 %>% select(V1,V3,V4); colnames(d_pos2)=c("chr","pos_19","snp")
d_pos=merge(d_pos1,d_pos2,by=c("chr","snp"))

df=merge(d_var,d_pos,by="snp")
df$chr_num=gsub("[a-z]","", df$chr)
df$pos_19_num=as.numeric(as.character(df$pos_19))
df$chr_num=as.numeric(as.character(df$chr_num))
df$chromosome=as.numeric(as.character(df$chromosome))
df <- df[(df$chr_num>=1 & df$chr_num<=22),]
df=df[df$chr_num==df$chromosome,]

cols <- c('chr_num','pos_19_num','ref','alt')
df <- unite(df,"allele1", cols, sep=":", na.rm = TRUE, remove = FALSE)

cols <- c('chr_num','pos_19_num','alt','ref')
df <- unite(df,"allele2", cols, sep=":", na.rm = TRUE, remove = FALSE)


dg1=df[df$allele1 %in% d_maf$V1,]
dg2=df[df$allele2 %in% d_maf$V1,]

dg1$SNP=dg1$allele1
dg2$SNP=dg2$allele2

dg=rbind(dg1,dg2)
dg=dg %>% select(-c(snp,variant,chromosome,position,allele1,allele2,pos_38,pos_19_num,chr_num))

write.table(dg,file=outfile,quote = FALSE,sep="\t",row.names = F)
