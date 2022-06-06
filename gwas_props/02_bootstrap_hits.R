
## Script for obtaining bootstrapped samples of GWAS hits by resampling traits and LD blocks

args <- commandArgs(TRUE)
iter <- args[1] #bootstrapping iteration; performed 1000 times
gwasfile="filter_indep_gwas.assoc"
blockfile="snp_annots/snps.LD_blocks.txt"
outfile=paste0("boot_",as.character(iter),".txt")

library(tidyverse)
library(data.table)
set.seed(as.numeric(iter))

d_gwas=fread(gwasfile)

d_blocks=fread(blockfile,header = F)
colnames(d_blocks)=c("SNP","Block")

bad_snps=d_blocks[duplicated(d_blocks$SNP),]$SNP #remove a few SNPs assigned to multiple blocks
d_blocks=d_blocks[!(d_blocks$SNP %in% bad_snps),]

df=d_gwas[(d_gwas$SNP %in% d_blocks$SNP),]
df=left_join(df,d_blocks,by="SNP")

#####

gwas_blocks=unique(df$Block)
gwas_traits=unique(df$pheno)

blocks_temp=sample(gwas_blocks,size = length(gwas_blocks),replace = T)
traits_temp=sample(gwas_traits,size = length(gwas_traits),replace = T)

#bootstrap traits
d_boot1=df[df$pheno=="XXXX",]
tx=aggregate(traits_temp, list(traits_temp), length); colnames(tx)=c("trait","count")
for (j in 1:nrow(tx)){
  dx=do.call("rbind", replicate(tx$count[j], df[df$pheno==tx$trait[j],], simplify = FALSE))
  d_boot1=rbind(d_boot1,dx)
}

#bootstrap LD blocks
d_boot2=df[df$pheno=="XXXX",]
tx=aggregate(blocks_temp, list(blocks_temp), length); colnames(tx)=c("block","count")
for (j in 1:nrow(tx)){
  dx=do.call("rbind", replicate(tx$count[j], d_boot1[d_boot1$Block==tx$block[j],], simplify = FALSE))
  d_boot2=rbind(d_boot2,dx)
}
  

#####

d_out=d_boot2 %>% select(-c(Block))
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)
