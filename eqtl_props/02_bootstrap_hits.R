
## Script for obtaining bootstrapped samples of eQTLs by resampling tissues and LD blocks

args <- commandArgs(TRUE)
iter <- args[1] #bootstrapping iteration; performed 1000 times
eqtlfile="filter_eqtls.pruned_tissues.assoc"
blockfile="snp_annots/snps.LD_blocks.txt"
outfile=paste0("boot_",as.character(iter),".txt")

library(tidyverse)
library(data.table)
set.seed(as.numeric(iter))

d_eqtl=fread(eqtlfile)

d_blocks=fread(blockfile,header = F)
colnames(d_blocks)=c("SNP","Block")

bad_snps=d_blocks[duplicated(d_blocks$SNP),]$SNP #remove a few SNPs assigned to multiple blocks
d_blocks=d_blocks[!(d_blocks$SNP %in% bad_snps),]

df=d_eqtl[(d_eqtl$SNP %in% d_blocks$SNP),]
df=left_join(df,d_blocks,by="SNP")

#####

eqtl_blocks=unique(df$Block)
eqtl_tissues=unique(df$Tissue)

blocks_temp=sample(eqtl_blocks,size = length(eqtl_blocks),replace = T)
tissues_temp=sample(eqtl_tissues,size = length(eqtl_tissues),replace = T)
  
#bootstrap tissues
d_boot1=df[df$Tissue=="XXXX",]
tx=aggregate(tissues_temp, list(tissues_temp), length); colnames(tx)=c("tissue","count")
for (j in 1:nrow(tx)){
  dx=do.call("rbind", replicate(tx$count[j], df[df$Tissue==tx$tissue[j],], simplify = FALSE))
  d_boot1=rbind(d_boot1,dx)
}

#bootstrap LD blocks
d_boot2=df[df$Tissue=="XXXX",]
tx=aggregate(blocks_temp, list(blocks_temp), length); colnames(tx)=c("block","count")
for (j in 1:nrow(tx)){
  dx=do.call("rbind", replicate(tx$count[j], d_boot1[d_boot1$Block==tx$block[j],], simplify = FALSE))
  d_boot2=rbind(d_boot2,dx)
}
  

#####

d_out=d_boot2 %>% select(-c(Block))
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)
