
## Script to filter eQTLs and prune brain tissues

#file concatenating the output of codes in "eqtl_process" across all tissue-eGene pairs
eqtlfile="eqtl_props/eqtls.all_tissues.clumped.txt"
snpfile="snp_annotations/filter_snps.txt"
genefile="gene_annotations/genes.protein_coding.v39.gtf"
outfile1="eqtl_props/filter_eqtls.all_tissues.assoc"
outfile2="eqtl_props/filter_eqtls.pruned_tissues.assoc"

library(tidyverse)
library(data.table)
set.seed(0)

d_pc_genes=fread(genefile)
d_snp=fread(snpfile)

d_eqtl=fread(eqtlfile,header=F)
colnames(d_eqtl)=c("Tissue","Gene","SNP","Pval","rank")
d_eqtl$SNP=as.character(d_eqtl$SNP)
d_eqtl <- (d_eqtl) %>% separate(Gene,c("GeneSymbol","nn"),sep="\\.",remove=F) %>% select(-nn)
d_eqtl <- d_eqtl %>% select(Tissue,GeneSymbol,SNP,Pval)

df=d_eqtl[(d_eqtl$SNP %in% d_snp$SNP),]
df=df[(df$GeneSymbol %in% d_pc_genes$GeneSymbol),]

#######

tissues_to_remove=unique(df$Tissue)[grepl("Brain", unique(df$Tissue), fixed = F)]
tissues_to_remove=tissues_to_remove[!(tissues_to_remove %in% c("Brain_Cortex","Brain_Cerebellum"))]
dg=df[!(df$Tissue %in% tissues_to_remove),]

######

write.table(df,file=outfile1,quote=F,sep="\t",row.names=F)
write.table(dg,file=outfile2,quote=F,sep="\t",row.names=F)
