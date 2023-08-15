
set.seed(0)

library(tidyverse)
library(data.table)

snpfile="snp_annotations/filter_snps.txt" #data on Zenodo
pc_file="gene_annotations/genes.protein_coding.v39.gtf"

d_snp=fread(snpfile)
d_snp$SNP=as.character(d_snp$SNP)

d_pc_genes=fread(pc_file)

####

eqtlfile="extended_eqtl_analyses/eqtlcatalouge/ge.clumped.assoc.txt"
outfile="supp_note_data/eqtl_catalogue/ge_eqtls.clumped.txt" #deposited on Zenodo

d_eqtl=fread(eqtlfile)
colnames(d_eqtl)[1]="GeneSymbol"

d_out=d_eqtl[(d_eqtl$GeneSymbol %in% d_pc_genes$GeneSymbol),]
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)

###

eqtlfile="extended_eqtl_analyses/eqtlcatalouge/exon.clumped.assoc.txt"
outfile="supp_note_data/eqtl_catalogue/exon_eqtls.clumped.txt" #deposited on Zenodo

d_eqtl=fread(eqtlfile)
d_eqtl <- (d_eqtl) %>% separate(molecular_trait_id,c("GeneSymbol","nn"),sep="\\.",remove=F) %>% select(-nn)

d_out=d_eqtl[(d_eqtl$GeneSymbol %in% d_pc_genes$GeneSymbol),]
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)

