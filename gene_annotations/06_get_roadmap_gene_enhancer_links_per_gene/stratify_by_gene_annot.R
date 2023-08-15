bedfile <- "Liu_et_al_roadmap_gene_enhancer_links/processed/all_regions.merged_biosamples.txt"
annotfile <- "all_annots_pc_genes.txt"
pcfile <- "genes.protein_coding.v39.gtf"
convfile_gene <- "ncbi_genes.hgnc_id.txt"
convfile_ensembl <- "ensembl_versions_lookup.txt"
outdir <- "Liu_et_al_roadmap_gene_enhancer_links/processed/by_gene_constraint/"

library(tidyverse)
library(dplyr)
library(data.table)

### to convert Ensembl gene IDs to HGNC 
convert_ensembl_to_hgnc <- function(gene_ids,d_pc,d_conv){
  
  dx=d_conv[ !is.na(d_conv$hgnc_id),]
  dx=dx[dx$hgnc_id!="",]
  dx=dx[dx$hgnc_id!=" ",]
  
  dx=dx[dx$ensembl_gene_id %in% gene_ids,]
  dx=dx[!duplicated(dx$ensembl_gene_id),]
  
  dx=dx[dx$hgnc_id %in% d_pc$hgnc_id,]
  dx=dx[!duplicated(dx$hgnc_id),]
  dx=dx %>% select(ensembl_gene_id,hgnc_id)
  colnames(dx)=c("GeneSymbol","hgnc_id")
  return(dx)
  
}

##
d_pc_genes=fread(pcfile)
d_conv_gene=fread(convfile_gene)
d_conv_ensm=fread(convfile_ensembl)

#process constraint
d_annot=fread(annotfile)
d_annot$LOEUF_rank=ntile(d_annot$LOEUF,5) #by quantile of LOEUF
d_annot_LOEUF1=d_annot[d_annot$LOEUF_rank==1,]
d_annot_LOEUF2=d_annot[d_annot$LOEUF_rank==2,]
d_annot_LOEUF3=d_annot[d_annot$LOEUF_rank==3,]
d_annot_LOEUF4=d_annot[d_annot$LOEUF_rank==4,]
d_annot_LOEUF5=d_annot[d_annot$LOEUF_rank==5,]

#
d_bed=fread(bedfile,header=F)
colnames(d_bed)=c("chr","start","end","GeneSymbol")

gene_ids=d_bed$GeneSymbol
d_genes=convert_ensembl_to_hgnc(gene_ids,d_pc_genes,d_conv)
d_bed=left_join(d_genes,d_bed)

d_bed=d_bed %>% separate(chr, c("xx", "chr"),"r") %>% select(-xx)
d_bed$chr=as.numeric(d_bed$chr)
d_bed$start=as.numeric(d_bed$start)
d_bed$end=as.numeric(d_bed$end)
d_bed=d_bed %>% arrange(chr,start,end) 
d_bed$chr=paste0("chr",as.character(d_bed$chr))

####

d_out=d_bed[d_bed$hgnc_id %in% d_annot_LOEUF1$hgnc_id,] %>% select(chr,start,end,hgnc_id)
outfile=paste0(outdir,"/all_regions.LOEUF1_genes.txt")
write.table(d_out,outfile,quote=F,sep="\t",row.names=F,col.names=F)

d_out=d_bed[d_bed$hgnc_id %in% d_annot_LOEUF2$hgnc_id,] %>% select(chr,start,end,hgnc_id)
outfile=paste0(outdir,"/all_regions.LOEUF2_genes.txt")
write.table(d_out,outfile,quote=F,sep="\t",row.names=F,col.names=F)

d_out=d_bed[d_bed$hgnc_id %in% d_annot_LOEUF3$hgnc_id,] %>% select(chr,start,end,hgnc_id)
outfile=paste0(outdir,"/all_regions.LOEUF3_genes.txt")
write.table(d_out,outfile,quote=F,sep="\t",row.names=F,col.names=F)

d_out=d_bed[d_bed$hgnc_id %in% d_annot_LOEUF4$hgnc_id,] %>% select(chr,start,end,hgnc_id)
outfile=paste0(outdir,"/all_regions.LOEUF4_genes.txt")
write.table(d_out,outfile,quote=F,sep="\t",row.names=F,col.names=F)

d_out=d_bed[d_bed$hgnc_id %in% d_annot_LOEUF5$hgnc_id,] %>% select(chr,start,end,hgnc_id)
outfile=paste0(outdir,"/all_regions.LOEUF5_genes.txt")
write.table(d_out,outfile,quote=F,sep="\t",row.names=F,col.names=F)


