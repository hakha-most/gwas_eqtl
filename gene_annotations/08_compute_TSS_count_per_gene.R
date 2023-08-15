########
##compute TSS count per gene from FANTOM5 data

library(data.table)
library(tidyverse)
library(dplyr)

##fantom5 directory is on Zenodo in dir "gene_annots"

### function to convert gene names to HGNC 
convert_name_to_hgnc <- function(gene_ids,d_pc,d_conv){
  
  dx=d_conv[d_conv$Symbol %in% gene_ids,]
  dx=dx[!duplicated(dx$Symbol),]
  
  dx=dx[dx$hgnc_id %in% d_pc$hgnc_id,]
  dx=dx[!duplicated(dx$hgnc_id),]
  
  colnames(dx)=c("gene","hgnc_id")
  return(dx)
  
}

###

d_pc_genes=fread("genes.protein_coding.v39.gtf")
d_conv=fread("ncbi_genes.hgnc_id.txt")

d_gene=fread("fantom5/raw_data/hg19.cage_peak_phase1and2combined_ann.txt.gz",skip = 7)
d_gene <- d_gene %>% select("00Annotation","short_description",hgnc_id)
colnames(d_gene)=c("name","SD","hgnc_id")
d_gene=d_gene %>% separate(SD, c("p", "gene"), sep="@" ,extra = "drop", fill = "right")
d_gene <- d_gene %>% select(name,gene,hgnc_id)

##
#genes with assigned HGNC ID in the data
d_gene1=d_gene[!is.na(d_gene$hgnc_id),]
d_gene1=d_gene1[d_gene1$hgnc_id %in% d_pc_genes$hgnc_id,] #filter protein-coding genes
d_count1=aggregate((d_gene1$hgnc_id), list(d_gene1$hgnc_id), length) %>% arrange(-x) #count regions per gene
colnames(d_count1)=c("hgnc_id","promoter_count")

##
#genes with no assigned HGNC ID in the data
d_gene2=d_gene[is.na(d_gene$hgnc_id),]
d_gene2=d_gene2[d_gene2$gene %in% d_conv$Symbol,]

d_count2=aggregate((d_gene2$gene), list(d_gene2$gene), length) %>% arrange(-x)
colnames(d_count2)=c("gene","promoter_count")

gene_ids=d_count2$gene 
d_genes2=convert_name_to_hgnc(gene_ids,d_pc_genes,d_conv) #conver gene symbol to HGNC ID
d_count2=left_join(d_genes2,d_count2) %>% select(-gene)

####
#combine all
d_count=rbind(d_count1,d_count2)
d_count=d_count[!duplicated(d_count$hgnc_id),]

#assign zero to genes with no annotation in the data
d_count=left_join(d_pc_genes %>% select(hgnc_id),d_count)
d_count[is.na(d_count$promoter_count),]$promoter_count=0

outfile="fantom5/processed/protein_coding_genes.TSS_count.txt"
write.table(d_count,file=outfile,quote=F,sep="\t",row.names=F)

