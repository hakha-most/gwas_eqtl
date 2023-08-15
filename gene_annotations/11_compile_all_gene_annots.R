
## complie all gene features into one file

library(tidyverse)
library(data.table)
library(dplyr)

set.seed(1)

####The following directories are on Zenode in dir "gene_annots"
#ABC
#fantom5
#saha_et_al
#Liu_et_al_roadmap_gene_enhancer_links


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


### to convert gene names to HGNC 
convert_name_to_hgnc <- function(gene_ids,d_pc,d_conv){
  
  dx=d_conv[d_conv$Symbol %in% gene_ids,]
  dx=dx[!duplicated(dx$Symbol),]
  
  dx=dx[dx$hgnc_id %in% d_pc$hgnc_id,]
  dx=dx[!duplicated(dx$hgnc_id),]
  
  colnames(dx)=c("gene","hgnc_id")
  return(dx)
  
}

#====

#conversion data, curated via NCBI gene info and BioMart

d_conv_gene=fread("ncbi_genes.hgnc_id.txt")
d_conv_ensm=fread("ensembl_versions_lookup.txt")

#====

#GENCODE protein-coding genes, v39 in GRCH37
d_pc_genes=fread("genes.protein_coding.v39.gtf")
d_pc_genes$length=abs(d_pc_genes$tss-d_pc_genes$tes)

#====

#length of coding sequence from Ensembl BioMart download (GRCH38) 
d_ensembl=fread("ensemble_biomart_cds_length_grch38.processed.txt")

#====

#TSS density around each gene
d_TSSD=fread("genic_density_around_TSS.txt")

#====

#pLI and LOEUF data from gnomAD
d_gnomad=fread("gnomAD.txt")

gene_ids=d_gnomad$gene
d_genes=convert_name_to_hgnc(gene_ids,d_pc_genes,d_conv_gene)

d_gnomad=left_join(d_genes,d_gnomad) %>% select(-gene)

##

#hs estimates from Agarwal et al. 2022
d_hs=fread("hs.txt")
d_hs$hs=10^(d_hs$log10_map)
d_hs=d_hs %>% select(Gene,hs)
colnames(d_hs)=c("gene","hs")

gene_ids=d_hs$gene
d_genes=convert_name_to_hgnc(gene_ids,d_pc_genes,d_conv_gene)

d_hs=left_join(d_genes,d_hs) %>% select(-gene)

##

#promoter counts per gene from FANTOM5
d_promoter=fread("fantom5/processed/protein_coding_genes.TSS_count.txt")

#enhancer features per gene from ABC paper
d_ABC=fread("ABC/processed/ABC_enh_counts_convert_ids.txt")
d_ABC$ABC_length_per_type=d_ABC$ABC_bp_per_type/1000 #bp to Kb
d_ABC=d_ABC %>% select(hgnc_id,ABC_count1,ABC_length_per_type)
colnames(d_ABC)=c("hgnc_id","ABC_count","ABC_length_per_type")

#enhancer features per gene from Liu et al.
d_Road=fread("Liu_et_al_roadmap_gene_enhancer_links/processed/Roadmap_enh_props.txt")
d_Road$Roadmap_length_per_type=d_Road$Roadmap_bp_per_type/1000 #bp to Kb
d_Road=d_Road %>% select(hgnc_id,Roadmap_count1,Roadmap_length_per_type)
colnames(d_Road)=c("hgnc_id","Roadmap_count","Roadmap_length_per_type")

#connectedness rank from co-expression networks (cross-tissue) from Saha et al.
d_connect=fread("saha_et_al/twns/processed/genes_by_combined_rank.txt")
d_connect$rank=1:nrow(d_connect)
d_connect$connect_decile=ntile(1/d_connect$rank,n = 10)
d_connect$connect_quantile=ntile(1/d_connect$rank,n = 5)
d_connect=d_connect %>% select(hgnc_id,connect_decile,connect_quantile)
d1=merge(d_connect,d_pc_genes,by="hgnc_id") %>% select(hgnc_id,connect_decile,connect_quantile) #include genes with no connections
d2=d_pc_genes[!(d_pc_genes$hgnc_id %in% d_connect$hgnc_id),] %>% select(hgnc_id)
d2$connect_decile=0; d2$connect_quantile=0
d1$connectedness=1; d2$connectedness=0
d_connect=rbind(d1,d2)

#Transcription factors
d_TF=fread("TF_genes.txt")
colnames(d_TF)=c("GeneSymbol")

gene_ids=d_TF$GeneSymbol
d_genes=convert_ensembl_to_hgnc(gene_ids,d_pc_genes,d_conv_ensm)

d_TF=left_join(d_genes,d_TF) %>% select(-GeneSymbol)

#connectedness in InWeb protein-protein interaction networks
d_PPI=fread("PPI_degree.txt"); d_PPI=d_PPI %>% select(-gene,-PPI_degree)


##
#aggregate all
d_genes=left_join(d_pc_genes,d_ensembl)
d_genes=d_genes %>% select(hgnc_id,GeneSymbol,gene,length,CDS_length)
d_genes=left_join(d_genes,d_TSSD,by="hgnc_id")
d_genes=left_join(d_genes,d_gnomad,by="hgnc_id")
d_genes=left_join(d_genes,d_ABC,by="hgnc_id")
d_genes=left_join(d_genes,d_Road,by="hgnc_id")
d_genes=left_join(d_genes,d_promoter,by="hgnc_id")
d_genes=left_join(d_genes,d_connect,by="hgnc_id")
d_genes=left_join(d_genes,d_PPI,by="hgnc_id")
d_genes$TF=0; d_genes[d_genes$hgnc_id %in% d_TF$hgnc_id,]$TF=1
d_genes=left_join(d_genes,d_hs,by="hgnc_id")

#====

d_genes$length=d_genes$length/1000
d_genes$CDS_length=d_genes$CDS_length/1000

d_genes$HI=d_genes$pLI
d_genes[!is.na(d_genes$pLI) & d_genes$pLI>0.9,]$HI=1
d_genes[!is.na(d_genes$pLI) & d_genes$pLI<=0.9,]$HI=0

#====

d_out=d_genes[!duplicated(d_genes$hgnc_id),]
outfile="all_annots_pc_genes.txt"
write.table(d_out,file=outfile,quote = FALSE,sep="\t",row.names = F)



