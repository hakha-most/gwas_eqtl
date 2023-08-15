
## Script to compute basic properties of eQTL SNPs similar to "04_calc_properties.R", but for each tissue separately

# This script, as provided, is run on eQTL SNPs. 
# The same script was run using 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

eqtlfile="eqtl_props/filter_eqtls.pruned_tissues.assoc" 
infofile="snp_annotations/filter_snps.txt"
genefile="gene_annotations/all_annots_pc_genes.txt"
outfile1="eqtl_props/eqtl.SNP_count.by_tissue"
outfile2="eqtl_props/eqtl.basic_props.by_tissue"

d_info=fread(infofile)
d_gene=fread(genefile)
d_gene$gene=d_gene$hgnc_id

d_gene$LOEUF_cat=ntile(d_gene$LOEUF,5)
d_gene$Road_length_cat=ntile(d_gene$Roadmap_length_per_type,5)
d_gene$Road_count_cat=ntile(d_gene$Roadmap_count,5)
d_gene$ABC_length_cat=ntile(d_gene$ABC_length_per_type,5)
d_gene$ABC_count_cat=ntile(d_gene$ABC_count,5)
d_gene$TSS_cat=ntile(d_gene$promoter_count,5)

d_gene$tRoad_length=d_gene$Road_length_cat; d_gene[d_gene$Road_length_cat==5,]$tRoad_length=1; d_gene[d_gene$Road_length_cat<5,]$tRoad_length=0
d_gene$tRoad_count=d_gene$Road_count_cat; d_gene[d_gene$Road_count_cat==5,]$tRoad_count=1; d_gene[d_gene$Road_count_cat<5,]$tRoad_count=0
d_gene$tABC_length=d_gene$ABC_length_cat; d_gene[d_gene$ABC_length_cat==5,]$tABC_length=1; d_gene[d_gene$ABC_length_cat<5,]$tABC_length=0
d_gene$tABC_count=d_gene$ABC_count_cat; d_gene[d_gene$ABC_count_cat==5,]$tABC_count=1; d_gene[d_gene$ABC_count_cat<5,]$tABC_count=0

d_gene$tLOEUF=d_gene$LOEUF_cat; d_gene[d_gene$LOEUF_cat==1,]$tLOEUF=1; d_gene[d_gene$LOEUF_cat>1,]$tLOEUF=0
d_gene$tTSS=d_gene$TSS_cat; d_gene[d_gene$TSS_cat==5,]$tTSS=1; d_gene[d_gene$TSS_cat<5,]$tTSS=0
d_gene$tconnect=d_gene$connect_quantile; d_gene[d_gene$connect_quantile>=4,]$tconnect=1; d_gene[d_gene$connect_quantile<4,]$tconnect=0
d_gene$tPPI=d_gene$PPI_degree_quantile; d_gene[d_gene$PPI_degree_quantile==5,]$tPPI=1; d_gene[d_gene$PPI_degree_quantile<5,]$tPPI=0

d_all=left_join(d_info,(d_gene %>% select(-TSSD)),by="gene") #merge SNP and genic features

d_eqtl=fread(eqtlfile)
d_eqtl=d_eqtl[d_eqtl$Pval<(5e-8),] #restrict to top eQTLs using a p-value cutoff matching GWAS

d_eqtl_annots=left_join(d_eqtl,d_all,by="SNP")

#####

tx_temp=aggregate((d_eqtl_annots$SNP), list(d_eqtl_annots$Tissue), length)
tx_temp=tx_temp %>% arrange(-x)
colnames(tx_temp)=c("tissue","count")

write.table(tx_temp,file=outfile1,quote=F,sep=",",row.names=F)

#####

props_HI=aggregate((d_eqtl_annots$HI), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)
props_TF=aggregate((d_eqtl_annots$TF), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)
props_tTSS=aggregate((d_eqtl_annots$tTSS), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)
props_tLOEUF=aggregate((d_eqtl_annots$tLOEUF), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)
props_tconnect=aggregate((d_eqtl_annots$tconnect), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)
props_tPPI=aggregate((d_eqtl_annots$tPPI), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)

props_tRoad_length=aggregate((d_eqtl_annots$tRoad_length), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)
props_tRoad_count=aggregate((d_eqtl_annots$tRoad_count), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)
props_tABC_length=aggregate((d_eqtl_annots$tABC_length), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)
props_tABC_count=aggregate((d_eqtl_annots$tABC_count), list(d_eqtl_annots$Tissue), mean,na.rm=TRUE) %>% arrange(Group.1)


props=data.frame(tissue=props_HI$Group.1,prop_HI=props_HI$x,prop_TF=props_TF$x,prop_tTSS=props_tTSS$x,prop_tLOEUF=props_tLOEUF$x,prop_tconnect=props_tconnect$x,prop_tPPI=props_tPPI$x,prop_tRoad_length=props_tRoad_length$x,prop_tRoad_count=props_tRoad_count$x,prop_ABC_length=props_tABC_length$x,prop_ABC_count=props_tABC_count$x)

d_out=data.frame(t(props %>% select(-tissue)))
colnames(d_out)=props$tissue
d_out$annot=rownames(d_out)

write.table(d_out,file=outfile2,quote=F,sep=",",row.names=F)











