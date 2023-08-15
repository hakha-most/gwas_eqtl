library(data.table)
library(tidyverse)
library(dplyr)

##########

set.seed(0)

#all SNP annots are on Zenodo
d_snp1=fread("snp_annotations/maf01.kgp.snps"); colnames(d_snp1)=c("SNP","rsID") #filtered SNPs, level 1
d_snp2=fread("snp_annotations/filter_snps.txt");  #filtered SNPs, level 2 (level 1 SNPs that are within 1MB of TSSs)

#data on Zenodo
d_gc=fread("supp_note_data/gwas_atlas/gwasATLAS_v20191115_GC.txt") #genetic correlation data
d_trait=fread("supp_note_data/gwas_atlas/gwasATLAS_v20191115.txt") #trait list
d_wata=fread("supp_note_data/gwas_atlas/watanabe_S3_table.txt") #filtered trait list by Watanabe et al.
d_loci=fread("supp_note_data/gwas_atlas/gwasATLAS_v20191115_riskloci.txt") #risk loci

d_disease=fread("extended_gwas_analyses/gwas_atlas/disease.ids",header = F); colnames(d_disease)="id" #disease list extracted from the trait list
d_disease_disorder=fread("extended_gwas_analyses/gwas_atlas/disease_disorder.ids",header = F); colnames(d_disease_disorder)="id" #disease or disorder list extracted from the trait list

#######

d_loci=d_loci[d_loci$rsID %in% d_snp1$rsID,] #filter risk loci

#######
#curate groups of traits

#filter GWAS with European-descent samples
all_pops=sort(unique(d_trait$Population))
eur_pops=all_pops[c(9,31:34)] 
d_trait=d_trait[d_trait$Population %in% eur_pops,]

#compute number of GWAS stduies per trait
tx_temp=aggregate((d_trait$uniqTrait), list((d_trait$uniqTrait)), length)
tx_temp=tx_temp %>% arrange(-x)
colnames(tx_temp)=c("uniqTrait","count")
trait_counts=tx_temp

#compute number of loci per study
tx_temp=aggregate((d_loci$rsID), list(d_loci$id), length)
tx_temp=tx_temp %>% arrange(-x)
colnames(tx_temp)=c("id","hit_count")
loci_counts=tx_temp

#for each trait select the GWAS study with highest number of hits 
dx=d_trait %>% select(id,uniqTrait)
dx=left_join(dx,trait_counts,by="uniqTrait")
dx=left_join(dx,loci_counts,by="id")
dx[is.na(dx$hit_count),]$hit_count=0
dx=dx %>% arrange(-hit_count)
dx_no_rep=dx[!duplicated(dx$uniqTrait),]

#grouping
d_trait_no_rep=d_trait[d_trait$id %in% dx_no_rep$id,]
d_disease_no_rep=d_disease[d_disease$id %in% dx_no_rep$id,]
d_disease_disorder_no_rep=d_disease_disorder[d_disease_disorder$id %in% dx_no_rep$id,]

#######
#filter trait choices in Watanabe et al.
d_wata_yes=d_wata[d_wata$Selected=="Yes",]

#grouping
d_trait_filter=d_trait[d_trait$id %in% d_wata_yes$`atlas ID`,]
d_trait_disease=d_trait_filter[d_trait_filter$id %in% d_disease$id,]
d_trait_disease_disorder=d_trait_filter[d_trait_filter$id %in% d_disease_disorder$id,]

########
#function to prune traits based on genetic correlation data such
get_indep_phenos <- function(d_target,d_loci,d_gc,rg_cutoff){
  
  d_gc_filter=d_gc[(d_gc$id1 %in% d_target$id) & (d_gc$id2 %in% d_target$id),]
  
  d_gc1=d_gc_filter[(d_gc_filter$rg>0) & ((d_gc_filter$rg - 2*d_gc_filter$se) > rg_cutoff),]
  d_gc2=d_gc_filter[(d_gc_filter$rg<0) & ((d_gc_filter$rg + 2*d_gc_filter$se) < (-rg_cutoff)),]
  
  d_gc_filter=rbind(d_gc1,d_gc2)
  
  ####
  
  tx_temp=aggregate((d_loci$rsID), list(d_loci$id), length)
  tx_temp=tx_temp %>% arrange(-x)
  colnames(tx_temp)=c("id","count")
  
  tx_temp_filter=tx_temp[tx_temp$id %in% d_target$id,]
  
  done_prune=0
  pruned_tx_temp=tx_temp_filter
  i=0
  while (done_prune==0){
    
    i=i+1
    print(i)
    pheno_temp=pruned_tx_temp$id[i]
    d_rg_temp=d_gc_filter[(d_gc_filter$id1==pheno_temp | d_gc_filter$id2==pheno_temp),]
    phenos_to_remove=unique(c(d_rg_temp$id1,d_rg_temp$id2)); phenos_to_remove=phenos_to_remove[phenos_to_remove!=pheno_temp]
    pruned_tx_temp=pruned_tx_temp[!(pruned_tx_temp$id %in% phenos_to_remove),]
    if (nrow(pruned_tx_temp)==i){done_prune=1}
    
  }
  
  return(pruned_tx_temp$id)
  
}

#prune traits such that no pair of traits have genetic correlation >0.5
indep_filter_phenos=get_indep_phenos(d_trait_filter,d_loci,d_gc,0.5) #ALL filtered traits
indep_filter_disease=get_indep_phenos(d_trait_disease,d_loci,d_gc,0.5) #ALL filtered traits
indep_filter_disease_disorder=get_indep_phenos(d_trait_disease_disorder,d_loci,d_gc,0.5) #ALL filtered traits

######
#output number of loci per trait group
d_test=d_loci[d_loci$rsID %in% d_snp2$rsID,]

nrow(d_test[(d_test$id %in% d_trait_no_rep$id),])
nrow(d_test[(d_test$id %in% d_disease_no_rep$id),])
nrow(d_test[(d_test$id %in% d_disease_disorder_no_rep$id),])

nrow(d_test[(d_test$id %in% indep_filter_phenos),])
nrow(d_test[(d_test$id %in% indep_filter_disease),])
nrow(d_test[(d_test$id %in% indep_filter_disease_disorder),])

#output number of traits per trait group

c=d_test[(d_test$id %in% d_trait_no_rep$id),]
length(unique(c$id))

c=d_test[(d_test$id %in% d_disease_no_rep$id),]
length(unique(c$id))

c=d_test[(d_test$id %in% d_disease_disorder_no_rep$id),]
length(unique(c$id))

c=d_test[(d_test$id %in% indep_filter_phenos),]
length(unique(c$id))

c=d_test[(d_test$id %in% indep_filter_disease),]
length(unique(c$id))

c=d_test[(d_test$id %in% indep_filter_disease_disorder),]
length(unique(c$id))

######

d_out=d_trait %>% select(id,Domain,ChapterLevel,uniqTrait)

d_out$no_rep_phenos=0
d_out$no_rep_disease=0
d_out$no_rep_disease_disorder=0
d_out$indep_phenos=0
d_out$indep_disease=0
d_out$indep_disease_disorder=0

d_out[d_out$id %in% d_trait_no_rep$id,]$no_rep_phenos=1
d_out[d_out$id %in% d_disease_no_rep$id,]$no_rep_disease=1
d_out[d_out$id %in% d_disease_disorder_no_rep$id,]$no_rep_disease_disorder=1
d_out[d_out$id %in% indep_filter_phenos,]$indep_phenos=1
d_out[d_out$id %in% indep_filter_disease,]$indep_disease=1
d_out[d_out$id %in% indep_filter_disease_disorder,]$indep_disease_disorder=1

write.table(d_out,file="extended_gwas_analyses/gwas_atlas/atlas.phenos.txt",quote=F,sep="\t",row.names=F)

######
#divide independent trait groups by type 

dx=d_out[d_out$indep_phenos==1, ]
dx=dx[!(dx$id %in% d_disease$id),]
dx=dx[!(dx$id %in% d_disease_disorder$id),]

dx_cog=dx[dx$Domain=="Cognitive" | dx$Domain=="Neurological" | dx$Domain=="Psychiatric" | dx$uniqTrait=="Educational attainment",]
dx_rep=dx[dx$Domain=="Reproduction" | dx$uniqTrait=="Birth weight",]
dx_envi=dx[dx$Domain=="Environment" & dx$uniqTrait!="Educational attainment",]
dx_phys=dx[(dx$Domain=="Skeletal" | dx$Domain=="Body Structures") & dx$uniqTrait!="Birth weight",]
dx_metal=dx[dx$Domain=="Metabolic" & dx$uniqTrait!="Birth weight",] 
dx_immun=dx[dx$Domain=="Immunological",] 
dx_disease=d_out[d_out$indep_disease==1, ]
dx_disease_disorder=d_out[d_out$indep_disease_disorder==1, ]

dx_cog$cat="Cognitive"
dx_rep$cat="Reproduction"
dx_metal$cat="Metabolic"
dx_phys$cat="Physical"
dx_immun$cat="Immunological"
dx_disease$cat="Disease"
dx_disease_disorder$cat="Disease, Disorder"

d_out_cats=rbind(dx_cog,dx_rep,dx_metal,dx_phys,dx_immun,dx_disease,dx_disease_disorder)

write.table(d_out_cats,file="extended_gwas_analyses/gwas_atlas/atlas.phenos.indep_cats.txt",quote=F,sep="\t",row.names=F)

##########
#curate tables to release on Zenodo

d_pheno_x=fread("extended_gwas_analyses/gwas_atlas/atlas.phenos.txt")
d_pheno_y=fread("extended_gwas_analyses/gwas_atlas/atlas.phenos.indep_cats.txt")
d_pheno_y=d_pheno_y[d_pheno_y$cat!="Disease, Disorder"] %>% select(id,cat)
colnames(d_pheno_y)=c("id","category")

d_pheno=d_pheno_x %>% select(-Domain,-ChapterLevel)
colnames(d_pheno)=c("id","trait","all_traits","all_diseases","all_diseases_disorders","indep_traits","indep_diseases","indep_diseases_disorders")
d_pheno=left_join(d_pheno,d_pheno_y)

write.table(d_pheno,file="extended_gwas_analyses/gwas_atlas/zenodo_phenos.txt",quote=F,sep="\t",row.names=F)

#
d_risk=d_loci[(d_loci$rsID %in% d_snp2$rsID),] %>% arrange(id,p) %>% dplyr::select(id,rsID,p)
write.table(d_risk,file="extended_gwas_analyses/gwas_atlas/zenodo_risk_loci.txt",quote=F,sep="\t",row.names=F)


