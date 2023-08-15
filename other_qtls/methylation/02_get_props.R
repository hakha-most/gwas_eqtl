set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

##########
#functions to normalize features for regression

#by variant
normalize_na_rm<-function(v){
  m=mean(v,na.rm=T)
  s=sd(v,na.rm=T)
  v=(v-m)/s
}

#by gene
normalize_na_rm_by_proxy<-function(v,vp){
  m=mean(vp,na.rm=T)
  s=sd(vp,na.rm=T)
  v=(v-m)/s
}

##########

infofile="snp_annotations/filter_snps.txt"
genefile="gene_annotations/all_annots_pc_genes.txt"
multifile="gene_annotations/12_GO_annotations/genes_multiGO.txt"

d_multi=fread(multifile) %>% select(gene,GO_BP_count_400_auc)
d_gene=fread(genefile)
d_gene$gene=d_gene$hgnc_id
d_gene=left_join(d_gene,d_multi)

#colnames(d_gene_annots)[2]="GeneName"

#######

d_info=fread(infofile)
d_all=left_join(d_info,(d_gene %>% select(-TSSD)),by="gene")

###

d_gene_norm=d_gene
d_gene_norm$connect_decile=d_gene_norm$connect_decile*1/1
d_gene_norm$PPI_degree_decile=d_gene_norm$PPI_degree_decile*1/1

d_gene_norm$GO_BP_count_400_auc=normalize_na_rm_by_proxy(d_gene_norm$GO_BP_count_400_auc,d_all$GO_BP_count_400_auc)
d_gene_norm$TF=normalize_na_rm_by_proxy(d_gene_norm$TF,d_all$TF)
d_gene_norm$pLI=normalize_na_rm_by_proxy(d_gene_norm$pLI,d_all$pLI)
d_gene_norm$Roadmap_length_per_type=normalize_na_rm_by_proxy(d_gene_norm$Roadmap_length_per_type,d_all$Roadmap_length_per_type)
d_gene_norm$Roadmap_count=normalize_na_rm_by_proxy(d_gene_norm$Roadmap_count,d_all$Roadmap_count)
d_gene_norm$ABC_count=normalize_na_rm_by_proxy(d_gene_norm$ABC_count,d_all$ABC_count)
d_gene_norm$ABC_length_per_type=normalize_na_rm_by_proxy(d_gene_norm$ABC_length_per_type,d_all$ABC_length_per_type)
d_gene_norm$promoter_count=normalize_na_rm_by_proxy(d_gene_norm$promoter_count,d_all$promoter_count)
d_gene_norm[d_gene_norm$connectedness==1,]$connect_decile=normalize_na_rm_by_proxy(d_gene_norm[d_gene_norm$connectedness==1,]$connect_decile,d_all[d_all$connectedness==1,]$connect_decile)
d_gene_norm[d_gene_norm$PPI_degree_cat==1,]$PPI_degree_decile=normalize_na_rm_by_proxy(d_gene_norm[d_gene_norm$PPI_degree_cat==1,]$PPI_degree_decile,d_all[d_all$PPI_degree_cat==1,]$PPI_degree_decile)


########

d_gwas=fread("other_qtls/methylation/cpgs.txt")
d_gwas$hgnc_id=d_gwas$gene
d_gwas$TSSD_cat=ntile(d_gwas$TSSD,n = 20)

d_gwas_annots_norm=left_join(d_gwas,(d_gene_norm %>% select(-TSSD,-gene)))

###########

d_cpgs=fread("other_qtls/methylation/st8_conditional_analysis.txt") #supp table from Hawe et al.
d_cpgs=d_cpgs[d_cpgs$cpg %in% d_gwas_annots_norm$cpg,]

###########
#logistic regression
dx=d_gwas_annots_norm 
dx$type=0
dx[dx$cpg %in% d_cpgs$cpg,]$type=1

dx=dx %>% select(-colnames(dx)[1:4])

base_count=23

mod1=glm(type~TSSD+factor(TSSD_cat)+length+CDS_length+TF,data=dx,family = "binomial")
TF_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
TF_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

mod1=glm(type~TSSD+factor(TSSD_cat)+length+CDS_length+connect_decile,data=dx[dx$connectedness==1,],family = "binomial")
connect_rank_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
connect_rank_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

mod1=glm(type~TSSD+factor(TSSD_cat)+length+CDS_length+pLI,data=dx,family = "binomial")
pLI_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
pLI_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

mod1=glm(type~TSSD+factor(TSSD_cat)+length+CDS_length+promoter_count,data=dx,family = "binomial")
promoter_count_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
promoter_count_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

mod1=glm(type~TSSD+factor(TSSD_cat)+length+CDS_length+PPI_degree_decile,data=dx[dx$PPI_degree_cat==1,],family = "binomial")
PPI_rank_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
PPI_rank_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

mod1=glm(type~TSSD+factor(TSSD_cat)+length+CDS_length+ABC_length_per_type+ABC_count,data=dx,family = "binomial")
ABC_length_coef_count=unname(summary(mod1)$coef[,1][base_count+1])
ABC_count_coef_length=unname(summary(mod1)$coef[,1][base_count+2])
ABC_length_coef_count_SE=unname(summary(mod1)$coef[,2][base_count+1])
ABC_count_coef_length_SE=unname(summary(mod1)$coef[,2][base_count+2])

mod1=glm(type~TSSD+factor(TSSD_cat)+length+CDS_length+Roadmap_length_per_type+Roadmap_count,data=dx,family = "binomial")
Road_length_coef_count=unname(summary(mod1)$coef[,1][base_count+1])
Road_count_coef_length=unname(summary(mod1)$coef[,1][base_count+2])
Road_length_coef_count_SE=unname(summary(mod1)$coef[,2][base_count+1])
Road_count_coef_length_SE=unname(summary(mod1)$coef[,2][base_count+2])

mod1=glm(type~TSSD+factor(TSSD_cat)+length+CDS_length+GO_BP_count_400_auc,data=dx,family = "binomial")
multiGO_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
multiGO_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

d_coef=data.frame(TF_coef_marginal,connect_rank_coef_marginal,pLI_coef_marginal,promoter_count_coef_marginal,PPI_rank_coef_marginal,ABC_length_coef_count,ABC_count_coef_length,Road_length_coef_count,Road_count_coef_length,multiGO_coef_marginal)
d_coef_SE=data.frame(TF_coef_marginal_SE,connect_rank_coef_marginal_SE,pLI_coef_marginal_SE,promoter_count_coef_marginal_SE,PPI_rank_coef_marginal_SE,ABC_length_coef_count_SE,ABC_count_coef_length_SE,Road_length_coef_count_SE,Road_count_coef_length_SE,multiGO_coef_marginal_SE)
colnames(d_coef_SE)=colnames(d_coef)

d_coef$estimate="coeff"
d_coef_SE$estimate="coeff_SE"

d_out=rbind(d_coef,d_coef_SE)

outfile="other_qtls/methylation/cpg.regression_results.txt"
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)


