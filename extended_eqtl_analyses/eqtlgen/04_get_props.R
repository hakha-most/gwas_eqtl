
set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="supp_note_data/eqtlgen/eqtlgen_eqtls.clumped.txt" #data on Zenodo
infofile="snp_annotations/filter_snps.txt" #data on Zenodo
genefile="gene_annotations/all_annots_pc_genes.txt"
multifile="gene_annotations/12_GO_annotations/genes_multiGO.txt"
outfile="extended_eqtl_analyses/eqtlgen/eqtlgen.regression_results.txt"

#####
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

#####

d_multi=fread(multifile) %>% select(gene,GO_BP_count_400_auc)
colnames(d_multi)[1]="hgnc_id"

d_gene=fread(genefile)
d_gene=left_join(d_gene,d_multi)
d_gene$gene=d_gene$hgnc_id

d_info=fread(infofile)

d_all=left_join(d_info,(d_gene %>% select(-TSSD)),by="gene")
d_all$d_tss_signed=(d_all$d_tss)
d_all$d_tss=abs(d_all$d_tss)
d_all$d_tss_cat=ntile(d_all$d_tss,n = 20)
d_all$MAF_cat=ntile(d_all$MAF,n = 20)
d_all$L2_cat=ntile(d_all$L2,n = 20)
d_all$TSSD_cat=ntile(d_all$TSSD,n = 20)

#########
#normalize features for logistic regression model
d_all_norm=d_all
d_all_norm$connect_decile=d_all_norm$connect_decile*1/1
d_all_norm$PPI_degree_decile=d_all_norm$PPI_degree_decile*1/1

d_all_norm$GO_BP_count_400_auc=normalize_na_rm(d_all_norm$GO_BP_count_400_auc)
d_all_norm$TF=normalize_na_rm(d_all_norm$TF)
d_all_norm$pLI=normalize_na_rm(d_all_norm$pLI)
d_all_norm$Roadmap_length_per_type=normalize_na_rm(d_all_norm$Roadmap_length_per_type)
d_all_norm$Roadmap_count=normalize_na_rm(d_all_norm$Roadmap_count)
d_all_norm$ABC_count=normalize_na_rm(d_all_norm$ABC_count)
d_all_norm$ABC_length_per_type=normalize_na_rm(d_all_norm$ABC_length_per_type)
d_all_norm$promoter_count=normalize_na_rm(d_all_norm$promoter_count)
d_all_norm[d_all_norm$connectedness==1,]$connect_decile=normalize_na_rm(d_all_norm[d_all_norm$connectedness==1,]$connect_decile)
d_all_norm[d_all_norm$PPI_degree_cat==1,]$PPI_degree_decile=normalize_na_rm(d_all_norm[d_all_norm$PPI_degree_cat==1,]$PPI_degree_decile)

#

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

#####

d_gwas=fread(gwasfile,header=F)
colnames(d_gwas)=c("eGeneSymbol","SNP","Pval","rank")
d_gwas=d_gwas[d_gwas$SNP %in% d_all$SNP,]
d_gwas=d_gwas %>% arrange(Pval)

#deciles of p-value
d_gwas$Pval_cat=ntile(d_gwas$Pval,10)

d_gwas_annots=left_join(d_gwas,d_all,by="SNP")
d_gwas_annots_norm=left_join(d_gwas,d_all_norm,by="SNP")

##
#link SNP-eGenes to gene properties
d1=d_all_norm %>% select(SNP,rsID,MAF,L2,TSSD,d_tss,gene)
d2=d_all_norm %>% select(SNP,d_tss_signed,d_tss_cat,MAF_cat,L2_cat,TSSD_cat)

dx=d_gwas_annots_norm %>% select(colnames(d_gwas_annots_norm)[1:5])
dx=left_join(dx,d1,by="SNP")
dx$GeneSymbol=dx$eGeneSymbol
dx=left_join(dx,(d_gene_norm %>% select(-gene,-TSSD)),by="GeneSymbol")
dx=left_join(dx,d2,by="SNP")
d_eqtl_eGene=dx %>% arrange(Pval)

######################
#function to run logistic regression models

compute_reg <- function(d_test_glm,d2){
	
	d_test_glm=d_test_glm[complete.cases(d_test_glm),]
	d1=d_test_glm %>% select(-eGeneSymbol,-Pval,-Pval_cat,-rank)
	d1$type=1
	dx=rbind(d1,d2)

	base_count=83

	mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+TF,data=dx,family = "binomial")
	TF_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
	TF_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

	mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+connect_decile,data=dx[dx$connectedness==1,],family = "binomial")
	connect_rank_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
	connect_rank_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

	mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+pLI,data=dx,family = "binomial")
	pLI_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
	pLI_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

	mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+promoter_count,data=dx,family = "binomial")
	promoter_count_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
	promoter_count_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

	mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+PPI_degree_decile,data=dx[dx$PPI_degree_cat==1,],family = "binomial")
	PPI_rank_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
	PPI_rank_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])

	mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+ABC_length_per_type+ABC_count,data=dx,family = "binomial")
	ABC_length_coef_count=unname(summary(mod1)$coef[,1][base_count+1])
	ABC_count_coef_length=unname(summary(mod1)$coef[,1][base_count+2])
	ABC_length_coef_count_SE=unname(summary(mod1)$coef[,2][base_count+1])
	ABC_count_coef_length_SE=unname(summary(mod1)$coef[,2][base_count+2])

	mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+Roadmap_length_per_type+Roadmap_count,data=dx,family = "binomial")
	Road_length_coef_count=unname(summary(mod1)$coef[,1][base_count+1])
	Road_count_coef_length=unname(summary(mod1)$coef[,1][base_count+2])
	Road_length_coef_count_SE=unname(summary(mod1)$coef[,2][base_count+1])
	Road_count_coef_length_SE=unname(summary(mod1)$coef[,2][base_count+2])

	mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+GO_BP_count_400_auc,data=dx,family = "binomial")
	multiGO_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
	multiGO_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])
	#

	d_coef=data.frame(TF_coef_marginal,connect_rank_coef_marginal,pLI_coef_marginal,promoter_count_coef_marginal,PPI_rank_coef_marginal,ABC_length_coef_count,ABC_count_coef_length,Road_length_coef_count,Road_count_coef_length,multiGO_coef_marginal)
	d_coef_SE=data.frame(TF_coef_marginal_SE,connect_rank_coef_marginal_SE,pLI_coef_marginal_SE,promoter_count_coef_marginal_SE,PPI_rank_coef_marginal_SE,ABC_length_coef_count_SE,ABC_count_coef_length_SE,Road_length_coef_count_SE,Road_count_coef_length_SE,multiGO_coef_marginal_SE)
	colnames(d_coef_SE)=colnames(d_coef)

	d_coef$estimate="coeff"
	d_coef_SE$estimate="coeff_SE"

	d_out=rbind(d_coef,d_coef_SE)
	return(d_out)
}

######################
#regression models by decile of p-value

d2=d_all_norm[complete.cases(d_all_norm),]
d2=d2[sample(1:nrow(d2),100000,replace=T),]
d2$type=0

#

pcat=1
d_test_glm=d_eqtl_eGene[d_eqtl_eGene$Pval_cat==pcat,]
d_out_temp=compute_reg(d_test_glm,d2)
d_out_temp$Pval_cat=pcat
d_out=d_out_temp

for (pcat in 2:10){
	print(pcat)
	d_test_glm=d_eqtl_eGene[d_eqtl_eGene$Pval_cat==pcat,]
	d_out_temp=compute_reg(d_test_glm,d2)
	d_out_temp$Pval_cat=pcat
	d_out=rbind(d_out,d_out_temp)	
}

write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)



