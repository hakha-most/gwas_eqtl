
set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="extended_gwas_analyses/ukbiobank/gwas_hits.filtered_phenos.clumped_sorted.txt" #GWAS hits for all filtered traits
bloodfile="colocalization/blood_traits.txt" #list of blood-related trats
filterfile="gwas_props/filter_indep_gwas.assoc" #GWAS results for the main 44 traits
infofile="snp_annotations/filter_snps.txt"
genefile="gene_annotations/all_annots_pc_genes.txt"
multifile="gene_annotations/12_GO_annotations/genes_multiGO.txt"

d_info=fread(infofile)
d_multi=fread(multifile) %>% select(gene,GO_BP_count_400_auc)
d_gene=fread(genefile)
d_gene$gene=d_gene$hgnc_id
d_gene=left_join(d_gene,d_multi)

d_all=left_join(d_info,(d_gene %>% select(-TSSD)),by="gene")
d_all$d_tss_signed=(d_all$d_tss)
d_all$d_tss=abs(d_all$d_tss)

#bin distance to TSS, MAF, LD score, and gene density for regression models
d_all$d_tss_cat=ntile(d_all$d_tss,n = 20)
d_all$MAF_cat=ntile(d_all$MAF,n = 20)
d_all$L2_cat=ntile(d_all$L2,n = 20)
d_all$TSSD_cat=ntile(d_all$TSSD,n = 20)

######################
#normalize features for logistic regression model

normalize_na_rm<-function(v){
  m=mean(v,na.rm=T)
  s=sd(v,na.rm=T)
  v=(v-m)/s
}

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

##########
#function to run logistic regression models
compute_reg <- function(d_test_glm,d2){
	
	d_test_glm=d_test_glm[complete.cases(d_test_glm),]
	d1=d_test_glm %>% select(-pheno,-Pval)
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

##########
#split GWAS data by trait category

d_gwas=fread(gwasfile)
d_blood=fread(bloodfile,sep=",")
d_filter=fread(filterfile)

d_gwas_annots=left_join(d_gwas,d_all,by="SNP")
d_gwas_annots_norm=left_join(d_gwas,d_all_norm,by="SNP")

d_gwas_filter=d_gwas_annots_norm[(d_gwas_annots_norm$pheno %in% d_filter$pheno),]
d_gwas_blood=d_gwas_filter[(d_gwas_filter$pheno %in% d_blood$pheno) ,]
d_gwas_not_blood=d_gwas_filter[!(d_gwas_filter$pheno %in% d_blood$pheno),]

######################
#regression models by trait category

d2=d_all_norm[complete.cases(d_all_norm),]
d2=d2[sample(1:nrow(d2),100000,replace=T),]
d2$type=0

d_test_glm=d_gwas_annots_norm
d_out1=compute_reg(d_test_glm,d2)
d_out1$trait_cat="all_traits"

d_test_glm=d_gwas_filter
d_out2=compute_reg(d_test_glm,d2)
d_out2$trait_cat="filter_traits"

d_test_glm=d_gwas_blood
d_out3=compute_reg(d_test_glm,d2)
d_out3$trait_cat="filter_blood_traits"

d_test_glm=d_gwas_not_blood
d_out4=compute_reg(d_test_glm,d2)
d_out4$trait_cat="filter_not_blood_traits"

d_out=rbind(d_out1,d_out2,d_out3,d_out4)


#######

outfile="extended_gwas_analyses/ukbiobank/ukb.regression_results.by_trait.txt"
write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)

