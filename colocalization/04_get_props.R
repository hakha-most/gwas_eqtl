set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="colocalization/coloc_status.txt"
infofile="snp_annotations/filter_snps.txt"
genefile="gene_annotations/all_annots_pc_genes.txt"
multifile="gene_annotations/12_GO_annotations/genes_multiGO.txt"
hitfile="colocalization/blood_trait.snps.txt"

d_multi=fread(multifile) %>% select(gene,GO_BP_count_400_auc)
d_gene=fread(genefile)
d_gene$gene=d_gene$hgnc_id
d_gene=left_join(d_gene,d_multi)

d_info=fread(infofile)
d_all=left_join(d_info,(d_gene %>% select(-TSSD)),by="gene")
d_all$d_tss_signed=(d_all$d_tss)
d_all$d_tss=abs(d_all$d_tss)
d_all$d_tss_cat=ntile(d_all$d_tss,n = 20)
d_all$MAF_cat=ntile(d_all$MAF,n = 20)
d_all$L2_cat=ntile(d_all$L2,n = 20)
d_all$TSSD_cat=ntile(d_all$TSSD,n = 20)

######################
#normalize features for logistic regression model

normalize_na_rm_by_proxy<-function(v,vp){
	m=mean(vp,na.rm=T)
	s=sd(vp,na.rm=T)
	v=(v-m)/s
}

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

##

d_hit=fread(hitfile)
colnames(d_hit)=c("SNP","d_tss","gene")

d_all=d_info
d_all$d_tss_signed=(d_all$d_tss)
d_all$d_tss=abs(d_all$d_tss)
d_all$d_tss_cat=ntile(d_all$d_tss,n = 20)
d_all$MAF_cat=ntile(d_all$MAF,n = 20)
d_all$L2_cat=ntile(d_all$L2,n = 20)
d_all$TSSD_cat=ntile(d_all$TSSD,n = 20)

d_all=d_all[d_all$SNP %in% d_hit$SNP,]
d_all=d_all %>% select(-d_tss,-gene)
d_all=left_join(d_all,d_hit,by="SNP")
d_all=left_join(d_all,(d_gene_norm %>% select(-TSSD)),by="gene")

d_all_norm=d_all
###

d_gwas=fread(gwasfile)
colnames(d_gwas)[1]="SNP"

##########
#function to run logistic regression models

compute_reg <- function(dx){
	
	dx=dx[complete.cases(dx),]
	
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
#GTEx regression

snp1=d_gwas[is.na(d_gwas$GTEx_P),]$SNP #GWAS hit without linked eQTL
snp2=d_gwas[!is.na(d_gwas$GTEx_P),]$SNP #GWAS hit with linked eQTL

d1=d_all_norm[d_all_norm$SNP %in% snp1,]
d2=d_all_norm[d_all_norm$SNP %in% snp2,]

d1$type=0
d2$type=1
dx=rbind(d1,d2)

d_out=compute_reg(dx)
d_out$cat="GTEx"

######################
#eQTLGen regression

gen_cols=paste0("eQTL_P",as.character(1:10))

for (k in 1:10){

print(k)

d_temp=d_gwas %>% select(SNP,gen_cols[k])
colnames(d_temp)=c("SNP","P")

snp1=d_temp[is.na(d_temp$P),]$SNP
snp2=d_temp[!is.na(d_temp$P),]$SNP

d1=d_all_norm[d_all_norm$SNP %in% snp1,]
d2=d_all_norm[d_all_norm$SNP %in% snp2,]

d1$type=0
d2$type=1
dx=rbind(d1,d2)

d_out_temp=compute_reg(dx)
d_out_temp$cat=gen_cols[k]
d_out=rbind(d_out,d_out_temp)

}


#######

outfile="colocalization/coloc_vs_no_coloc.regression_results.txt"
write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)


