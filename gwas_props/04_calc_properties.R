
## Script to compute basic properties of GWAS SNPs.

# This script, as provided, is run on GWAS SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

gwasfile="gwas_props/filter_indep_gwas.assoc" 
infofile="snp_annotations/filter_snps.txt"
genefile="gene_annotations/all_annots_pc_genes.txt"
outfile="gwas_props/gwas.basic_props"

d_info=fread(infofile) #SNP features
d_gene=fread(genefile) #genic features
d_gene$gene=d_gene$hgnc_id

#bin continuous genic features
d_gene$LOEUF_cat=ntile(d_gene$LOEUF,5)
d_gene$Road_length_cat=ntile(d_gene$Roadmap_length_per_type,5)
d_gene$Road_count_cat=ntile(d_gene$Roadmap_count,5)
d_gene$ABC_length_cat=ntile(d_gene$ABC_length_per_type,5)
d_gene$ABC_count_cat=ntile(d_gene$ABC_count,5)
d_gene$TSS_cat=ntile(d_gene$promoter_count,5)

d_all=left_join(d_info,(d_gene %>% select(-TSSD)),by="gene")
d_all$d_tss_signed=(d_all$d_tss)
d_all$d_tss=abs(d_all$d_tss)

#bin distance to TSS, MAF, LD score, and gene density for regression models
d_all$d_tss_cat=ntile(d_all$d_tss,n = 20)
d_all$MAF_cat=ntile(d_all$MAF,n = 20)
d_all$L2_cat=ntile(d_all$L2,n = 20)
d_all$TSSD_cat=ntile(d_all$TSSD,n = 20)

d_gwas=fread(gwasfile)
d_gwas_annots=left_join(d_gwas,d_all,by="SNP")

########

#compute props
d_test=d_gwas_annots 

prop_HI=mean(d_test$HI,na.rm=T)
prop_TF=mean(d_test$TF,na.rm=T)
prop_connect_top_decile=nrow(d_test[d_test$connect_decile==10,])/nrow(d_test)
prop_connect_top_half=nrow(d_test[d_test$connect_decile>=5,])/nrow(d_test)
mean_promoter_count=mean(d_test$promoter_count,na.rm=T)
mean_connectedness=mean(d_test$connectedness,na.rm=T)

d_x1=d_test[d_test$connectedness==1,] %>% group_by(connect_quantile) %>% count()
prop_connect_qs=t(d_x1$n/nrow(d_test))
colnames(prop_connect_qs)=paste0("prop_connect_q",as.character(1:5))

d_x1=d_test[d_test$PPI_degree_cat==1,] %>% group_by(PPI_degree_quantile) %>% count()
prop_PPI_qs=t(d_x1$n/nrow(d_test))
colnames(prop_PPI_qs)=paste0("prop_PPI_q",as.character(1:5))

d_x1=d_test[!is.na(d_test$LOEUF_cat),] %>% group_by(LOEUF_cat) %>% count()
prop_LOEUF_qs=t(d_x1$n/nrow(d_test))
colnames(prop_LOEUF_qs)=paste0("prop_LOEUF_q",as.character(1:5))

d_x1=d_test %>% group_by(ABC_length_cat) %>% count()
prop_ABC_length_qs=t(d_x1$n/nrow(d_test))
colnames(prop_ABC_length_qs)=paste0("prop_ABC_length_qs",as.character(1:5))

d_x1=d_test %>% group_by(ABC_count_cat) %>% count()
prop_ABC_count_qs=t(d_x1$n/nrow(d_test))
colnames(prop_ABC_count_qs)=paste0("prop_ABC_count_qs",as.character(1:5))

d_x1=d_test %>% group_by(Road_length_cat) %>% count()
prop_Roadmap_length_qs=t(d_x1$n/nrow(d_test))
colnames(prop_Roadmap_length_qs)=paste0("prop_Roadmap_length_qs",as.character(1:5))

d_x1=d_test %>% group_by(Road_count_cat) %>% count()
prop_Roadmap_count_qs=t(d_x1$n/nrow(d_test))
colnames(prop_Roadmap_count_qs)=paste0("prop_Roadmap_count_qs",as.character(1:5))

d_out1=data.frame(prop_HI,prop_TF,prop_connect_top_decile,prop_connect_top_half,mean_promoter_count,mean_connectedness,prop_connect_qs,prop_PPI_qs,prop_LOEUF_qs,prop_Roadmap_length_qs,prop_Roadmap_count_qs,prop_ABC_length_qs,prop_ABC_count_qs)

######################
######################
######################

#logistic regression models

#normalize all features
normalize_na_rm<-function(v){
	m=mean(v,na.rm=T)
	s=sd(v,na.rm=T)
	v=(v-m)/s
}

d_all_norm=d_all
d_all_norm$connect_decile=d_all_norm$connect_decile*1/1
d_all_norm$PPI_degree_decile=d_all_norm$PPI_degree_decile*1/1

d_all_norm$TF=normalize_na_rm(d_all_norm$TF)
d_all_norm$pLI=normalize_na_rm(d_all_norm$pLI)
d_all_norm$Roadmap_length_per_type=normalize_na_rm(d_all_norm$Roadmap_length_per_type)
d_all_norm$Roadmap_count=normalize_na_rm(d_all_norm$Roadmap_count)
d_all_norm$ABC_count=normalize_na_rm(d_all_norm$ABC_count)
d_all_norm$ABC_length_per_type=normalize_na_rm(d_all_norm$ABC_length_per_type)
d_all_norm$promoter_count=normalize_na_rm(d_all_norm$promoter_count)
d_all_norm[d_all_norm$connectedness==1,]$connect_decile=normalize_na_rm(d_all_norm[d_all_norm$connectedness==1,]$connect_decile)
d_all_norm[d_all_norm$PPI_degree_cat==1,]$PPI_degree_decile=normalize_na_rm(d_all_norm[d_all_norm$PPI_degree_cat==1,]$PPI_degree_decile)

d_gwas_annots_norm=left_join(d_gwas,d_all_norm,by="SNP")

#

d_test_glm=d_gwas_annots_norm[complete.cases(d_gwas_annots_norm),] # GWAS SNPs
d1=d_test_glm %>% select(-pheno,-Pval)
d1$type=1

d2=d_all_norm[complete.cases(d_all_norm),] # random SNPs
d2=d2[sample(1:nrow(d2),100000,replace=T),]
d2$type=0

dx=rbind(d1,d2)

#joint model
mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+TF+connectedness+connect_decile+Roadmap_count+Roadmap_length_per_type+pLI+promoter_count+PPI_degree_cat+PPI_degree_decile,data=dx,family = "binomial")

#extract coeffs

base_count=83

TF_coef=unname(summary(mod1)$coef[,1][base_count+1])
connectedness_coef=unname(summary(mod1)$coef[,1][base_count+2])*sd(d_all_norm$connectedness,na.rm=T)
connect_rank_coef=unname(summary(mod1)$coef[,1][base_count+3])
Road_coef1=unname(summary(mod1)$coef[,1][base_count+4])
Road_coef2=unname(summary(mod1)$coef[,1][base_count+5])
pLI_coef=unname(summary(mod1)$coef[,1][base_count+6])
promoter_count_coef=unname(summary(mod1)$coef[,1][base_count+7])
PPI_cat_coef=unname(summary(mod1)$coef[,1][base_count+8])*sd(d_all_norm$PPI_degree_cat,na.rm=T)
PPI_rank_coef=unname(summary(mod1)$coef[,1][base_count+9])

#################
##marginal models

mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+TF,data=dx,family = "binomial")
TF_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])

mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+connect_decile,data=dx[dx$connectedness==1,],family = "binomial")
connect_rank_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])

mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+pLI,data=dx,family = "binomial")
pLI_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])

mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+promoter_count,data=dx,family = "binomial")
promoter_count_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])

mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+PPI_degree_decile,data=dx[dx$PPI_degree_cat==1,],family = "binomial")
PPI_rank_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])

mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+ABC_length_per_type+ABC_count,data=dx,family = "binomial")
ABC_length_coef_count=unname(summary(mod1)$coef[,1][base_count+1])
ABC_count_coef_length=unname(summary(mod1)$coef[,1][base_count+2])

mod1=glm(type~d_tss+MAF+L2+TSSD+factor(TSSD_cat)+factor(d_tss_cat)+factor(MAF_cat)+factor(L2_cat)+length+CDS_length+Roadmap_length_per_type+Roadmap_count,data=dx,family = "binomial")
Road_length_coef_count=unname(summary(mod1)$coef[,1][base_count+1])
Road_count_coef_length=unname(summary(mod1)$coef[,1][base_count+2])


#

d_out2=data.frame(TF_coef,connectedness_coef,connect_rank_coef,Road_coef1,Road_coef2,pLI_coef,promoter_count_coef,PPI_cat_coef,PPI_rank_coef)
d_out3=data.frame(TF_coef_marginal,connect_rank_coef_marginal,pLI_coef_marginal,promoter_count_coef_marginal,PPI_rank_coef_marginal,ABC_length_coef_count,ABC_count_coef_length,Road_length_coef_count,Road_count_coef_length)

######################
######################
######################

d_out=cbind(d_out1,d_out2,d_out3)
write.table(d_out,file=outfile,quote=F,sep=",",row.names=F)

