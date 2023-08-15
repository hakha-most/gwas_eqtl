library(data.table)
library(tidyverse)
library(dplyr)

##########
#code to generate data for Supplementary Fig. 4

set.seed(0)

multifile="gene_annotations/12_GO_annotations/genes_multiGO.txt"
genefile="gene_annotations/all_annots_pc_genes.txt"

d_multi=fread(multifile) %>% select(gene,GO_BP_count_400_auc)
colnames(d_multi)[1]="hgnc_id"
d_gene_annots=fread(genefile)
d_gene_annots=left_join(d_gene_annots,d_multi)


#######
#normalize features for logistic regression model

normalize_na_rm<-function(v){
  m=mean(v,na.rm=T)
  s=sd(v,na.rm=T)
  v=(v-m)/s
}

d_gene_annots_norm=d_gene_annots
d_gene_annots_norm$connect_decile=d_gene_annots_norm$connect_decile*1/1
d_gene_annots_norm$PPI_degree_decile=d_gene_annots_norm$PPI_degree_decile*1/1

d_gene_annots_norm$GO_BP_count_400_auc=normalize_na_rm(d_gene_annots_norm$GO_BP_count_400_auc)
d_gene_annots_norm$TF=normalize_na_rm(d_gene_annots_norm$TF)
d_gene_annots_norm$pLI=normalize_na_rm(d_gene_annots_norm$pLI)
d_gene_annots_norm$Roadmap_length_per_type=normalize_na_rm(d_gene_annots_norm$Roadmap_length_per_type)
d_gene_annots_norm$Roadmap_count=normalize_na_rm(d_gene_annots_norm$Roadmap_count)
d_gene_annots_norm$ABC_count=normalize_na_rm(d_gene_annots_norm$ABC_count)
d_gene_annots_norm$ABC_length_per_type=normalize_na_rm(d_gene_annots_norm$ABC_length_per_type)
d_gene_annots_norm$promoter_count=normalize_na_rm(d_gene_annots_norm$promoter_count)
d_gene_annots_norm[d_gene_annots_norm$connectedness==1,]$connect_decile=normalize_na_rm(d_gene_annots_norm[d_gene_annots_norm$connectedness==1,]$connect_decile)
d_gene_annots_norm[d_gene_annots_norm$PPI_degree_cat==1,]$PPI_degree_decile=normalize_na_rm(d_gene_annots_norm[d_gene_annots_norm$PPI_degree_cat==1,]$PPI_degree_decile)


#######
#function to run logistic regression models
calc_reg <- function(dx){
  
  dx=dx %>% select(type,length,CDS_length,TSSD,pLI,ABC_count,ABC_length_per_type,Roadmap_count,Roadmap_length_per_type,promoter_count,connect_decile,connectedness,PPI_degree_decile,PPI_degree_cat,TF,GO_BP_count_400_auc)
  dx=dx[complete.cases(dx),]
  
  base_count=4
  
  mod1=glm(type~TSSD+length+CDS_length+TF,data=dx,family = "binomial")
  TF_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
  TF_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])
  
  mod1=glm(type~TSSD+length+CDS_length+connect_decile,data=dx[dx$connectedness==1,],family = "binomial")
  connect_rank_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
  connect_rank_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])
  
  mod1=glm(type~TSSD+length+CDS_length+pLI,data=dx,family = "binomial")
  pLI_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
  pLI_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])
  
  mod1=glm(type~TSSD+length+CDS_length+promoter_count,data=dx,family = "binomial")
  promoter_count_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
  promoter_count_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])
  
  mod1=glm(type~TSSD+length+CDS_length+PPI_degree_decile,data=dx[dx$PPI_degree_cat==1,],family = "binomial")
  PPI_rank_coef_marginal=unname(summary(mod1)$coef[,1][base_count+1])
  PPI_rank_coef_marginal_SE=unname(summary(mod1)$coef[,2][base_count+1])
  
  mod1=glm(type~TSSD+length+CDS_length+ABC_length_per_type+ABC_count,data=dx,family = "binomial")
  ABC_length_coef_count=unname(summary(mod1)$coef[,1][base_count+1])
  ABC_count_coef_length=unname(summary(mod1)$coef[,1][base_count+2])
  ABC_length_coef_count_SE=unname(summary(mod1)$coef[,2][base_count+1])
  ABC_count_coef_length_SE=unname(summary(mod1)$coef[,2][base_count+2])
  
  mod1=glm(type~TSSD+length+CDS_length+Roadmap_length_per_type+Roadmap_count,data=dx,family = "binomial")
  Road_length_coef_count=unname(summary(mod1)$coef[,1][base_count+1])
  Road_count_coef_length=unname(summary(mod1)$coef[,1][base_count+2])
  Road_length_coef_count_SE=unname(summary(mod1)$coef[,2][base_count+1])
  Road_count_coef_length_SE=unname(summary(mod1)$coef[,2][base_count+2])
  
  mod1=glm(type~TSSD+length+CDS_length+GO_BP_count_400_auc,data=dx,family = "binomial")
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

#######

d_data1=fread("extended_gwas_analyses/PoPS/UKB_AllMethods_GenePrioritization.txt") #data from Weeks et al. 
d_data2=fread("extended_gwas_analyses/PoPS/PASS_AllMethods_GenePrioritization.txt")  

d_data1$trait=paste0(d_data1$trait,"_UKB")
d_data2$trait=paste0(d_data2$trait,"_PASS")

d_data=rbind((d_data1 %>% select(trait,lead_variant,ensgid,pops_rank)), (d_data2 %>% select(trait,lead_variant,ensgid,pops_rank)))
d_data=d_data[d_data$ensgid %in% d_gene$GeneSymbol,] #filter protein-coding genes

#######

traits=unique(d_data$trait)

#loop over traits
kk=1
trait_temp=traits[kk]
print(paste0(as.character(kk),":",trait_temp))
d_data_temp=d_data[d_data$trait==trait_temp,]
d_data_temp=d_data_temp %>% arrange(pops_rank)
d_data_temp=d_data_temp[!duplicated(d_data_temp$lead_variant),] #get top pops_rank gene per lead variant
d2=d_gene_annots_norm[!(d_gene_annots_norm$GeneSymbol %in% d_data_temp$ensgid),]; d2$type=0
d1=d_gene_annots_norm[(d_gene_annots_norm$GeneSymbol %in% d_data_temp$ensgid),]; d1$type=1
dx=rbind(d2,d1)
d_out_temp=calc_reg(dx)
d_out_temp$trait=trait_temp
d_out_temp$count=nrow(d1)
d_out=d_out_temp
d_data_cumulative=d_data_temp

for (kk in 2:length(traits)){
  
  trait_temp=traits[kk]
  print(paste0(as.character(kk),":",trait_temp))
  d_data_temp=d_data[d_data$trait==trait_temp,]
  d_data_temp=d_data_temp %>% arrange(pops_rank) #get top pops_rank gene per lead variant
  d_data_temp=d_data_temp[!duplicated(d_data_temp$lead_variant),] 
  
  #regression
  d2=d_gene_annots_norm[!(d_gene_annots_norm$GeneSymbol %in% d_data_temp$ensgid),]; d2$type=0
  d1=d_gene_annots_norm[(d_gene_annots_norm$GeneSymbol %in% d_data_temp$ensgid),]; d1$type=1
  dx=rbind(d2,d1)
  d_out_temp=calc_reg(dx)
  d_out_temp$trait=trait_temp
  d_out_temp$count=nrow(d1)
  d_out=rbind(d_out,d_out_temp)
  d_data_cumulative=rbind(d_data_cumulative,d_data_temp)
  
}

#print stats
length(unique(d_data_cumulative$trait))
nrow(d_data_cumulative)

#save

outfile="extended_gwas_analyses/PoPS/PoPS.regression_results.txt"
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)

outfile="extended_gwas_analyses/PoPS/PoPS.SNP_gene_links.txt"
write.table(d_data_cumulative,file=outfile,quote=F,sep="\t",row.names=F)


