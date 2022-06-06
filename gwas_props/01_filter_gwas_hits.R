
## Script to filter and prune the trait list such that no pairs of traits have rho_g > 0.5.

gwasfile="gwas_hits.all_phenos.clumped_sorted.txt"
snpfile="snp_annots/filter_snps.txt"
rg_file="ukb_neale_lab/rg.pairwise_min_0.5.txt"
h2_file="ukb_neale_lab/h2.2Oct2019.tsv"
phenos_file="ukb_neale_lab/phenos_list.txt"
outfile="filter_indep_gwas.assoc"

library(tidyverse)
library(data.table)
set.seed(0)

##########

#read files

d_snp=fread(snpfile)

d_gwas=fread(gwasfile,header=F) 
colnames(d_gwas)=c("pheno","SNP","Pval")
d_gwas=d_gwas[(d_gwas$SNP %in% d_snp$SNP),] #restrict to filtered SNPs
  
d_rg=fread(rg_file)

d_h2=fread(h2_file)

good_phenos=fread(phenos_file,header = F)
colnames(good_phenos)=c("p","p_conv") # p: trait label; p_conv: trait label modified for categorical traits

##########

#filter & get independent traits

#filter by h2

df1=d_h2[d_h2$isBinary=="TRUE",]
df1_top=df1[df1$n_cases/(df1$n_cases+df1$n_controls)>0.05,] #restrict to binary traits with prevalence > 5%
df2=d_h2[d_h2$isBinary=="FALSE",]
df=rbind(df1_top,df2)

df <- df %>% select(phenotype,h2_liability,h2_liability_se)
colnames(df)=c("p","h2","h2_se")

df=df[complete.cases(df),]
df$h2_lower=df$h2-2*df$h2_se
df$h2_upper=df$h2+2*df$h2_se

df_top=df[df$h2_lower>0.05,] #restrict to traits with h2 significantly > 0.05
df_top=df_top[(df_top$p %in% good_phenos$p),] #restrict to "meaningful" traits 

dx=good_phenos
colnames(dx)=c("p","p_conv")
df_top=left_join(df_top,dx,by="p")

d_h2_top=df_top

#prune traits to keep pariwise independent traits, defined as rho_g < 0.5

tx_temp=aggregate((d_gwas$SNP), list(d_gwas$pheno), length)
tx_temp=tx_temp %>% arrange(-x)
colnames(tx_temp)=c("pheno","count")

tx_temp_good=tx_temp[tx_temp$pheno %in% good_phenos$p_conv,]
tx_temp_good=tx_temp_good[tx_temp_good$pheno %in% d_h2_top$p_conv,]

done_prune=0
pruned_tx_temp=tx_temp_good
i=0
while (done_prune==0){
  
  i=i+1
  print(i)
  pheno_temp=pruned_tx_temp$pheno[i]
  d_rg_temp=d_rg[(d_rg$p1==pheno_temp | d_rg$p2==pheno_temp),]
  phenos_to_remove=unique(c(d_rg_temp$p1,d_rg_temp$p2)); phenos_to_remove=phenos_to_remove[phenos_to_remove!=pheno_temp]
  pruned_tx_temp=pruned_tx_temp[!(pruned_tx_temp$pheno %in% phenos_to_remove),]
  if (nrow(pruned_tx_temp)==i){done_prune=1}
  
}

pruned_phenos=pruned_tx_temp$pheno
top_phenos=tx_temp_good[tx_temp_good$count>=50,]$pheno #restrict to traits with > 50 hits

d_gwas_filter_indep=d_gwas[d_gwas$pheno %in% intersect(pruned_phenos,top_phenos),]
d_gwas_filter_indep=d_gwas_filter_indep %>% distinct(`pheno`, `SNP`, .keep_all = TRUE)

write.table(d_gwas_filter_indep,file=outfile,quote=F,sep="\t",row.names=F)
  



