library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(sfsmisc)
library(cowplot)
library(RColorBrewer)
library(MASS)

theme_set(theme_cowplot(font_size = 9))

##########

outdir="model_simulations/"

##########

N=10000000 #number of variants
w_share=0.5 #fraction of shared variants across two traits
N_gene=20000 #number of genes

compute_stats <- function(N,w_share,w1,N_gene,selection){
  
  w2=1-w1
  b_distr=rnorm(N,0,1)
  
  #########
  #trait-specific gene effects
  a1_not_share=rnorm(N_gene*(1-w_share),0,1)
  a2_not_share=rnorm(N_gene*(1-w_share),0,1)
  
  #shared gene effects
  a1=rnorm(N_gene*w_share,0,1)
  a2=rnorm(N_gene*w_share,0,1)
  
  a1_index=sort(a1, index.return=TRUE)$ix
  a2_index=sort(a2, index.return=TRUE)$ix
  
  a1 <- a1[a1_index]
  a2 <- a2[a2_index]
  
  #rearrange to induce correlation of effect for shared SNPs
  rho=0.75
  
  mus=c(mean(a1),mean(a2))
  s1=sd(a1); s2=sd(a2); sigma <- matrix(c(s1^2, s1*s2*rho,s2*s1*rho,s2^2),2)  
  
  mvdat <- mvrnorm(n = length(a1), mu = mus, Sigma = sigma, empirical = TRUE)
  rx <- rank(mvdat[ , 1], ties.method = "first")
  ry <- rank(mvdat[ , 2], ties.method = "first")
  
  a1_share <- a1[rx]
  a2_share <- a2[ry]
  
  #########
  
  a1_gene_distr=c(a1_share,a1_not_share)
  a2_gene_distr=c(a2_share,a2_not_share)
  d_eff=data.frame(effect1=a1_gene_distr,effect2=a2_gene_distr,gene=1:N_gene)
  
  d_genes=data.frame(gene=sample(1:N_gene,size=N,replace = T))
  d_genes=left_join(d_genes,d_eff)
  
  a1_distr=d_genes$effect1
  a2_distr=d_genes$effect2
  a1_distr=a1_distr/sd(a1_distr)
  a2_distr=a2_distr/sd(a2_distr)
  genes=d_genes$gene
  
  #compute aggregate effects as defined in the text
  net_effect= w1*(a1_distr*b_distr)^2 + w2*(a2_distr*b_distr)^2 
  
  #find selection parameter, K, corresponding to 10% reduction in phenotypic variance
  matched=0
  K=2.8
  while (matched==0){
    K=K+0.01
    v_plat_distr=K*(1-exp(-1*net_effect/K))
    v_neutral_distr=net_effect
    
    d=mean(v_plat_distr/v_neutral_distr)
    if (abs(d-0.9)<0.001){matched=1}
  }
  
  #variance terms
  v_p_neutral_distr=v_neutral_distr/net_effect
  v_p_plat_distr=v_plat_distr/v_neutral_distr
  
  v_e_neutral_distr=v_p_neutral_distr*b_distr^2
  v1_neutral_distr=v_p_neutral_distr * a1_distr^2*b_distr^2
  v2_neutral_distr=v_p_neutral_distr * a2_distr^2*b_distr^2
  
  v_e_plat_distr=v_p_plat_distr*b_distr^2
  v1_plat_distr=v_p_plat_distr * a1_distr^2*b_distr^2
  v2_plat_distr=v_p_plat_distr * a2_distr^2*b_distr^2
  
  ########
  
  if (selection==1){v_e_distr=v_e_plat_distr; v1_distr=v1_plat_distr; v2_distr=v2_plat_distr}
  if (selection==0){v_e_distr=v_e_neutral_distr; v1_distr=v1_neutral_distr; v2_distr=v2_neutral_distr}
  
  ###
  
  d_genes=data.frame(gene=genes,effect1=a1_distr,effect2=a2_distr,b_distr)
  d_genes$net_effect=w1*(d_genes$effect1*b_distr)^2 + w2*(d_genes$effect2*b_distr)^2
  
  ##########
  
  gwas1_gamma1=1:10
  gwas1_gamma2=1:10
  gwas1_beta2=1:10
  gwas2_gamma1=1:10
  gwas2_gamma2=1:10
  gwas2_beta2=1:10
  gwas_gamma1=1:10
  gwas_gamma2=1:10
  gwas_beta2=1:10
  eqtl_gamma1=1:10
  eqtl_gamma2=1:10
  eqtl_beta2=1:10
  
  gwas1_frac_shared=1:10
  gwas2_frac_shared=1:10
  gwas_frac_shared=1:10
  eqtl_frac_shared=1:10
  
  for (kk in 1:10){
    print(kk)
    
    ########effect########
    
    x=quantile(v_e_distr,probs = seq(0,1,0.1))
    c_eqtl=x[11-kk]
    
    x=quantile(v1_distr,probs = seq(0,1,0.1))
    c1_gwas=x[11-kk]
    
    x=quantile(v2_distr,probs = seq(0,1,0.1))
    c2_gwas=x[11-kk]
    
    gwas1_hits=which(v1_distr>c1_gwas)
    gwas2_hits=which(v2_distr>c2_gwas)
    gwas_hits=unique(c(gwas1_hits,gwas2_hits))
    eqtl_hits=which(v_e_distr>c_eqtl)
    
    gwas1_gamma1[kk]=mean(a1_distr[gwas1_hits]^2)
    gwas1_gamma2[kk]=mean(a2_distr[gwas1_hits]^2)
    gwas1_beta2[kk]=mean(b_distr[gwas1_hits]^2)
    
    gwas2_gamma1[kk]=mean(a1_distr[gwas2_hits]^2)
    gwas2_gamma2[kk]=mean(a2_distr[gwas2_hits]^2)
    gwas2_beta2[kk]=mean(b_distr[gwas2_hits]^2)
    
    gwas_gamma1[kk]=mean(a1_distr[gwas_hits]^2)
    gwas_gamma2[kk]=mean(a2_distr[gwas_hits]^2)
    gwas_beta2[kk]=mean(b_distr[gwas_hits]^2)
    
    eqtl_gamma1[kk]=mean(a1_distr[eqtl_hits]^2)
    eqtl_gamma2[kk]=mean(a2_distr[eqtl_hits]^2)
    eqtl_beta2[kk]=mean(b_distr[eqtl_hits]^2)
    
    ##pleiotropy

    d_genes_temp=d_genes[gwas1_hits,]
    gwas1_frac_shared[kk]=nrow(d_genes_temp[d_genes_temp$gene<=(N_gene*w_share),])/nrow(d_genes_temp)
    
    d_genes_temp=d_genes[gwas2_hits,]
    gwas2_frac_shared[kk]=nrow(d_genes_temp[d_genes_temp$gene<=(N_gene*w_share),])/nrow(d_genes_temp)
    
    d_genes_temp=d_genes[gwas_hits,]
    gwas_frac_shared[kk]=nrow(d_genes_temp[d_genes_temp$gene<=(N_gene*w_share),])/nrow(d_genes_temp)
    
    d_genes_temp=d_genes[eqtl_hits,]
    eqtl_frac_shared[kk]=nrow(d_genes_temp[d_genes_temp$gene<=(N_gene*w_share),])/nrow(d_genes_temp)
    
    #
    
  }
  
  dx1=data.frame(gwas1_gamma1,gwas1_gamma2,gwas1_beta2)
  dx2=data.frame(gwas2_gamma1,gwas2_gamma2,gwas2_beta2)
  dx3=data.frame(gwas_gamma1,gwas_gamma2,gwas_beta2)
  dx4=data.frame(eqtl_gamma1,eqtl_gamma2,eqtl_beta2)
  dy=cbind(dx1,dx2,dx3,dx4); dy$power=seq(0.1,1,0.1)
  
  dx1=data.frame(gwas1_frac_shared)
  dx2=data.frame(gwas2_frac_shared)
  dx3=data.frame(gwas_frac_shared)
  dx4=data.frame(eqtl_frac_shared)
  dz=cbind(dx1,dx2,dx3,dx4); dz$power=seq(0.1,1,0.1)
  
  d_out=list(dy,dz,K)
  return(d_out)
  
}

#different selection scenarios
w1=0.5; selection=1
results_selection=compute_stats(N,w_share,w1,N_gene,selection)

w1=1; selection=1
results_selection_1trait=compute_stats(N,w_share,w1,N_gene,selection)

w1=0.5; w2=1-w1; selection=0
results_neutral=compute_stats(N,w_share,w1,N_gene,selection)


######

brewer.pal(n = 5, name = "Set1")
display.brewer.pal(n = 5, name = 'Set1')
col_scale=c("Neutral"="#FF7F00","Selection (both traits)"="#4DAF4A","Selection (one trait)"="#984EA3","Selection (Trait 1)"="#E41A1C","Selection (Trait 2)"="#377EB8")

######

lab_x="GWAS study power (Trait 1)"
lab_y="Mean squared gene effect\n(Trait 1)"

d1=results_selection[[1]] %>% dplyr::select(power,gwas1_gamma1); colnames(d1)=c("power","effect2"); d1$Model="Selection (both traits)"
d2=results_selection_1trait[[1]]  %>% dplyr::select(power,gwas1_gamma1); colnames(d2)=c("power","effect2"); d2$Model="Selection (Trait 1)"
d2p=results_selection_1trait[[1]]  %>% dplyr::select(power,gwas2_gamma2); colnames(d2p)=c("power","effect2"); d2p$Model="Selection (Trait 2)"
d3=results_neutral[[1]] %>% dplyr::select(power,gwas1_gamma1); colnames(d3)=c("power","effect2"); d3$Model="Neutral"
d_data=rbind(d1,d2,d2p,d3)

p1<-ggplot(d_data,aes(x=power, y = effect2, col=Model)) +
  geom_point(size=1.5)  + 
  geom_line(size=0.5)+
  scale_color_manual(values = col_scale)+
  geom_hline(yintercept=1, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0,3.2))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_multi_pheno_GWAS_gamma2_by_power.pdf")
save_plot(outfile, p1,base_width = 3.5, base_height = 2)

#

######

lab_x="eQTL study power"
lab_y="Mean squared gene effect\n(Trait 1)"

d1=results_selection[[1]] %>% dplyr::select(power,eqtl_gamma1); colnames(d1)=c("power","effect2"); d1$Model="Selection (both traits)"
d2=results_selection_1trait[[1]]  %>% dplyr::select(power,eqtl_gamma1); colnames(d2)=c("power","effect2"); d2$Model="Selection (Trait 1)"
d2p=results_selection_1trait[[1]]  %>% dplyr::select(power,eqtl_gamma2); colnames(d2p)=c("power","effect2"); d2p$Model="Selection (Trait 2)"
d3=results_neutral[[1]] %>% dplyr::select(power,eqtl_gamma1); colnames(d3)=c("power","effect2"); d3$Model="Neutral"
d_data=rbind(d1,d2,d2p,d3)

p1<-ggplot(d_data,aes(x=power, y = effect2, col=Model)) +
  geom_point(size=1.5)  + 
  geom_line(size=0.5)+
  scale_color_manual(values = col_scale)+
  geom_hline(yintercept=1, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0,3.2))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_multi_pheno_eQTL_gamma2_by_power.pdf")
save_plot(outfile, p1,base_width = 3.5, base_height = 2)

#


######

lab_x="GWAS study power (Trait 1)"
lab_y="Mean squared SNP effect\non expression"

d1=results_selection[[1]] %>% dplyr::select(power,gwas1_beta2); colnames(d1)=c("power","beta2"); d1$Model="Selection (both traits)"
d2=results_selection_1trait[[1]]  %>% dplyr::select(power,gwas1_beta2); colnames(d2)=c("power","beta2"); d2$Model="Selection (Trait 1)"
d2p=results_selection_1trait[[1]]  %>% dplyr::select(power,gwas2_beta2); colnames(d2p)=c("power","beta2"); d2p$Model="Selection (Trait 2)"
d3=results_neutral[[1]] %>% dplyr::select(power,gwas1_beta2); colnames(d3)=c("power","beta2"); d3$Model="Neutral"
d_data=rbind(d1,d2,d2p,d3)

p1<-ggplot(d_data,aes(x=power, y = beta2, col=Model)) +
  geom_point(size=1.5)  + 
  geom_line(size=0.5)+
  scale_color_manual(values = col_scale)+
  geom_hline(yintercept=1, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0.5,4.5))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_multi_pheno_GWAS_beta2_by_power.pdf")
save_plot(outfile, p1,base_width = 3.5, base_height = 2)

#

######

lab_x="eQTL study power"
lab_y="Mean squared SNP effect\non expression"

d1=results_selection[[1]] %>% dplyr::select(power,eqtl_beta2); colnames(d1)=c("power","beta2"); d1$Model="Selection (both traits)"
d2=results_selection_1trait[[1]]  %>% dplyr::select(power,eqtl_beta2); colnames(d2)=c("power","beta2"); d2$Model="Selection (one trait)"
d3=results_neutral[[1]] %>% dplyr::select(power,eqtl_beta2); colnames(d3)=c("power","beta2"); d3$Model="Neutral"
d_data=rbind(d1,d2,d3)

p1<-ggplot(d_data,aes(x=power, y = beta2, col=Model)) +
  geom_point(size=1.5)  + 
  geom_line(size=0.5)+
  scale_color_manual(values = col_scale)+
  geom_hline(yintercept=1, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0.5,4.5))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_multi_pheno_eQTL_beta2_by_power.pdf")
save_plot(outfile, p1,base_width = 3.5, base_height = 2)


######

lab_x="GWAS study power (Trait 1)"
lab_y="Fraction of genes\nwith correlated effects"

d1=results_selection[[2]] %>% dplyr::select(power,gwas1_frac_shared); colnames(d1)=c("power","fract"); d1$Model="Selection (both traits)"
d2=results_selection_1trait[[2]]  %>% dplyr::select(power,gwas1_frac_shared); colnames(d2)=c("power","fract"); d2$Model="Selection (Trait 1)"
d2p=results_selection_1trait[[2]]  %>% dplyr::select(power,gwas2_frac_shared); colnames(d2p)=c("power","fract"); d2p$Model="Selection (Trait 2)"
d3=results_neutral[[2]] %>% dplyr::select(power,gwas1_frac_shared); colnames(d3)=c("power","fract"); d3$Model="Neutral"
d_data=rbind(d1,d2,d2p,d3)

p1<-ggplot(d_data,aes(x=power, y = fract, col=Model)) +
  geom_point(size=1.5)  + 
  geom_line(size=0.5)+
  scale_color_manual(values = col_scale)+
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0.45,0.52))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_multi_pheno_GWAS_fract_shared_by_power.pdf")
save_plot(outfile, p1,base_width = 3.5, base_height = 2)

#
lab_x="eQTL study power"
lab_y="Fraction of genes\nwith correlated effects"

d1=results_selection[[2]] %>% dplyr::select(power,eqtl_frac_shared); colnames(d1)=c("power","fract"); d1$Model="Selection (both traits)"
d2=results_selection_1trait[[2]]  %>% dplyr::select(power,eqtl_frac_shared); colnames(d2)=c("power","fract"); d2$Model="Selection (one trait)"
d3=results_neutral[[2]] %>% dplyr::select(power,eqtl_frac_shared); colnames(d3)=c("power","fract"); d3$Model="Neutral"
d_data=rbind(d1,d2,d3)

p1<-ggplot(d_data,aes(x=power, y = fract, col=Model)) +
  geom_point(size=1.5)  + 
  geom_line(size=0.5)+
  scale_color_manual(values = col_scale)+
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0.45,0.52))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_multi_pheno_eqt_fract_shared_by_power.pdf")
save_plot(outfile, p1,base_width = 3.5, base_height = 2)
