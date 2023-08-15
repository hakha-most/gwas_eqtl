library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(sfsmisc)
library(cowplot)
library(MASS)

theme_set(theme_cowplot(font_size = 9))

##########

outdir="model_simulations/"

################

K=2.986 #selection parameter
rho=0 #correlation between beta and gamma in the model
cor_str="no_correlation"

N=10000000
a_distr=rnorm(N,0,1)
b_distr=rnorm(N,0,1)

########
#rearrange to obtain a given correlation between effects

a_index=sort(a_distr^2, index.return=TRUE)$ix
b_index=sort(b_distr^2, index.return=TRUE)$ix

a_distr <- a_distr[a_index]
b_distr <- b_distr[b_index]

mus=c(mean(a_distr^2),mean(b_distr^2))
s1=sd(a_distr^2); s2=sd(b_distr^2); sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2)  

mvdat <- mvrnorm(n = length(a_distr), mu = mus, Sigma = sigma, empirical = TRUE)
rx <- rank(mvdat[ , 1], ties.method = "first")
ry <- rank(mvdat[ , 2], ties.method = "first")

a_distr <- a_distr[rx]
b_distr <- b_distr[ry]

########
#variance terms

v_plat_distr=K*(1-exp(-1*a_distr^2*b_distr^2/K))
v_neutral_distr=a_distr^2*b_distr^2

v_p_neutral_distr=v_neutral_distr/(a_distr^2*b_distr^2)
v_p_plat_distr=v_plat_distr/v_neutral_distr

v_e_neutral_distr=v_p_neutral_distr*b_distr^2
v_e_plat_distr=v_p_plat_distr*b_distr^2

########
#trends under selection

gwas_gamma2=1:10
gwas_beta2=1:10
eqtl_gamma2=1:10
eqtl_beta2=1:10
no_coloc_gamma2=1:10
no_coloc_beta2=1:10
coloc_gamma2=1:10
coloc_beta2=1:10
coloc_rate=1:10

for (kk in 1:10){
  print(kk)
  
  ########coloc#########
  
  jj=850+1 #power percentile
  x=quantile(v_plat_distr,probs = seq(0,1,0.001))
  c_gwas=x[jj]
  
  x=quantile(v_e_plat_distr,probs = seq(0,1,0.1))
  c_eqtl=x[11-kk]
  
  gwas_hits=which(v_plat_distr>c_gwas)
  eqtl_hits=which(v_e_plat_distr>c_eqtl)
  
  coloc_hits=intersect(gwas_hits,eqtl_hits)
  no_coloc_hits=gwas_hits[!(gwas_hits %in% coloc_hits)]
  coloc_rate[kk]=length(coloc_hits)/length(gwas_hits)
  
  coloc_gamma2[kk]=mean(a_distr[coloc_hits]^2)
  coloc_beta2[kk]=mean(b_distr[coloc_hits]^2)
  no_coloc_gamma2[kk]=mean(a_distr[no_coloc_hits]^2)
  no_coloc_beta2[kk]=mean(b_distr[no_coloc_hits]^2)
  
  ########effect########
  
  x=quantile(v_e_plat_distr,probs = seq(0,1,0.1))
  c_eqtl=x[11-kk]
  
  x=quantile(v_plat_distr,probs = seq(0,1,0.1))
  c_gwas=x[11-kk]
  
  gwas_hits=which(v_plat_distr>c_gwas)
  eqtl_hits=which(v_e_plat_distr>c_eqtl)
  
  gwas_gamma2[kk]=mean(a_distr[gwas_hits]^2)
  gwas_beta2[kk]=mean(b_distr[gwas_hits]^2)
  
  eqtl_gamma2[kk]=mean(a_distr[eqtl_hits]^2)
  eqtl_beta2[kk]=mean(b_distr[eqtl_hits]^2)
  
}

d_data=data.frame(power=seq(0.1,1,0.1),coloc_rate,coloc_gamma2,coloc_beta2,no_coloc_beta2,no_coloc_gamma2,gwas_gamma2,gwas_beta2,eqtl_gamma2,eqtl_beta2)

#########

lab_x="eQTL study power"
lab_y="Mean squared effect\nof co-discovered SNPs"

d1=d_data %>% dplyr::select(power,coloc_gamma2); colnames(d1)=c("power","effect2"); d1$cat="gamma"
d2=d_data %>% dplyr::select(power,coloc_beta2); colnames(d2)=c("power","effect2"); d2$cat="beta"
d_coloc=rbind(d1,d2)

p1<-ggplot(d_coloc,aes(x=power, y = effect2, shape=cat)) +
  geom_point(size=1.5)  + 
  geom_line(size=0.5)  + 
  scale_shape_manual(name = "",labels = c(expression(beta^2), expression(gamma^2)),values = c(16, 17))+
  geom_segment(aes(x=0.1,xend=1,y=effect2[10],yend=effect2[10]) , linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0,5.1))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_codisc_effect_by_power_",cor_str,"_cor.pdf")
save_plot(outfile, p1,base_width = 2.5, base_height = 2.3)


#########

col_scale=c("GWAS"="#1E56A0","eQTL"="#FF5959")

d1=d_data %>% dplyr::select(power,gwas_gamma2); colnames(d1)=c("power","effect2"); d1$cat="gamma"
d2=d_data %>% dplyr::select(power,gwas_beta2); colnames(d2)=c("power","effect2"); d2$cat="beta"
d_gwas=rbind(d1,d2)
d_gwas$study="GWAS"

d1=d_data %>% dplyr::select(power,eqtl_gamma2); colnames(d1)=c("power","effect2"); d1$cat="gamma"
d2=d_data %>% dplyr::select(power,eqtl_beta2); colnames(d2)=c("power","effect2"); d2$cat="beta"
d_eqtl=rbind(d1,d2)
d_eqtl$study="eQTL"

d_plot=rbind(d_gwas,d_eqtl)

lab_x="Power"
lab_y="Mean squared effect"

p1<-ggplot(d_plot,aes(x=power, y = effect2, shape=cat,col=study)) +
  geom_point(size=1.5)  +
  geom_line(size=0.5)  +
  scale_color_manual(name="",values=col_scale)+
  scale_shape_manual(name = "",labels = c(expression(beta^2), expression(gamma^2)),values = c(16, 17))+
  geom_hline(yintercept=1, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0,4.7))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_effect_by_power_",cor_str,"_cor.pdf")
save_plot(outfile, p1,base_width = 2.5, base_height = 2.5)

#

lab_x="eQTL study power"
lab_y="Colocalization rate"

p1<-ggplot(d_data,aes(x=power, y = coloc_rate)) +
  geom_point(size=1)  + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_coloc_rate_by_power_",cor_str,"_cor.pdf")
save_plot(outfile, p1,base_width = 2.5, base_height = 2.5)



########
#trends under neutrality

gwas_gamma2=1:10
gwas_beta2=1:10
eqtl_gamma2=1:10
eqtl_beta2=1:10
coloc_gamma2=1:10
coloc_beta2=1:10
coloc_rate=1:10

for (kk in 1:10){
  print(kk)
  
  ########coloc#########
  
  jj=850+1 #power percentile
  x=quantile(v_neutral_distr,probs = seq(0,1,0.001))
  c_gwas=x[jj]
  
  x=quantile(v_e_neutral_distr,probs = seq(0,1,0.1))
  c_eqtl=x[11-kk]
  
  gwas_hits=which(v_neutral_distr>c_gwas)
  eqtl_hits=which(v_e_neutral_distr>c_eqtl)
  
  coloc_hits=intersect(gwas_hits,eqtl_hits)
  coloc_rate[kk]=length(coloc_hits)/length(gwas_hits)
  
  coloc_gamma2[kk]=mean(a_distr[coloc_hits]^2)
  coloc_beta2[kk]=mean(b_distr[coloc_hits]^2)
  
  ########effect########
  
  x=quantile(v_e_neutral_distr,probs = seq(0,1,0.1))
  c_eqtl=x[11-kk]
  
  x=quantile(v_neutral_distr,probs = seq(0,1,0.1))
  c_gwas=x[11-kk]
  
  gwas_hits=which(v_neutral_distr>c_gwas)
  eqtl_hits=which(v_e_neutral_distr>c_eqtl)
  
  ###
  
  gwas_gamma2[kk]=mean(a_distr[gwas_hits]^2)
  gwas_beta2[kk]=mean(b_distr[gwas_hits]^2)
  
  eqtl_gamma2[kk]=mean(a_distr[eqtl_hits]^2)
  eqtl_beta2[kk]=mean(b_distr[eqtl_hits]^2)
  
}

d_data=data.frame(power=seq(0.1,1,0.1),coloc_rate,coloc_gamma2,coloc_beta2,gwas_gamma2,gwas_beta2,eqtl_gamma2,eqtl_beta2)

#########

lab_x="eQTL study power"
lab_y="Mean squared effect\nof colocalized SNPs"

d1=d_data %>% dplyr::select(power,coloc_gamma2); colnames(d1)=c("power","effect2"); d1$cat="gamma"
d2=d_data %>% dplyr::select(power,coloc_beta2); colnames(d2)=c("power","effect2"); d2$cat="beta"
d_coloc=rbind(d1,d2)

p1<-ggplot(d_coloc,aes(x=power, y = effect2, shape=cat)) +
  geom_point(size=1.5)  + 
  geom_line(size=0.5)  + 
  scale_shape_manual(name = "",labels = c(expression(beta^2), expression(gamma^2)),values = c(16, 17))+
  geom_hline(yintercept=d_coloc$effect2[10], linetype="dashed", size=0.5)+
  geom_hline(yintercept=d_coloc$effect2[20], linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  xlim(c(0,1))+
  ylim(c(0,5.1))+
  theme_cowplot(9)


outfile=paste0(outdir,"neutral_coloc_effect_by_power_",cor_str,"_cor.pdf")
save_plot(outfile, p1,base_width = 2.5, base_height = 2.5)

#########

col_scale=c("GWAS"="#1E56A0","eQTL"="#FF5959")

d1=d_data %>% dplyr::select(power,gwas_gamma2); colnames(d1)=c("power","effect2"); d1$cat="gamma"
d2=d_data %>% dplyr::select(power,gwas_beta2); colnames(d2)=c("power","effect2"); d2$cat="beta"
d_gwas=rbind(d1,d2)
d_gwas$study="GWAS"

d1=d_data %>% dplyr::select(power,eqtl_gamma2); colnames(d1)=c("power","effect2"); d1$cat="gamma"
d2=d_data %>% dplyr::select(power,eqtl_beta2); colnames(d2)=c("power","effect2"); d2$cat="beta"
d_eqtl=rbind(d1,d2)
d_eqtl$study="eQTL"

d_plot=rbind(d_gwas,d_eqtl)

lab_x="Power"
lab_y="Mean squared effect"

p1<-ggplot(d_plot,aes(x=power, y = effect2, shape=cat,col=study)) +
  geom_point(size=1.5)  +
  geom_line(size=0.5)  + 
  scale_color_manual(name="",values=col_scale)+
  scale_shape_manual(name = "",labels = c(expression(beta^2), expression(gamma^2)),values = c(16, 17))+
  geom_hline(yintercept=1, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0,4.7))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"neutral_effect_by_power_",cor_str,"_cor.pdf")
save_plot(outfile, p1,base_width = 2.5, base_height = 2.5)

#

lab_x="eQTL study power"
lab_y="Colocalization rate"

p1<-ggplot(d_data,aes(x=power, y = coloc_rate)) +
  geom_point(size=1)  + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed", size=0.5)+
  xlab(lab_x)+
  ylab(lab_y)+
  ylim(c(0,1))+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"neutral_coloc_rate_by_power_",cor_str,"_cor.pdf")
save_plot(outfile, p1,base_width = 2.5, base_height = 2.5)

