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

################

N=10000000 #number of variants
w_share=0.1 #proportion of shared vs cell type specific regulatory variants
w1=0.05; w2=0.2; w3=1-w1-w2 #proportions of different cell types
N_gene=20000 #number of genes

#gene effects per cell type
a1_gene_distr=rnorm(N_gene,0,4)
a2_gene_distr=rnorm(N_gene,0,2)
a3_gene_distr=rnorm(N_gene,0,1)

#assign gene effects to variants
d_gene_distr=data.frame(a1_gene_distr,a2_gene_distr,a3_gene_distr,gene=1:N_gene)
d_a_distr=data.frame(gene=sample(1:N_gene,size=N,replace = T))
d_a_distr=left_join(d_a_distr,d_gene_distr)

a1_distr=d_a_distr$a1_gene_distr
a2_distr=d_a_distr$a2_gene_distr
a3_distr=d_a_distr$a3_gene_distr

#########
#assign expression effects to variants
b1_not_share=rnorm(N*(1-w_share)/3,0,1);b1_not_share=c(b1_not_share,rep(0,2*N*(1-w_share)/3))
b2_not_share=rnorm(N*(1-w_share)/3,0,1);b2_not_share=c(rep(0,N*(1-w_share)/3),b2_not_share,rep(0,N*(1-w_share)/3))
b3_not_share=rnorm(N*(1-w_share)/3,0,1);b3_not_share=c(rep(0,2*N*(1-w_share)/3),b3_not_share)

b1=rnorm(N*w_share,0,1)
b2=rnorm(N*w_share,0,1)
b3=rnorm(N*w_share,0,1)

b1_index=sort(b1, index.return=TRUE)$ix
b2_index=sort(b2, index.return=TRUE)$ix
b3_index=sort(b3, index.return=TRUE)$ix

b1 <- b1[b1_index]
b2 <- b2[b2_index]
b3 <- b3[b3_index]

#rearrange to induce correlation of effect for shared SNPs
rho=0.75

mus=c(mean(b1),mean(b2),mean(b3))
s1=sd(b1); s2=sd(b2); s3=sd(b3); sigma <- matrix(c(s1^2, s1*s2*rho, s1*s3*rho,s1*s2*rho, s2^2,s2*s3*rho,s1*s3*rho,s2*s3*rho,s3^2),3)  

mvdat <- mvrnorm(n = length(b1), mu = mus, Sigma = sigma, empirical = TRUE)
rx <- rank(mvdat[ , 1], ties.method = "first")
ry <- rank(mvdat[ , 2], ties.method = "first")
rz <- rank(mvdat[ , 3], ties.method = "first")

b1_share <- b1[rx]
b2_share <- b2[ry]
b3_share <- b3[rz]

#combine shared and specific SNPs
b1_distr=c(b1_share,b1_not_share)
b2_distr=c(b2_share,b2_not_share)
b3_distr=c(b3_share,b3_not_share)

########
#compute aggregate effects as defined in the text

b_distr=w1*b1_distr + w2*b2_distr + w3*b3_distr
a_distr=(w1*b1_distr*a1_distr + w2*b2_distr*a2_distr + w3*b3_distr*a3_distr)/(w1*b1_distr + w2*b2_distr + w3*b3_distr)

########
#find selection parameter, K, corresponding to 10% reduction in phenotypic variance

matched=0
K=0.5
while (matched==0){
  K=K+0.01
  v_plat_distr=K*(1-exp(-1*a_distr^2*b_distr^2/K))
  v_neutral_distr=a_distr^2*b_distr^2
  
  d=mean(v_plat_distr/v_neutral_distr)
  if (abs(d-0.9)<0.001){matched=1}
}

#variance terms
v_p_neutral_distr=v_neutral_distr/(a_distr^2*b_distr^2)
v_p_plat_distr=v_plat_distr/v_neutral_distr

v_e_neutral_distr=v_p_neutral_distr*b_distr^2
v_e_plat_distr=v_p_plat_distr*b_distr^2

##########SELECTION

coloc_rate=1:10

coloc_shared_enrich=1:10
coloc_c1_enrich=1:10
coloc_c2_enrich=1:10
coloc_c3_enrich=1:10

eqtl_shared_enrich=1:10
eqtl_c1_enrich=1:10
eqtl_c2_enrich=1:10
eqtl_c3_enrich=1:10

gwas_shared_enrich=1:10
gwas_c1_enrich=1:10
gwas_c2_enrich=1:10
gwas_c3_enrich=1:10

coloc_gamma2_shared=1:10
coloc_gamma2_c1=1:10
coloc_gamma2_c2=1:10
coloc_gamma2_c3=1:10

gwas_gamma2_shared=1:10
gwas_gamma2_c1=1:10
gwas_gamma2_c2=1:10
gwas_gamma2_c3=1:10

eqtl_gamma2_shared=1:10
eqtl_gamma2_c1=1:10
eqtl_gamma2_c2=1:10
eqtl_gamma2_c3=1:10


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
  
  gwas_shared_hits=gwas_hits[which(gwas_hits<(N*w_share))]
  gwas_c1_hits=gwas_hits[which(gwas_hits>(N*w_share) & gwas_hits<(N*w_share + N*(1-w_share)/3))]
  gwas_c2_hits=gwas_hits[which(gwas_hits>(N*w_share + N*(1-w_share)/3) & gwas_hits<(N*w_share + 2*N*(1-w_share)/3))]
  gwas_c3_hits=gwas_hits[which(gwas_hits>(N*w_share + 2*N*(1-w_share)/3) & gwas_hits<(N*w_share + 3*N*(1-w_share)/3))]
  
  coloc_hits=intersect(gwas_hits,eqtl_hits)
  coloc_rate[kk]=length(coloc_hits)/length(gwas_hits)
  
  coloc_shared_hits=coloc_hits[which(coloc_hits<(N*w_share))]
  coloc_c1_hits=coloc_hits[which(coloc_hits>(N*w_share) & coloc_hits<(N*w_share + N*(1-w_share)/3))]
  coloc_c2_hits=coloc_hits[which(coloc_hits>(N*w_share + N*(1-w_share)/3) & coloc_hits<(N*w_share + 2*N*(1-w_share)/3))]
  coloc_c3_hits=coloc_hits[which(coloc_hits>(N*w_share + 2*N*(1-w_share)/3) & coloc_hits<(N*w_share + 3*N*(1-w_share)/3))]
  
  coloc_shared_enrich[kk]=(length(coloc_shared_hits) / length(coloc_hits)) / (length(gwas_shared_hits) / length(gwas_hits))
  coloc_c1_enrich[kk]=(length(coloc_c1_hits) / length(coloc_hits)) / (length(gwas_c1_hits) / length(gwas_hits))
  coloc_c2_enrich[kk]=(length(coloc_c2_hits) / length(coloc_hits)) / (length(gwas_c2_hits) / length(gwas_hits))
  coloc_c3_enrich[kk]=(length(coloc_c3_hits) / length(coloc_hits)) / (length(gwas_c3_hits) / length(gwas_hits))
  
  coloc_gamma2_shared[kk]=mean(rowSums(cbind(w1*a1_distr[coloc_shared_hits],w2*a2_distr[coloc_shared_hits],w3*a3_distr[coloc_shared_hits]))^2) 
  coloc_gamma2_c1[kk]=mean(w1^2*a1_distr[coloc_c1_hits]^2) 
  coloc_gamma2_c2[kk]=mean(w2^2*a2_distr[coloc_c2_hits]^2) 
  coloc_gamma2_c3[kk]=mean(w3^2*a3_distr[coloc_c3_hits]^2) 
  
  ########effect########
  
  x=quantile(v_e_plat_distr,probs = seq(0,1,0.1))
  c_eqtl=x[11-kk]
  
  x=quantile(v_plat_distr,probs = seq(0,1,0.1))
  c_gwas=x[11-kk]
  
  gwas_hits=which(v_plat_distr>c_gwas)
  eqtl_hits=which(v_e_plat_distr>c_eqtl)
  
  ###cell type
  
  eqtl_shared_hits=eqtl_hits[which(eqtl_hits<(N*w_share))]
  eqtl_c1_hits=eqtl_hits[which(eqtl_hits>(N*w_share) & eqtl_hits<(N*w_share + N*(1-w_share)/3))]
  eqtl_c2_hits=eqtl_hits[which(eqtl_hits>(N*w_share + N*(1-w_share)/3) & eqtl_hits<(N*w_share + 2*N*(1-w_share)/3))]
  eqtl_c3_hits=eqtl_hits[which(eqtl_hits>(N*w_share + 2*N*(1-w_share)/3) & eqtl_hits<(N*w_share + 3*N*(1-w_share)/3))]
  
  eqtl_shared_enrich[kk]=length(eqtl_shared_hits) / length(eqtl_hits)/ w_share
  eqtl_c1_enrich[kk]=length(eqtl_c1_hits) / length(eqtl_hits)  / ((1-w_share)/3)
  eqtl_c2_enrich[kk]=length(eqtl_c2_hits) / length(eqtl_hits) / ((1-w_share)/3)
  eqtl_c3_enrich[kk]=length(eqtl_c3_hits) / length(eqtl_hits) / ((1-w_share)/3)
  
  eqtl_gamma2_shared[kk]=mean(rowSums(cbind(w1*a1_distr[eqtl_shared_hits],w2*a2_distr[eqtl_shared_hits],w3*a3_distr[eqtl_shared_hits]))^2) 
  eqtl_gamma2_c1[kk]=mean(w1^2*a1_distr[eqtl_c1_hits]^2) 
  eqtl_gamma2_c2[kk]=mean(w2^2*a2_distr[eqtl_c2_hits]^2) 
  eqtl_gamma2_c3[kk]=mean(w3^2*a3_distr[eqtl_c3_hits]^2)  
  
  #

  gwas_shared_hits=gwas_hits[which(gwas_hits<(N*w_share))]
  gwas_c1_hits=gwas_hits[which(gwas_hits>(N*w_share) & gwas_hits<(N*w_share + N*(1-w_share)/3))]
  gwas_c2_hits=gwas_hits[which(gwas_hits>(N*w_share + N*(1-w_share)/3) & gwas_hits<(N*w_share + 2*N*(1-w_share)/3))]
  gwas_c3_hits=gwas_hits[which(gwas_hits>(N*w_share + 2*N*(1-w_share)/3) & gwas_hits<(N*w_share + 3*N*(1-w_share)/3))]
  
  
  gwas_shared_enrich[kk]=length(gwas_shared_hits) / length(gwas_hits)/ w_share
  gwas_c1_enrich[kk]=length(gwas_c1_hits) / length(gwas_hits)  / ((1-w_share)/3)
  gwas_c2_enrich[kk]=length(gwas_c2_hits) / length(gwas_hits) / ((1-w_share)/3)
  gwas_c3_enrich[kk]=length(gwas_c3_hits) / length(gwas_hits) / ((1-w_share)/3)
  
  gwas_gamma2_shared[kk]=mean(rowSums(cbind(w1*a1_distr[gwas_shared_hits],w2*a2_distr[gwas_shared_hits],w3*a3_distr[gwas_shared_hits]))^2) 
  gwas_gamma2_c1[kk]=mean(w1^2*a1_distr[gwas_c1_hits]^2) 
  gwas_gamma2_c2[kk]=mean(w2^2*a2_distr[gwas_c2_hits]^2) 
  gwas_gamma2_c3[kk]=mean(w3^2*a3_distr[gwas_c3_hits]^2) 
  
}

d_data_coloc=data.frame(power=seq(0.1,1,0.1),coloc_rate,coloc_shared_enrich,coloc_c1_enrich,coloc_c2_enrich,coloc_c3_enrich,coloc_gamma2_shared,coloc_gamma2_c1,coloc_gamma2_c2,coloc_gamma2_c3)
d_data_gwas=data.frame(gwas_shared_enrich,gwas_c1_enrich,gwas_c2_enrich,gwas_c3_enrich,gwas_gamma2_shared,gwas_gamma2_c1,gwas_gamma2_c2,gwas_gamma2_c3)
d_data_eqtl=data.frame(eqtl_shared_enrich,eqtl_c1_enrich,eqtl_c2_enrich,eqtl_c3_enrich,eqtl_gamma2_shared,eqtl_gamma2_c1,eqtl_gamma2_c2,eqtl_gamma2_c3)

d_data=cbind(d_data_coloc,d_data_gwas,d_data_eqtl)
d_data[is.na(d_data)]=0

#########

brewer.pal(n = 4, name = "Set1")
display.brewer.pal(n = 4, name = 'Set1')
col_scale=c("cell type 1"="#E41A1C","cell type 2"="#377EB8","cell type 3"="#4DAF4A","shared"="#984EA3")

########

lab_x="eQTL study power"
lab_y="Mean squared weighted genic effect\nof co-discovered SNPs"

d_share=d_data %>% dplyr::select(power,coloc_gamma2_shared); colnames(d_share)=c("power","effect2"); d_share$cat="shared"
d1=d_data %>% dplyr::select(power,coloc_gamma2_c1); colnames(d1)=c("power","effect2"); d1$cat="cell type 1"
d2=d_data %>% dplyr::select(power,coloc_gamma2_c2); colnames(d2)=c("power","effect2"); d2$cat="cell type 2"
d3=d_data %>% dplyr::select(power,coloc_gamma2_c3); colnames(d3)=c("power","effect2"); d3$cat="cell type 3"
d_plot=rbind(d_share,d1,d2,d3)

p1<-ggplot(d_plot,aes(x=power, y = effect2,col=cat)) +
  geom_point(size=1)  +
  geom_line(size=0.5)  + 
  geom_segment(aes(x=0.1,xend=1,y=effect2[10],yend=effect2[10]) , linetype="dashed", size=0.5, color="#984EA3")+
  geom_segment(aes(x=0.1,xend=1,y=effect2[20],yend=effect2[20]) , linetype="dashed", size=0.5, color="#E41A1C")+
  geom_segment(aes(x=0.1,xend=1,y=effect2[30],yend=effect2[30]) , linetype="dashed", size=0.5, color="#377EB8")+
  geom_segment(aes(x=0.1,xend=1,y=effect2[40],yend=effect2[40]) , linetype="dashed", size=0.5, color="#4DAF4A")+
  scale_color_manual(name = "context",values = col_scale)+
  xlab(lab_x)+
  ylab(lab_y)+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_cell_type_mixture_coloc_effect_by_power.pdf")
save_plot(outfile, p1,base_width = 3, base_height = 2.5)

########

lab_x="eQTL study power"
lab_y="Cell type enrichment\nof co-discovered SNPs"

d0=d_data %>% dplyr::select(power,coloc_shared_enrich); colnames(d0)=c("power","enrich"); d0$cat="shared"
d1=d_data %>% dplyr::select(power,coloc_c1_enrich); colnames(d1)=c("power","enrich"); d1$cat="cell type 1"
d2=d_data %>% dplyr::select(power,coloc_c2_enrich); colnames(d2)=c("power","enrich"); d2$cat="cell type 2"
d3=d_data %>% dplyr::select(power,coloc_c3_enrich); colnames(d3)=c("power","enrich"); d3$cat="cell type 3"
d_plot=rbind(d0,d1,d2,d3)

p1<-ggplot(d_plot,aes(x=power, y = enrich,col=cat)) +
  geom_point(size=1)  +
  geom_line(size=0.5)  + 
  geom_segment(aes(x=0.1,xend=1,y=1,yend=1) , linetype="dashed", size=0.5, color="black")+
  scale_color_manual(name = "context",values = col_scale)+
  xlab(lab_x)+
  ylab(lab_y)+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_cell_type_mixture_coloc_rate_by_power.pdf")
save_plot(outfile, p1,base_width = 3, base_height = 2.5)


#########

col_scale=c("GWAS"="#1E56A0","eQTL"="#FF5959")
col_scale=c("cell type 1"="#E41A1C","cell type 2"="#377EB8","cell type 3"="#4DAF4A","shared"="#984EA3")

d0=d_data %>% dplyr::select(power,gwas_shared_enrich); colnames(d0)=c("power","enrich"); d0$cat="shared"
d1=d_data %>% dplyr::select(power,gwas_c1_enrich); colnames(d1)=c("power","enrich"); d1$cat="cell type 1"
d2=d_data %>% dplyr::select(power,gwas_c2_enrich); colnames(d2)=c("power","enrich"); d2$cat="cell type 2"
d3=d_data %>% dplyr::select(power,gwas_c3_enrich); colnames(d3)=c("power","enrich"); d3$cat="cell type 3"
d_gwas=rbind(d0,d1,d2,d3)
d_gwas$study="GWAS"

d0=d_data %>% dplyr::select(power,eqtl_shared_enrich); colnames(d0)=c("power","enrich"); d0$cat="shared"
d1=d_data %>% dplyr::select(power,eqtl_c1_enrich); colnames(d1)=c("power","enrich"); d1$cat="cell type 1"
d2=d_data %>% dplyr::select(power,eqtl_c2_enrich); colnames(d2)=c("power","enrich"); d2$cat="cell type 2"
d3=d_data %>% dplyr::select(power,eqtl_c3_enrich); colnames(d3)=c("power","enrich"); d3$cat="cell type 3"
d_eqtl=rbind(d0,d1,d2,d3)
d_eqtl$study="eQTL"

d_plot=rbind(d_gwas,d_eqtl)

lab_x="Power"
lab_y="Cell type enrichment\nof discovered SNPs"

p1<-ggplot(d_plot,aes(x=power, y = enrich,col=cat)) +
  geom_point(size=1)  +
  geom_line(size=0.5,aes(linetype=study))  + 
  scale_color_manual(name="context",values=col_scale)+
  geom_segment(aes(x=0.1,xend=1,y=1,yend=1) , linetype="dashed", size=0.5, color="black")+
  xlab(lab_x)+
  ylab(lab_y)+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_cell_type_mixture_discovery_by_power.pdf")
save_plot(outfile, p1,base_width = 3, base_height = 2.5)

#

col_scale=c("cell type 1"="#E41A1C","cell type 2"="#377EB8","cell type 3"="#4DAF4A","shared"="#984EA3")

lab_x="GWAS study power"
lab_y="Cell type enrichment\nof discovered SNPs"

p1<-ggplot(d_gwas,aes(x=power, y = enrich,col=cat)) +
  geom_point(size=1)  +
  geom_line(size=0.5)  + 
  scale_color_manual(name="context",values=col_scale)+
  geom_segment(aes(x=0.1,xend=1,y=1,yend=1) , linetype="dashed", size=0.5, color="black")+
  xlab(lab_x)+
  ylab(lab_y)+
  xlim(c(0,1))+
  ylim(c(0,3.5))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_cell_type_mixture_discovery_by_power_GWAS.pdf")
save_plot(outfile, p1,base_width = 3, base_height = 2.5)

#

col_scale=c("cell type 1"="#E41A1C","cell type 2"="#377EB8","cell type 3"="#4DAF4A","shared"="#984EA3")

lab_x="eQTL study power"
lab_y="Cell type enrichment\nof discovered SNPs"

p1<-ggplot(d_eqtl,aes(x=power, y = enrich,col=cat)) +
  geom_point(size=1)  +
  geom_line(size=0.5)  + 
  scale_color_manual(name="context",values=col_scale)+
  geom_segment(aes(x=0.1,xend=1,y=1,yend=1) , linetype="dashed", size=0.5, color="black")+
  xlab(lab_x)+
  ylab(lab_y)+
  xlim(c(0,1))+
  ylim(c(0,3.5))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_cell_type_mixture_discovery_by_power_eQTL.pdf")
save_plot(outfile, p1,base_width = 3, base_height = 2.5)

####

d0=d_data %>% dplyr::select(power,gwas_gamma2_shared); colnames(d0)=c("power","effect2"); d0$cat="shared"
d1=d_data %>% dplyr::select(power,gwas_gamma2_c1); colnames(d1)=c("power","effect2"); d1$cat="cell type 1"
d2=d_data %>% dplyr::select(power,gwas_gamma2_c2); colnames(d2)=c("power","effect2"); d2$cat="cell type 2"
d3=d_data %>% dplyr::select(power,gwas_gamma2_c3); colnames(d3)=c("power","effect2"); d3$cat="cell type 3"
d_gwas=rbind(d0,d1,d2,d3)
d_gwas$study="GWAS"

d0=d_data %>% dplyr::select(power,eqtl_gamma2_shared); colnames(d0)=c("power","effect2"); d0$cat="shared"
d1=d_data %>% dplyr::select(power,eqtl_gamma2_c1); colnames(d1)=c("power","effect2"); d1$cat="cell type 1"
d2=d_data %>% dplyr::select(power,eqtl_gamma2_c2); colnames(d2)=c("power","effect2"); d2$cat="cell type 2"
d3=d_data %>% dplyr::select(power,eqtl_gamma2_c3); colnames(d3)=c("power","effect2"); d3$cat="cell type 3"
d_eqtl=rbind(d0,d1,d2,d3)
d_eqtl$study="eQTL"

d_plot=rbind(d_gwas,d_eqtl)

lab_x="Power"
lab_y="Mean squared weighted genic effect\nof discovered SNPs"

col_scale=c("cell type 1"="#E41A1C","cell type 2"="#377EB8","cell type 3"="#4DAF4A","shared"="#984EA3")

p1<-ggplot(d_plot,aes(x=power, y = effect2,col=cat)) +
  geom_point(size=1,alpha=0.8)  +
  geom_line(size=0.5,alpha=0.6,aes(linetype=study))  + 
  scale_color_manual(name="context",values=col_scale)+
  xlab(lab_x)+
  ylab(lab_y)+
  xlim(c(0,1))+
  theme_cowplot(9)


outfile=paste0(outdir,"selection_model_cell_type_mixture_effect_by_power.pdf")
save_plot(outfile, p1,base_width = 3, base_height = 2.5)

############

plot_size=2.3

N=10000000
x1=0.05
x2=0.2
x3=0.75


a1=rnorm(N*x1,0,4);a2=rnorm(N*x2,0,2);a3=rnorm(N*x3,0,1)

dx1=data.frame(x=a1,cat="cell type 1")
dx2=data.frame(x=a2,cat="cell type 2")
dx3=data.frame(x=a3,cat="cell type 3")

d_plot=rbind(dx1,dx2,dx3)

lab_x=("Genic effect size (per cell)")

col_scale=c("cell type 1"="#E41A1C","cell type 2"="#377EB8","cell type 3"="#4DAF4A","shared"="#984EA3")


p0<-ggplot(d_plot, aes(x=x, color=cat)) +
  geom_density(show.legend=FALSE)+
  stat_density(geom="line",position="identity")+
  scale_color_manual(name="context",values=col_scale)+
  xlab(lab_x)+
  theme_cowplot(9)+
  theme(legend.position = "none")

outfile=paste0(outdir,"selection__model_cell_type_mixture_gamma_distr.pdf")
save_plot(outfile, p0,base_width = plot_size, base_height = 0.75*plot_size)
