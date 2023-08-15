library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(sfsmisc)
library(cowplot)
theme_set(theme_cowplot(font_size = 9))

##############
#simulate variants and effect sizes

N=10000000 #number of variants
a_distr=rnorm(N,0,1) #gene effect distr
b_distr=rnorm(N,0,1) #variant effect distr

K=2.986 #selection parameter; see Supp Note

v_neutral_distr=a_distr^2*b_distr^2 #phenotypic variance under neutrality
v_plat_distr=K*(1-exp(-1*a_distr^2*b_distr^2/K)) #phenotypic variance under selection

v_p_neutral_distr=v_neutral_distr/(a_distr^2*b_distr^2) #allele frequency variance under neutrality
v_p_plat_distr=v_plat_distr/v_neutral_distr #allele frequency variance under selection

v_e_neutral_distr=v_p_neutral_distr*b_distr^2 #expression variance under neutrality
v_e_plat_distr=v_p_plat_distr*b_distr^2 #expression variance under selection

max_xy=qchisq(0.95, df=1,lower.tail = T) #plot region

##############
#grid over the parameter space

a_sq_max=8 
b_sq_max=8
diff=0.0025

a_sq=seq(diff,a_sq_max,diff)
b_sq=seq(diff,b_sq_max,diff)
abg <- xy.grid(b_sq,a_sq)

#compute variance terms for each point on the grid
v_plat=K*(1-exp(-1*abg[,1]*abg[,2]/K))
v_neutral=abg[,1]*abg[,2]

v_p_neutral=v_neutral/(abg[,1]*abg[,2])
v_p_plat=v_plat/(abg[,1]*abg[,2])

abg=cbind(abg,v_p_neutral)
abg=cbind(abg,v_p_plat)

v_e_neutral=v_p_neutral*abg[,1]
v_e_plat=v_p_plat*abg[,1]
v_e_mild=v_p_mild*abg[,1]

abg=cbind(abg,v_e_neutral)
abg=cbind(abg,v_e_plat)

v_y_neutral=v_p_neutral*abg[,1]*abg[,2]
v_y_plat=v_p_plat*abg[,1]*abg[,2]

abg=cbind(abg,v_y_neutral)
abg=cbind(abg,v_y_plat)

################
#plot parameters

lab_x=expression("Variant effect on expression,"~paste(beta^2))
lab_y=expression("Gene effect on trait,"~paste(gamma^2))

col_eqtl="#FF5959"
col_gwas="#1E56A0"

plot_size=2.3
jj=850+1 #power percentile corresponding to 15% discovery

########SELECTION########
#########################

###find grid plots near eQTL discovery line at 15%
selection_type=6 #column 6; eQTL discovery under selection
x=quantile(v_e_plat_distr,probs = seq(0,1,0.001))

c=x[jj]
c1=x[jj-1]
c2=x[jj+1]

d_data=data.frame(b=as.numeric(),a=as.numeric(),c=as.numeric()) 
d_fit=data.frame(b=as.numeric(),a=as.numeric(),c=as.numeric()) #fit to smoothen the discovery line

index=which(abg[,selection_type]>(c1) & abg[,selection_type]<(c2))
abg_plot=abg[index,]
abg_plot=abg_plot[abg_plot[,1]<max_xy & abg_plot[,2]<max_xy,]

d_data_temp=data.frame(b=abg_plot[,1],a=abg_plot[,2])
d_data_temp$c=c
d_data=rbind(d_data,d_data_temp)

lo <- loess(abg_plot[,2]~abg_plot[,1],span=0.10)
fit_plot=cbind(abg_plot[,1],predict(lo))
fit_plot=fit_plot[order(fit_plot[,1],decreasing=FALSE),]
fit_plot=data.frame(b=fit_plot[,1],a=fit_plot[,2])
fit_plot$c=c
d_fit=rbind(d_fit,fit_plot)

d_data_eqtl=d_data
d_fit_eqtl=d_fit
d_fit_eqtl$cat="eqtl"
d_fit_eqtl$lim=min(d_fit_eqtl$a)

###find grid plots near GWAS discovery line at 15%
selection_type=8 #column 8; GWAS discovery under selection
x=quantile(v_plat_distr,probs = seq(0,1,0.001))

c=x[jj]
c1=x[jj-1]
c2=x[jj+1]

d_data=data.frame(b=as.numeric(),a=as.numeric(),c=as.numeric())
d_fit=data.frame(b=as.numeric(),a=as.numeric(),c=as.numeric())

index=which(abg[,selection_type]>(c1) & abg[,selection_type]<(c2))
abg_plot=abg[index,]
abg_plot=abg_plot[abg_plot[,1]<max_xy & abg_plot[,2]<max_xy,]

d_data_temp=data.frame(b=abg_plot[,1],a=abg_plot[,2])
d_data_temp$c=c
d_data=rbind(d_data,d_data_temp)

lo <- loess(abg_plot[,2]~abg_plot[,1],span=0.10)
fit_plot=cbind(abg_plot[,1],predict(lo))
fit_plot=fit_plot[order(fit_plot[,1],decreasing=FALSE),]
fit_plot=data.frame(b=fit_plot[,1],a=fit_plot[,2])
fit_plot$c=c
d_fit=rbind(d_fit,fit_plot)

d_data_gwas=d_data
d_fit_gwas=d_fit
d_fit_gwas$cat="gwas"
d_fit_gwas$lim=max(d_fit_gwas$a)

#####

d_plot=rbind(d_fit_gwas,d_fit_eqtl)

pA<-ggplot(NULL,aes(x=b, y = a)) +
  geom_line(data=d_fit_gwas,size=0.5,color=col_gwas)  + 
  geom_ribbon(data=d_fit_gwas,aes(ymin = a, ymax = lim), fill = col_gwas, alpha = 0.3)+
  geom_line(data=d_fit_eqtl,size=0.5,color=col_eqtl)  + 
  geom_ribbon(data=d_fit_eqtl,aes(ymin = lim, ymax = a), fill = col_eqtl, alpha = 0.3)+
  xlim(min(d_plot$lim),max(d_plot$lim))+
  ylim(min(d_plot$lim),max(d_plot$lim))+
  xlab(lab_x)+
  ylab(lab_y)+
  theme(axis.title.x=element_text(size=10),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none",
        axis.line = element_line(colour = "black"))


outfile="model_simulations/selection_model_discovery.pdf"
save_plot(outfile, pA,base_width = plot_size, base_height = plot_size)


########NEUTRAL########
#######################

###find grid plots near eQTL discovery line at 15%
selection_type=5 #column 5; eQTL discovery under neutrality
x=quantile(v_e_neutral_distr,probs = seq(0,1,0.001))

c=x[jj]
c1=x[jj-1]
c2=x[jj+1]

d_data=data.frame(b=as.numeric(),a=as.numeric(),c=as.numeric())
d_fit=data.frame(b=as.numeric(),a=as.numeric(),c=as.numeric())

index=which(abg[,selection_type]>(c1) & abg[,selection_type]<(c2))
abg_plot=abg[index,]
abg_plot=abg_plot[abg_plot[,1]<max_xy & abg_plot[,2]<max_xy,]

d_data_temp=data.frame(b=abg_plot[,1],a=abg_plot[,2])
d_data_temp$c=c
d_data=rbind(d_data,d_data_temp)

lo <- loess(abg_plot[,1]~abg_plot[,2],span=0.10)
fit_plot=cbind(predict(lo),abg_plot[,2])
fit_plot=fit_plot[order(fit_plot[,1],decreasing=FALSE),]
fit_plot=fit_plot[which(fit_plot[,2]<=max_xy & fit_plot[,2]>=0),]
fit_plot=data.frame(b=fit_plot[,1],a=fit_plot[,2])
fit_plot$c=c
d_fit=rbind(d_fit,fit_plot)

d_data_eqtl=d_data
d_fit_eqtl=d_fit
d_fit_eqtl$cat="eqtl"

###find grid plots near GWAS discovery line at 15%
selection_type=7 #column 5; GWAS discovery under neutrality
x=quantile(v_neutral_distr,probs = seq(0,1,0.001))

c=x[jj]
c1=x[jj-1]
c2=x[jj+1]

d_data=data.frame(b=as.numeric(),a=as.numeric(),c=as.numeric())
d_fit=data.frame(b=as.numeric(),a=as.numeric(),c=as.numeric())

index=which(abg[,selection_type]>(c1) & abg[,selection_type]<(c2))
abg_plot=abg[index,]
abg_plot=abg_plot[abg_plot[,1]<=max_xy & abg_plot[,2]<=max_xy,]

d_data_temp=data.frame(b=abg_plot[,1],a=abg_plot[,2])
d_data_temp$c=c
d_data=rbind(d_data,d_data_temp)

lo <- loess(abg_plot[,2]~abg_plot[,1],span=0.10)
fit_plot=cbind(abg_plot[,1],predict(lo))
fit_plot=fit_plot[order(fit_plot[,1],decreasing=FALSE),]
fit_plot=fit_plot[which(fit_plot[,2]<=max_xy & fit_plot[,2]>=0),]
fit_plot=data.frame(b=fit_plot[,1],a=fit_plot[,2])
fit_plot$c=c
d_fit=rbind(d_fit,fit_plot)

d_data_gwas=d_data
d_fit_gwas=d_fit
d_fit_gwas$cat="gwas"
d_fit_gwas$lim=max(d_fit_gwas$a)

#####

B_star=mean(d_fit_eqtl$b) #straight vertical line for eQTL discovery

pB<-ggplot(NULL,aes(x=b, y = a)) +
  geom_line(data=d_fit_gwas,size=0.5,color=col_gwas)  + 
  geom_ribbon(data=d_fit_gwas,aes(ymin = a, ymax = lim), fill = col_gwas, alpha = 0.3)+
  geom_segment(aes(x = B_star, xend = B_star, y = 0, yend = max(d_fit_gwas$a)),color=col_eqtl,size=0.5)+
  annotate("rect", xmin = B_star, xmax = max(d_fit_gwas$b), ymin = 0, ymax = max(d_fit_gwas$a), fill=col_eqtl, alpha = 0.3)+
  xlim(min(d_plot$lim),max(d_fit_gwas$b))+
  ylim(0,max(d_fit_gwas$a))+
  xlab(lab_x)+
  ylab(lab_y)+
  theme(axis.title.x=element_text(size=10),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none",
        axis.line = element_line(colour = "black"))


outfile="model_simulations/neutral_model_discovery.pdf"
save_plot(outfile, pB,base_width = plot_size, base_height = plot_size)

