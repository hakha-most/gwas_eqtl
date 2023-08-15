library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(sfsmisc)
library(cowplot)
library(MASS)
library(segmented)

theme_set(theme_cowplot(font_size = 9))

##########

outdir="model_simulations/"

##########

S=-0.25 #selection exponent

N=20000000
a_distr=rnorm(N,0,1)
b_distr=rnorm(N,0,1)

###
#exponential allele freq. distr.
u=runif(N,min=0,max=1); a=20; p_min=0.001; p=-1/a*(log(exp(-a*p_min)-u*(exp(-a*p_min)-exp(-a*0.5)))); p1=p   # minor allele frequency distribution (a: exponential width chosen arbitrarily)

#scale effect sizes
effect_unadj=a_distr^2*b_distr^2
effect_adj=effect_unadj * (2*p*(1-p))^(S) 
f_root=sqrt((2*p*(1-p))^(S) )

a_distr=a_distr*sqrt(f_root)
b_distr=b_distr*sqrt(f_root)

a_distr=a_distr/sd(a_distr)
b_distr=b_distr/sd(b_distr)

###
#variance terms
effect_distr=a_distr^2*b_distr^2
v_p_distr=2*p*(1-p)
v_distr=v_p_distr * a_distr^2*b_distr^2
v_e_distr=v_p_distr * b_distr^2

#########
#alpha model describes E[effect^2 | p]
#here we find E[p | effect^2]

v_diff=0.02
v_grid=seq(v_diff,16,v_diff)
p_expected=rep(0,length(v_grid))
for (iter in 1:length(v_grid)){
  print(iter)
  v_index=which(effect_distr>(v_grid[iter]-v_diff) & effect_distr<(v_grid[iter]+v_diff))
  p_expected[iter]=mean(p[v_index])
  
}

V_p_expected=2*p_expected*(1-p_expected)

#piecewise linear fit
lin.mod <- lm(V_p_expected~v_grid) 
segmented.mod <- segmented(lin.mod, seg.Z = ~v_grid, psi=NA,control=seg.control(display=TRUE, K=3, fix.npsi = T))

plot(v_grid,V_p_expected)
plot(segmented.mod, add=T)

##########
#plot discovery lines

max_xy=max_xy=qchisq(0.95, df=1,lower.tail = T)
lab_x=expression("Variant effect on expression,"~paste(beta^2))
lab_y=expression("Gene effect on trait,"~paste(gamma^2))

col_eqtl="#FF5959"
col_gwas="#1E56A0"

plot_size=2.3
jj=850+1 #power percentile

#grid
a_sq_max=4
b_sq_max=4
diff=0.0025

a_sq=seq(diff,a_sq_max,diff)
b_sq=seq(diff,b_sq_max,diff)
abg <- xy.grid(b_sq,a_sq)

#allele freq variance for grid points under the learned model given effect^2
if (S!=0){v_p_plat=predict(segmented.mod,newdata=data.frame(v_grid=abg[,1]*abg[,2]))}
if (S==0){v_p_plat=predict(lin.mod,newdata=data.frame(v_grid=abg[,1]*abg[,2]))}

v_plat=abg[,1]*abg[,2]*v_p_plat
v_e_plat=abg[,1]*v_p_plat

abg=cbind(abg,v_p_plat,v_e_plat,v_plat)

###find grid plots near eQTL discovery line at 15%
selection_type=4
x=quantile(v_e_distr,probs = seq(0,1,0.001))

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

lo <- loess(abg_plot[,2]~abg_plot[,1],span=1)
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
selection_type=5
x=quantile(v_distr,probs = seq(0,1,0.001))

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

################

d_plot=rbind(d_fit_gwas,d_fit_eqtl)

pA<-ggplot(NULL,aes(x=b, y = a)) +
  geom_line(data=d_fit_gwas,size=0.5,color=col_gwas)  + 
  geom_ribbon(data=d_fit_gwas,aes(ymin = a, ymax = lim), fill = col_gwas, alpha = 0.3)+
  geom_line(data=d_fit_eqtl,size=0.5,color=col_eqtl)  + 
  geom_ribbon(data=d_fit_eqtl,aes(xmin = b, xmax = max_xy), fill = col_eqtl, alpha = 0.3)+
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


outfile=paste0(outdir,"selection_model_discovery_alpha_model_",as.character(abs(S)),".pdf")
save_plot(outfile, pA,base_width = plot_size, base_height = plot_size)

