library(ggplot2)
library(viridis)
library(reshape)
library(metR)
library(akima)
library(interp)
exper = read.csv("~/Documents/GitHub/publication_model/hetero_cost/het_vs_dorm 1 50 2025-12-08 140823.854472",header = F)
exper = exper[-which(exper$V4 < 1),]
exper2 = read.csv("~/Documents/GitHub/publication_model/hetero_cost/het_vs_dorm 2 50 2025-12-08 14:07:50.552769", header = F)
exper3 = read.csv("~/Documents/GitHub/publication_model/hetero_cost/het_vs_dorm 3 50 2025-12-08 14:07:06" ,header = F)
exper4 = read.csv("~/Documents/GitHub/publication_model/hetero_cost/het_vs_dorm 4 50 2025-12-05 160703.319167" ,header = F)
exper2 = exper2[-1,]
exper3 = exper3[-1,]
exper4 = exper4[-1,]
cn = as.vector(exper[1,])
cn[5]= 'rep'
exper = exper[-1,]
exper=rbind(exper,exper2,exper3,exper4)
colnames(exper) = cn
exper[,2] = as.factor(exper[,2])
exper[,3] = as.numeric(exper[,3])
exper[,1] = as.numeric(exper[,1])
exper = aggregate(exper[,c(1)],by = list(exper[,2],exper[,3],exper[,4]),FUN = sum)
colnames(exper) = c('Lifestyle','Degree_of_heterogeneity','Cost','Percent')
exper[,3] = as.numeric(exper[,3])
exper$Cost = 1/exper$Cost
exper = exper[which(exper$Cost<.5),]
# exper = exper[which(exper$Degree_of_heterogeneity<.025),]
exper$Percent = exper$Percent
het_only = exper[which(exper[,1]=="Heterogeneity Only"),]
het_dorm = exper[which(exper[,1]=="Dormancy and Heterogeneity"),]  
intrpl = 1000
het_only2 = interp(as.numeric(het_only$Cost),as.numeric(het_only$Degree_of_heterogeneity),as.numeric(het_only$Percent),
                  xo = seq(min(as.numeric(het_only$Cost)),max(as.numeric(het_only$Cost)),length = intrpl),
                  yo=seq(min(as.numeric(het_only$Degree_of_heterogeneity)),max(as.numeric(het_only$Degree_of_heterogeneity)),length = intrpl)
                  )
het_dorm2 = interp(as.numeric(het_dorm$Cost),as.numeric(het_dorm$Degree_of_heterogeneity),as.numeric(het_dorm$Percent),
                   xo = seq(min(as.numeric(het_dorm$Cost)),max(as.numeric(het_dorm$Cost)),length = intrpl),
                   yo=seq(min(as.numeric(het_dorm$Degree_of_heterogeneity)),max(as.numeric(het_dorm$Degree_of_heterogeneity)),length= intrpl)
                   )
ho = data.frame(interp2xyz(het_only2))
hd = data.frame(interp2xyz(het_dorm2))
colnames(ho) = c('Cost','Degree_of_heterogeneity','Percent')
colnames(hd) = c('Cost','Degree_of_heterogeneity','Percent')
# ho = ho[which(ho$Cost<=7),]
# hd = hd[which(hd$Cost<=7),]
ho_plot = ggplot()+
  geom_contour_filled(data = ho,aes(x = as.numeric(Cost),y=as.numeric(Degree_of_heterogeneity),z = as.numeric(Percent/200)),bins = 20)+
  # geom_vline(xintercept  = .2)+
  # geom_point(data = exper,aes(x = Cost,y=Degree_of_heterogeneity), size = .25)+
  # ggtitle('Only Heterogeneity')+
  xlab("Cost")+
  ylab("Amount of Constitutive Expression")+
  # scale_fill_viridis()+
  # scale_x_discrete(label=abbreviate)+
  # scale_x_log10(trans = 'log10')+
  theme_classic()+
  theme(text = element_text(size = 35),legend.position = 'none')
plot(ho_plot)
k = seq(1,100,1)
r = matrix(rep(k,100),nrow = 100)
r = cbind(melt(data.frame(r)),seq(1,1000))
colnames(r) = c('Cost','Degree_of_heterogeneity','Percent')
hd_plot = ggplot()+
  geom_contour_filled(data = hd,aes(x = Cost,y=Degree_of_heterogeneity/280.4095 * 30,z = Percent/200),bins = 20)+
  # geom_point(data = exper,aes(x = as.numeric(Cost),y=as.numeric(Degree_of_heterogeneity)),size = .25)+
  # ggtitle('Heterogeneity and Dormancy')+
  # geom_vline(xintercept = 0.2)+
  ylab("Amount of Constitutive Expression")+
  xlab("Cost")+
  # scale_x_log10(trans = 'log10')+
  # scale_fill_viridis()+
  # scale_x_discrete(label=abbreviate)+
  theme_classic()+
  theme(text = element_text(size = 35),legend.position = 'none')
plot(hd_plot) 
leg_data = data.frame(seq(1,100,1),seq(1,100,1))
colnames(leg_data) = c("Value","Percent")
legend = ggplot()+
  geom_point(data = leg_data,aes(x=Value,y=Percent,color = Value))+
  # scale_color_continuous(palette('iridis'))
  geom_point(data = exper,aes(x = Cost,y=Degree_of_heterogeneity/280.4095 * 30))+
  # ggtitle('Only Heterogeneity')+
  # xlab("Cost")+
  # ylab("Amount of heterogeneity")+
  scale_fill_brewer(palette = viridis)+
  # scale_fill_viridis()+
  # scale_x_discrete(label=abbreviate)+
  # scale_x_log10(trans = 'log10')+
  theme_classic()+
  theme(text = element_text(size = 25))
  