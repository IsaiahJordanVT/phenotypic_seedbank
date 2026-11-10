library(ggplot2)
library(viridis)
library(reshape)
library(metR)
exper = read.csv("~/Documents/GitHub/Phenotypic_seedbank/hetero_dorm/het_vs_dorm 1 500 2026-08-14 13:14:26.737561",header = F)
exper = exper[-which(exper$V4 < 1),]
exper[,2] = as.factor(exper[,2])
exper[,3] = as.numeric(exper[,3])
exper[,1] = as.numeric(exper[,1])
exper = aggregate(exper[,c(1)],by = list(exper[,2],exper[,3],exper[,4]),FUN = mean)
colnames(exper) = c('Phenotype','Degree_of_heterogeneity','Cost','Percent')
exper[,3] = as.numeric(exper[,3])
exper$Cost = 1/exper$Cost
exper = exper[which(exper$Cost<.5),]
het_only = exper[which(exper[,1]=="Heterogeneity Only"),]
het_dorm = exper[which(exper[,1]=="Dormancy and Heterogeneity"),]  
ho = het_only
hd = het_dorm
colnames(ho) = c('Phenotype','Degree_of_heterogeneity','Cost','Percent')
colnames(hd) = c('Phenotype','Degree_of_heterogeneity','Cost','Percent')
ho = ho[which(ho$Cost<=7),]
hd = hd[which(hd$Cost<=7),]
ho_plot = ggplot()+
  geom_contour_filled(data = ho,aes(x = as.numeric(Cost),y=as.numeric(Degree_of_heterogeneity),z = as.numeric(Percent)),breaks = c(seq(0, 1.1, .075)))+
  xlab("Cost")+
  ylab("Amount of Constitutive Expression")+
  scale_fill_viridis_d(begin = 0,
                       end = 1,
                       name = "Proportion\nSurviving",
                       direction = 1,
                       guide = guide_colorsteps(reverse = F,
                                                barheight = 15))+
  scale_y_continuous(
    # labels = function(x) {
    #   # Convert numbers to standard e-notation
    #   labs <- scales::scientific(x)
    #   # Reformat e-notation to plotmath style
    #   labs <- gsub("0e\\+0", "0", labs)
    #   labs <- gsub("e\\+0", " %*% 10^", labs)
    #   labs <- gsub("e\\+", " %*% 10^", labs)
    #   labs <- gsub("e-0", " %*% 10^-", labs)
    #   labs <- gsub("e-", " %*% 10^-", labs)
    #   # # Clean up 0 and clean up 1 * 10^y to just 10^y if preferred
    #   labs <- gsub("^0 %\\*% .*$", "0", labs)
    #   # labs <- gsub("^1 %\\*% ", "", labs) 
    #   parse(text = labs)
    # },
    breaks = seq(0,.0012,.0003)
  )+
  theme_classic()+
  theme(text = element_text(size = 28),legend.position = 'none')
plot(ho_plot)
hd_plot = ggplot()+
  geom_contour_filled(data = hd,aes(x = as.numeric(Cost),y=as.numeric(Degree_of_heterogeneity),z = as.numeric(Percent)),breaks = c(seq(0, 1.1, .075)),bins = 50)+
  xlab("Cost")+
  ylab("Amount of Constitutive Expression")+
  scale_fill_viridis_d(begin = 0,
                       end = 1,
                       name = "Proportion\nSurviving",
                       direction = 1,
                       guide = guide_colorsteps(reverse = F,
                                                barheight = 15))+
  scale_y_continuous(
    # labels = function(x) {
    #   # Convert numbers to standard e-notation
    #   labs <- scales::scientific(x)
    #   # Reformat e-notation to plotmath style
    #   labs <- gsub("0e\\+0", "0", labs)
    #   labs <- gsub("e\\+0", " %*% 10^", labs)
    #   labs <- gsub("e\\+", " %*% 10^", labs)
    #   labs <- gsub("e-0", " %*% 10^-", labs)
    #   labs <- gsub("e-", " %*% 10^-", labs)
    #   # # Clean up 0 and clean up 1 * 10^y to just 10^y if preferred
    #   labs <- gsub("^0 %\\*% .*$", "0", labs)
    #   # labs <- gsub("^1 %\\*% ", "", labs)
    # 
    #   parse(text = labs)
    # },
    breaks = seq(0,.0012,.0003)
  )+
  theme_classic()+
  theme(text = element_text(size = 28),legend.position = 'none')
plot(hd_plot) 
# leg_data = data.frame(seq(1,100,1),seq(1,100,1))
# colnames(leg_data) = c("Value","Percent")
# legend = ggplot()+
#   geom_point(data = leg_data,aes(x=Value,y=Percent,color = Value))+
#   # scale_color_continuous(palette('iridis'))
#   geom_point(data = exper,aes(x = Cost,y=Degree_of_heterogeneity/280.4095 * 30))+
#   # ggtitle('Only Heterogeneity')+
#   # xlab("Cost")+
#   # ylab("Amount of heterogeneity")+
#   scale_fill_brewer(palette = viridis)+
#   # scale_fill_viridis()+
#   # scale_x_discrete(label=abbreviate)+
#   # scale_x_log10(trans = 'log10')+
#   theme_classic()+
#   theme(text = element_text(size = 25))
ggsave("~/Documents/GitHub/Phenotypic_seedbank//plots/fig5/fig5B.png",hd_plot,height = 6.64, width = 7.62)
ggsave("~/Documents/GitHub/Phenotypic_seedbank//plots/fig5/fig5A.png",ho_plot,height = 6.64, width = 7.62)
