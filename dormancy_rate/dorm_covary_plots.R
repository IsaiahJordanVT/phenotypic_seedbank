library(ggplot2)
library(scales)
library(geomtextpath)

file_path <- '~/Documents/GitHub/publication_model/dormancy_rate/dormancy_rate 200 1 2026-10-14 15:23:37.425772'
results_1 <- read.csv(file = file_path, header = FALSE)
results = aggregate(results_1[,1],by=results_1[,c(2,3)],FUN = mean)
colnames(results) = c("Dormancy","Resuscitation",'Survival')

frac_pops = data.frame(matrix(nrow = 0,ncol=3))
param_space = 10^seq(-14,-1,.5)
res_param_space = 10^seq(-13,-.1,.75)
colnames(frac_pops)= c("Dormancy",'Resuscitation','Fraction')
for(i in param_space){
  for(k in res_param_space){
    pop = vbmax_equibpop_determination(i,k,spont_death)[1]
    frac_pop = sum(pop[[1]][[2]])/sum(pop[[1]][[1]]+pop[[1]][[2]])
    if(frac_pop > 1){
      frac_pop = NA
    } else if(frac_pop < 0){
      frac_pop = NA
    }
    new_row = c(i,k,frac_pop)
    frac_pops = rbind(frac_pops,new_row)
  }
}
colnames(frac_pops)= c("Dormancy",'Resuscitation','Fraction')
ggplot(data = results, aes(x = as.numeric(Dormancy), y = as.numeric(Resuscitation), z = as.numeric(Survival))) +
  geom_contour_filled(breaks = c(seq(0, 1.1, .1))) +
  xlab("Dormancy") +
  ylab("Resuscitation") +
  geom_textcontour(data = subset(frac_pops, Fraction > 0),
                   aes(
                     x = Dormancy, 
                     y = Resuscitation, 
                     z = Fraction,
                     label = after_stat(paste0("10^", log10(level)))
                   ),
                   breaks = c(.001, .0001, .00001, .000001), 
                   color = 'white',
                   size = 8,
                   linewidth = 1,
                   gap = TRUE,      
                   parse = TRUE     
  ) +
  scale_fill_viridis_d(begin = 0,
                       end = .8,
                       name = "Proportion\nSurviving",
                       direction = 1,
                       guide = guide_colorsteps(reverse = F,
                                                barheight = 15))+
  scale_x_log10(
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(
                labels = trans_format("log10", math_format(10^.x))) +
  
  scale_size_binned(n.breaks = 25) + 
  geom_point(aes(x = dormancy_rate, y = resuscitation_rate), color = 'orange', size = 4,pch = 8) +
  theme_classic()+
  theme(text = element_text(size = 25))
mixed = results
mixed$Fraction = NA
for(i in 1:nrow(frac_pops)){
  mixed$Fraction[which(mixed$Dormancy == frac_pops[i,1]&mixed$Resuscitation == frac_pops[i,2])] = frac_pops$Fraction[i]
}
