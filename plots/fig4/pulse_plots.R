library(Hmisc)
library(ggplot2)
results = read.csv("~/Documents/GitHub/publication_model/pulse/pulse_intensity 200 2025-11-21 16:11:04.116009 .txt",header = F)
# res2 = read.csv("~/Documents/GitHub/publication_model/results/pulse_intensity 10 2025-11-13 13:48:07.989277 .txt",header = F)
# results = rbind(results,res2)
colnames(results) = c('survived',"pulse_intensity",'rep','final_time','dormancy rate','seed')
results_sums = aggregate(results[,c(1)],by = list(results[,2],results[,5]),FUN = sum)
print(results_sums)
num_reps = c()
for (i in unique(results$pulse_intensity)){
num_reps = c(length(which(results$pulse_intensity == i)),num_reps)
}
num_reps = 200
results_sums[,3]=results_sums[,3]/num_reps
colnames(results_sums) = c('pulse_intensity','dormancy_rate','survived')
results_sums$dormancy_rate = as.factor(results_sums$dormancy_rate)
results_sums = as.data.frame(results_sums)
# results_sums = results_sums[-which(results_sums$dormancy_rate == 0),]
# results_sums = results_sums[,-3]
results_sums$lower = binconf(results_sums$survived*num_reps,num_reps)[,2]
results_sums$upper = binconf(results_sums$survived*num_reps,num_reps)[,3]
results_sums = results_sums[which(results_sums$pulse_intensity >= 2),]
pulse_plot = ggplot(data = results_sums, aes(x=pulse_intensity,y=survived,color = dormancy_rate))+
  # geom_point(size = 3)+
  # ylim(c(0,1))+
  ylab("Proportion Surviving\n Second Pulse")+
  xlab("First Pulse Intensity")+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, size = .8)+
  guides(color = guide_legend(title = paste("Dormancy Rate (","\u03b4",')',sep = '')),text = element_text(size = 15),hjust = -.5)+
  scale_color_hue(labels = c("0", "0.0001"))+
  scale_color_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  geom_line(size = 1.5)+
  theme_classic()+
  theme(text =element_text(size = 30),legend.position = 'none')
plot(pulse_plot)
res = read.csv("~/Documents/GitHub/publication_model/pulse/kill_intensity 200 2025-11-20 23:09:44.509584 .txt")
res_means = aggregate(res[1],by = res[,c(2,3)],FUN = mean)
res_se = aggregate(res[1],by = res[,c(2,3)],FUN = sd)
res_se = res_se[,3]/sqrt(reps)
colnames(res_means) = c('pulse','dormancy','popsize')
colnames(res) = c('popsize','pulse','dormancy')
res_means$dormancy = as.factor(res_means$dormancy)
res$dormancy = as.factor(res$dormancy)

res_means$lower = res_means$popsize - res_se
res_means$upper = res_means$popsize + res_se
# res_means[which(res_means[,3] == 0),3] = NA
ggplot(res_means,aes(x = pulse,y=popsize,group=dormancy,color = dormancy))+
  geom_line(size = 1.5)+
  geom_point(size = 4)+
  # geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, size = .8)+
  scale_color_manual(labels = c("0","0.0000001"),values = c('#F8766D','#000000'))+
  scale_y_log10(labels = c(expression(10^2),expression(10^6),expression(10^10)),breaks = c(10^2,10^6,10^10))+
  ylab('Population Size\nBefore Second Pulse')+
  xlab('Pulse Intensity')+
  theme_classic()+
  theme(legend.position = 'none',text =element_text(size = 30))
# ggsave('pulse_plot.png',pulse_plot,device = png)# ggsavepulse('pulse_plot.png',pulse_plot,device = png)
