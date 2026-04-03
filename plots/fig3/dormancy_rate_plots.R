library(ggplot2)
library(gridExtra)
# scaleFUN <- function(x) sprintf("%.4f", x)
reps = 200
  results = read.csv("~/Documents/GitHub/publication_model/dormancy_rate/dormancy_rate 200 2026-04-01 15:31:30.054541",header = F)
  results[,3] = results[,3]
results_sum = aggregate(results[,1],list(results[,2]),FUN = mean)
# results_sum[,2] = results_sum[,2]/200
colnames(results_sum) = c('dormancy_rate','survived')
results_sum$lower = binconf(results_sum$survived*reps,reps)[,2]
results_sum$upper = binconf(results_sum$survived*reps,reps)[,3]
dorm_rate_plot = ggplot(results_sum,aes(x= dormancy_rate/10^-7,y = survived))+
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  # xlab(expression(paste('Dormancy Rate (Resuscitation Rate X', 10^-4,')')))+
  xlab('Multiplier Of Both Rates')+
  # geom_vline(xintercept = .1,color = 'green',linetype = 2,lwd = 1.25)+
  ylab('Proportion Surviving')+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.25, size = 1)+
  # scale_x_log10(breaks = c(10^-7,10^-4,10^-1),labels = c(expression(10^-7),expression(10^-4),expression(10^-1)))+
  scale_x_log10(breaks = c(.01,.1,1,10,100,100000),labels = c(expression(10^-2),expression(10^-1),1,expression(10^1),expression(10^2),expression(10^6)))+
  ggtitle('')+
  theme_classic()+
  theme(legend.position="none",text = element_text(size = 30))
plot(dorm_rate_plot)
res = read.csv("~/Documents/GitHub/publication_model/dormancy_rate/dormancy_rate_res 200 2025-12-03 15:26:55.399393",header =F)
res_results_sum = aggregate(res[,c(1)],by = list(res[,(2)]),FUN = mean)
colnames(res_results_sum) = c('dormancy_rate','survived')
res_results_sum$lower = binconf(res_results_sum$survived*reps,reps)[,2]
res_results_sum$upper = binconf(res_results_sum$survived*reps,reps)[,3]
colnames(res_results_sum) = c('dormancy_rate','survived','lower','upper')
res_plot=ggplot(res_results_sum,aes(x = dormancy_rate/280.4095 * 30,y=survived))+
  geom_line(size = 1.5)+
  geom_point(size =3)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.25, size = 1)+
  geom_vline(xintercept = dormancy_rate/280.4095 * 30, linetype = 'dashed')+
  ylab('Proportion surviving')+
  xlab('Resuscitation rate')+
  scale_x_log10(labels = c(expression(10^-11),expression(10^-8),expression(10^-5),expression(10^-2)),breaks = c(10^-11,10^-8,10^-5,10^-2))+
  theme_classic()+
  theme(text = element_text(size = 30))
plot(res_plot)
dor = read.csv("~/Documents/GitHub/publication_model/dormancy_rate/dormancy_rate_dorm 200 2025-12-03 15:26:55.3986",header =F)
dor_results_sum = aggregate(dor[,c(1)],by = list(dor[,(2)]),FUN = mean)
colnames(dor_results_sum) = c('dormancy_rate','survived')
dor_results_sum$lower = binconf(dor_results_sum$survived*reps,reps)[,2]
dor_results_sum$upper = binconf(dor_results_sum$survived*reps,reps)[,3]
colnames(dor_results_sum) = c('dormancy_rate','survived','lower','upper')
dorm_plot =ggplot(dor_results_sum,aes(x = dormancy_rate/280.4095 * 30,y=survived))+
  geom_line(size = 1.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.25, size = 1)+
  geom_point(size = 3)+
  geom_vline(xintercept = resuscitation_rate/280.4095 * 30, linetype = 'dashed')+
  ylab('Proportion surviving')+
  xlab('Dormancy rate')+
  scale_x_log10(labels = c(expression(10^-8),expression(10^-6),expression(10^-4)),breaks = c(10^-8,10^-6,10^-4))+
  theme_classic()+
  theme(text = element_text(size = 30))
plot(dorm_plot)
