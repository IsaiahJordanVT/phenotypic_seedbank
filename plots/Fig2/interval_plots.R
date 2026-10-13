  library(ggplot2)
  library(cowplot)
  library(Hmisc)
  results = read.csv("~/Documents/GitHub/publication_model/interval/interval 1000 2026-08-04 09:02:49.553642 .csv",header = F)
  # results = results[,-1]
  # results =  rbind(results,results2,results3,results4,short_results1,short_results2,short_results3,short_results4)#,intermed_results6,intermed_results5,intermed_results4,intermed_results3,intermed_results2,intermed_results1)
  colnames(results)= c('survived','interval','rep','final_time','genotype','max_time',
                       'death_rate','seed')
  reps = 1000
  colnames(results)= c('survived','interval','rep','final_time','genotype','max_time',
                       'death_rate','final_size')
  results = results[-which(results$interval==1), ]
  
  results_sum = aggregate(results[,c(1)],results[,c(2,5,7)],FUN = mean)
  results_size = aggregate(results[,c(8)],results[,c(2,5,7)],FUN = mean)
  results_sum = results_sum[-which(results_sum$interval == 251),]
  colnames(results_sum) = c('interval','genotype','death_rate','survival_rate')
  results_num = aggregate(results[,1],results[,c(2,6,3)],FUN = NROW)
  colnames(results_size) = c('interval','genotype','death_rate','population_size')
  colnames(results_num) = c('interval','genotype','death_rate','num_reps')
  results_size[,c(2)] = factor(results_size[,c(2)])
  results_size[,c(3)] = factor(results_size[,c(3)])
  results_size[,1]=results_size[,1]
  lwd=c(2,1)
  sizes = ggplot(results_size,aes(y=population_size,x=interval/280.4095 * 30,linetype = death_rate,color = genotype,size=genotype))+
    geom_line(aes())+
    scale_color_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
    ylab("Population Size At\nSecond Pulse")+
    xlab("Time Between Pulses (Minutes)")+
    # ggtitle("B.")+
    geom_line(aes(y = 10^12),color = 'blue',linetype = 'dotted')+
    scale_y_log10(labels = c(expression(10^8,10^8.5,10^9,10^9.5,10^10,10^10.5,10^11,10^11.5,10^12,10^12.5)),breaks = c(10^seq(8,12.5,.5)))+
    scale_size_manual(values = lwd)+
    theme_classic()+
    theme(text=element_text(size = 20))+
    theme(axis.text.y = element_text(angle = 0),plot.title = element_text(hjust = -0.07))+
  
    guides(fill = "none", color = "none", linetype = "none", shape = "none",size='none')
  results_sum$death_rate = factor(results_sum$death_rate)
  results_sum$genotype = factor(results_sum$genotype)
  results_sum$lower = binconf(results_sum$survival_rate*reps,reps)[,2]
  results_sum$upper = binconf(results_sum$survival_rate*reps,reps)[,3]
  results_sum = results_sum[-which(results_sum$death_rate == 1e-4),] ##remove intermediate death rate replicates, not sufficiently different to be visual clear
  survival = ggplot(data = results_sum, aes(x=interval/280.4095 * 30,y=survival_rate,linetype = death_rate,color = genotype))+
    geom_line(size = 1.5)+
    scale_color_manual(values = c('#F8766D','#000000'))+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=15, size = 1.1,linetype = 1)+
    geom_vline(xintercept = 3401/280.4095 * 30)+
    scale_shape(labels = c('0.00005','0.0001'))+
    # geom_point(aes(shape = death_rate),size = 5)+
    xlab("Minutes")+
    ylab('Proportion Surviving')+
    # ggtitle("A.")+
    theme_classic()+
    theme(text=element_text(size = 25),plot.title = element_text(hjust = -0.10))+
    guides(color="none",shape ='none',linetype = 'none')
  plot(survival)
  # legend = ggplot(data = results_sum, aes(x=interval,y=survival_rate,linetype = death_rate,color = genotype))+
  #   geom_errorbar(aes(ymin=lower, ymax=upper), width=80, size = 1.1)+
  #   scale_color_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  #   # geom_point(aes(shape = death_rate),size = 5)+
  #   xlab("Time Between Pulses (Minutes)")+
  #   ylab('Proportion Surviving')+
  #   geom_line(size = 4)+
  #   labs(color = 'Dormancy\nRate (\u03b4)',linetype = expression('Spontaneous\nDeath Rate (D'[c]*')'))+
  #   stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.2) +
  #   theme_classic()+
  #   theme(text=element_text(size = 20))
  dft = data.frame(i = as.integer(c(1,2,3,4)))
  legend2 = ggplot(data = dft, aes(x=i,y=i,linetype = (c('a','b','c','d')),color = (c('a','b','c','d'))))+
    geom_line(size = 1)+
    scale_color_manual(labels = c(paste("\u03b4 = 0, ", expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.000005",'\n',sep = ""),
                                    paste("\u03b4 = 0.0001, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.000005",'\n',sep = ""),
                                    paste("\u03b4 = 0, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.001",'\n',sep = ""),
                                    paste("\u03b4 = 0.0001, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.001",'\n',sep = "")),
    values = c('#F8766D','#000000','#F8766D','#000000'))+
    scale_linetype_manual(labels = c(paste("\u03b4 = 0, ", expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.000005",'\n',sep = ""),
                                     paste("\u03b4 = 0.0001, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.000005",'\n',sep = ""),
                                     paste("\u03b4 = 0, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.001",'\n',sep = ""),
                                     paste("\u03b4 = 0.0001, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.001",'\n',sep = "")),
    values = c(1,2,2,2)
    )+
    guides(linetype = guide_legend(override.aes = list(linetype = c(1,1,2,2))),
           color = guide_legend(override.aes = list(c('#F8766D','#000000','#F8766D','#000000'))),
           labels = guide_legend(override.aes
                                 = list(c(paste("\u03b4 = 0, ", expression('Spontaneous\nDeath Rate '),'(','(',expression(D[n]),')',')'," = 0.000005",'\n',sep = ""),
                                          paste("\u03b4 = 0.0001, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.000005",'\n',sep = ""),
                                          paste("\u03b4 = 0, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.001",'\n',sep = ""),
                                          paste("\u03b4 = 0.0001, ",expression('Spontaneous\nDeath Rate '),'(',expression(D[n]),')'," = 0.001",'\n',sep = "")))),
           title = guide_legend(''))+
    # xlab("Time Between Pulses (Minutes)")+
    # ylab('Proportion Surviving')+
    # scale_size_manual(values = c(.25,.25),guide= 'legend')+
    labs('')+
    # labs(color = 'Dormancy\nRate (\u03b4)',linetype = expression('Spontaneous\nDeath Rate (D'[c]*')'))+
    #guides(color = guide_legend(expression('Spontaneous\nDeath Rate (D'[c]*')')),linetype = guide_legend(expression('Spontaneous\nDeath Rate (D'[c]*')')))+
    theme_classic()+
    theme(text=element_text(size = 20))
  
  plot(legend2)
  # fig2 = grid.arrange(survival,sizes,get_legend(legend), layout_matrix = rbind(c(1,1,1,1,3),c(1,1,1,1,3),c(1,1,1,1,3),c(2,2,2,2,NA),c(2,2,2,2,NA)))
  # plot(fig2)
  ggsave('survival.png',survival,height =4.5,width = 7.85)
  # ggsave('sizes.png',sizes,height =3.5,width = 5.25)
  # ggsave('legend.png',legend,height =3.5,width = 5.25)
  # ggsave2('fig2',fig2)
  plot(survival)
