library(gridExtra)
library(grid)
library(ggplot2)
library(gtable)
library(cowplot)
set.seed(25)
scaleFUN <- function(x) sprintf("%.0f", x)
popul = vbmax_equibpop_determination(dormancy_rate,resuscitation_rate,spont_death)[[1]]
dorm = simulation(popul,cost = 5, hetero,vbmax,init_carb = c_eq,max_time = (3701), toxin_times = c(1,3401), toxin_strength = c(4,20),
init_ts = 0.083,spont_death = spont_death,dormancy_rate = dormancy_rate, d_r = d_r)
popul[[1]] = popul[[1]] + popul[[2]]
popul[[2]] = popul[[3]]
no_dorm = simulation(popul,cost = 5, hetero,vbmax,init_carb = c_eq,max_time = (3701),toxin_times = c(1,3401), toxin_strength = c(4,20),
init_ts = 0.083,spont_death = spont_death,dormancy_rate = 0,d_r= d_r)
dorm_init = as.data.frame(cbind(rowSums(dorm$pops[,,min(which(dorm$tox_list>3.8))+500]),c(1:15)))
dorm_pre = as.data.frame(cbind(rowSums(dorm$pops[,,min(which(dorm$tox_list>19.8))]),c(1:15)))
dorm_post = as.data.frame(cbind(rowSums(dorm$pops[,,min(which(dorm$tox_list>19.8)+200)]),c(1:15)))
colnames(dorm_init) = c('Number','Tolerance')
colnames(dorm_pre) = c('Number','Tolerance')
colnames(dorm_post) = c('Number','Tolerance')
dorm_init$Tolerance = as.factor(dorm_init$Tolerance)
dorm_pre$Tolerance = as.factor(dorm_pre$Tolerance)
dorm_post$Tolerance = as.factor(dorm_post$Tolerance)
dorm_init$Geno = as.factor(rep('Does Form Seed Bank',length(dorm_init$Number)))
dorm_pre$Geno = as.factor(rep('Does Form Seed Bank',length(dorm_init$Number)))
dorm_post$Geno = as.factor(rep('Does Form Seed Bank',length(dorm_init$Number)))
no_dorm_init = as.data.frame(cbind(rowSums(no_dorm$pops[,,min(which(no_dorm$tox_list>3.8))+500]),c(1:15)))
no_dorm_pre = as.data.frame(cbind(rowSums(no_dorm$pops[,,min(which(no_dorm$tox_list>19.8))]),c(1:15)))
no_dorm_post = as.data.frame(cbind(rowSums(no_dorm$pops[,,length(no_dorm$times)]),c(1:15)))
colnames(no_dorm_init) = c('Number','Tolerance')
colnames(no_dorm_pre) = c('Number','Tolerance')
colnames(no_dorm_post) = c('Number','Tolerance')
no_dorm_pre$Tolerance = as.factor(no_dorm_pre$Tolerance)
no_dorm_post$Tolerance = as.factor(no_dorm_post$Tolerance)
no_dorm_init$Tolerance = as.factor(no_dorm_init$Tolerance)
no_dorm_init$Geno = as.factor(rep('Does Not Form Seed Bank',length(no_dorm_init$Number)))
no_dorm_pre$Geno = as.factor(rep('Does Not Form Seed Bank',length(no_dorm_pre$Number)))
no_dorm_post$Geno = as.factor(rep('Does Not Form Seed Bank',length(no_dorm_post$Number))) 
init = rbind(no_dorm_init,dorm_init)
init[which(init$Number==0),1] = NA
pre = rbind(no_dorm_pre,dorm_pre)
pre[which(pre$Number==0),1] = NA
post = rbind(no_dorm_post,dorm_post)
post[which(post$Number==0),1] = NA
init_plt = ggplot(as.data.frame(init),aes(x=Tolerance,y=Number,color = Geno,fill = Geno))+
  geom_bar(stat = 'identity',show.legend = F,position = "dodge")+
  scale_color_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  scale_fill_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  scale_size_discrete(range = c(3.5,1))+
  # scale_color_brewer(palette = "Dark2")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(10^seq(0,12,2)),labels = c(1,  expression(10^2), expression(10^4),expression(10^6),expression(10^8),expression(10^10),expression(10^12)))+
  ylab("Population Size")+
  scale_x_discrete(breaks = seq(1,15,2))+
  # ggtitle('F.')+
  theme_classic()+
  theme(plot.title = element_text(hjust = -0.2,size = 60),text = element_text(size = 20))
pre_plt = ggplot(as.data.frame(pre),aes(x=Tolerance,y=Number,color = Geno,fill = Geno),alpha = c(.2))+
  geom_bar(stat = 'identity',show.legend = F,position = "dodge")+
  scale_size_discrete(range = c(3.5,1))+
  scale_color_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  scale_fill_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  ylab("Population Size")+
  # scale_color_brewer(palette = "Dark2")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(10^seq(0,12,2)),labels = c(1,  expression(10^2), expression(10^4),expression(10^6),expression(10^8),expression(10^10),expression(10^12)))+
  scale_x_discrete(breaks = seq(1,15,2))+
  # ggtitle('G.')+
  theme_classic()+
  theme(plot.title = element_text(hjust = -0.2,size = 60),text = element_text(size = 20))
legend = ggplot(as.data.frame(post),aes(x=Tolerance,y=Number,color = Geno,group = Geno,),position = "dodge",size = 4)+
  geom_line(stat = 'identity',show.legend = T)+
  guides(color = guide_legend(title='Dormancy\nRate (\u03b4)',override.aes = list(size=c(2.5,1))))+
  scale_color_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  scale_fill_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  scale_size_discrete(range = c(3.5,1))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(10^seq(0,12,2)),labels = c(1,  10^2, 10^4,10^6,10^8,10^10,10^12))+
  scale_x_continuous(breaks = seq(1,15,2))+
  ylab("Population Size")+
  theme_classic()+
  theme(plot.title = element_text(hjust = -0.2,size = 60),text = element_text(size = 27))

post_plt = ggplot(as.data.frame(post),aes(x=Tolerance,y=Number,color = Geno,fill = Geno),alpha = 1)+
  geom_bar(stat = 'identity',show.legend = F,position = "dodge")+
  scale_size_discrete(range = c(3.5,1))+
  scale_color_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  scale_fill_manual(labels = c("0","0.0001"),values = c('#F8766D','#000000'))+
  # scale_color_brewer(palette = "Dark2")+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(1,5),labels = c(1,5))+
  scale_x_discrete(breaks = seq(1,15,2))+
  ylab("Population Size")+
  theme_classic()+
  theme(plot.title = element_text(hjust = -0.2,size = 60),text = element_text(size = 20))

tox = cbind(dorm$tox_list,(dorm$times))
colnames(tox) = c("ToxinConcentration",'Time')
tox = data.frame(tox)
toxins = ggplot(tox,aes(x=Time/280.4095 * 30, y=ToxinConcentration))+
  geom_line(size = 2)+
  # ggtitle('A.')+
  ylab('')+
  xlab('Time (Minutes)')+
  # scale_x_continuous(expand = c(0.01, 0.5))+
  scale_y_continuous(breaks = c())+
  theme_classic()+
  theme(plot.title = element_text(hjust = -0.,size = 30),text = element_text(size = 25))
plot(init_plt)
plot(pre_plt)
plot(post_plt)
legend = get_legend(legend)
ggsave('/Users/ijordan/Documents/GitHub/Phenotypic_seedbank/plots/fig1/post.png',post_plt,height =3.5,width = 4)
ggsave('/Users/ijordan/Documents/GitHub/Phenotypic_seedbank/plots/fig1/pre.png',pre_plt,height =3.5,width =4)
ggsave('/Users/ijordan/Documents/GitHub/Phenotypic_seedbank/plots/fig1/init.png',init_plt,height =3.5,width = 4)
ggsave('/Users/ijordan/Documents/GitHub/Phenotypic_seedbank/plots/fig1/toxins.png',toxins,height = 2,width = 4*3)
ggsave('/Users/ijordan/Documents/GitHub/Phenotypic_seedbank/plots/fig1/lgnd.png',legend,height = 3.5,width = 4)
