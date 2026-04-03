  library(doParallel)
  library(flock)
  setwd("~/Documents/GitHub/publication_model/dormancy_rate/")
  registerDoParallel(cores = 5)
  
  max_time = 3501
  reps = 200
  fn = paste('dormancy_rate',reps, Sys.time())
  param_fn = paste('dormancy_rate_params',reps, Sys.time())
  dorm_fn = paste('dormancy_rate_dorm',reps, Sys.time())
  dorm_param_fn = paste('dormancy_rate_params',reps, Sys.time())
  res_fn = paste('dormancy_rate_res',reps, Sys.time())
  res_param_fn = paste('dormancy_rate_res_params',reps, Sys.time())
  tmp_fn = paste('tmp')
  first_pulse = 1
  ##prep seeds
  set.seed(1)
  param_space = 10^seq(-9,-1,.5)
  res_param_space = 10^seq(-13,-.1,1)
  dorm_param_space =  10^seq(-9,-3,.5)
  dorm_seeds = sample.int(10^9,reps*length(dorm_param_space)*4)
  res_seeds = sample.int(10^9,reps*length(res_param_space)*4)
  seeds = sample.int(10^9,reps*length(param_space)*4)
  

  
  
  globe_pars = c()
  for(u in ls()){
    if(is.numeric(get(u)))
      globe_pars = c(globe_pars,(paste(u,get(u))))
  }
  globe_pars = Filter(function(x) !any(grepl("seeds", x)), globe_pars)
  globe_pars = Filter(function(x) !any(grepl("dorm", x)), globe_pars)
  globe_pars = Filter(function(x) !any(grepl("spont", x)), globe_pars)
  
  t = Sys.time()
foreach(l = 1:length(param_space)) %do%{
    foreach(i = 1:reps)%dopar%{
      set.seed(seeds[(l-1)*reps+i])
      test_dorm = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = c(4,20),
      toxin_times = c(1,max_time-100),max_time = max_time,
      init_ts = 0.01,spont_death = spont_death,dormancy_rate = param_space[l],resuscitation_rate = 10^-4*param_space[l],d_r = d_r)
      env_pars = c()
      for(u in ls()){
        if(is.numeric(get(u)))
          env_pars = c(env_pars,(paste(u,get(u))))
      }
      env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i])
      env_pars = t(env_pars)
      if(tail(test_dorm$times,1)>(max_time - 1)){
        rep_row = c(1,param_space[l],i,tail(test_dorm$times,1),seeds[(l-1)*reps+i])
        rep_row = t(rep_row)
        locked = flock::lock(fn)
        locked_param = flock::lock(param_fn)
        write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
        write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
        flock::unlock(locked)
        flock::unlock(locked_param)
      } else {
        rep_row = c(0,param_space[l],i,tail(test_dorm$times,1),seeds[(l-1)*reps+i])
        rep_row = t(rep_row)
        locked = flock::lock(fn)
        locked_param = flock::lock(param_fn)
        write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
        write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
        flock::unlock(locked)
        flock::unlock(locked_param)
      }
    }
}
foreach(l = 1:length(dorm_param_space)) %do%{
  foreach(i = 1:reps)%dopar%{
    set.seed(dorm_seeds[(l-1)*reps+i])
    test_dorm = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = c(4,20),
                           toxin_times = c(1,max_time-100),max_time = max_time,
                           init_ts = 0.01,spont_death = spont_death,dormancy_rate = dorm_param_space[l],resuscitation_rate = resuscitation_rate,d_r = d_r)
    env_pars = c()
    for(u in ls()){
      if(is.numeric(get(u)))
        env_pars = c(env_pars,(paste(u,get(u))))
    }
    env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i])
    env_pars = t(env_pars)
    if(tail(test_dorm$times,1)>(max_time - 1)){
      rep_row = c(1,dorm_param_space[l],i,tail(test_dorm$times,1),dorm_seeds[(l-1)*reps+i])
      rep_row = t(rep_row)
      locked = flock::lock(dorm_fn)
      locked_param = flock::lock(dorm_param_fn)
      write.table(env_pars,file = dorm_param_fn,append = T,col.names = F,row.names = F,sep = ',')
      write.table(rep_row,file = dorm_fn,append = T,col.names = F,row.names = F,sep = ',')
      flock::unlock(locked)
      flock::unlock(locked_param)
    } else {
      rep_row = c(0,dorm_param_space[l],i,tail(test_dorm$times,1),dorm_seeds[(l-1)*reps+i])
      rep_row = t(rep_row)
      locked = flock::lock(dorm_fn)
      locked_param = flock::lock(dorm_param_fn)
      write.table(env_pars,file = dorm_param_fn,append = T,col.names = F,row.names = F,sep = ',')
      write.table(rep_row,file = dorm_fn,append = T,col.names = F,row.names = F,sep = ',')
      flock::unlock(locked)
      flock::unlock(locked_param)
    }
  }
}
foreach(l = 1:length(res_param_space)) %do%{
  foreach(i = 1:reps)%dopar%{
    set.seed(res_seeds[(l-1)*reps+i])
    test_dorm = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = c(4,20),
                           toxin_times = c(1,max_time-100),max_time = max_time,
                           init_ts = 0.01,spont_death = spont_death,dormancy_rate = dormancy_rate,resuscitation_rate = res_param_space[l],d_r = d_r)
    env_pars = c()
    for(u in ls()){
      if(is.numeric(get(u)))
        env_pars = c(env_pars,(paste(u,get(u))))
    }
    env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i])
    env_pars = t(env_pars)
    if(tail(test_dorm$times,1)>(max_time - 1)){
      rep_row = c(1,res_param_space[l],i,tail(test_dorm$times,1),res_seeds[(l-1)*reps+i])
      rep_row = t(rep_row)
      locked = flock::lock(res_fn)
      locked_param = flock::lock(res_param_fn)
      write.table(env_pars,file = res_param_fn,append = T,col.names = F,row.names = F,sep = ',')
      write.table(rep_row,file = res_fn,append = T,col.names = F,row.names = F,sep = ',')
      flock::unlock(locked)
      flock::unlock(locked_param)
    } else {
      rep_row = c(0,res_param_space[l],i,tail(test_dorm$times,1),res_seeds[(l-1)*reps+i])
      rep_row = t(rep_row)
      locked = flock::lock(res_fn)
      locked_param = flock::lock(res_param_fn)
      write.table(env_pars,file = res_param_fn,append = T,col.names = F,row.names = F,sep = ',')
      write.table(rep_row,file = res_fn,append = T,col.names = F,row.names = F,sep = ',')
      flock::unlock(locked)
      flock::unlock(locked_param)
    }
  }
}

dorm_pop_fracs = c()
res_pop_fracs = c()
param_space = 10^seq(-25,-.1,.5)
res_times = 1/ param_space
for(i in param_space){
  popul = vbmax_equibpop_determination(i,resuscitation_rate,spont_death)[[1]]
  dorm_frac = sum(popul[[2]])/sum(popul[[1]]+popul[[2]])
  popul = vbmax_equibpop_determination(dormancy_rate,i,spont_death)[[1]]
  res_frac = sum(popul[[2]])/sum(popul[[1]]+popul[[2]])
  dorm_pop_fracs = c(dorm_pop_fracs,dorm_frac)
  res_pop_fracs = c(res_pop_fracs,res_frac)
}
plot(y=dorm_pop_fracs,x= (param_space),type = 'l')
abline(v = (10^-4))
plot(y=res_pop_fracs,x= (param_space),type = 'l')
# abline(v = (10^-4))
print(Sys.time()-t)

res = read.csv("~/Documents/GitHub/publication_model/dormancy_rate/dormancy_rate_res 10 2025-11-26 13:00:59.590548",header =F)
results_sum = aggregate(res[,c(1)],by = list(res[,(2)]),FUN = mean)
ggplot(results_sum,aes(x = Group.1,y=x))+
  geom_line()+
  geom_point()+
  ylab('Proportion surviving')+
  xlab('Resuscitation rate')+
  scale_x_log10()+
  theme_classic()+
  theme(text = element_text(size = 30))
dor = read.csv("~/Documents/GitHub/publication_model/dormancy_rate/dormancy_rate_dorm 10 2025-11-26 13:00:59.585746",header =F)
results_sum = aggregate(dor[,c(1)],by = list(dor[,(2)]),FUN = mean)
ggplot(results_sum,aes(x = Group.1,y=x))+
  geom_line()+
  geom_point()+
  ylab('Proportion surviving')+
  xlab('Dormancy rate')+
  scale_x_log10()+
  theme_classic()+
  theme(text = element_text(size = 30))
# test_dorm = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = c(4.5,20),
#                        toxin_times = c(1,max_time-100),max_time = max_time,
#                        init_ts = 0.01,spont_death = spont_death,dormancy_rate = dormancy_rate,resuscitation_rate = 10e-1*10^-4,d_r = d_r)

