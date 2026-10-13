  library(doParallel)
  library(flock)
  setwd("~/Documents/GitHub/publication_model/dormancy_rate/")
  registerDoParallel(cores = 5)
  
  max_time = 3501
  reps = 200
  seed = 1
  tmp_fn = paste('tmp')
  first_pulse = 1
  ##prep seeds
  set.seed(seed)
  fn = paste('dormancy_rate',reps,seed, Sys.time())
  param_fn = paste('dormancy_rate_params',reps,seed, Sys.time())
  dorm_fn = paste('dormancy_rate_dorm',reps, Sys.time())
  dorm_param_fn = paste('dormancy_rate_params',reps, Sys.time())
  res_fn = paste('dormancy_rate_res',reps, Sys.time())
  res_param_fn = paste('dormancy_rate_res_params',reps,seed, Sys.time())
  param_space = 10^seq(-14,-9,.5) 
  param_space = c(param_space,0)
  res_param_space = 10^seq(-13,-.1,.75)
  seeds = sample.int(10^9,reps*length(param_space)*length(res_param_space))
  
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
    foreach(k = 1:length(res_param_space))%do%{
      foreach(i = 1:reps)%dopar%{
        set.seed(seeds[i+ reps*(k-1)+reps*length(res_param_space)*(l-1)])
      test_dorm = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = c(4,20),
      toxin_times = c(1,max_time-100),max_time = max_time,
      init_ts = 0.01,spont_death = spont_death,dormancy_rate = param_space[l],resuscitation_rate = res_param_space[k],d_r = d_r)
      env_pars = c()
      for(u in ls()){
        if(is.numeric(get(u)))
          env_pars = c(env_pars,(paste(u,get(u))))
      }
      env_pars = c(env_pars,globe_pars,seeds[i+ reps*(k-1)+reps*length(res_param_space)*(l-1)])
      env_pars = t(env_pars)
      if(tail(test_dorm$times,1)>(max_time - 1)){
        rep_row = c(1,param_space[l],res_param_space[k],i,tail(test_dorm$times,1),seeds[i+ reps*(k-1)+reps*length(res_param_space)*(l-1)])
        rep_row = t(rep_row)
        locked = flock::lock(fn)
        locked_param = flock::lock(param_fn)
        write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
        write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
        flock::unlock(locked)
        flock::unlock(locked_param)
      } else {
        rep_row = c(0,param_space[l],res_param_space[k],i,tail(test_dorm$times,1),seeds[i+ reps*(k-1)+reps*length(res_param_space)*(l-1)])
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
  }

print(Sys.time()-t)

res = read.csv(fn,header =F)

t= Sys.time()

