    library(doParallel)
    library(flock)
    # setwd("Users/isaiahjordan/Documents/Phenotypic_dormancy_output/pulse_intensity")
    setwd('~/Documents/GitHub/publication_model/results/')
    registerDoParallel(cores = 5)
    # param_space = seq()
    param_space = seq(2,6,.5)
    d_r = .01
    max_time = 3501
    first_pulse = 1
    reps = 200
    fn = paste('pulse_intensity',reps, Sys.time(),'.txt')
    param_fn = paste('pulse_intesity_params',reps, Sys.time(),'.txt')
    first_pulse = 1
    set.seed(1)
    seeds = sample.int(10^9,reps*length(param_space))
    globe_pars = c()
    t = Sys.time()
    for(u in ls()){
      if(is.numeric(get(u)))
        globe_pars = c(globe_pars,(paste(u,get(u))))
    }
    globe_pars = Filter (function(x) !any(grepl("seeds", x)), globe_pars)
    globe_pars = Filter(function(x) !any(grepl("dorm", x)), globe_pars)
    globe_pars = Filter(function(x) !any(grepl("spont", x)), globe_pars)
    globe_pars = Filter(function(x) !any(grepl("param_space", x)), globe_pars)
    dormancy_rate = 10^-7
  resuscitation_rate = dormancy_rate*10^4
    
  foreach(l = 1:length(param_space)) %do% {
    ###dorm case
    
      foreach(i = 1:reps)%dopar%{
        dormancy_rate = dormancy_rate
        resuscitation_rate = resuscitation_rate
        set.seed(seeds [(l-1)*reps+i])
        pulse = simulation(popul,cost = cost,hetero,vbmax,init_carb = c_eq,toxin_strength = c(param_space[l],20),
                toxin_times = c(first_pulse,(max_time - 100)),max_time = ((max_time)),
                init_ts = .083,spont_death = spont_death,
                dormancy_rate = dormancy_rate,resuscitation_rate = resuscitation_rate,d_r = d_r)
        env_pars = c()
        for(u in ls()){
          if(is.numeric(get(u)))
            env_pars = c(env_pars,(paste(u,get(u))))
        }
        env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i])
        env_pars = t(env_pars)
        if(tail(pulse$times,1)>(max_time-1)){
          rep_row = c(1,param_space[l],i,tail(pulse$times,1),dormancy_rate,
                      seeds[(l-1)*reps+i])
          rep_row = t(rep_row)
          locked = flock::lock(fn)
          locked = flock::lock(param_fn)
          write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
          write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
          flock::unlock(locked)
        } else {
          rep_row = c(0,param_space[l],i,tail(pulse$times,1),dormancy_rate,
                      seeds[(l-1)*reps+i])
          rep_row = t(rep_row)
          locked = flock::lock(fn)
          locked = flock::lock(param_fn)
          write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
          write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
          flock::unlock(locked)
        }
        }
        
      ####no dorm case  
      foreach(i = 1:reps)%dopar%{
        dormancy_rate = 0
        resuscitation_rate = 0
        set.seed(seeds[(l-1)*reps+i])
        pulse = simulation(popul,cost = cost,hetero,vbmax,init_carb = c_eq,toxin_strength = c(param_space[l],20),
                           toxin_times = c(first_pulse,(max_time - 100)),max_time = ((max_time)),
                           init_ts = .083,spont_death = spont_death,
                           dormancy_rate = dormancy_rate,resuscitation_rate = resuscitation_rate,d_r = d_r)
        env_pars = c()
        for(u in ls()){
          if(is.numeric(get(u)))
            env_pars = c(env_pars,(paste(u,get(u))))
        }
        env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i])
        env_pars = t(env_pars)
        if(tail(pulse$times,1)>(max_time-1)){
          rep_row = c(1,param_space[l],i,tail(pulse$times,1),dormancy_rate,
                      seeds[(l-1)*reps+i])
          rep_row = t(rep_row)
          locked = flock::lock(fn)
          locked = flock::lock(param_fn)
          write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
          write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
          flock::unlock(locked)
        } else {
          rep_row = c(0,param_space[l],i,tail(pulse$times,1),dormancy_rate,
                      seeds[(l-1)*reps+i])
          rep_row = t(rep_row)
          locked = flock::lock(fn)
          param_locked = flock::lock(param_fn)
          write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
          write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
          flock::unlock(locked)
          flock::unlock(param_locked)
        }
      }
  }
  tms = .01
  time_to_decay = c()
  for(i in param_space){
    decays = c()
    tox = i
    time = 0
    while(tox > .01){
      tox = tox - d_r *tox *tms
      time = time + tms
      
    }
    time_to_decay = c(time_to_decay,time)
  kills = data.frame()
  }
  kill_fn = paste('kill_intensity',reps, Sys.time(),'.txt')
  foreach(l = 1:length(param_space))%do%{
  foreach(i = 1:reps)%dopar%{
    killing = simulation(popul,cost = cost,hetero,vbmax,init_carb = c_eq,toxin_strength = c(param_space[l],0),
                       toxin_times = c(first_pulse,(max_time)),max_time = ((max_time-1)),
                       init_ts = .083,spont_death = spont_death,
                       dormancy_rate = dormancy_rate,resuscitation_rate = resuscitation_rate,d_r = d_r)
    rep_row = c(killing$pop_list[length(killing$pop_list)],param_space[l],dormancy_rate)
    rep_row = t(rep_row)
    locked = flock::lock(fn)
   write.table(rep_row,file = kill_fn,append = T,col.names = F,row.names = F,sep = ',')
   flock::unlock(locked)
  }
  }
  foreach(l = 1:length(param_space))%do%{
    foreach(i = 1:reps)%dopar%{
      killing = simulation(popul,cost = cost,hetero,vbmax,init_carb = c_eq,toxin_strength = c(param_space[l],0),
                           toxin_times = c(first_pulse,(max_time)),max_time = ((max_time-1)),
                           init_ts = .083,spont_death = spont_death,
                           dormancy_rate = 0,resuscitation_rate = resuscitation_rate,d_r = d_r)
      rep_row = c(killing$pop_list[length(killing$pop_list)],param_space[l],0)
      rep_row = t(rep_row)
       locked = flock::lock(fn)
      write.table(rep_row,file = kill_fn,append = T,col.names = F,row.names = F,sep = ',')
      flock::unlock(locked)
    }
  }

print(Sys.time()-t)
      