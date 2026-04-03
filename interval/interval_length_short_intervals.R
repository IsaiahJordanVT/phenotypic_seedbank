      library(doParallel)
      registerDoParallel(cores = 5)
      reps = 200
      setwd("~/Documents/GitHub/publication_model/interval")
      d_r = 0.01
      tox_pulses = c(4,20)
      ##interval parameters
      first_pulse = 1
      param_space = c()
      param_space = seq(501,4501,1000)
      param_space = c(param_space,seq(251,4251,1000))
      param_space = c(param_space,seq(751,4751,1000))
 
    
      ##prep seeds
      seed = 1
      set.seed(seed)
      seeds = sample.int(10^9,reps*length(param_space)*4)
      fn = paste('interval_short',seed,reps, Sys.time(),'.csv')
      param_fn = paste('interval_short_params',seed,reps, Sys.time(),'.csv')
    
      globe_pars = c()
      for(u in ls()){
        if(is.numeric(get(u)))
          globe_pars = c(globe_pars,(paste(u,get(u))))
      }
      globe_pars = Filter(function(x) !any(grepl("seeds", x)), globe_pars)
      globe_pars = Filter(function(x) !any(grepl("dorm", x)), globe_pars)
      globe_pars = Filter(function(x) !any(grepl("spont", x)), globe_pars)
      globe_pars = Filter(function(x) !any(grepl("param_space", x)), globe_pars)
         t = Sys.time()
        foreach(l = 1:length(param_space))%do%{
          max_time = param_space[l] +100
          foreach(i = 1:reps)%dopar%{
          dormancy_rate = dormancy_rate
          spont_death = spont_death
          set.seed(seeds[(l-1)*reps+i])
          int_sim = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = tox_pulses,
                               toxin_times = c(first_pulse,param_space[l]),max_time = (max_time),
                               init_ts = .083,dormancy_rate = dormancy_rate,spont_death =spont_death,d_r = d_r)
          env_pars = c()
          for(u in ls()){
            if(is.numeric(get(u)))
              env_pars = c(env_pars,(paste(u,get(u))))
          }
          env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i],param_space[i])
          env_pars = t(env_pars)
          if(tail(int_sim$times,1)>(max_time-1)){
            rep_row = c(1,param_space[l],i,tail(int_sim$times,1),dormancy_rate,(max_time),spont_death,seeds[(l-1)*reps+i])
            rep_row = t(rep_row)
            locked = flock::lock(fn)
            locked_param = flock::lock(param_fn)
            write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
            write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
            flock::unlock(locked)
            flock::unlock(locked_param)
          } else {
            rep_row = c(0,param_space[l],i,tail(int_sim$times,1),dormancy_rate,(max_time),spont_death,seeds[(l-1)*reps+i])
            rep_row = t(rep_row)
            locked = flock::lock(fn)
            locked_param = flock::lock(param_fn)
            write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
            write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
            flock::unlock(locked)
            flock::unlock(locked_param)
          }
          }
    
    
    
        foreach(i = 1:reps)%dopar%{
          # for(i in 1:reps){
          dormancy_rate = 0
          spont_death = spont_death
          set.seed(seeds[(l-1)*reps+i])
          int_sim = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = tox_pulses,
                                 toxin_times = c(first_pulse,param_space[l]),max_time = (max_time),
                                 init_ts = .083,dormancy_rate = 0,spont_death =spont_death,d_r = d_r)
          env_pars = c()
          for(u in ls()){
            if(is.numeric(get(u)))
              env_pars = c(env_pars,(paste(u,get(u))))
          }
          env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i],param_space[i])
          env_pars = t(env_pars)
          if(tail(int_sim$times,1)>(max_time-1)){
            rep_row = c(1,param_space[l],i,tail(int_sim$times,1),dormancy_rate,(max_time),spont_death,seeds[(l-1)*reps+i])
            rep_row = t(rep_row)
            locked = flock::lock(fn)
            locked_param = flock::lock(param_fn)
            write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
            write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
            flock::unlock(locked)
            flock::unlock(locked_param)
          } else {
            rep_row = c(0,param_space[l],i,tail(int_sim$times,1),dormancy_rate,(max_time),spont_death,seeds[(l-1)*reps+i])
            rep_row = t(rep_row)
            locked = flock::lock(fn)
            locked_param = flock::lock(param_fn)
            write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
            write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
            flock::unlock(locked)
            flock::unlock(locked_param)
          }
        }
    
    
    
            foreach(i = 1:reps)%dopar%{
          dormancy_rate = dormancy_rate
          spont_death = 0.001
          # popul = vbmax_equibpop_determination(dormancy_rate,dormancy_rate*50,spont_death)[[1]]
          set.seed(seeds[(l-1)*reps+i])
          int_sim = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = tox_pulses,
                               toxin_times = c(first_pulse,param_space[l]),max_time = (max_time),
                               init_ts = .083,dormancy_rate = dormancy_rate,spont_death =0.001,d_r = d_r)
          env_pars = c()
          for(u in ls()){
            if(is.numeric(get(u)))
              env_pars = c(env_pars,(paste(u,get(u))))
          }
          env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i])
          env_pars = t(env_pars)
          if(tail(int_sim$times,1)>(max_time-1)){
            rep_row = c(1,param_space[l],i,tail(int_sim$times,1),dormancy_rate,(max_time),spont_death,seeds[(l-1)*reps+i])
            rep_row = t(rep_row)
            locked = flock::lock(fn)
            locked_param = flock::lock(param_fn)
            write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
            write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
            flock::unlock(locked)
            flock::unlock(locked_param)
          } else {
            rep_row = c(0,param_space[l],i,tail(int_sim$times,1),dormancy_rate,(max_time),spont_death,seeds[(l-1)*reps+i])
            rep_row = t(rep_row)
            locked = flock::lock(fn)
            locked_param = flock::lock(param_fn)
            write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
            write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
            flock::unlock(locked)
            flock::unlock(locked_param)
          }
            }
    
    
    
    
            foreach(i = 1:reps)%dopar%{
          dormancy_rate = 0
          spont_death = 0.001
          set.seed(seeds[(l-1)*reps+i])
          int_sim = simulation(popul,cost = cost, hetero,vbmax,init_carb = c_eq,toxin_strength = tox_pulses,
                               toxin_times = c(first_pulse,param_space[l]),max_time = (max_time),
                               init_ts = .083,dormancy_rate = 0,spont_death =0.001,d_r = d_r)
          env_pars = c()
          for(u in ls()){
            if(is.numeric(get(u)))
              env_pars = c(env_pars,(paste(u,get(u))))
          }
          env_pars = c(env_pars,globe_pars,seeds[(l-1)*reps+i])
          env_pars = t(env_pars)
          if(tail(int_sim$times,1)>(max_time-1)){
            rep_row = c(1,param_space[l],i,tail(int_sim$times,1),dormancy_rate,(max_time),spont_death,seeds[(l-1)*reps+i])
            rep_row = t(rep_row)
            locked = flock::lock(fn)
            locked_param = flock::lock(param_fn)
            write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
            write.table(rep_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
            flock::unlock(locked)
            flock::unlock(locked_param)
          } else {
            rep_row = c(0,param_space[l],i,tail(int_sim$times,1),dormancy_rate,(max_time),spont_death,seeds[(l-1)*reps+i])
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
        print(Sys.time()-t)
       