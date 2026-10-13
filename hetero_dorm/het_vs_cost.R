  library(doParallel)
  library(flock)
  reps = 500
  seed  = 1
  setwd("~/Documents/GitHub/publication_model/hetero_cost/")
  registerDoParallel(cores = 5)
  
  ##set filenames
  fn = paste('het_vs_dorm',seed,reps, Sys.time())
  param_fn = paste('het_dorm_params',seed,reps, Sys.time(),'.csv')
  
  ##other parameters
  exp_run_time = Sys.time()
  first_pulse = 1
  max_time = 3501
  # costs = c(1,10)
  costs = seq(.001,15,.25)
  # hetero_rates = c(9.5e-4,2e-3)
  hetero_rates = seq(0,.0015,0.00015625)

  
  ##prep seeds
  set.seed(seed) ##different seeds used to distribute task, seeds used:1-4, with 50 reps each
  seeds = sample.int(10^9,reps*length(costs)*length(hetero_rates))
  globe_pars = c()
  for(u in ls()){
    if(is.numeric(get(u)))
      globe_pars = c(globe_pars,(paste(u,get(u))))
  }
  globe_pars = Filter(function(x) !any(grepl("seeds", x)), globe_pars)
  globe_pars = Filter(function(x) !any(grepl("cost", x)), globe_pars)
  globe_pars = Filter(function(x) !any(grepl("hetero", x)), globe_pars)
  
  first_time = Sys.time()
  ##ONLY HETEROGENEITY CASE##
  het_survive = c()
  popul = vbmax_equibpop_determination(dormancy_rate,resuscitation_rate,spont_death)[[1]]
  popul[[1]] = popul[[1]] + popul[[2]]
  popul[[2]] = popul[[3]]
  foreach(j = 1:length(costs)) %do% {
  for (k in 1:length(hetero_rates)) {
    is_alive = 0
    foreach (l = 1:reps) %dopar% {
      print(paste(k, ',', l, ',', j))
      set.seed(seeds[((j-1)*reps*length(hetero_rates)+(k-1)*reps)+l])
      pops = simulation(popul,cost = costs[j], hetero_rates[k],vbmax,init_carb = c_eq,toxin_strength = c(4,20),
                        toxin_times = c(first_pulse,(max_time - 100)),max_time = (max_time),
                        init_ts = 0.01,spont_death = spont_death,dormancy_rate = 0,d_r=0.01)
      env_pars = c()
      for(u in ls()){
        if(is.numeric(get(u)))
          env_pars = c(env_pars,(paste(u,get(u))))
      }
      env_pars = c(env_pars,globe_pars,seeds[((j-1)*reps*length(hetero_rates)+(k-1)*reps)+l])
      env_pars = t(env_pars)
      num_cells = sum(pops[[1]][, , dim(pops$pops)[3]])
      if (num_cells > 0 & dim(pops$pops)[3] > (max_time - 1)) {
        het_sur_row = c(1, "Heterogeneity Only", hetero_rates[k], costs[j],l,seeds[((j-1)*reps*length(hetero_rates)+(k-1)*reps)+l])
        het_sur_row = t(as.data.frame(het_sur_row))
        locked = flock::lock(fn)
        locked_param = flock::lock(param_fn)
        write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
        write.table(het_sur_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
        flock::unlock(locked)
        flock::unlock(locked_param)
      } else {
        het_sur_row = c(0, "Heterogeneity Only", hetero_rates[k], costs[j],l,seeds[((j-1)*reps*length(hetero_rates)+(k-1)*reps)+l])
        het_sur_row = t(as.data.frame(het_sur_row))
        locked = flock::lock(fn)
        locked_param = flock::lock(param_fn)
        write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
        write.table(het_sur_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
        flock::unlock(locked)
        flock::unlock(locked_param)
      }
    }
  }
  }
  last_time = Sys.time()
  print(last_time - first_time)

  ##ONLY DORMANCY CASE##
  popul = vbmax_equibpop_determination(dormancy_rate,resuscitation_rate,spont_death)[[1]]
  dorm_survive = c()
  hetero = 0
  is_alive = 0
  foreach(j = 1:length(costs)) %do% {
  foreach (l = 1:reps) %dopar% {
    print(paste(l, ',', j))
    set.seed(seeds[((j-1)*reps*length(hetero_rates)+l)])
    pops = simulation(popul,cost = costs[j], 0,vbmax,init_carb = c_eq,toxin_strength = c(4,20),
                      toxin_times = c(first_pulse,(max_time - 100)),max_time = max_time,
                      init_ts = 0.01,spont_death = spont_death,dormancy_rate = dormancy_rate,resuscitation_rate = resuscitation_rate,d_r = 0.01)
    env_pars = c()
    for(u in ls()){
      if(is.numeric(get(u)))
        env_pars = c(env_pars,(paste(u,get(u))))
    }
    env_pars = c(env_pars,globe_pars,seeds[((j-1)*reps*length(hetero_rates)+l)])
    env_pars = t(env_pars)
    num_cells = sum(pops[[1]][, , dim(pops$pops)[3]])
    if (num_cells > 0 & dim(pops$pops)[3] > (max_time - 1)) {
      dorm_sur_row = c(1, "Dormancy Only", 0, costs[j],l,seeds[((j-1)*reps*length(hetero_rates)+l)])
      dorm_sur_row = t(as.data.frame(dorm_sur_row))
      locked = flock::lock(fn)
      locked_param = flock::lock(param_fn)
      write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
      write.table(dorm_sur_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
      flock::unlock(locked)
      flock::unlock(locked_param)
    } else {
      dorm_sur_row = c(0, "Dormancy Only", 0, costs[j],l,seeds[((j-1)*reps*length(hetero_rates)+l)])
      dorm_sur_row = t(as.data.frame(dorm_sur_row))
      locked = flock::lock(fn)
      locked_param = flock::lock(param_fn)
      write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
      write.table(dorm_sur_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
      flock::unlock(locked)
      flock::unlock(locked_param)
    }
  }
  }

  ## DORMANCY AND HETEROGENEITY CASE##
  dorm_het_survive = c()
  is_alive = 0
  # first_time = Sys.time()
  foreach(j = 1:length(costs)) %do% {
  for (k in 1:length(hetero_rates)) {
    is_alive = 0
    foreach (l = 1:reps) %dopar%{
      print(paste(k, ',', l, ',', j))
      set.seed(seeds[((j-1)*reps*length(hetero_rates)+(k-1)*reps)+l])
      pops = simulation(popul,cost = costs[j], hetero_rates[k],vbmax,init_carb = c_eq,toxin_strength = c(4,20),
                        toxin_times = c(first_pulse,(max_time-100)),max_time = (max_time),
                        init_ts = 0.01,spont_death = spont_death,dormancy_rate = dormancy_rate,d_r= .01)
      env_pars = c()
      for(u in ls()){
        if(is.numeric(get(u)))
          env_pars = c(env_pars,(paste(u,get(u))))
      }
      env_pars = c(env_pars,globe_pars,seeds[((j-1)*reps*length(hetero_rates)+(k-1)*reps)+l])
      env_pars = t(env_pars)
      num_cells = sum(pops[[1]][, , dim(pops$pops)[3]])
      if (num_cells > 0 & dim(pops$pops)[3] > (max_time - 1)) {
        dorm_het_sur_row = c(1, "Dormancy and Heterogeneity", hetero_rates[k], costs[j],l,seeds[((j-1)*reps*length(hetero_rates)+(k-1)*reps)+l])
        dorm_het_sur_row = t(as.data.frame(dorm_het_sur_row))
        locked = flock::lock(fn)
        locked_param = flock::lock(param_fn)
        write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
        write.table(dorm_het_sur_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
        flock::unlock(locked)
        flock::unlock(locked_param)
      } else {
        dorm_het_sur_row = c(0, "Dormancy and Heterogeneity", hetero_rates[k], costs[j],l,seeds[((j-1)*reps*length(hetero_rates)+(k-1)*reps)+l])
        dorm_het_sur_row = t(as.data.frame(dorm_het_sur_row))
        locked = flock::lock(fn)
        locked_param = flock::lock(param_fn)
        write.table(env_pars,file = param_fn,append = T,col.names = F,row.names = F,sep = ',')
        write.table(dorm_het_sur_row,file = fn,append = T,col.names = F,row.names = F,sep = ',')
        flock::unlock(locked)
        flock::unlock(locked_param)
      }
    }

  }
  }
  time_to_complete = Sys.time() - first_time
  print(time_to_complete)