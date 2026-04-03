 # setwd("~/Documents/methods and code/")
  ##life history traits - the number of stages of growth before division and of phenotypes  
  pheno_bins = 15
  stages = 3


  ##rate at which phenotype spontaneously declines due to normal protein degradation
  active_pheno_loss = 0.00

  ##rates at which individuals become dormant and are resuscitated from dormancy
  dormancy_rate = 10^-7
  resuscitation_rate = dormancy_rate*10^4

  ##rate at which spontaneous death occurs
  spont_death = 2.777778e-06
  
  ##amount of random expression of resistance genes
  hetero = 0 ##
  
  ##the cost of expressing the phenotype, ie how much longer you expect to stay in the penalty box than you would expect to wait to grow outside of it
  cost = 5 ##
  allowed_error = 0.005
  
  ##Michealis - Menten constants for death caused by toxin. hill_coef gives steepness of killing curve
  hill_coef = 15 ##
  km_death = 1.7 ##
  vmax_death = 4
  
  ##M-M constants for induction of phenotype caused by toxicity 
  km_toxin = .01 ##
  vmax_toxin = 3 ## 
  
  ##default timestep (0.1)
  tmst = 0.1
  #default end time of simulation
  max_time = 1250
  
  ##km for birth rate
  km_carb = 10^-5 ## smaller = higher growth at lower carbon
  
  ##decay rate of the toxin
  d_r = .01
  
  ##desired carbon equilibrium
  c_eq =.00000001
  ##rate at which carbon is added to the system
  source = c_eq##
  
  ##population size at equilibrium
  N = 10^12
  ##initial v max for birth rate - will change

  vbmax_equibpop_determination = function(dormancy_rate,resuscitation_rate,spont_death){
  ##initialize variable to numerically determine the value of vbmax at equilibrium under the specified parameters
  vbmax = 1
  eigen_diff = 1
  old_eigen_diff = 1
  step = 1e-1
  leslie_rows = stages * 2
  
  while (
    abs(eigen_diff) > 1e-15){
    
    ##add some number to current vbmax
    vbmax = vbmax + step
    
    ##calculate new leslie matrix
    leslie = matrix(nrow = leslie_rows,ncol = leslie_rows,0)
    leslie[leslie_rows/2+1,leslie_rows] = 2 * growth_rate(c_eq,vbmax,1) 
    inc = 1
    for (i in (leslie_rows/2+1):leslie_rows){
      if (i < leslie_rows){
        leslie[i+1,i] =  growth_rate(c_eq,vbmax,1)
        leslie[i,i] = leslie[i,i] + (1 - dormancy_rate) * (1 - spont_death) * (1 - growth_rate(c_eq,vbmax,1))
      }else{
        leslie[i,i] = leslie[i,i] + (1 - dormancy_rate) * (1 - spont_death) * (1 - growth_rate(c_eq,vbmax,1))
      }
      leslie[inc,i] =leslie[inc,i] + dormancy_rate
      inc = inc + 1
    }
    for (i in 1:(leslie_rows/2)){
      leslie[i,i] =(1 - resuscitation_rate) * (1 - spont_death)
      leslie[inc,i] = leslie[inc,i] + resuscitation_rate
      inc = inc + 1
    }
    
    ##find difference between 1 and first eigenvalue of the leslie matrix
    eigen_diff =as.numeric( 1 - abs(eigen(leslie)$values[1]))
    
    ##adjust step size 
    if (
      abs(eigen_diff) < abs(old_eigen_diff)){
    } else {
      step = -step * .9
    }
    if (c_eq < 0){
      c_eq = 0
      step = abs(step)
    }
   old_eigen_diff = eigen_diff
   # print(eigen_diff)
  }
  


  tmp_mat = matrix(0, nrow = pheno_bins, ncol = stages)
  popul = list("active" = tmp_mat,"dormant" = tmp_mat, "induced" = tmp_mat, "non-dormant" = tmp_mat,"non-dormant-induced" = tmp_mat)
  ##confusingly, dormants are the first half of the leslie matrix and actives are the second half
  if(all(as.numeric(eigen(leslie)$vector[1:(leslie_rows/2),1]) == 0)){
    popul[[2]][1,] = rep(0,ncol(popul[[2]]))
  } else {
  popul[[2]][1,] = round(N * as.numeric(eigen(leslie)$vector[1:(leslie_rows/2),1]/sum(eigen(leslie)$vector[,1])))}
  popul[[1]][1,] = round(N * as.numeric(eigen(leslie)$vector[(leslie_rows/2 + 1):(leslie_rows),1]/sum(eigen(leslie)$vector[,1])))
  pop_vec = matrix(c(popul[[2]][1:3],popul[[1]][1:3]))
  pv= pop_vec
  for (i in 1:1000000){
    # pop_vec = leslie %*% pop_vec
  }
  pop_vec = round(pop_vec)
  print(pop_vec[4:6])
  vmax_carb_active <<- source / (N * (sum(popul[[1]])/sum(popul[[1]]+popul[[2]])) * (c_eq /(km_carb + c_eq)))

  return(list(popul=popul,vbmax=vbmax))
  }
  popul = vbmax_equibpop_determination(dormancy_rate,resuscitation_rate,spont_death)[[1]]
  vbmax = vbmax_equibpop_determination(dormancy_rate,resuscitation_rate,spont_death)[[2]]
  sum(popul[[1]])/sum(popul[[1]]+popul[[2]])
  