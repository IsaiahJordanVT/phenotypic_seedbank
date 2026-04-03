library(reshape2)
library(reshape)
library(ggplot2)
library(vegan)

## The following block (to the next note) contains functions relevant for determining the rates at which individuals
## transition, as well as functions which actually determine which transistions will actually happen for a given timestep
## and process

##simple function to calculate michealis menten growth rate
growth_rate = function(carb, vbmax, tmst) {
  birth =  vbmax * carb / (km_carb + carb) * tmst
  return(birth)
}

##same as above, but with an added coefficient to make the rate higher or lower
cost_penalty = function(carb, vbmax, cost, tmst) {
  birth = cost * vbmax * carb / (km_carb + carb) * tmst
  return(birth)
}


##calculate where 'offspring' end up - last column results in births that get send to the middle + 1th column
births = function(carb, vbmax, tmst, population) {
  # print(population)
  births = 0
  ##calculate current rate
  curr_bmax = growth_rate(carb, vbmax, tmst)
  ##matrix to store results in
  grown_individuals_matrix = matrix(rep(0, length(population)), ncol = ncol(population))
  
  ##these loops determine how many offspring each element produces
  for (j in 1:(stages)) {
    for (i in seq(1, length(grown_individuals_matrix[, 1]))) {
      if (population[i, j] > 5e20) {
        indiv_sampled_to_grow = rpois(1, population[i, j] * curr_bmax)
        if (indiv_sampled_to_grow > population[i,j]){
          indiv_sampled_to_grow = population[i,j]
        }
      } else if (population[i, j] == 0) {
        indiv_sampled_to_grow= 0
      } else {
        indiv_sampled_to_grow= rbinom(1, (population[i, j]), 1 - exp(-curr_bmax))
      }
      grown_individuals_matrix[i, j] = indiv_sampled_to_grow
    }
  }
  pop_change = matrix(ncol = ncol(grown_individuals_matrix), nrow = nrow(grown_individuals_matrix), 0)
  for (i in 1:nrow(grown_individuals_matrix)) {
    ##is this individual in the last column? If so, it divides.
    for (j in 1:ncol(grown_individuals_matrix)) {
      next_col = j + 1
      if (next_col > ncol(grown_individuals_matrix)) {
        next_col = 1
        divided = TRUE
      } else {
        divided = FALSE
      }
      
      ##depleted_phenotype is used to determine which row the new individual should appear in. 
      depleted_phenotype = i - i / (2 * ncol(grown_individuals_matrix))
      # print(c('1',depleted_phenotype))
      ##If depleted_phenotype = 1, phenotype cannot decline, and all offspring simply move to the next column
      if (any(depleted_phenotype < 1)) {
        depleted_phenotype = 1
      }
      
      ##if depleted_phenotype isn't an integer, it is rounded up and down, and the remainder and its complement stored
      if (depleted_phenotype !=  as.integer(depleted_phenotype)) {
        remainder = c(1 - (depleted_phenotype - floor(depleted_phenotype)), (depleted_phenotype - floor(depleted_phenotype)))
        depleted_phenotype = c(floor(depleted_phenotype), ceiling(depleted_phenotype))
      }
      
      ##if depleted_phenotype contains numbers greater than one, determine where offspring should go
      if (length(depleted_phenotype) > 1) {
        if (divided == T) {
          ## if, after rounding, the sum of multiple of the remainder vector times the element of the
          ##growth matrix is less than one, add an individual to the lower phenotypic bin
          if (sum(round(remainder * grown_individuals_matrix[i, j])) < grown_individuals_matrix[i, j]) {
            pop_change[depleted_phenotype, next_col] = pop_change[depleted_phenotype, next_col] + (round(remainder  * grown_individuals_matrix[i, j]) + c(1, 0)) * 2
            births = births + grown_individuals_matrix[i,j]
            ##otherwise, don't.
          } else if (sum(round(remainder * grown_individuals_matrix[i, j])) > grown_individuals_matrix[i, j]) {
            pop_change[depleted_phenotype, next_col] = pop_change[depleted_phenotype, next_col] + (round(remainder * grown_individuals_matrix[i, j]) + c(0, -1)) * 2
            births = births + grown_individuals_matrix[i,j]
          } else {
             pop_change[depleted_phenotype, next_col] = pop_change[depleted_phenotype, next_col] + round(remainder * 2 * grown_individuals_matrix[i, j])
            births = births + grown_individuals_matrix[i,j]
          }

          ## if we aren't in the last column, we don't multiply by two.
        } else {
          if (sum(round(remainder * grown_individuals_matrix[i, j])) < grown_individuals_matrix[i, j]) {
            pop_change[depleted_phenotype, next_col] = pop_change[depleted_phenotype, next_col] + round(remainder * grown_individuals_matrix[i, j]) + c(1, 0)
          } else if (sum(round(remainder * grown_individuals_matrix[i, j])) > grown_individuals_matrix[i, j]) {
            pop_change[depleted_phenotype, next_col] = pop_change[depleted_phenotype, next_col] + round(remainder * grown_individuals_matrix[i, j]) + c(0, -1)
          } else {
            pop_change[depleted_phenotype, next_col] = pop_change[depleted_phenotype, next_col] + round(remainder * grown_individuals_matrix[i, j])
          }
        }
      } else {
        #if depleted_phenotype happens to be integer-valued (e.g., 1) ignore the above and just add the individuals
        if (divided == T) {
          pop_change[depleted_phenotype, next_col] = pop_change[depleted_phenotype, next_col] + 2 * grown_individuals_matrix[i, j]
          births = births + grown_individuals_matrix[i,j]
        } else {
          pop_change[depleted_phenotype, next_col] = pop_change[depleted_phenotype, next_col] + grown_individuals_matrix[i, j]
        }
      }
    }
  }
  
  ## remove the individuals who grow from their elements.
  pop_change = pop_change - grown_individuals_matrix
  return(list(pop_change,births))
}

##This function samples individuals in the "penalty box," which have recently induced a more tolerant phenotype
births_penalty = function(cost, carb, vbmax, tmst, population){
  curr_bmax = cost_penalty(carb, vbmax,cost, tmst)

  ##matrix to store results in
  grown_individuals_matrix = matrix(rep(0, length(population)), ncol = ncol(population))
  
  ##these loops determine how many offspring each element produces
  for (j in 1:(ncol(population))) {
    for (i in seq(1, length(grown_individuals_matrix[, 1]))) {
      if (population[i, j] > 5e20) {
        indiv_sampled_to_grow = rpois(1, population[i, j] * curr_bmax)
        if (indiv_sampled_to_grow > population[i,j]){
          indiv_sampled_to_grow = population[i,j]
        }
      } else {
        indiv_sampled_to_grow= rbinom(1, (population[i, j]), 1 - exp(-curr_bmax))
      }
      grown_individuals_matrix[i, j] = indiv_sampled_to_grow
    }
  }
  return(list("grown_individuals_matrix" = grown_individuals_matrix))
}


##death rates. In the absence of toxin, death occurs at a constant rate
death = function(toxin, selected_population,spont_death, tmst) {
  ##set up a matrix to store the number of deaths in each element
  death_mat = matrix(ncol = ncol(selected_population),
                     nrow = nrow(selected_population),
                     0)
  for (i in 1:nrow(selected_population)) {
    ##calculate the rate of death for each row - that is, phenotype
    pheno = i
    rate_death = vmax_death * (toxin ^ (hill_coef) / ((km_death * pheno) ^ hill_coef + toxin ^
                                           (hill_coef)))
    for (j in 1:ncol(selected_population)){
      if (selected_population[i, j] > 5e20) {
        ##sum the two contributors of death - spontaneous and toxin induced
        death_mat[i, j] = rpois(1,
                                rate_death * selected_population[i, j] * tmst + spont_death * (selected_population[i, j]) *
                                  tmst)
      } else {
        death_mat[i, j] = rbinom(1, (selected_population[i, j]), (1 - exp(-rate_death * tmst)) +
                                   (1 - exp(-spont_death * tmst)))
      }
    }
  }
  return(death_mat)
}

##dormancy rates - should be  fairly self explanatory at this point. from the first half of the
##population matrix, sample individuals to resuscitate. From the second half, sample individuals
##to go dormant.
dorm_trans = function(tmst,dormancy_rate,resuscitation_rate, population) {
  dormancy_rate_scaled =  dormancy_rate * tmst
  resuscitation_rate_scaled =  resuscitation_rate * tmst
  new_dormant_matrix = matrix(nrow = pheno_bins, ncol = ncol(population[[1]]), 0)
  for (i in seq(1, pheno_bins)) {
    for (j in (1):ncol(population[[1]])) {
      if (population[[1]][i, j] > 5e20) {
        new_dormants = sum(rpois(1, population[[1]][i, j] * dormancy_rate_scaled))
      } else if (population[[1]][i, j] == 0) {
        new_dormants = 0
      } else {
        new_dormants = sum(rbinom(1, (population[[1]][i, j]), 1 - exp(-dormancy_rate_scaled)))
      }
      new_dormant_matrix[i, j] = new_dormants
    }
  }
  resuscitated_matrix = matrix(nrow = pheno_bins, ncol = ncol(population[[2]]), 0)
  for (i in seq(1, pheno_bins)) {
    for (j in 1:(ncol(population[[1]]))) {
      if (population[[2]][i, j] > 5e20) {
        new_resuscitated = sum(rpois(1, population[[2]][i, j] * resuscitation_rate_scaled))
        } else if (population[[2]][i, j] == 0) {
        new_resuscitated = 0
      } else {
        new_resuscitated = sum(rbinom(1, (population[[2]][i, j]), 1 - exp(-resuscitation_rate_scaled)))
      }
      resuscitated_matrix[i, j] = new_resuscitated
    }
  }
  return(cbind(resuscitated_matrix, new_dormant_matrix))
}

##same as above, but for indiviudals in the 'penalty box'
dorm_trans_pb = function(tmst,dormancy_rate,resuscitation_rate, population) {
  dormancy_rate_scaled =  dormancy_rate * tmst
  resuscitation_rate_scaled =  resuscitation_rate * tmst
  new_dormant_matrix = matrix(nrow = pheno_bins, ncol = (ncol(population[[1]])), 0)
  for (i in seq(1, pheno_bins)) {
    for (j in (1:ncol(population[[1]]))) {
      if (population[[3]][i, j] > 5e20) {
        new_dormants = sum(rpois(1, population[[3]][i, j] * dormancy_rate_scaled))
      } else if (population[[3]][i, j] == 0) {
        new_dormants = 0
      } else {
        new_dormants = sum(rbinom(1, (population[[3]][i, j]), 1 - exp(-dormancy_rate_scaled)))
      }
      new_dormant_matrix[i, j] = new_dormants
    }
  }
  return(new_dormant_matrix)
}

## each individual transitions through the compartments at a constant rate dependent on toxin levels. They are 
## however first consigned to a 'penalty box,' which they must be sampled for growth to escape
induction_rates = function(toxin, population, tmst, hetero) {
  induced_matrix = matrix(nrow = (pheno_bins - 1),
                        ncol = ncol(population),
                        0)
  induction = (vmax_toxin * (toxin / (km_toxin + toxin)) + hetero) * tmst
  for (i in seq(1, pheno_bins - 1)) {
    for (j in 1:ncol(population)) {
      if (population[i, j] > 5e20) {
        indiv_sampled_to_induce = rpois(1, population[i, j] * induction)
      } else if (population[i, j] == 0) {
        indiv_sampled_to_induce = 0
      } else {
        indiv_sampled_to_induce = sum(rbinom(1, (population[i, j]), 1 - exp(-induction)))
      }
      induced_matrix[i, j] = indiv_sampled_to_induce
    }
  }
  penalty_box = rbind(rep(0, ncol(induced_matrix)), induced_matrix)
  growers = - rbind(induced_matrix, rep(0, ncol(induced_matrix)))
  returns = list("penalty_box" = penalty_box, "growers" = growers)
  return(returns)
}

##this is a function to determine intrinsic loss of phenotype
phenotype_loss = function(sub_population,tmst){
  pheno_loss = matrix(0,nrow = pheno_bins,ncol = stages)
  for (j in 1:ncol(sub_population)){
    for (i in 2:pheno_bins){
      if (sub_population[i,j] > 5e20){
        num_lost = rpois(1,sub_population[i,j] * active_pheno_loss*tmst)
        pheno_loss[(i-1),j] = pheno_loss[(i-1),j] + num_lost
        pheno_loss[i,j] = pheno_loss[i,j] - num_lost
      } else {
      num_lost  = rbinom(1,sub_population[i,j],1-exp(-active_pheno_loss*tmst))
      pheno_loss[(i-1),j] = pheno_loss[(i-1),j] + num_lost
      pheno_loss[i,j] = pheno_loss[i,j] - num_lost
      } 
      }
  }
  return(pheno_loss)
  }

## carbon utilization rate - MM rate * number of cells actively growing
carbon_util = function(carb, population,tmst) {
  logcarb = log(carb)
  active = (sum(population[[1]])+sum(population[[3]])+sum(population[[4]]+sum(population[[5]]))) * vmax_carb_active * tmst * ((carb) / (km_carb + (carb)))
  return(exp((logcarb+(-active/carb + ((source * tmst)/carb)))))
  # return(carb)
}


## simulate the system
##arguments:
##population - matrix with the initial population
##cost - the cost modifier for the birth rate, which determines the rate of leaving the penalty box
##hetero - the amount of heterogeneity
##vbmax - maximum growth rate, solved for in leslie_matrix_calculations.
##init_carb - initial carbon in the system
##toxin_strength - vector with toxin concentrations
##toxin_times - vector with times that toxin pulses occur
##max_time - if the population doesn't go extinct how long to simulate for
##init_ts - the initial timestep
##mutation_rate - probability that a mutant arises. for metapopulation - invasion model
##migrant_frequency - rate at which indiviudals colonize the environment. in the metapopulation model will be 
## rate of migration between populations. Used in regular model
##toxin_mode - the type of toxin pulse received - constant pulse providing constant toxin concentration (T), or typical
##pulse and decay type (F)

##The below block runs the simulation with costs, and returns the results.
simulation = function(population,
                      cost,
                      hetero,
                      vbmax,
                      init_carb = c_eq,
                      toxin_strength,
                      toxin_times = 20,
                      max_time = 500,
                      init_ts = 0.01,
                      dormancy_rate = .0001,
                      resuscitation_rate = dormancy_rate * 10000,
                      spont_death = 0.001,
                      d_r = .1
){
  ##re-tries simulating the system when an error occurs, shouldn't be needed
  done = F
  attempts = 0
  while(done == F & attempts < 10){
    try({
      attempts = attempts +1
      
      ##preallocation of the array that holds population data for each timestep. saves a lot of time
      pops = array(0, c(pheno_bins, stages*2, max_time * 10^3))
      
      ##for determining runtime of a simulation
      start_time = Sys.time()
      
      ##initializing tmst (timestep)
      tmst = init_ts
      
      ##initializing population
      pops[, , 1] = cbind(rowSums(population[[2]]),rowSums(population[[1]]+population[[3]]))
      
      ##initializing a list of carbon concentrations
      carbs = c()
      
      ##initializing time
      time = 0
      
      ##initailizing list of timesteps, useless if tmst is fixed
      tss = c()
      
      ##initializing list of times - identical to seq(0,max_time,tmst) if tmst is fixed
      times = c(0)
      
      ##initializing toxin
      toxin = 0
      tox_flag = c(0,0)
      
      ##initializing a list of toxin concentrations over time
      tox_list = c()
      
      ##initializing a list of total population size at each time
      pop_list = c(sum(population[[1]]+population[[2]]+population[[3]]))
      
      ##below, functions defining the different events that can occur are defined, and applied to the population
      
      ##induction
      induction_routine = function(tmp_pop){
        induced = induction_rates(toxin, tmp_pop[[1]], tmst, hetero)
        tmp_pop[[1]] <- tmp_pop[[1]] + induced$growers
        tmp_pop[[3]] <- tmp_pop[[3]] + induced$penalty_box
        return(tmp_pop)
      }
      
      ## loss of phenotype
      loss_routine = function(tmp_pop){
        p_l = phenotype_loss(tmp_pop[[1]],tmst)
        tmp_pop[[1]] <- tmp_pop[[1]] + p_l
        return(tmp_pop)
      }
      
      ## transition between dormant and non-dormant states
      dorm_routine = function(tmp_pop){
        dorm_transition <- dorm_trans(tmst,dormancy_rate,resuscitation_rate, tmp_pop)
        pb_dorm_transition <- dorm_trans_pb(tmst,dormancy_rate,resuscitation_rate,tmp_pop)
        tmp_pop[[1]]<- tmp_pop[[1]] - dorm_transition[, (ncol(dorm_transition) /
                                                           2 + 1):ncol(dorm_transition)]
        tmp_pop[[1]] <- tmp_pop[[1]] + dorm_transition[, 1:(ncol(population[[1]]))]
        tmp_pop$tmp_pop_2 <- tmp_pop$tmp_pop_2 - dorm_transition[, 1:(ncol(population[[1]]))]
        tmp_pop$tmp_pop_2 <- tmp_pop$tmp_pop_2 + dorm_transition[, (ncol(dorm_transition) /
                                                                      2 + 1):ncol(dorm_transition)]
        tmp_pop$tmp_pop_2 <- tmp_pop$tmp_pop_2 + pb_dorm_transition
        tmp_pop$tmp_pop_3 <- tmp_pop$tmp_pop_3 - pb_dorm_transition
        return(tmp_pop)
      }
      
      
      
      ##birth
      birth_routine = function(tmp_pop){
        pop_change <- births(carb, vbmax, tmst, tmp_pop[[1]])
        pop_change_cost <- births_penalty(cost,carb, vbmax, tmst, tmp_pop[[3]])
        tmp_pop[[1]] <- tmp_pop[[1]] + pop_change[[1]]
        tmp_pop[[1]] <- tmp_pop[[1]] + pop_change_cost$grown_individuals_matrix
        tmp_pop[[3]] <- tmp_pop[[3]] - pop_change_cost$grown_individuals_matrix
        return(list(tmp_pop,sum(pop_change[[2]])))
      }
      
      ##death
      death_routine = function(tmp_pop){
        death_pop_1 = death(toxin, tmp_pop[[1]],spont_death, tmst)
        death_pop_2 = death(toxin,tmp_pop[[2]],spont_death,tmst)
        death_pop_3 = death(toxin, tmp_pop[[3]],spont_death, tmst)
        tmp_pop[[1]] <- tmp_pop[[1]] - death_pop_1
        tmp_pop[[2]] <- tmp_pop[[2]] - death_pop_2
        tmp_pop[[3]] <- tmp_pop[[3]] - death_pop_3
        return(list(tmp_pop,sum(death_pop_1,death_pop_2,death_pop_3)))
      }
      
      ##list of processes that occur in each timestep, used for randomization
      order_ops = c('induction_routine','loss_routine','dorm_routine','birth_routine','death_routine')
      
      ##current carbon levels
      carb = init_carb
      tmst = init_ts
      r = 1
      
      
      ##vector to store births and deaths in
      bir_dea = rep(NA, max_time * 10^3)
      number = 0
      ##beginning of the simulation loop
      while (time < max_time & sum(population[[1]]+population[[2]]+population[[3]]) > 0)
      {
        carbs = c(carbs,carb)
        number = number +1
        
        ##Is it time for a carbon pulse? Did the pulse already occur?
          if(time > toxin_times[1] & tox_flag[1] == 0) {
            toxin = toxin + toxin_strength[1]
            tox_flag[1] = 1
          }else{
            if(time > toxin_times[2] & tox_flag[2] == 0 & length(toxin_strength) > 1) {
              toxin = toxin + toxin_strength[2]
              tox_flag[2] = 1
            }
          }
          toxin = toxin -  toxin * tmst * d_r
          if(toxin < 1.2e-11){
            toxin = 0
          }
        
        
        ##record toxin level
        tox_list = c(tox_list, toxin)
        tmp_pop = list("tmp_pop_1" = population[[1]],
                       "tmp_pop_2" = population[[2]],
                       "tmp_pop_3" = population[[3]])
        
        ##randomize the order of operations
        rand_ops = sample(1:5,5)
        for(i in rand_ops){
          func = get(order_ops[i])
          if(i == 4){
            birth_op = func(tmp_pop)
            tmp_pop = birth_op[[1]]
            num_births = birth_op[[2]]
          }else if(i == 5){
            death_op = func(tmp_pop)
            tmp_pop = death_op[[1]]
            num_deaths = death_op[[2]]
            scavangeable_carb = death_op[[2]]
          }else{
            tmp_pop = func(tmp_pop)
          }
        }
        
        tmp_pop_1 = tmp_pop[[1]]
        tmp_pop_2 = tmp_pop[[2]]
        tmp_pop_3 = tmp_pop[[3]]
        
        total_actives <- (tmp_pop_1 + tmp_pop_3)
        total_pop <- cbind(rowSums(tmp_pop_2),rowSums(total_actives))
        
        ##adjusting timestep according to change of population in previous step
        r <- max(.01,abs(sum(total_pop) - sum(pops[,,length(times)]))/(sum(pops[,,length(times)])))
        tmst <- allowed_error/r
        bir_rate = vbmax * carb / (km_carb + carb) * tmst
        dea_rate = vmax_death * (toxin ^ (hill_coef) / ((km_death * 1) ^ hill_coef + toxin ^(hill_coef)))* tmst
        ind_rate = (vmax_toxin * (toxin / (km_toxin + toxin)) + hetero) * tmst
        probabilities = 1-exp(-c(bir_rate,dea_rate,ind_rate))
        while(any(probabilities == 1)){
        tmst = tmst * .5
        bir_rate = vbmax * carb / (km_carb + carb) * tmst
        dea_rate = vmax_death * (toxin ^ (hill_coef) / ((km_death * 1) ^ hill_coef + toxin ^(hill_coef)))* tmst
        ind_rate = (vmax_toxin * (toxin / (km_toxin + toxin)) + hetero) * tmst
        probabilities = 1-exp(-c(bir_rate,dea_rate,ind_rate))
        }
        if(tmst > 10*init_ts){
          tmst = init_ts*10 
        }
        
        ##store the outcome of this timestep
        population[[1]] <- tmp_pop_1
        population[[2]] <- tmp_pop_2
        population[[3]] <- tmp_pop_3
        
        carb = carbon_util(carb, population, tmst)
        bir_dea[length(times)] = num_births - num_deaths
        tss <- c(tss, tmst)
        times <- c(times, time)
        time <- time + tmst
        pops[, , length(times)] <- total_pop
        pop_list <- c(pop_list, sum(total_pop))
        
      }
      ##trim unused portion of the array. to initialize the population matrix storing data for each timestep, a much bigger
      ##empty matrix is generated
      pops = pops[, , 1:length(times)]
      bir_dea = bir_dea[1:length(times)]
      return(
        list(
          'pops' = pops,
          'times' = times,
          'tss' = tss,
          'pop_list' = pop_list,
          'carbs' = carbs,
          'tox_list' = tox_list,
          'last_pop' = population, 
          'bir_dea' = bir_dea
        )
      )
      
      ##did this simulation complete without errors? if so, don't try again. belongs to the try loop, probably unnecessary now.
      
    })
    done = T
    
  }
}
###This block generates the "rainbow" plot of subpopulations of phenotypes and dormancy state'''

## shaping the data for plotting. also generates a plot. Useful if you want to see changes in phenotype over time
plotting = function(population, pops) {
  times = pops$times
  time_list = rep(0, length(times))
  inc = 0
  pop = melt.array(pops$pops[,,seq(1,dim(pops$pops)[[3]],50)])
  colnames(pop) = c('Tolerance','State','Time','Population')
  pop$Time = times[as.factor(pop$Time)]
  dorm_act = c("Dormant","Active")
  pop$State = dorm_act[as.factor(pop$State)]
  # print(tail(pop))
  pop$Population = as.numeric(pop$Population)
  pop$State = as.character(pop$State)
  pop$Tolerance = factor(pop$Tolerance,levels = as.character(seq(1,pheno_bins)))
  pop$Time = as.numeric(pop$Time)*50/280.4095 * 30
  # print(unique(pop$Time))
  #pop = pop[pop$Tolerance == 1,]
  #pop = pop[pop$State == 'Active',]
  pop_no_dorm = subset(pop,select = -State)
  plt = ggplot(pop,
               aes(
                 x = Time,
                 y = Population,
                 color = Tolerance,
                 linetype = State
               )) +
    geom_line() +
    geom_vline(xintercept = 11)+
    xlim(c(0,pop$Time[length(pop$Time)]))+
    geom_hline(yintercept = 0) +
    #scale_alpha_discrete(c(0.4,1)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(10 ^ seq(0, 12, 3))) +
    # theme_classic()+
    theme(plot.title = element_text(hjust = 0.5,size = 30),text = element_text(size = 17),axis.text = element_text(size =20))
  plot(plt)
}
