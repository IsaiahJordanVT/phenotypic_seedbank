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



births <- function(carb, vbmax, tmst, population) {
  # 1. Initialization and dimension checks
  num_rows <- nrow(population)
  num_cols <- ncol(population)
  
  births <- 0
  curr_bmax <- growth_rate(carb, vbmax, tmst)
  
  # Matrix to store the number of individuals transitioning
  grown_indiv_mat <- matrix(0, nrow = num_rows, ncol = num_cols)
  
  # 2. Vectorized sampling for transitions
  # Use rbinom for safe values; switch to rpois only for extremely large values
  p_trans <- 1 - exp(-curr_bmax)
  
  # Vectorized calculation over the whole population matrix
  # If population is 0, rbinom handles it naturally (returns 0)
  grown_indiv_mat[] <- rbinom(length(population), population, p_trans)
  
  # Correct for extreme values (> 5e20) where rbinom fails due to integer limits
  extreme_idx <- which(population > 5e20)
  if (length(extreme_idx) > 0) {
    poisson_samples <- rpois(length(extreme_idx), population[extreme_idx] * curr_bmax)
    grown_indiv_mat[extreme_idx] <- pmin(poisson_samples, population[extreme_idx])
  }
  
  # Matrix tracking net population movement
  pop_change <- matrix(0, nrow = num_rows, ncol = num_cols)
  
  # 3. Process transitions and phenotype depletion
  for (j in 1:num_cols) {
    # Determine the destination column (with wrap-around)
    next_col <- if (j == num_cols) 1 else j + 1
    divided  <- (j == num_cols)
    multiplier <- if (divided) 2 else 1
    
    for (i in 1:num_rows) {
      count <- grown_indiv_mat[i, j]
      if (count == 0) next
      
      # Track global birth metrics if dividing
      if (divided) {
        births <- births + count
      }
      
      # Calculate continuous phenotype drop mapping
      depleted_phenotype <- max(1, i - i / (2 * num_cols))
      floor_p <- floor(depleted_phenotype)
      
      # Handle integer vs non-integer phenotypes cleanly
      if (depleted_phenotype == floor_p) {
        pop_change[floor_p, next_col] <- pop_change[floor_p, next_col] + (multiplier * count)
      } else {
        ceil_p <- ceiling(depleted_phenotype)
        
        # Calculate remainder weights
        rem_ceil <- depleted_phenotype - floor_p
        rem_floor <- 1 - rem_ceil
        
        # Distribute counts using proportional rounding
        count_floor <- round(rem_floor * count)
        count_ceil  <- round(rem_ceil * count)
        total_round <- count_floor + count_ceil
        
        # Adjust rounding errors to match total count perfectly
        if (total_round < count) {
          count_floor <- count_floor + 1
        } else if (total_round > count) {
          count_ceil <- count_ceil - 1
        }
        
        # Apply changes to the grid
        pop_change[floor_p, next_col] <- pop_change[floor_p, next_col] + (multiplier * count_floor)
        pop_change[ceil_p, next_col]  <- pop_change[ceil_p, next_col]  + (multiplier * count_ceil)
      }
    }
  }
  
  # 4. Subtract individuals leaving their original bins
  pop_change <- pop_change - grown_indiv_mat
  
  return(list(pop_change = pop_change, births = births))
}

births_penalty = function(cost, carb, vbmax, tmst, population){
  curr_bmax = cost_penalty(carb, vbmax, cost, tmst)
  
  # Logical mask for the condition population > 5e20
  high_pop = population > 5e20
  
  # Initialize the result matrix
  grown_individuals_matrix = matrix(0, nrow = nrow(population), ncol = ncol(population))
  
  # Vectorized Poisson branch for high populations
  if (any(high_pop)) {
    lambda = population[high_pop] * curr_bmax
    # rpois can handle vector inputs for lambda
    sampled_poisson = rpois(sum(high_pop), lambda)
    # Ensure sampled offspring do not exceed population size
    grown_individuals_matrix[high_pop] = pmin(sampled_poisson, population[high_pop])
  }
  
  # Vectorized Binomial branch for low/normal populations
  low_pop = !high_pop
  if (any(low_pop)) {
    size_vals = population[low_pop]
    prob_vals = 1 - exp(-curr_bmax)
    # rbinom vectorized over size parameter
    grown_individuals_matrix[low_pop] = rbinom(sum(low_pop), size = size_vals, prob = prob_vals)
  }
  
  return(list("grown_individuals_matrix" = grown_individuals_matrix))
}

##death rates. In the absence of toxin, death occurs at a constant rate
death = function(toxin, selected_population, spont_death, tmst) {
  # Calculate vector of rates per phenotype (row)
  phenos <- 1:nrow(selected_population)
  rate_death <- vmax_death * (toxin^hill_coef / ((km_death * phenos)^hill_coef + toxin^hill_coef))
  
  # Match rates to matrix dimensions
  rate_matrix <- matrix(rate_death, nrow = nrow(selected_population), ncol = ncol(selected_population), byrow = FALSE)
  
  # Correct joint probability: P(Death) = 1 - exp(-(Rate_Toxin + Rate_Spont) * tmst)
  # This mathematically prevents probabilities from exceeding 1.0 completely
  prob_bin <- 1 - exp(-(rate_matrix + spont_death) * tmst)
  # Execute stochastic draw for the entire matrix at once
  death_mat <- matrix(rbinom(length(selected_population), selected_population, prob_bin), 
                      nrow = nrow(selected_population))
  return(death_mat)
}

##dormancy rates - should be  fairly self explanatory at this point. from the first half of the
##population matrix, sample individuals to resuscitate. From the second half, sample individuals
##to go dormant.
dorm_trans = function(tmst, dormancy_rate, resuscitation_rate, population) {
  # Calculate transition probabilities
  prob_dorm  <- 1 - exp(-dormancy_rate * tmst)
  prob_resus <- 1 - exp(-resuscitation_rate * tmst)
  
  # Extract active and dormant sub-populations from the list
  active_pop  <- population[[1]]
  dormant_pop <- population[[2]]
  
  # Sample dormancy for all elements of the active matrix at once
  new_dormant_matrix <- matrix(
    rbinom(length(active_pop), active_pop, prob_dorm), 
    nrow = nrow(active_pop)
  )
  
  # Sample resuscitation for all elements of the dormant matrix at once
  resuscitated_matrix <- matrix(
    rbinom(length(dormant_pop), dormant_pop, prob_resus), 
    nrow = nrow(dormant_pop)
  )
  
  # Bind matrices matching your original layout: [Resuscitated | New Dormant]
  return(cbind(resuscitated_matrix, new_dormant_matrix))
}

##same as above, but for indiviudals in the 'penalty box'
dorm_trans_pb = function(tmst, dormancy_rate, resuscitation_rate, population) {
  # Calculate transition probability
  prob_dorm <- 1 - exp(-dormancy_rate * tmst)
  
  # Extract the penalty box sub-population matrix
  pb_pop <- population[[3]]
  
  # Sample dormancy for all elements of the penalty box matrix simultaneously
  new_dormant_matrix <- matrix(
    rbinom(length(pb_pop), pb_pop, prob_dorm), 
    nrow = nrow(pb_pop)
  )
  
  return(new_dormant_matrix)
}

## each individual transitions through the compartments at a constant rate dependent on toxin levels. They are 
## however first consigned to a 'penalty box,' which they must be sampled for growth to escape
induction_rates = function(toxin, population, tmst, hetero) {
  # Calculate induction probability
  induction <- (vmax_toxin * (toxin / (km_toxin + toxin)) + hetero) * tmst
  prob_induce <- 1 - exp(-induction)
  
  # Isolate all rows except the final phenotype bin
  # (Matching your seq(1, pheno_bins - 1) logic)
  sub_pop <- population[1:(nrow(population) - 1), , drop = FALSE]
  
  # Draw transitions for all isolated elements simultaneously
  induced_draws <- matrix(
    rbinom(length(sub_pop), sub_pop, prob_induce), 
    nrow = nrow(sub_pop)
  )
  
  # Construct the output matrices using row-binding offsets
  # No need for rep(0, ncol(...)), R pads a single 0 across the row automatically
  penalty_box <- rbind(0, induced_draws)
  growers     <- -rbind(induced_draws, 0)
  
  return(list("penalty_box" = penalty_box, "growers" = growers))
}

##this is a function to determine intrinsic loss of phenotype
phenotype_loss = function(sub_population, tmst){
  # Calculate loss probability
  prob_loss <- 1 - exp(-active_pheno_loss * tmst)
  
  # Isolate all rows that can lose phenotype (Rows 2 to pheno_bins)
  loss_eligible_pop <- sub_population[2:nrow(sub_population), , drop = FALSE]
  
  # Draw the number of individuals losing their phenotype simultaneously
  num_lost <- matrix(
    rbinom(length(loss_eligible_pop), loss_eligible_pop, prob_loss), 
    nrow = nrow(loss_eligible_pop)
  )
  
  # Initialize the output matrix matching sub_population dimensions
  pheno_loss <- matrix(0, nrow = nrow(sub_population), ncol = ncol(sub_population))
  
  # Row offset assignment: 
  # Cells that lose phenotype subtract from their current row (Rows 2 to pheno_bins)
  # and add to the row below them (Rows 1 to pheno_bins-1)
  pheno_loss[1:(nrow(pheno_loss) - 1), ] <- pheno_loss[1:(nrow(pheno_loss) - 1), ] + num_lost
  pheno_loss[2:nrow(pheno_loss), ]       <- pheno_loss[2:nrow(pheno_loss), ] - num_lost
  
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

simulation = function(population, cost, hetero, vbmax, init_carb = c_eq, toxin_strength, 
                      toxin_times = 20, max_time = 500, init_ts = 0.01, dormancy_rate = .0001, 
                      resuscitation_rate = dormancy_rate * 10000, spont_death = 0.001, d_r = .1 ){
  
  done = FALSE
  attempts = 0
  
  while(done == FALSE && attempts < 10){
    try({
      attempts = attempts + 1
      
      # 1. Structural Pre-allocations (Protects RAM from exploding)
      max_steps <- max_time * 10^3
      pops <- array(0, c(pheno_bins, stages * 2, max_steps))
      
      # Pre-allocate tracking vectors instead of using slow c() growth
      carbs    <- numeric(max_steps)
      tss      <- numeric(max_steps)
      times    <- numeric(max_steps)
      tox_list <- numeric(max_steps)
      pop_list <- numeric(max_steps)
      bir_dea  <- numeric(max_steps)
      
      # Initial conditions
      tmst <- init_ts
      time <- 0
      toxin <- 0
      tox_flag <- c(0, 0)
      carb <- init_carb
      step_idx <- 1
      
      # Track starting totals
      pops[, , 1] <- cbind(rowSums(population[[2]]), rowSums(population[[1]] + population[[3]]))
      carbs[1]    <- carb
      times[1]    <- time
      tox_list[1] <- toxin
      pop_list[1] <- sum(population[[1]] + population[[2]] + population[[3]])
      
      # ------------------------------------------------------------
      # 2. Re-defined Sub-routines using Vectorized Engines
      # ------------------------------------------------------------
      induction_routine = function(tmp_pop){
        induced <- induction_rates(toxin, tmp_pop[[1]], tmst, hetero)
        tmp_pop[[1]] <- tmp_pop[[1]] + induced$growers
        tmp_pop[[3]] <- tmp_pop[[3]] + induced$penalty_box
        return(tmp_pop)
      }
      
      loss_routine = function(tmp_pop){
        tmp_pop[[1]] <- tmp_pop[[1]] + phenotype_loss(tmp_pop[[1]], tmst)
        return(tmp_pop)
      }
      
      dorm_routine = function(tmp_pop){
        dorm_transition    <- dorm_trans(tmst, dormancy_rate, resuscitation_rate, tmp_pop)
        pb_dorm_transition <- dorm_trans_pb(tmst, dormancy_rate, resuscitation_rate, tmp_pop)
        
        # Pull exact column widths matching matrix splits
        half_cols <- ncol(dorm_transition) / 2
        
        # Active Grower updates
        tmp_pop[[1]] <- tmp_pop[[1]] - dorm_transition[, (half_cols + 1):ncol(dorm_transition)]
        tmp_pop[[1]] <- tmp_pop[[1]] + dorm_transition[, 1:half_cols]
        
        # Dormant updates
        tmp_pop[[2]] <- tmp_pop[[2]] - dorm_transition[, 1:half_cols]
        tmp_pop[[2]] <- tmp_pop[[2]] + dorm_transition[, (half_cols + 1):ncol(dorm_transition)]
        tmp_pop[[2]] <- tmp_pop[[2]] + pb_dorm_transition
        
        # Penalty Box updates
        tmp_pop[[3]] <- tmp_pop[[3]] - pb_dorm_transition
        return(tmp_pop)
      }
      
      birth_routine = function(tmp_pop){
        pop_change      <- births(carb, vbmax, tmst, tmp_pop[[1]])
        pop_change_cost <- births_penalty(cost, carb, vbmax, tmst, tmp_pop[[3]])
        
        tmp_pop[[1]] <- tmp_pop[[1]] + pop_change$pop_change
        tmp_pop[[1]] <- tmp_pop[[1]] + pop_change_cost$grown_individuals_matrix
        tmp_pop[[3]] <- tmp_pop[[3]] - pop_change_cost$grown_individuals_matrix
        return(list(tmp_pop, pop_change$births))
      }
      
      death_routine = function(tmp_pop){
        death_pop_1 <- death(toxin, tmp_pop[[1]], spont_death, tmst)
        death_pop_2 <- death(toxin, tmp_pop[[2]], spont_death, tmst)
        death_pop_3 <- death(toxin, tmp_pop[[3]], spont_death, tmst)
        
        tmp_pop[[1]] <- tmp_pop[[1]] - death_pop_1
        tmp_pop[[2]] <- tmp_pop[[2]] - death_pop_2
        tmp_pop[[3]] <- tmp_pop[[3]] - death_pop_3
        return(list(tmp_pop, sum(death_pop_1, death_pop_2, death_pop_3)))
      }
      
      order_ops = c('induction_routine', 'loss_routine', 'dorm_routine', 'birth_routine', 'death_routine')
      
      # ------------------------------------------------------------
      # 3. Main Operational Loop
      # ------------------------------------------------------------
      while (time < max_time && sum(population[[1]] + population[[2]] + population[[3]]) > 0) {
        step_idx <- step_idx + 1
        
        # Expand index tracking allocations if simulation outlives guess limits safely
        if(step_idx > length(times)) {
          length(carbs)    <- step_idx * 2
          length(tss)      <- step_idx * 2
          length(times)    <- step_idx * 2
          length(tox_list) <- step_idx * 2
          length(pop_list) <- step_idx * 2
          length(bir_dea)  <- step_idx * 2
        }
        
        carbs[step_idx] <- carb
        
        # Toxin Pulsing Logic
        if(time > toxin_times[1] && tox_flag[1] == 0) {
          toxin <- toxin + toxin_strength[1]
          tox_flag[1] <- 1
        } else if(length(toxin_strength) > 1 && time > toxin_times[2] && tox_flag[2] == 0) {
          toxin <- toxin + toxin_strength[2]
          tox_flag[2] <- 1
        }
        
        toxin <- toxin - toxin * tmst * d_r
        if(toxin < 1.2e-11) toxin <- 0
        tox_list[step_idx] <- toxin
        
        # Setup operation list states
        tmp_pop <- list(population[[1]], population[[2]], population[[3]])
        rand_ops <- 1:5#sample(1:5, 5)
        num_births <- 0
        num_deaths <- 0
        
        for(i in rand_ops){
          func <- get(order_ops[i])
          if(i == 4){
            birth_op   <- func(tmp_pop)
            tmp_pop    <- birth_op[[1]]
            num_births <- birth_op[[2]]
          } else if(i == 5){
            death_op   <- func(tmp_pop)
            tmp_pop    <- death_op[[1]]
            num_deaths <- death_op[[2]]
          } else {
            tmp_pop    <- func(tmp_pop)
          }
        }
        
        population[[1]] <- tmp_pop[[1]]
        population[[2]] <- tmp_pop[[2]]
        population[[3]] <- tmp_pop[[3]]
        
        total_actives <- (population[[1]] + population[[3]])
        total_pop     <- cbind(rowSums(population[[2]]), rowSums(total_actives))
        
        # Adaptive Timestep Calculations
        prev_pop_total <- pop_list[step_idx - 1]
        r <- max(.01, abs(sum(total_pop) - prev_pop_total) / prev_pop_total)
        tmst <- allowed_error / r
        
        bir_rate <- vbmax * carb / (km_carb + carb) * tmst
        # Adjusted to track row 1 phenotype threshold natively
        dea_rate <- vmax_death * (toxin^(hill_coef) / ((km_death * 1)^hill_coef + toxin^(hill_coef))) * tmst
        ind_rate <- (vmax_toxin * (toxin / (km_toxin + toxin)) + hetero) * tmst
        probabilities <- 1 - exp(-c(bir_rate, dea_rate, ind_rate))
        
        while(any(probabilities >= 1)){
          tmst <- tmst * .5
          bir_rate <- vbmax * carb / (km_carb + carb) * tmst
          dea_rate <- vmax_death * (toxin^(hill_coef) / ((km_death * 1)^hill_coef + toxin^(hill_coef))) * tmst
          ind_rate <- (vmax_toxin * (toxin / (km_toxin + toxin)) + hetero) * tmst
          probabilities <- 1 - exp(-c(bir_rate, dea_rate, ind_rate))
        }
        
        if(tmst > 10 * init_ts) tmst <- init_ts * 10
        
        # Update Environments and Outputs
        carb <- carbon_util(carb, population, tmst)
        bir_dea[step_idx - 1] <- num_births - num_deaths
        tss[step_idx - 1]     <- tmst
        times[step_idx]       <- time + tmst
        time                  <- time + tmst
        pops[, , step_idx]    <- total_pop
        pop_list[step_idx]    <- sum(total_pop)
      }
      
      # ------------------------------------------------------------
      # 4. Clean Vector Outliers & Trim Traces
      # ------------------------------------------------------------
      pops     <- pops[, , 1:step_idx, drop = FALSE]
      times    <- times[1:step_idx]
      tss      <- tss[1:(step_idx - 1)]
      pop_list <- pop_list[1:step_idx]
      carbs    <- carbs[1:step_idx]
      tox_list <- tox_list[1:step_idx]
      bir_dea  <- bir_dea[1:(step_idx - 1)]
      
      return(list(
        'pops'     = pops,
        'times'    = times,
        'tss'      = tss,
        'pop_list' = pop_list,
        'carbs'    = carbs,
        'tox_list' = tox_list,
        'last_pop' = population,
        'bir_dea'  = bir_dea
      ))
    })
    done = TRUE
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
