library(ggplot2)
library(metR)
library(geomtextpath)
# 1. Load data
file_path <- '~/Documents/GitHub/publication_model/dormancy_rate/dormancy_rate 200 1 2026-08-24 12:28:29.694037'
results_1 <- read.csv(file = file_path, header = FALSE)

# 2. Setup seed space
set.seed(1)
reps <- 200
total_seeds <- reps * length(param_space) * length(res_param_space)
seeds <- sample.int(10^9, total_seeds)

# 3. Extract columns explicitly (Ensuring numeric conversion)
dorm_param <- as.numeric(results_1[, 2])
rep_val    <- as.numeric(results_1[, 4])
seed_val   <- as.integer(round(as.numeric(results_1[, 6]))) # Double check if column 6 holds the seeds

# 4. Vectorized match to find 'ind' for all rows at once
# match() returns NA instead of crashing if a seed isn't found
inds <- match(seed_val, seeds)

# 5. Vectorized match to find 'l' for all rows
l <- match(round(dorm_param,9), round(param_space,9))

# 6. Calculate 'k' for all rows instantly
ks <- as.integer(round(((inds - rep_val - (reps * length(res_param_space) * (l - 1))) / reps) + 1))

# 7. Safety check: Check if any seeds failed to map
if (any(is.na(ks))) {
  warning(paste(sum(is.na(ks)), "rows could not find a matching seed! Check if column 6 is correct."))
}

# 8. Directly map values back to column 3
results_1[, 3] <- res_param_space[ks]
results = aggregate(results_1[,1],by=results_1[,c(2,3)],FUN = mean)
colnames(results) = c("Dormancy","Resuscitation",'Survival')

frac_pops = data.frame(matrix(nrow = 0,ncol=3))
frac_pops_param_space = 10^seq(-9,-1,.15)
frac_pops_res_param_space = 10^seq(-13,-.1,.15)
colnames(frac_pops)= c("Dormancy",'Resuscitation','Fraction')
for(i in param_space){
  for(k in res_param_space){
    pop = vbmax_equibpop_determination(i,k,spont_death)[1]
    frac_pop = sum(pop[[1]][[2]])/sum(pop[[1]][[1]]+pop[[1]][[2]])
    if(frac_pop > 1){
      frac_pop = NA
    } else if(frac_pop < 0){
      frac_pop = NA
    }
    new_row = c(i,k,frac_pop)
    frac_pops = rbind(frac_pops,new_row)
  }
}
colnames(frac_pops)= c("Dormancy",'Resuscitation','Fraction')
ggplot(data = results, aes(x = as.numeric(Dormancy), y = as.numeric(Resuscitation), z = as.numeric(Survival)), bins = 10) +
  geom_contour_filled(breaks = c(seq(0, 1.1, .1))) +
  xlab("Dormancy") +
  ylab("Resuscitation") +
  geom_textcontour(data = subset(frac_pops, Fraction > 0),
                   aes(
                     x = Dormancy, 
                     y = Resuscitation, 
                     z = Fraction,
                     label = after_stat(paste0("10^", log10(level)))
                   ),
                   breaks = c(.001, .0001, .00001, .000001), 
                   color = 'black',
                   size = 5,           
                   gap = TRUE,      
                   parse = TRUE     
  ) +
  
  scale_x_log10() +
  scale_y_log10() +
  scale_size_binned(n.breaks = 25) + 
  geom_point(aes(x = dormancy_rate, y = resuscitation_rate), color = 'black', size = 3) + 
  theme_classic()+
  theme(text = element_text(size = 25))
mixed = results
mixed$Fraction = NA
for(i in 1:nrow(frac_pops)){
  mixed$Fraction[which(mixed$Dormancy == frac_pops[i,1]&mixed$Resuscitation == frac_pops[i,2])] = frac_pops$Fraction[i]
}
