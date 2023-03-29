#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity #####
##### Methods: Computational simulation #############
#### PREPARATION IN CONFIG - SPECIES AND SYSTEM ####

config$gen3sis$general$start_time <- 20

config$gen3sis$general$end_time <- 1

timesteps_total <- length(config$gen3sis$general$start_time:config$gen3sis$general$end_time)

config$gen3sis$general$max_number_of_species <- 50000
config$gen3sis$general$max_number_of_coexisting_species <- 200000

rep <- 1
quantitity_rep <- 1

config$gen3sis$general$end_of_timestep_observer <- function(data, vars, config){
  save_traits()
}

config$gen3sis$speciation$divergence_threshold <- 10


plasti <- c(0, 0.25, 1)

pos <- 0
pos2 <- 0

finalresult <- data.frame(plasticity = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          replications = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          speciation = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          extinction = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          diversif = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          traitevolution = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          timesimulation = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0))
