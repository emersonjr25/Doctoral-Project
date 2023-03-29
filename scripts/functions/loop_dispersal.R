#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity #####
##### Methods: Computational simulation #############
##### FUNCTION: DISPERSION WITH LIMITS PER SPECIES #####

loop_dispersal2 <- function (config, data, vars) {
  if (config$gen3sis$general$verbose >= 3) {
    cat(paste("entering dispersal module \n"))
  }
  data$all_species <- lapply(data$all_species, disperse2, data$landscape, 
                             data$distance_matrix, config)
  if (config$gen3sis$general$verbose >= 3) {
    cat(paste("exiting dispersal module \n"))
  }
  return(list(config = config, data = data, vars = vars))
}