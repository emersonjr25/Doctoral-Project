#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity #####
##### Methods: Computational simulation #############
##### FUNCTION: APPLY EVOLUTION WITHOUT CLUSTERS #####

# THIS FUNCTION WAS ONLY USE TO TEST MODEL PREMISES #

# config$gen3sis$mutation$apply_evolution <- function(species, cluster_indices, landscape, config) {
# 
#   trait_evolutionary_power <- 0.001
#   traits <- species[["traits"]]
#   #mutations
#   mutation_deltas <-rnorm(length(traits[, "temp"]), mean=0, sd=trait_evolutionary_power)
#   traits[, "temp"] <- traits[, "temp"] + mutation_deltas
#   return(traits)
# }