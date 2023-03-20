#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Computational simulation #############
##### FUNCTION: ALTERING INITIAL TRAIT VALUES OF THE SPECIES #####

# config$gen3sis$initialization$create_ancestor_species <- function(landscape, config) {
# 
#   range <- c(-180, 180, -90, 90)
#   co <- landscape$coordinates
#   selection <- co[, "x"] >= range[1] &
#     co[, "x"] <= range[2] &
#     co[, "y"] >= range[3] &
#     co[, "y"] <= range[4]
#   initial_cells <- rownames(co)[selection]
#   new_species <- create_species(initial_cells, config)
#   #set local adaptation to max optimal temp equals local temp
#   new_species$traits[ , "temp"] <- rnorm(length(landscape$environment[,"temp"]), 0.5, 0.5)
#   new_species$traits[ , "temp"][new_species$traits[ , "temp"] <= 0]  <- runif(1, 0, 0.3)
#   new_species$traits[ , "temp"][new_species$traits[ , "temp"] >= 1]  <- runif(1, 0.7, 1)
#   new_species$traits[ , "dispersal"] <- 1
# 
#   return(list(new_species))
# }
# 