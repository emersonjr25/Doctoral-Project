#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity #####
##### Methods: Computational simulation #############
##### FUNCTION: APPLY ECOLOGY WITH PLASTICITY #####

config$gen3sis$ecology$apply_ecology <- function(abundance, traits, landscape, 
                                                 config, plast_value) {
  abundance_scale = 1.12
  abundance_threshold = 1
  #abundance threshold
  survive <- abundance>=abundance_threshold
  abundance[!survive] <- 0
  
  # Modification for plasticity:
  # Traits are now subtracted from a variance of optima
  # Such optima are defined by the parameter plast
  # Higher plast equals more plasticity
  # plast = 0, means no plasticity
  
  plasticity <- function(x, plast) {
    return(seq(x - (x * plast), x + (x * plast), 0.01))
  }
  
  plasticity2 <- function(x, land) {
    min(abs(x - land))
  }
  
  traits_sub <- lapply(traits[, 'temp'], plasticity, plast_value)
  traits_sub2 <- mapply(plasticity2, traits_sub, landscape[,'temp'])
  abundance <- ((1 - traits_sub2)*abundance_scale)*as.numeric(survive)
  
  #abundance threshold
  abundance[abundance<abundance_threshold] <- 0
  return(abundance)
}