#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Computational simulation #############
### FUNCTION: MODIFYING INPUT TEMPERATURE - + VARIATION ###

modify_input_temperature <- function(config, data, vars){
  temp_ancient <- data[["inputs"]][["environments"]][['temp']]
  new_temp <- matrix(NA, 
                     nrow = nrow(temp_ancient),
                     ncol = ncol(temp_ancient),
                     dimnames = dimnames(temp_ancient))
  temp <- 0
  for(i in 1:ncol(temp_ancient)){
    temp <- temp_ancient[, i]
    temp[!is.na(temp)] <- sapply(temp[!is.na(temp)], function(x){
      x <- abs(rnorm(1, x, sd = 0.5))
    })
    temp[!is.na(temp) & temp >= 1] <- 1
    temp[!is.na(temp) & temp <= 0] <- 0.1
    new_temp[, i] <- temp
    temp <- 0
  }
  data[["inputs"]][["environments"]][['temp']] <- new_temp
  return(list(config = config, data = data, vars = vars))
}