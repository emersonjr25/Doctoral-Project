#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity #####
##### Methods: Computational simulation #############
### FUNCTION: MODIFYING INPUT TEMPERATURE - + VARIATION ###

modify_input_temperature <- function(config, data, vars, type_envir = 'random'){
  #browser()
  temp_ancient <- data[["inputs"]][["environments"]][['temp']]
  #temp_ancient <- val[['data']][["inputs"]][["environments"]][['temp']]
  new_temp <- matrix(NA, 
                     nrow = nrow(temp_ancient),
                     ncol = ncol(temp_ancient),
                     dimnames = dimnames(temp_ancient))
  temp <- 0
  if(type_envir == 'random'){
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
  } else if (type_envir == 'stable_low') {
    temp_sequence <- seq(mean(temp_ancient[!is.na(temp_ancient)]), 0, -0.001)
    temp_sequence2 <- seq(0, 1, 0.001)
    temp_sequence3 <- seq(1, mean(temp_ancient[!is.na(temp_ancient)]), -0.001)
    temp_sequence_final <- round(c(temp_sequence, temp_sequence2, temp_sequence3), 3)
    temp_sequence_final <- rep_len(temp_sequence_final, ncol(temp_ancient))
    for(i in 1:ncol(temp_ancient)){
      temp <- temp_ancient[, i]
      temp[!is.na(temp)] <- sapply(temp[!is.na(temp)], function(x){
        x <- temp_sequence_final[i]
      })
      temp[!is.na(temp) & temp >= 1] <- temp_sequence_final[i]
      temp[!is.na(temp) & temp <= 0] <- temp_sequence_final[i]
      new_temp[, i] <- temp
      temp <- 0
    } 
    data[["inputs"]][["environments"]][['temp']] <- new_temp
    return(list(config = config, data = data, vars = vars))
  } else if (type_envir == 'stable_fast'){
    temp_sequence <- seq(mean(temp_ancient[!is.na(temp_ancient)]), 0, -0.01)
    temp_sequence2 <- seq(0, 1, 0.01)
    temp_sequence3 <- seq(1, mean(temp_ancient[!is.na(temp_ancient)]), -0.01)
    temp_sequence_final <- round(c(temp_sequence, temp_sequence2, temp_sequence3), 3)
    temp_sequence_final <- rep_len(temp_sequence_final, ncol(temp_ancient))
    for(i in 1:ncol(temp_ancient)){
      temp <- temp_ancient[, i]
      temp[!is.na(temp)] <- sapply(temp[!is.na(temp)], function(x){
        x <-  temp_sequence_final[i]
      })
      temp[!is.na(temp) & temp >= 1] <- temp_sequence_final[i]
      temp[!is.na(temp) & temp <= 0] <- temp_sequence_final[i]
      new_temp[, i] <- temp
      temp <- 0
    } 
    data[["inputs"]][["environments"]][['temp']] <- new_temp
    return(list(config = config, data = data, vars = vars))
  } else {
    message(paste0('Error in input ', type_envir, 
                   ' : temperature does not change, chose a correct option - random, stable_low, or stable_fast'))
  }
}