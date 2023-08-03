#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity #####
##### Methods: Computational simulation #############
### FUNCTION: MODIFYING INPUT TEMPERATURE - + VARIATION ###

modify_input_temperature <- function(config, data, vars, type_envir = 'random'){
  temp_ancient <- data[["inputs"]][["environments"]][['temp']]
  X <- data[["inputs"]][["coordinates"]][, 1]
  
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
    
    ### creating the temperature temporal sequence from average ###
    mean_temp <- round(mean(temp_ancient[!is.na(temp_ancient)]), 3)
    temp_sequence1 <- seq(mean_temp, 0, -0.001)
    temp_sequence2 <- seq(0, 1, 0.001)
    temp_sequence3 <- seq(1, mean_temp, -0.001)
    temp_without_right_len <- round(c(temp_sequence1, temp_sequence2, temp_sequence3), 3)
    temp_sequence_ended <- rep_len(temp_without_right_len, ncol(temp_ancient))
    
    ### organizing in way that function setup_landscape read ###
    #good temperature begin in timestep that begin simulation #
    temp_vector_final <- rep(0, length(temp_sequence_ended))
    pos_mean <- which(temp_sequence_ended == mean_temp)[1]
    first_sequence_part <- temp_sequence_ended[pos_mean:config$gen3sis$general$start_time]
    second_sequence_part <- rep(mean_temp, (length(temp_sequence_ended) - config$gen3sis$general$start_time))
    temp_vector_final[config$gen3sis$general$start_time:pos_mean] <- first_sequence_part
    temp_vector_final[(config$gen3sis$general$start_time + 1):length(temp_sequence_ended)] <- second_sequence_part
    ##### longitudinal variation ####
    names_longitude <- seq(1, 800, 1)
    for(i in 1:ncol(temp_ancient)){
      temp <- temp_ancient[, i]
      longi_below <- seq(temp_vector_final[i], temp_vector_final[i] - 0.35, -0.000877)
      longi_below[longi_below < 0 ] <- 0
      longi_below_increasing <- sort(longi_below, decreasing = FALSE)
      longitudinal_variation <- c(longi_below_increasing, longi_below)
      if(length(longitudinal_variation) < 800){
        difference <- abs(length(longitudinal_variation) - length(names_longitude))
        values_to_add_vector <- rep(longitudinal_variation[length(longitudinal_variation)], difference)
        longitudinal_variation <- c(longitudinal_variation, values_to_add_vector)
      } 
      if(length(longitudinal_variation) > 800){
        longitudinal_variation <- longitudinal_variation[1:800]
      }
      names(longitudinal_variation) <- names_longitude
      temp[!is.na(temp)] <- longitudinal_variation[names(longitudinal_variation) %in% names(temp[!is.na(temp)])]
      ### latitudinal variation ###
      unique_X <- unique(X)
      for(j in seq_along(unique_X)){
        latitudinal_variation <- names(X[unique_X[j] == X])
        positions_to_change <- which(latitudinal_variation %in% names(temp[!is.na(temp)]))
        latitudinal_position <- latitudinal_variation[positions_to_change]
        temp[latitudinal_position] <- sapply(temp[latitudinal_position], FUN = function(x){
          x <- round(abs(rnorm(1, mean=x, sd=0.01)), 4)
        })
      }
      temp <- unlist(temp)
      new_temp[, i] <- temp
    } 
    data[["inputs"]][["environments"]][['temp']] <- new_temp
    return(list(config = config, data = data, vars = vars))
  } else if (type_envir == 'stable_fast'){
    
    ### creating the temperature temporal sequence from average ###
    mean_temp <- round(mean(temp_ancient[!is.na(temp_ancient)]), 3)
    temp_sequence1 <- seq(mean_temp, 0, -0.002)
    temp_sequence2 <- seq(0, 1, 0.002)
    temp_sequence3 <- seq(1, mean_temp, -0.002)
    temp_without_right_len <- round(c(temp_sequence1, temp_sequence2, temp_sequence3), 3)
    temp_sequence_ended <- rep_len(temp_without_right_len, ncol(temp_ancient))
    
    ### organizing in way that function setup_landscape read ###
    #good temperature begin in timestep that begin simulation #
    temp_vector_final <- rep(0, length(temp_sequence_ended))
    pos_mean <- which(temp_sequence_ended == mean_temp)[1]
    first_sequence_part <- temp_sequence_ended[pos_mean:config$gen3sis$general$start_time]
    second_sequence_part <- rep(mean_temp, (length(temp_sequence_ended) - config$gen3sis$general$start_time))
    temp_vector_final[config$gen3sis$general$start_time:pos_mean] <- first_sequence_part
    temp_vector_final[(config$gen3sis$general$start_time + 1):length(temp_sequence_ended)] <- second_sequence_part
    
    ##### longitudinal variation ####
    names_longitude <- seq(1, 800, 1)
    for(i in 1:ncol(temp_ancient)){
      temp <- temp_ancient[, i]
      longi_below <- seq(temp_vector_final[i], temp_vector_final[i] - 0.35, -0.000877)
      longi_below[longi_below < 0 ] <- 0
      longi_below_increasing <- sort(longi_below, decreasing = FALSE)
      longitudinal_variation <- c(longi_below_increasing, longi_below)
      if(length(longitudinal_variation) < 800){
        difference <- abs(length(longitudinal_variation) - length(names_longitude))
        values_to_add_vector <- rep(longitudinal_variation[length(longitudinal_variation)], difference)
        longitudinal_variation <- c(longitudinal_variation, values_to_add_vector)
      } 
      if(length(longitudinal_variation) > 800){
        longitudinal_variation <- longitudinal_variation[1:800]
      }
      names(longitudinal_variation) <- names_longitude
      temp[!is.na(temp)] <- longitudinal_variation[names(longitudinal_variation) %in% names(temp[!is.na(temp)])]
      
      ### latitudinal variation ###
      unique_X <- unique(X)
      for(j in seq_along(unique_X)){
        latitudinal_variation <- names(X[unique_X[j] == X])
        positions_to_change <- which(latitudinal_variation %in% names(temp[!is.na(temp)]))
        latitudinal_position <- latitudinal_variation[positions_to_change]
        temp[latitudinal_position] <- sapply(temp[latitudinal_position], FUN = function(x){
          x <- round(abs(rnorm(1, mean=x, sd=0.01)), 4)
        })
      }
      temp <- unlist(temp)
      new_temp[, i] <- temp
    } 
    data[["inputs"]][["environments"]][['temp']] <- new_temp
    return(list(config = config, data = data, vars = vars))
  } else {
    message(paste0('Error in input ', type_envir, 
                   ' : temperature does not change, chose a correct option - random, stable_low, or stable_fast'))
  }
}