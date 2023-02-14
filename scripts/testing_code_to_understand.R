library(gen3sis)

#######################################################
##### FIRST TEST ######
#######################################################
#without loop 
plasti <- 0.1
trait <- list(a = 1:5)
result <- vector("list", length = length(plasti))

testing <- function(data, plasti){
  plasticity <- function(x, plast) {
    return(seq(x - (x * plast), x + (x * plast), 0.01))
  }
  traits_sub <- lapply(data[[1]], plasticity, plasti)
  traits_sub
}
result <- testing(trait, plasti)

#######################################################
##### with loop ######

plasti <- seq(0.1, 1, 0.1)
trait <- list(a = 1:5)
result <- vector("list", length = length(plasti))

for (i in 1:length(plasti)){
  testing <- function(data, plasti){
    plasticity <- function(x, plast) {
      return(seq(x - (x * plast), x + (x * plast), 0.01))
    }
    traits_sub <- lapply(data[[1]], plasticity, plasti[[i]])
    traits_sub
  }
  result[[i]] <- testing(trait, plasti)
}

#######################################################
#### SECOND TEST ####
#######################################################

plasti <- 0.1
trait <- list(a = 1:5)
trait_sub <- vector("list", length = length(plasti))
local_environment <- list(6:10)
#abundance <- rep(100, 10)

testing <- function(data, plasti){
  plasticity <- function(x, plast) {
    return(seq(x - (x * plast), x + (x * plast), 0.01))
  }
  traits_sub <- lapply(data[[1]], plasticity, plasti)
  traits_sub
}
trait_sub <- testing(trait, plasti)

plasticity2 <- function(x, land) {
  min(abs(x - land))
}

plasticity2(trait_sub[[1]], local_environment[[1]][[1]])

traits_sub2 <- mapply(plasticity2, trait_sub, land = local_environment[[1]])




config <- function(abundance, trait, environment) {
  
  plasticity <- function(x, plast) {
    return(seq(x - (x * plast), x + (x * plast), 0.01))
  }
  
  plasticity2 <- function(x, land) {
    min(abs(x - land))
  }
  
  traits_sub <- lapply(trait[[1]], plasticity, 0.1)
  traits_sub2 <- sapply(traits_sub, plasticity2, land = environment[[1]])
  abundance <- ((1 - traits_sub2)*10)
  return(abundance)
}  

NEW_abd <- config(abundance, trait, local_environment)

config$gen3sis$ecology$apply_ecology <- function(abundance, traits, landscape, config, plasti) {
  abundance_scale = 10
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
  traits_sub <- lapply(traits[, 'temp'], plasticity, 0.1)
  traits_sub2 <- sapply(traits_sub, plasticity2, land = landscape[,'temp'])
  abundance <- ((1 - traits_sub2)*abundance_scale)*as.numeric(survive)
  
  
  #abundance threshold
  abundance[abundance<abundance_threshold] <- 0
  k <- ((landscape[,'area']*(landscape[,'arid']+0.1)*(landscape[,'temp']+0.1))
        *abundance_scale^2)
  total_ab <- sum(abundance)
  subtract <- total_ab-k
  if (subtract > 0) {
    # print(paste("should:", k, "is:", total_ab, "DIFF:", round(subtract,0) ))
    while (total_ab>k){
      alive <- abundance>0
      loose <- sample(1:length(abundance[alive]),1)
      abundance[alive][loose] <- abundance[alive][loose]-1
      total_ab <- sum(abundance)
    }
    #set negative abundances to zero
    abundance[!alive] <- 0
  }
  return(abundance)
}  

loop_ecology2 <- function (config, data, vars, plasti) {
  if (config$gen3sis$general$verbose >= 3) {
    cat(paste("entering ecology module @ time", vars$ti, 
              "\n"))
  }
  all_cells <- rownames(data$landscape$environment)
  all_species_presence <- do.call(cbind, lapply(data$all_species, 
                                                FUN = function(sp) {
                                                  all_cells %in% names(sp$abundance)
                                                }))
  rownames(all_species_presence) <- all_cells
  occupied_cells <- rownames(all_species_presence)[rowSums(all_species_presence) > 
                                                     0]
  for (cell in occupied_cells) {
    local_environment = data$landscape[["environment"]][cell, 
                                                        , drop = FALSE]
    coo_sp <- which(all_species_presence[cell, ])
    traits <- matrix(nrow = length(coo_sp), ncol = length(config$gen3sis$general$trait_names))
    abundance <- numeric(length(coo_sp))
    colnames(traits) <- config$gen3sis$general$trait_names
    i <- 1
    for (spi in coo_sp) {
      traits[i, ] <- data$all_species[[spi]][["traits"]][cell, 
                                                         config$gen3sis$general$trait_names]
      abundance[i] <- data$all_species[[spi]][["abundance"]][cell]
      i <- i + 1
    }
    max_n_sp_idi <- config$gen3sis$general$max_number_of_coexisting_species
    if (length(coo_sp) > max_n_sp_idi) {
      vars$flag <- "max_number_coexisting_species"
      paste0("Maximum number of species per cell (i.e. max_n_sp_idi) reached. Specifically ", 
             length(coo_sp), "(>", max_n_sp_idi, ") species @ t", 
             vars$ti, " idi", cell)
      return(list(config = config, data = data, vars = vars))
    }
    rownames(traits) <- coo_sp
    names(abundance) <- coo_sp
    NEW_abd <- config$gen3sis$ecology$apply_ecology(abundance, 
                                                    traits, local_environment, config, plasti)
    names(NEW_abd) <- coo_sp
    shalldie <- NEW_abd == 0
    for (spi in coo_sp) {
      data$all_species[[spi]][["abundance"]][cell] <- NEW_abd[[toString(spi)]]
    }
    die_sure <- as.integer(names(NEW_abd)[NEW_abd == 0])
    if (length(die_sure) > 0) {
      chars <- as.character(die_sure)
    }
  }
  species_list <- list()
  for (species in data$all_species) {
    cells <- names(species[["abundance"]])[species[["abundance"]] != 
                                             0]
    updated_species <- limit_species_to_cells(species, cells)
    species_list <- append(species_list, list(updated_species))
  }
  data$all_species <- species_list
  if (config$gen3sis$general$verbose >= 3) {
    cat(paste("exiting ecology module @ time", vars$ti, 
              "\n"))
  }
  return(list(config = config, data = data, vars = vars))
}


##### ANCIENT WAY #####
plasti <- seq(0.1, 1, 0.1)

for(p in 1:length(plasti)){
  
  config$gen3sis$ecology$apply_ecology <- function(abundance, traits, landscape, config) {
    abundance_scale = 10
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
    traits_sub <- lapply(traits[, 'temp'], plasticity, plast = plasti[p])
    traits_sub2 <- sapply(traits_sub, plasticity2, land = landscape[,'temp'])
    abundance <- ((1 - traits_sub2)*abundance_scale)*as.numeric(survive)
    
    
    #abundance threshold
    abundance[abundance<abundance_threshold] <- 0
    k <- ((landscape[,'area']*(landscape[,'arid']+0.1)*(landscape[,'temp']+0.1))
          *abundance_scale^2)
    total_ab <- sum(abundance)
    subtract <- total_ab-k
    if (subtract > 0) {
      # print(paste("should:", k, "is:", total_ab, "DIFF:", round(subtract,0) ))
      while (total_ab>k){
        alive <- abundance>0
        loose <- sample(1:length(abundance[alive]),1)
        abundance[alive][loose] <- abundance[alive][loose]-1
        total_ab <- sum(abundance)
      }
      #set negative abundances to zero
      abundance[!alive] <- 0
    }
    return(abundance)
  }  
}