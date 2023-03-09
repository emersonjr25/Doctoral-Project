#####################################################
########### PLASTICITY AND EVOLUTION STUDY ##########
#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Computational simulation #############
##### Script to input, modify and run model #########

#### PACKAGES ####
library(gen3sis)
library(here)
library(dplyr)

################## SIMULATION #########################
#### CARRYING CONFIGURATIONS AND PATHS ####
datapath <- here("data/raw/WorldCenter")
attach(loadNamespace('gen3sis'), name = 'gen3sis_all')
config = file.path(datapath, "config/config_worldcenter.R")
landscape = file.path(datapath, "landscape_new")
output_directory = NA
timestep_restart = NA
save_state = NA
call_observer = "all"
enable_gc = FALSE
verbose = 1

system_time_start <- Sys.time()
directories <- prepare_directories(
  config_file = config,
  input_directory = landscape,
  output_directory = output_directory
)

if (is.na(config)[1]) {
  stop("please provide either a config file or a config object")
} else if (class(config) == "gen3sis_config") {
  config[["directories"]] <- directories
} else if (class(config) == "character") {
  file.copy(config, directories$output)
  config <- create_input_config(config_file = config)
  config[["directories"]] <- directories
} else {
  stop("this is not a known config, please provide either a config file or a config object")
}
if (!verify_config(config)) {
  stop("config verification failed")
}

#### MODIFICATIONS IN CONFIG - SPECIES AND SYSTEM ####
config$gen3sis$general$start_time <- 100

config$gen3sis$general$end_time <- 1

timesteps_total <- length(config$gen3sis$general$start_time:config$gen3sis$general$end_time)

config$gen3sis$general$max_number_of_species <- 10000

config$gen3sis$general$end_of_timestep_observer <- function(data, vars, config){
  save_traits()
}

rep <- 1

plasti <- c(0.25)

pos <- 0
pos2 <- 0

finalresult <- data.frame(plasticidade = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          replications = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          speciation = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          extinction = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          diversif = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          traitevolution = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          timestep = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          timesimulation = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          speciestotal = runif(rep * length(plasti) * timesteps_total, 0, 0))

# ALTERING INITIAL TRAIT VALUES #

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
# config$gen3sis$mutation$apply_evolution <- function(species, cluster_indices, landscape, config) {
# 
#   trait_evolutionary_power <- 0.001
#   traits <- species[["traits"]]
#   cells <- rownames(traits)
#   #mutations
#   mutation_deltas <-rnorm(length(traits[, "temp"]), mean=0, sd=trait_evolutionary_power)
#   traits[, "temp"] <- traits[, "temp"] + mutation_deltas
#   return(traits)
# }

# ADD PLASTICITY AND MODIFICATIONS IN ECOLOGY #

config$gen3sis$ecology$apply_ecology <- function(abundance, traits, landscape, config, plasticidade) {
  #browser()
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

  traits_sub <- lapply(traits[, 'temp'], plasticity, plasticidade)
  traits_sub2 <- mapply(plasticity2, traits_sub, landscape[,'temp'])
  abundance <- ((1 - traits_sub2)*abundance_scale)*as.numeric(survive)

  #abundance threshold
  abundance[abundance<abundance_threshold] <- 0
  # k <- ((landscape[,'area']*(landscape[,'arid']+0.1)*(landscape[,'temp']+0.1))
  #       *abundance_scale^2)
  # total_ab <- sum(abundance)
  # subtract <- total_ab-k
  # if (subtract > 0) {
  #   # print(paste("should:", k, "is:", total_ab, "DIFF:", round(subtract,0) ))
  #   while (total_ab>k){
  #     alive <- abundance>0
  #     loose <- sample(1:length(abundance[alive]),1)
  #     abundance[alive][loose] <- abundance[alive][loose]-1
  #     total_ab <- sum(abundance)
  #   }
  #   #set negative abundances to zero
  #   abundance[!alive] <- 0
  # }
  return(abundance)
}

loop_ecology2 <- function (config, data, vars, plasticidade) {
   #if(ti == 45){
  # browser()
  #     }
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
    NEW_abd <- config$gen3sis$ecology$apply_ecology(abundance, traits, local_environment, config, plasticidade)
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

# MODIFICATION IN INPUT TEMPERATURE - MORE VARIATION #

modify_input_temperature <- function(config, data, vars, seed = 1){
  set.seed(seed)
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

################### FOR TESTS ####################################

disperse2 <- function (species, landscape, distance_matrix, config) 
{
 # if (ti == 97){
 #   browser()
 # }
  if (!length(species[["abundance"]])) {
    return(species)
  }
  presence_spi_ti <- names(species[["abundance"]])
  all_cells <- rownames(landscape$coordinates)
  free_cells <- all_cells[!(all_cells %in% presence_spi_ti)]
  if(length(free_cells) >= 1){
    if(length(free_cells) >= 100){
      free_cells <- free_cells %>%
                    sample(100) %>% as.numeric() %>%
                    sort() %>% as.character()
    } else {
      half <- round(length(free_cells) - (length(free_cells) * 0.5))
      free_cells <- free_cells %>%
                   sample(half) %>% as.numeric() %>%
                   sort() %>% as.character()
    }
  }
  num_draws <- length(free_cells) * length(presence_spi_ti)
  r_disp <- config$gen3sis$dispersal$get_dispersal_values(num_draws, 
                                                          species, landscape, config)
  geo_disp <- distance_matrix[presence_spi_ti, free_cells, 
                              drop = FALSE]
  geo_disp <- geo_disp <= r_disp
  colonized <- rep(FALSE, length(all_cells))
  names(colonized) <- all_cells
  colonized[free_cells] <- apply(geo_disp, 2, any)
  tep_occ_id <- all_cells[colonized]
  if (length(tep_occ_id) > 0) {
    dest <- which(colonized == TRUE)
   # if(length(dest) >= 5){
    #  dest <- sample(dest, 3)
   # }
    if (length(presence_spi_ti) == 1) {
      orig <- rep(1, length(dest))
    }
    else {
      orig <- apply(geo_disp[, colonized[free_cells], 
                             drop = FALSE], 2, function(x) {
                               a <- which(x)
                               ifelse(length(a) > 1, sample(a, 1), a)
                             })
    }
    orig <- as.numeric(presence_spi_ti[orig])
    dest <- tep_occ_id
    species <- disperse_species(species, as.character(orig), 
                                dest, config)
  }
   #if (config$gen3sis$general$verbose >= 1){
  #   cat(paste('free cells', length(free_cells)))
  # }
  return(species)
}

loop_dispersal2 <- function (config, data, vars) {
#  if (ti == 97){
#    browser()
#  }
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


#### REMOVING TRAITS OF ANTERIOR SIMULATIONS ####
caminho <- here("data", "raw", "WorldCenter", "output", "config_worldcenter", "traits")
listfiles <- list.files("data/raw/WorldCenter/output/config_worldcenter/traits")
filestoread <- length(list.files("data/raw/WorldCenter/output/config_worldcenter/traits"))
cam <- 0
for(l in 1:filestoread){
  cam[l] <- caminho
}
camatualizado <- 0
for(k in 1:length(cam)){
  camatualizado[[k]] <- paste(cam[k], listfiles[k], sep = "/", collapse = "--")
}
file.remove(camatualizado)

for(p in 1:length(plasti)){
  
  for(r in 1:rep){
    val <- list(data = list(),
                vars = list(),
                config = config)
    val$config <- complete_config(val$config)
    val$config$gen3sis$general$verbose <- verbose
    val <- setup_inputs(val$config, val$data, val$vars)
    val <- setup_variables(val$config, val$data, val$vars)
    val <- setup_landscape(val$config, val$data, val$vars)
    val$data$landscape$id <- val$data$landscape$id + 1
    val <- init_attribute_ancestor_distribution(val$config,
                                                val$data, val$vars)
    
    val <- init_simulation(val$config, val$data, val$vars)
    val <- init_summary_statistics(val$data, val$vars, val$config)
    if (is.na(call_observer)) {
      save_steps <- c(val$config$gen3sis$general$start_time,
                      val$config$gen3sis$general$end_time)
    } else if (call_observer == "all") {
      save_steps <-
        val$config$gen3sis$general$start_time:val$config$gen3sis$general$end_time
    } else {
      steps <- as.integer(call_observer) + 2
      save_steps <-
        ceiling(
          seq(
            val$config$gen3sis$general$start_time,
            val$config$gen3sis$general$end_time,
            length.out = steps
          )
        )
    }
    val$vars$save_steps <- save_steps
    val$vars$steps <-
      val$config$gen3sis$general$start_time:val$config$gen3sis$general$end_time
    if (!is.na(timestep_restart)) {
      val <- restore_state(val, timestep_restart)
    }
    pos2 <- 0
    val <- modify_input_temperature(val$config, val$data, val$vars)
    
    for (ti in val$vars$steps) {
       #if (ti == 48) {
       #  break 
      # }
      val$vars$n_new_sp_ti <- 0
      val$vars$n_ext_sp_ti <- 0
      val$vars$n_sp_added_ti <- 0
      val$vars$ti <- ti
      if (verbose >= 2) {
        cat("loop setup \n")
      }
      
      val <- setup_landscape(val$config, val$data, val$vars)
      val <- restrict_species(val$config, val$data, val$vars)
      val <- setup_distance_matrix(val$config, val$data, val$vars)
      if (verbose >= 2) {
        cat("speciation \n")
      }
      val <- loop_speciation(val$config, val$data, val$vars)
      val <- update1.n_sp.all_geo_sp_ti(val$config, val$data,
                                        val$vars)
      val <- update2.n_sp_alive.geo_sp_ti(val$config, val$data,
                                          val$vars)
      if (verbose >= 2) {
        cat("dispersal \n")
      }
      val <- loop_dispersal2(val$config, val$data, val$vars)
      if (verbose >= 2) {
        cat("evolution \n")
      }
      #val <- loop_evolution(val$config, val$data, val$vars)
      if (verbose >= 2) {
        cat("ecology \n")
      }
      if (verbose >= 2) {
        cat("end of loop updates \n")
      }
      val$vars$n_sp_alive <- sum(sapply(val$data$all_species,
                                        function(sp) {
                                          ifelse(length(sp[["abundance"]]), 1, 0)
                                        }))
      val$vars$n_sp <- length(val$data$all_species)
      val <- update_loop_steps_variable(val$config, val$data,
                                        val$vars)
      if (val$vars$ti %in% val$vars$save_steps) {
        call_main_observer(val$data, val$vars, val$config)
      }
      val <- update_summary_statistics(val$data, val$vars,
                                       val$config)
      save_val(val, save_state)
      if (verbose >= 1) {
        cat(
          "step =",
          ti,
          ", species alive =",
          val$vars$n_sp_alive,
          ", species total =",
          val$vars$n_sp,
          "\n"
        )
      }
      if (val$vars$n_sp_alive >= val$config$gen3sis$general$max_number_of_species) {
        val$vars$flag <- "max_number_species"
        print("max number of species reached, breaking loop")
        break
      }
      val <- loop_ecology2(val$config, val$data, val$vars, plasti[p])
      
      if (verbose >= 0 & val$vars$flag == "OK") {
        cat("Simulation finished. All OK \n")
      } else if (verbose >= 0 & val$vars$flag == "max_number_species") {
        cat("Simulation finished. Early abort due to exceeding max number of species")
      } else if (verbose >= 0 &
                 val$vars$flag == "max_number_coexisting_species") {
        cat("Simulation finished. Early abort due to exceeding max number of co-occuring species")
      }
      val <- update.phylo(val$config, val$data, val$vars)
      
      system_time_stop <- Sys.time()
      total_runtime <- difftime(system_time_stop, system_time_start,
                                units = "hours")[[1]]
      #write_runtime_statisitics(val$data, val$vars, val$config,
      #                          total_runtime)
      sgen3sis <- make_summary(val$config, val$data, val$vars,
                               total_runtime, save_file = FALSE)
      #plot_summary(sgen3sis)
      
      ################## TRAIT EVOLUTION #####################
      
      ##### PATHWAY TO TRAITS DATA ####
      
      caminho <- here("data", "raw", "WorldCenter", "output", "config_worldcenter", "traits")
      
      listfiles <- list.files("data/raw/WorldCenter/output/config_worldcenter/traits")
      
      filestoread <- length(list.files("data/raw/WorldCenter/output/config_worldcenter/traits"))
      
      
      #### Organizing selection of trait files ####
      cam <- 0
      for(l in 1:filestoread){
        cam[l] <- caminho
      }
      
      camatualizado <- 0
      for(k in 1:length(cam)){
        camatualizado[[k]] <- paste(cam[k], listfiles[k], sep = "/", collapse = "--")
      }
      
      # SORT #
      camatualizado <- camatualizado[order(as.numeric(gsub("[^0-9]+", "", camatualizado)), decreasing = TRUE)]
      cam_length <- length(camatualizado)
      if(cam_length >= 2) {
        camatualizado <- camatualizado[c(cam_length - 1, cam_length)]
      }
      
      # READ FILES #
      datafinal <- list()
      for(d in 1:length(camatualizado)){
        datafinal[[d]] <- readRDS(camatualizado[[d]])
      }
      
      #### CALCULATING MEAN TRAIT PER TIME STEP #####
      
      ### change name of list ###
      for(i in 1:length(datafinal)){
        names(datafinal) <- paste0("timestep", length(datafinal):1)
      }
      
      #### function to return column ###
      colunatemp <- function(dados){
        return(dados[, 1])
      }
      
      ### only values of species temperature ###
      for(i in 1:length(datafinal)){
        datafinal[[i]] <- lapply(datafinal[[i]], colunatemp)
      }
      
      for(i in 1:length(datafinal)){
        datafinal[[i]] <- lapply(datafinal[[i]], mean)
      }
      
      for(i in 1:length(datafinal)){
        result <- lapply(datafinal[[i]], is.nan)
        datafinal[[i]] <- datafinal[[i]][which(result == FALSE)]
      }
      
      ####### FINAL MEAN PER TIME STEP #######
      if(ti <= (val$vars$steps[1] - 1)){
        datafinal_less_last <- datafinal[-length(datafinal)]
        datafinal_less_first <- datafinal[-1]
        list_difference <- vector("list", sum(lengths(datafinal_less_last)))
        time <- 0
        for(i in 1:length(datafinal_less_last)){
          for(k in seq_along(datafinal_less_last[[i]])){
            time <- time + 1
            posicao <- which(names(datafinal_less_last[[i]][k]) == names(datafinal_less_first[[i]]))
            list_difference[[time]] <- abs(as.numeric(datafinal_less_last[[i]][k]) - as.numeric(datafinal_less_first[[i]][posicao]))
          }
        }
        time <- 0
    
        list_difference2 <- list()
        if(length(list_difference) >= 1){
          for(i in 1:length(list_difference)){
            if(length(list_difference[[i]]) == 1) {
              list_difference2[[i]] <- list_difference[[i]]
            } else {
              list_difference2[[i]] <- NULL
            }
          } 
        } else {
          list_difference <- 0
          list_difference2 <- 0
          break
        }
        
        datafinal_result <- unlist(list_difference2)
        
        traitevolution <- abs((mean(datafinal_result)) / sum(length(datafinal_less_last) + 1))
        
      } else {
        traitevolution <- 0
      }
      
      ##### SPECIATION AND EXTINCTION ####
      pos <- pos + 1
      pos2 <- pos2 + 1
      
      ratespeciation <- round(sum(sgen3sis$summary$phylo_summary[, 3]) / pos2, digits = 2)
      
      rateextinction <- round(sum(sgen3sis$summary$phylo_summary[, 4]) / pos2, digits = 2)
      
      diversification <- ratespeciation - rateextinction
      
      
      if(ti %% 2 == 1) {
        finalresult$plasticidade[pos] <- plasti[[p]]
        finalresult$replications[pos] <- r
        finalresult$speciation[pos] <- ratespeciation
        finalresult$extinction[pos] <- rateextinction
        finalresult$diversif[pos] <- diversification
        finalresult$traitevolution[pos] <- traitevolution
        finalresult$timestep[pos] <- ti
        finalresult$timesimulation[pos] <- pos2
        finalresult$speciestotal[pos] <- sgen3sis$summary$phylo_summary[, 1][length(sgen3sis$summary$phylo_summary[, 1])]
      }
      #if(pos2 == 50){
     #   saveRDS(val, file = paste0(plasti, '_', pos2, 'simulation', '.rds'))
     #   break
      # }
    } 
  }

  pos2 <- 0
  caminho <- here("data", "raw", "WorldCenter", "output", "config_worldcenter", "traits")
  listfiles <- list.files("data/raw/WorldCenter/output/config_worldcenter/traits")
  filestoread <- length(list.files("data/raw/WorldCenter/output/config_worldcenter/traits"))
  cam <- 0
  for(l in 1:filestoread){
    cam[l] <- caminho
  }
  camatualizado <- 0
  for(k in 1:length(cam)){
    camatualizado[[k]] <- paste(cam[k], listfiles[k], sep = "/", collapse = "--")
  }
  file.remove(camatualizado)  
# rm(val, sgen3sis, rateextinction, ratespeciation, diversification, traitevolution, result, datafinal, datafinal_less_last, datafinal_less_first, list_difference, list_difference2)
}

path <- here("output")

write.csv2(finalresult, file.path(path, "finalresult.csv"), row.names = FALSE)
