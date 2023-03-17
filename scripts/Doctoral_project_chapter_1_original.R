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

#### PREPARATION IN CONFIG - SPECIES AND SYSTEM ####
config$gen3sis$general$start_time <- 50

config$gen3sis$general$end_time <- 1

timesteps_total <- length(config$gen3sis$general$start_time:config$gen3sis$general$end_time)

config$gen3sis$general$max_number_of_species <- 50000

config$gen3sis$general$end_of_timestep_observer <- function(data, vars, config){
  save_traits()
}

rep <- 1

plasti <- c(0)

pos <- 0
pos2 <- 0

finalresult <- data.frame(plasticidade = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          replications = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          speciation = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          extinction = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          diversif = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          traitevolution = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          timestep = runif(rep * length(plasti) * timesteps_total, 0, 0),
                          timesimulation = runif(rep * length(plasti) * timesteps_total, 0, 0))
#### CARRYING FUNCTIONS #####
folder <- here('scripts', 'functions')
functions_path <- paste0(folder, '/', list.files(folder))
for(i in functions_path){
  source(i)
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

#### EXECUTION SIMULATION ####

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
       # if (ti == 192) {
       #   break 
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
      
      current_time <- length(sgen3sis$summary$phylo_summary[, 3])
      current_time_less_one <- current_time - 1
      
      ratespeciation <- abs(as.numeric(sgen3sis$summary$phylo_summary[, 3][current_time] / sgen3sis$summary$phylo_summary[, 2][current_time_less_one]))
      
      rateextinction <- abs(as.numeric(sgen3sis$summary$phylo_summary[, 4][current_time] / sgen3sis$summary$phylo_summary[, 2][current_time_less_one]))
  
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
      }
     # if(pos2 == 16){
     #   saveRDS(val, file = paste0(plasti, '_', pos2, 'simulation', '.rds'))
       # break
       #}
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
 #rm(val, sgen3sis, rateextinction, ratespeciation, diversification, traitevolution, result, datafinal, datafinal_less_last, datafinal_less_first, list_difference, list_difference2)
}

path <- here("output")
write.csv2(finalresult, file.path(path, "finalresult.csv"), row.names = FALSE)
