#####################################################
########### PLASTICITY AND EVOLUTION STUDY ##########
#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity ####
##### Methods: Computational simulation #############
##### Script to input, modify and run model #########

#### PACKAGES
library(gen3sis)
library(here)
library(dplyr)
#library(raster)

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

#### CARRYING NEW FUNCTIONS AND CONFIG CHANGES #####
folder <- here('scripts', 'functions')
functions_path <- paste0(folder, '/', list.files(folder))
for(i in functions_path){
  source(i)
}

config$gen3sis$general$start_time <- 1000

config$gen3sis$general$end_time <- 1

timesteps_total <- length(config$gen3sis$general$start_time:config$gen3sis$general$end_time)

config$gen3sis$general$max_number_of_species <- 50000
config$gen3sis$general$max_number_of_coexisting_species <- 200000

rep <- 1
quantitity_rep <- 1

config$gen3sis$general$end_of_timestep_observer <- function(data, vars, config){
  save_traits()
}

config$gen3sis$speciation$divergence_threshold <- 10

environment_type <- c('random', 'stable_low', 'stable_fast')
environment_type_chose <- environment_type[2]

plasti <- c(0, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1)

pos <- 0
pos2 <- 0

finalresult <- data.frame(plasticity = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          replications = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          enviroment_type = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          speciation = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          extinction = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          diversif = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          traitevolution = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          timesimulation = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          abundance = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          occupancy = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0),
                          alive_spec = runif(quantitity_rep * length(plasti) * timesteps_total, 0, 0))


#### EXECUTION SIMULATION ####
for(p in 1:length(plasti)){
  
  for(r in rep:rep){
    
    #### REMOVING TRAITS OF ANTERIOR SIMULATIONS ####
    output_files <- here("data", "raw", "WorldCenter", "output", "config_worldcenter")
    general_files <- list.files(output_files)
    traits_path <- paste0('traits', r)
    if(sum(traits_path == general_files) == 1){
      general_path <- here("data", "raw", "WorldCenter", "output", "config_worldcenter", traits_path)
      listfiles <- list.files(paste0("data/raw/WorldCenter/output/config_worldcenter/", traits_path))
      filestoread <- length(listfiles)
      path_temporary <- 0
      for(l in 1:filestoread){
        path_temporary[l] <- general_path
      }
      path_update <- 0
      for(k in 1:length(path_temporary)){
        path_update[[k]] <- paste(path_temporary[k], listfiles[k], sep = "/", collapse = "--")
      }
      file.remove(path_update, showWarnings = FALSE)
    } 
    
    
    ####### CREATING INITIAL WORLD ######
    
    val <- list(data = list(),
                vars = list(),
                config = config)
    val$config$gen3sis$general$random_seed <- r
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
    val <- modify_input_temperature(val$config, val$data, val$vars, environment_type_chose)
   
     for (ti in val$vars$steps) {
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
      val <- loop_evolution(val$config, val$data, val$vars)
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

      sgen3sis <- make_summary(val$config, val$data, val$vars,
                                total_runtime, save_file = FALSE)
      
      # raster_data <- cbind(val[["data"]][["landscape"]][["coordinates"]],
      #       val[["data"]][["landscape"]][["environment"]][, 1])
      # colnames(raster_data) <- c('x', 'y', 'temp')
      # ras <- rasterFromXYZ(raster_data)
      # max_ras <- 1
      # min_ras <- 0
      # #sequencia <- seq(0.1, 1, 0.1)
      # rc <- c('#4f75e8', '#0A2F51', '#ffa600', '#fe9700', '#fc8700', '#f97600',
      #   '#f66504', '#f2520e', '#ed3c16', '#e81f1c')
      # image(ras, col=rc, bty = "o", xlab = "", ylab = "", las=1, asp = 1)
      # mtext(4, text="Temperature", line=1, cex=1.2)
      # raster::plot(rasterFromXYZ(raster_data), legend.only=TRUE, add=TRUE,col=rc)
      
      #  if(ti <= 999){
      #   ras <- rasterFromXYZ(sgen3sis$summary$`richness-final`)
      #   max_ras <- max(ras@data@values, na.rm=TRUE)
      #   min_ras <- min(ras@data@values, na.rm=TRUE)
      #   # rc <- color_richness(max(ras@data@values, na.rm=TRUE) + 1)
      #   #terrain color
      #   zerorichness_col <- "navajowhite3"
      #   if (max_ras==0){ #if all extinct
      #     rc <-  zerorichness_col
      #   } else {
      #     rc <- color_richness(max_ras)
      #     if (min_ras==0){ #if there is zero-richness (i.e. inhabited sites)
      #       rc <- c(zerorichness_col, rc)
      #     }
      #   }
      #   image(ras, col=rc, bty = "o", xlab = "", ylab = "", las=1, asp = 1)
      #   mtext(4, text="Final \u03B1 richness", line=1, cex=1.2)
      #   raster::plot(rasterFromXYZ(sgen3sis$summary$`richness-final`), legend.only=TRUE, add=TRUE,col=rc)
      # }
      ################## TRAIT EVOLUTION #####################
      
      ##### PATHWAY TO TRAITS DATA ####
      traits_path <- paste0('traits', r)
      general_path <- here("data", "raw", "WorldCenter", "output", "config_worldcenter", traits_path)
      listfiles <- list.files(paste0("data/raw/WorldCenter/output/config_worldcenter/", traits_path))
      filestoread <- length(listfiles)
      
      #### Organizing selection of trait files ####
      path_temporary <- 0
      for(l in 1:filestoread){
        path_temporary[l] <- general_path
      }
      
      path_update <- 0
      for(k in 1:length(path_temporary)){
        path_update[[k]] <- paste(path_temporary[k], listfiles[k], sep = "/", collapse = "--")
      }
      
      # SORT #
      path_update <- path_update[order(as.numeric(gsub("[^0-9]+", "", path_update)), decreasing = TRUE)]
      path_length <- length(path_update)
      if(path_length >= 2) {
        path_update <- path_update[c(path_length - 1, path_length)]
      }
      
      # READ FILES #
      datafinal <- list()
      for(d in 1:length(path_update)){
        datafinal[[d]] <- readRDS(path_update[[d]])
      }
      
      #### CALCULATING MEAN TRAIT PER TIME STEP #####
      
      ### change name of list ###
      for(i in 1:length(datafinal)){
        names(datafinal) <- paste0("timestep", length(datafinal):1)
      }
      
      #### function to return column ###
      col_temp <- function(x){
        return(x[, 1])
      }
      
      ### only values of species temperature ###
      for(i in 1:length(datafinal)){
        datafinal[[i]] <- lapply(datafinal[[i]], col_temp)
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
        trait_ancient_less_last <- datafinal[-length(datafinal)]
        trait_new_less_first <- datafinal[-1]
        list_difference <- vector("list", sum(lengths(trait_ancient_less_last)))
        position_list <- vector("list", sum(lengths(trait_ancient_less_last)))
        time <- 0
        for(i in 1:length(trait_ancient_less_last)){
          for(k in seq_along(trait_ancient_less_last[[i]])){
            time <- time + 1
            position <- which(names(trait_ancient_less_last[[i]][k]) == names(trait_new_less_first[[i]]))
            position_list[[time]] <- as.numeric(trait_ancient_less_last[[i]][position])
            list_difference[[time]] <- abs(as.numeric(trait_ancient_less_last[[i]][k]) - as.numeric(trait_new_less_first[[i]][position]))
          }
        }
        time <- 0
        
        if(length(list_difference) < 1){
          list_difference <- 0
          break
        }
        
        datafinal_result <- unlist(list_difference)
        trait_previous <- unlist(position_list)
        
        traitevolution <- abs(mean(datafinal_result / trait_previous))
        
      } else {
        traitevolution <- 0
      }
      ### extra results ###
      count_loop <- 1
      abundance_species <- vector(mode = 'integer', length(val[['data']][['all_species']]))
      
      for(i in seq_along(val[['data']][['all_species']])){
        abundance_species[count_loop] <- length(val[['data']][['all_species']][[count_loop]]$abundance)
        count_loop <- count_loop + 1
      }
 
      total_abundance <- (sum(abundance_species) / as.numeric(sgen3sis[["summary"]][["phylo_summary"]][, 2][length(sgen3sis[["summary"]][["phylo_summary"]][, 2])]))
      
      occupancy <- as.numeric(sgen3sis[["summary"]][["occupancy"]][length(sgen3sis[["summary"]][["occupancy"]])])
      
      alive_species <- as.numeric(sgen3sis[["summary"]][["phylo_summary"]][, 2][length(sgen3sis[["summary"]][["phylo_summary"]][, 2])])
      
      ##### SPECIATION AND EXTINCTION ####
      pos <- pos + 1
      pos2 <- pos2 + 1
      
      current_time <- length(sgen3sis$summary$phylo_summary[, 3])
      current_time_less_one <- current_time - 1
      
      ratespeciation <- abs(as.numeric(sgen3sis$summary$phylo_summary[, 3][current_time] / sgen3sis$summary$phylo_summary[, 2][current_time_less_one]))
      
      rateextinction <- abs(as.numeric(sgen3sis$summary$phylo_summary[, 4][current_time] / sgen3sis$summary$phylo_summary[, 2][current_time_less_one]))
      
      diversification <- ratespeciation - rateextinction
      
      finalresult$plasticity[pos] <- plasti[[p]]
      finalresult$replications[pos] <- r
      finalresult$enviroment_type[pos] <- environment_type_chose
      finalresult$speciation[pos] <- ratespeciation
      finalresult$extinction[pos] <- rateextinction
      finalresult$diversif[pos] <- diversification
      finalresult$traitevolution[pos] <- traitevolution
      finalresult$timesimulation[pos] <- pos2
      finalresult$abundance[pos] <- total_abundance
      finalresult$occupancy[pos] <- occupancy
      finalresult$alive_spec[pos] <- alive_species
    } 
  }
  pos2 <- 0
  traits_path <- paste0('traits', rep)
  general_path <- here("data", "raw", "WorldCenter", "output", "config_worldcenter", traits_path)
  listfiles <- list.files(paste0("data/raw/WorldCenter/output/config_worldcenter/", traits_path))
  filestoread <- length(listfiles)
  path_temporary <- 0
  for(l in 1:filestoread){
    path_temporary[l] <- general_path
  }
  path_update <- 0
  for(k in 1:length(path_temporary)){
    path_update[[k]] <- paste(path_temporary[k], listfiles[k], sep = "/", collapse = "--")
  }
  file.remove(path_update, showWarnings = FALSE)
  rm(val, sgen3sis, rateextinction, ratespeciation, diversification, traitevolution, result, datafinal, trait_ancient_less_last, trait_new_less_first, list_difference, position_list)
}

path <- here("output")
write.csv2(finalresult, file.path(path, paste0('rep_', rep, "_finalresult.csv")), row.names = FALSE)