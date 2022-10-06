#####################################################
########### PLASTICITY AND EVOLUTION STUDY ##########
#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Computational simulation #############
##### Script to input, modify and run model #########

#### PACKAGES ####
library(igraph)
library(stringi)
library(gen3sis)
library(raster)
library(truncnorm)
library(here)


#### SIMULATION ####

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

#pathway <- paste(datapath, "output/config_worldcenter/phy_res", sep = "/", collapse = "--")
#dir.create(pathway, showWarnings = FALSE)
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



### MODIFICATIONS IN CONFIG ###

config$gen3sis$general$start_time <- 20

config$gen3sis$general$end_time <- 1

timesteps_total <- length(config$gen3sis$general$start_time:config$gen3sis$general$end_time)

config$gen3sis$general$max_number_of_species <- 5000

config$gen3sis$general$end_of_timestep_observer <- function(data, vars, config){
  save_traits()
  #save_species()
  #plot_richness(data$all_species, data$landscape)
  
}

#########################################################
rep <- 1

plasti <- seq(0.1, 0.1, 0.1)

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
      seq((x * plast) - plast, (x * plast) + plast, 0.01)
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
    
    for (ti in val$vars$steps) {
      
      ##########################################################  
      
      # state <- "hightemp"
      # if(val$vars$steps[ti] == 10){
      #   if(state == "hightemp"){
      #     val$data$landscape$environment[, 1] <-  val$data$landscape$environment[, 1] + 0.9
      #     val$data$inputs$environments$temp  <- val$data$inputs$environments$temp + 0.9 }
      # }
      
      #val$data$inputs$environments$temp <- sapply(val$data$inputs$environments$temp, FUN = function(dados){
      #return(dados + abs(rnorm(1, 2, 0.1)))
      #})
      
      #val$data$landscape$environment[, 1] <- sapply(val$data$landscape$environment[, 1], FUN = function(dados){
      # return(dados + abs(rnorm(1, 0.3, 0.1)))
      #})
      
      #######################################################
      
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
      val <- loop_dispersal(val$config, val$data, val$vars)
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
      val <- loop_ecology(val$config, val$data, val$vars)
      #if (val$vars$flag == "max_number_coexisting_species") {
      # print("max number of coexisting species reached, breaking loop")
      # break
      
      if (verbose >= 0 & val$vars$flag == "OK") {
        cat("Simulation finished. All OK \n")
      } else if (verbose >= 0 & val$vars$flag == "max_number_species") {
        cat("Simulation finished. Early abort due to exceeding max number of species")
      } else if (verbose >= 0 &
                 val$vars$flag == "max_number_coexisting_species") {
        cat("Simulation finished. Early abort due to exceeding max number of co-occuring species")
      }
      val <- update.phylo(val$config, val$data, val$vars)
      # write.table(
      #  val$data$phy,
      # file = file.path(pathway,
      #                paste0("phy", "plast", p, "rep", r, ".txt")),
      # sep = "\t"
      # )
      #write_nex(
      #  phy = val$data$phy,
      #  label = "species",
      #  file.path(output_location = pathway,
      #     paste0("phy", "plast", p, "rep", r, ".nex"))
      # )
      system_time_stop <- Sys.time()
      total_runtime <- difftime(system_time_stop, system_time_start,
                                units = "hours")[[1]]
      #write_runtime_statisitics(val$data, val$vars, val$config,
      #                          total_runtime)
      sgen3sis <- make_summary(val$config, val$data, val$vars,
                               total_runtime, save_file = FALSE)
      #plot_summary(sgen3sis)
      
      # if (verbose >= 1) {
      #  cat("Simulation runtime:", total_runtime, "hours\n")
      # }
      ################## TRAIT EVOLUTION #####################
      
      ##### PATHWAY TO TRAITS DATA####
      
      caminho <- here("data", "raw", "WorldCenter", "output", "config_worldcenter", "traits")
      
      listfiles <- list.files("data/raw/WorldCenter/output/config_worldcenter/traits")
      
      filestoread <- length(list.files("data/raw/WorldCenter/output/config_worldcenter/traits"))
      
      
      #### Organizing selection of trait files ####
      cam <- 0
      for(l in 1:filestoread){
        cam[l] <- caminho
      }
      cam
      
      camatualizado <- 0
      for(k in 1:length(cam)){
        camatualizado[[k]] <- paste(cam[k], listfiles[k], sep = "/", collapse = "--")
      }
      
      # SORT #
      camatualizado <- camatualizado[order(as.numeric(gsub("[^0-9]+", "", camatualizado)), decreasing = TRUE)]
      
      
      # READ FILES #
      datafinal <- list()
      for(d in 1:filestoread){
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
       datafinal2 <- datafinal
       datafinal3 <- datafinal
       datafinal2 <- datafinal[-length(datafinal)]
       datafinal3 <- datafinal3[-1]
       list_difference <- list()
       list_difference <- vector("list", sum(lengths(datafinal2)))
       time <- 0
       for(i in 1:length(datafinal2)){
         for(k in seq_along(datafinal2[[i]])){
           time <- time + 1
           posicao <- which(names(datafinal2[[i]][k]) == names(datafinal3[[i]]))
           list_difference[[time]] <- abs(as.numeric(datafinal2[[i]][k]) - as.numeric(datafinal3[[i]][posicao]))
         }
       }
       time <- 0
       
       list_difference2 <- list()
       for(i in 1:length(list_difference)){
         if(length(list_difference[[i]]) == 1) {
           list_difference2[[i]] <- list_difference[[i]]
         } else {
           list_difference2[[i]] <- NULL
         }
       }
       final <- unlist(list_difference2)
       
       traitevolution <- mean(final) / sum(length(datafinal) + 1)
       
       #pos <- which(names(datafinal[[12]][21]) == names(datafinal2[[12]]))
       #list_difference[[1]] <- abs(as.numeric(datafinal[[12]][21]) - as.numeric(datafinal[[12]][pos]))
       #list_difference[[2]] <- abs(as.numeric(datafinal[[2]][1]) - as.numeric(datafinal[[3]][pos]))
      
     } else {
       traitevolution <- 0
     }
       
       ##### SPECIATION AND EXTINCTION ####
      pos <- pos + 1
      pos2 <- pos2 + 1
      
      ratespeciation <- round(sum(sgen3sis$summary$phylo_summary[, 3]) / pos2, digits = 2)
      
      rateextinction <- round(sum(sgen3sis$summary$phylo_summary[, 4]) / pos2, digits = 2)
      
      diversification <- ratespeciation - rateextinction
      
      
      # sequence <- sort(seq(from = config$gen3sis$general$end_time, to = config$gen3sis$general$start_time, by = 5), decreasing = TRUE)
      # 
      # #### OPTION 1 ####
      #  for(seque in 1:(length(sequence))){
      #   if(val$vars$steps[ti] == sequence[seque]) {
      #     finalresult$plasticidade[pos] <- plasti[p]
      #     finalresult$replications[pos] <- r
      #     finalresult$speciation[pos] <- ratespeciation
      #     finalresult$extinction[pos] <- rateextinction
      #     finalresult$diversif[pos] <- diversification
      #     finalresult$traitevolution[pos] <- traitevolution
      #     finalresult$timesimulation[pos] <- ti
      #   }
      #  }
      #### OPTION 2 ###
      if(ti %% 2 == 1) {
        finalresult$plasticidade[pos] <- plasti[p]
        finalresult$replications[pos] <- r
        finalresult$speciation[pos] <- ratespeciation
        finalresult$extinction[pos] <- rateextinction
        finalresult$diversif[pos] <- diversification
        finalresult$traitevolution[pos] <- traitevolution
        finalresult$timestep[pos] <- ti
        finalresult$timesimulation[pos] <- pos2
      }
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
  rm(val, sgen3sis, rateextinction, ratespeciation, diversification, traitevolution, result, datafinal, datafinal2, datafinal3, list_difference, list_difference2)
}

path <- here("output")
write.csv2(finalresult, file.path(path, "finalresult.csv"), row.names = FALSE)
#write.csv2(finalresult, file = "finalresult.csv", row.names = FALSE)
#saveRDS(finalresult, file.path(path, "finalresult.RDS"))
