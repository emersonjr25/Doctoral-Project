#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity #####
##### Methods: Computational simulation #############
##### FUNCTION: LOOP ECOLOGY WITH PLASTICITY #####

loop_ecology2 <- function (config, data, vars, plast_value) {
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
    NEW_abd <- config$gen3sis$ecology$apply_ecology(abundance, traits, local_environment, config, plast_value)
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