# ARQUIVO LIXO #
create_ancestor_species <- function(landscape, config) {
  
  all_species <- list()    
  suitablecells<- c(449, 450, 480, 481)
  for(cellID in 1: length(suitablecells)){ 
    cell <- as.character(suitablecells[cellID])
    new_species <- create_species(as.character(cell), config)  # use this function to create species, one can provide directly the initial cells
    new_species$traits[ , "dispersal"] <- 1
    new_species$traits[ , "temp"] <-  rnorm(1,20,0.5) 
    new_species$traits[ , "prec"] <- rnorm(1,500,50) 
    all_species <- append(all_species, list(new_species))
  }
  
  return(all_species)}


create_ancestor_species <- function(landscape, config) {
  range <- c(-95, -24, -68, 13)
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] &
    co[, "x"] <= range[2] &
    co[, "y"] >= range[3] &
    co[, "y"] <= range[4]
  
  new_species <- list()
  for(i in 1:10){
    initial_cells <- rownames(co)[selection]
    initial_cells <- sample(initial_cells, 1)
    new_species[[i]] <- create_species(initial_cells, config)
    #set local adaptation to max optimal temp equals local temp
    new_species[[i]]$traits[ , "temp"] <- landscape$environment[initial_cells,"temp"]
    new_species[[i]]$traits[ , "dispersal"] <- 1 
    #plot_species_presence(landscape, species=new_species[[i]])
  }
  
  return(new_species)
}

