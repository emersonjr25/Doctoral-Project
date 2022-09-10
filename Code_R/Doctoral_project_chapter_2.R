############ DOCTORAL CHAPTER 2 #################
##### EMERSON CAMPOS BARBOSA J?NIOR ####

### PACKAGES ####
library(dplyr)
library(stringr)
library(effectsize)
library(purrr)

### READING DATA ###

dados <- read.csv("Database.csv", header = TRUE, sep = ",")

subdata <- select(dados, paper_no, first_author_surname, pub_year, genus, species, population, source, trait_cat, simp_trait, T, mean)

subdata2 <- paste0(subdata$genus," ", subdata$species)

subdata$species_complete <- subdata2

### LIST SEPARATING STUDIES ####

unicos <- unique(subdata$paper_no)

separated_list <- list()

for(i in 1:length(unicos)){
  pos <- which(unicos[i] == subdata$paper_no)
  separated_list[[i]] <- subdata[pos, ]
  pos <- 0
}

### DISCOVERY OF POSITION OF THE TRAIT ###

name_traits <- list()

for(i in 1:length(separated_list)){
  name_traits[[i]] <- unique(separated_list[[i]]$simp_trait)
}

# SUB LIST #

only_traits <- list()

only_trait <- lapply(separated_list, FUN = function(dados){
  logic <- names(dados) == "simp_trait"
  posi <- which(logic == TRUE)
  return(return(dados[posi]))
})

### TRAITS POSITION  ###

out <- list()

for(i in 1:length(name_traits)){
  out[[i]] <- vector("list", length(name_traits[[i]]))
}


for(i in 1:length(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    out[[i]][[j]] <- which(name_traits[[i]][[j]] == only_trait[[i]])
  }
}

### SPECIES POSITION ###

pos_species <- list()

for(i in 1:length(species_per_study)){
  pos_species[[i]] <- vector("list", length(species_per_study[[i]]))
}

for(i in 1:length(species_per_study)){
  for(j in seq_along(species_per_study[[i]])){
    pos_species[[i]][[j]] <- which(species_per_study[[i]][[j]] == separated_list[[i]][["species_complete"]])
  }
}

### FILTERING TO FIND HEDGE's per species in studies ###

species_per_study <- list()

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    species_per_study[[i]] <- unique(separated_list[[i]]$species_complete[out[[i]][[j]]])
  }
}

number_species <- list()
for(i in 1:length(species_per_study)){
  if(species_per_study[[i]] == 1){
    number_species[[i]] <- 1
  } else {
    number_species[[i]] <- length(species_per_study[[i]])
  }
}
only_species_with_more_1 <- 0

for(i in 1:length(number_species)){
  if(number_species[[i]] == 1){
    only_species_with_more_1[[i]] <- 0
  } else {
    only_species_with_more_1[i] <- number_species[[i]]
  }
}

only_species_with_more_1 <- only_species_with_more_1 != 0

for(i in 1:length(only_species_with_more_1)){
  if(only_species_with_more_1[i] == TRUE){
    only_species_with_more_1[i] <- 1
  } else {
    only_species_with_more_1[i] <- 0
  }
}

#for(c in seq_along(number_species)){
  
  #if(number_species[[c]] == 1){
    
    ### MINIMUM ###
    
    minimo <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        minimo[[i]] <- min(separated_list[[i]]$T[out[[i]][[j]]])
      }
    }
    
    ### MAXIMUM ###
    
    maximo <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        maximo[[i]] <- max(separated_list[[i]]$T[out[[i]][[j]]])
      }
    }
    
    ### MIN POSITION ###
    
    posmin <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        posmin[[i]] <- which(separated_list[[i]]$T == minimo[[i]])
      }
    }
    
    ### MAX POSITION ###
    
    posmax <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        posmax[[i]] <- which(separated_list[[i]]$T == maximo[[i]])
      }
    }
    
    ### MIN MEAN ###
    
    mean_min <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        mean_min[[i]] <- mean(as.numeric(separated_list[[i]]$mean[posmin[[i]]]))
      }
    }
    
    ### MAX MEAN ###
    
    mean_max <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        mean_max[[i]] <- mean(as.numeric(separated_list[[i]]$mean[posmax[[i]]]))
      }
    }
    
    ### GROUP 1 ###
    
    group_1 <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        group_1[[i]] <- as.numeric(separated_list[[i]]$mean[posmin[[i]]])
      }
    }
    
    ### GROUP 2 ###
    
    group_2 <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        group_2 [[i]] <- as.numeric(separated_list[[i]]$mean[posmax[[i]]])
      }
    }
    
    ### HEDGES's g EFFECT###
    
    hedge_effect <- list()
    
    for(i in seq_along(name_traits)){
      for(j in seq_along(name_traits[[i]])){
        hedge_effect[[i]] <- (mean_min[[i]] - mean_max[[i]]) / sd_pooled(group_1[[i]], group_2[[i]])
      }
    }
 # } else {
    
  #}
#}

