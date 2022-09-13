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

pos_trait <- list()

for(i in 1:length(name_traits)){
  pos_trait[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in 1:length(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    pos_trait[[i]][[j]] <- which(name_traits[[i]][[j]] == only_trait[[i]])
  }
}

### SPECIES POSITION ###

species_per_study <- list()

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    species_per_study[[i]] <- unique(separated_list[[i]]$species_complete[pos_trait[[i]][[j]]])
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
#  only_species_with_more_1 <- 0
# 
#  for(i in 1:length(number_species)){
#    if(number_species[[i]] == 1){
#      only_species_with_more_1[[i]] <- 0
#    } else {
#      only_species_with_more_1[i] <- number_species[[i]]
#    }
#  }
#  
#  only_species_with_more_1 <- only_species_with_more_1 != 0
# 
#  for(i in 1:length(only_species_with_more_1)){
#    if(only_species_with_more_1[i] == TRUE){
#      only_species_with_more_1[i] <- 1
#   } else {
#      only_species_with_more_1[i] <- 0
#    }
#  }
# 
# pos_species <- list()
# 
# for(i in 1:length(species_per_study)){
#   pos_species[[i]] <- vector("list", length(species_per_study[[i]]))
# }
# 
# for(i in 1:length(species_per_study)){
#   for(j in seq_along(species_per_study[[i]])){
#     pos_species[[i]][[j]] <- which(species_per_study[[i]][[j]] == separated_list[[i]][["species_complete"]])
#   }
# }
# 
# pos_species_per_trait <- list()
# 
# for(i in 1:length(species_per_study)){
#   pos_species_per_trait[[i]] <- vector("list", length(species_per_study[[i]]))
# }


### FILTERING TO FIND HEDGE's per species in studies ###

### MINIMUM ###

minimo <- list()

for(i in 1:length(name_traits)){
  minimo[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    minimo[[i]][[j]] <- min(separated_list[[i]]$T[pos_trait[[i]][[j]]])
  }
}

### MAXIMUM ###

maximo <- list()

for(i in 1:length(name_traits)){
  maximo[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    maximo[[i]][[j]] <- max(separated_list[[i]]$T[pos_trait[[i]][[j]]])
  }
}

### MIN POSITION ###

posmin <- list()

for(i in 1:length(name_traits)){
  posmin[[i]] <- vector("list", length(name_traits[[i]]))
}

separada <- 0

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    separada <- separated_list[[i]]$T
    separada[pos_trait[[i]][[j]]] <- 1
    separada[separada != 1] <- 0
    separada[separada == 1] <- separated_list[[i]]$T[pos_trait[[i]][[j]]]
    posmin[[i]][[j]] <- which(separada == minimo[[i]][[j]])
  }
}

### MAX POSITION ###

posmax <- list()

for(i in 1:length(name_traits)){
  posmax[[i]] <- vector("list", length(name_traits[[i]]))
}

separada <- 0

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    separada <- separated_list[[i]]$T
    separada[pos_trait[[i]][[j]]] <- 1
    separada[separada != 1] <- 0
    separada[separada == 1] <- separated_list[[i]]$T[pos_trait[[i]][[j]]]
    posmax[[i]][[j]] <- which(separada == maximo[[i]][[j]])
  }
}

### MIN MEAN ###

mean_min <- list()

for(i in 1:length(name_traits)){
  mean_min[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    mean_min[[i]][[j]] <- mean(as.numeric(separated_list[[i]]$mean[posmin[[i]][[j]]]))
  }
}

### MAX MEAN ###

mean_max <- list()

for(i in 1:length(name_traits)){
  mean_max[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    mean_max[[i]][[j]] <- mean(as.numeric(separated_list[[i]]$mean[posmax[[i]][[j]]]))
  }
}

### GROUP 1 ###

group_1 <- list()

for(i in 1:length(name_traits)){
  group_1[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    group_1[[i]][[j]] <- as.numeric(separated_list[[i]]$mean[posmin[[i]][[j]]])
  }
}

### GROUP 2 ###

group_2 <- list()

for(i in 1:length(name_traits)){
  group_2[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    group_2[[i]][[j]] <- as.numeric(separated_list[[i]]$mean[posmax[[i]][[j]]])
  }
}

### HEDGES's g EFFECT###

hedge_effect <- list()

for(i in 1:length(name_traits)){
  hedge_effect[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    hedge_effect[[i]][[j]] <- abs((mean_min[[i]][[j]] - mean_max[[i]][[j]]) / sd_pooled(group_1[[i]][[j]], group_2[[i]][[j]]))
  }
}

hedge_effect_mean <- list()

for(i in 1:length(name_traits)){
  hedge_effect_mean[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(hedge_effect)){
    hedge_effect_mean[[i]] <- lapply(hedge_effect[[i]], is.na)
}

hedge_effect_out_na <- list()

for(i in 1:length(name_traits)){
  hedge_effect_out_na[[i]] <- vector("list", length(name_traits[[i]]))
}

for(i in seq_along(hedge_effect)){
  pos <- which(hedge_effect_mean[[i]] == FALSE)
  hedge_effect_out_na[[i]] <- as.numeric(hedge_effect[[i]][pos])
  pos <- 0
}

hedge_effect_mean_mean <- list()

for(i in seq_along(hedge_effect)){
  hedge_effect_mean_mean <- lapply(hedge_effect_out_na, mean)
}

#### ORGANIZING TO USE HEDGE's G ###

data_hedges <- do.call(rbind.data.frame, hedge_effect_mean_mean)


#pos <- which(number_species == 1)
#hedge_effect_final <- hedge_effect[pos]
