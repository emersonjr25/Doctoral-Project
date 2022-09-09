############ DOCTORAL CHAPTER 2 #################
##### EMERSON CAMPOS BARBOSA J?NIOR ####

### PACKAGES ####
library(dplyr)
library(stringr)
library(effectsize)
library(purrr)

dados <- read.csv("Database.csv", header = TRUE, sep = ",")

subdata <- select(dados, paper_no, first_author_surname, pub_year, genus, species, population, source, trait_cat, simp_trait, T, mean)

unicos <- unique(subdata$paper_no)

#pos <- which(unicos[1] == subdata$paper_no)
#subdata$paper_no[pos]

separated_list <- list()

for(i in 1:length(unicos)){
    pos <- which(unicos[i] == subdata$paper_no)
    separated_list[[i]] <- subdata[pos, ]
    pos <- 0
}

name_traits <- list()

for(i in 1:length(separated_list)){
  name_traits[[i]] <- unique(separated_list[[i]]$simp_trait)
}

only_traits <- list()

only_trait <- lapply(separated_list, FUN = function(dados){
  logic <- names(dados) == "simp_trait"
  posi <- which(logic == TRUE)
  return(return(dados[posi]))
})

out <- list()

for(i in 1:length(name_traits)){
  out[[i]] <- vector("list", length(name_traits[[i]]))
}

valor <- 0

for(i in 1:length(name_traits)){
  for(j in seq_along(name_traits[[i]])){
    if(valor < length(name_traits[[i]])){
      valor <- valor + 1
      out[[i]][[valor]] <- which(name_traits[[i]][[j]] == only_trait[[1]])
    } else {
      valor <- 0
      valor <- valor + 1
      out[[i]][[valor]] <- which(name_traits[[i]][[j]] == only_trait[[1]])
    }
  }
}

######################################################

pos_2 <- which(name_traits[[1]][[1]] == separated_list[[1]]$simp_trait)

unique(separated_list[[1]]$population[pos_2])
unique(separated_list[[1]]$source[pos_2])
unique(separated_list[[1]]$species[pos_2])

separated_list[[1]]$T[pos_2]

minimo <- min(separated_list[[1]]$T[pos_2])
maximo <- max(separated_list[[1]]$T[pos_2])

posmin <- which(separated_list[[1]]$T == minimo)
posmax <- which(separated_list[[1]]$T == maximo)

mean_min <- mean(as.numeric(separated_list[[1]]$mean[posmin]))
mean_max <- mean(as.numeric(separated_list[[1]]$mean[posmax]))

group_1 <- as.numeric(separated_list[[1]]$mean[posmin])

group_2 <- as.numeric(separated_list[[1]]$mean[posmax])

hedge_effect <- (mean_min - mean_max) / sd_pooled(group_1, group_2)

