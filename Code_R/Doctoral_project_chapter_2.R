############ DOCTORAL CHAPTER 2 #################
##### EMERSON CAMPOS BARBOSA JÚNIOR ####

### PACKAGES ####
library(dplyr)
library(stringr)

dados <- read.csv("Database.csv", header = TRUE, sep = ",")

subdata <- select(dados, paper_no, first_author_surname, pub_year, genus, species, population, source, trait_cat, simp_trait, T, mean)

unicos <- unique(subdata$paper_no)

pos <- list()

# for(i in 1:length(unicos)){
#   for(j in 1:length(subdata$paper_no)){
#     pos <- which(unicos[i] == subdata$paper_no)
#     separated_list[i] <- pos
#   }
# }

subdata$paper_no[which(unicos[1] == subdata$paper_no)]

pos <- which(unicos[1] == subdata$paper_no)
subdata$paper_no[pos]
subdata[pos, ]


separated_list <- list()

for(i in 1:length(unicos)){
    pos <- which(unicos[i] == subdata$paper_no)
    separated_list[[i]] <- subdata[pos, ]
    pos <- 0
}

nova <- list()

for(i in 1:length(separated_list)){
  nova[[i]] <- unique(separated_list[[i]]$simp_trait)
}

pos <- which(nova[[1]][[1]] == separated_list[[1]]$simp_trait)

separated_list[[1]]$T[pos]

### TESTE ###
unique(separated_list[[1]]$population[pos])
unique(separated_list[[1]]$source[pos])

minimo <- min(separated_list[[1]]$T[pos])
maximo <- max(separated_list[[1]]$T[pos])

posmin <- which(separated_list[[1]]$T == minimo)
posmax <- which(separated_list[[1]]$T == maximo)


as.numeric(separated_list[[1]]$mean)
mean(as.numeric(separated_list[[1]]$mean[posmin]))
mean(as.numeric(separated_list[[1]]$mean[posmax]))



for(i in 1:length(unicos)){
  for(j in 1:length(subdata$paper_no)){
    separated_list[i] <- subdata$paper_no[j][subdata$paper_no == unicos[i]] 
  }
}
