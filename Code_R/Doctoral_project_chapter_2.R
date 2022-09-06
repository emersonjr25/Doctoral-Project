############ DOCTORAL CHAPTER 2 #################
##### EMERSON CAMPOS BARBOSA JÚNIOR ####

### PACKAGES ####
library(dplyr)
library(stringr)

dados <- read.csv("Database.csv", header = TRUE, sep = ",")

subdata <- select(dados, paper_no, first_author_surname, pub_year, genus, species, population, source, trait_cat, simp_trait, T, mean)

new <- 0

for(i in 1:length(subdata$paper_no)){
  if(unique(subdata$paper_no[i] == subdata$paper_no))
   new[i] <- unique(subdata$paper_no)
}

data.frame()
unique(subdata$paper_no)
filter()
