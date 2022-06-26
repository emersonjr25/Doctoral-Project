#### DOCTORAL ANALYSIS ###
dados <- read.csv("C:/Users/Emerson Júnior/Google Drive/Doutorado/Base de Dados/Database.csv", header = TRUE, sep = ",")
metadado <- read.csv("C:/Users/Emerson Júnior/Google Drive/Doutorado/Base de Dados/Metadata.csv", header = TRUE, sep = ",")
metadado2 <- read.csv2("C:/Users/Emerson Júnior/Google Drive/Doutorado/Base de Dados/Metadata.csv")
write.csv2(metadado, file = "metadado.csv")
getwd()
View(dados)
View(metadado)
View(metadado2)
unique(dados$country)

install.packages("stringr")
library("stringr")

plot(dados$mean ~ dados$T)
unique(dados$state_province)
unique(dados$species)
count(dados$species)

str_count(dados$species, "nubila")
sum(str_count(dados$species, "nubila"))
str_length(dados$species)
