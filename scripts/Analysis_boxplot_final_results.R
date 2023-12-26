#####################################################
########### PLASTICITY AND EVOLUTION STUDY ##########
#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Computational simulation #############
##### Script to analysis results from simulation ####

#### packages ####
library(tidyverse)
library(ggplot2)
library(here)

#### data to graph ####
path_general <- here('output/files_to_results/stable_slow_envi_without_cost')
path_files <- list.files(path_general, pattern = 'csv')
path_complete <- paste0(path_general, "/", path_files)
data_general_one <- vector('list', length(path_complete))

for(i in seq_along(path_complete)){
  data_general_one[[i]] <- read.csv2(path_complete[i])
  assign(paste0('data', i), data_general_one[[i]])
}

all_data <- unlist(eapply(.GlobalEnv, is.data.frame))
data_general_one <- do.call(rbind, mget(names(all_data)[all_data]))
rownames(data_general_one) <- NULL
data_general_one <- data_general_one %>% filter(replications != 0)

rm(list = ls()[(ls() == 'data_general_one') == FALSE])
colnames(data_general_one)[1] <- c('plasticity')

########################################################
path_general2 <- here('output/files_to_results/stable_slow_envi_with_cost')
path_files <- list.files(path_general2, pattern = 'csv')
path_complete <- paste0(path_general2, "/", path_files)
data_general_two <- vector('list', length(path_complete))

for(i in seq_along(path_complete)){
  data_general_two[[i]] <- read.csv2(path_complete[i])
  assign(paste0('data', i), data_general_two[[i]])
}

all_data <- unlist(eapply(.GlobalEnv, is.data.frame))
names_test <- names(all_data)[all_data]
data_general_two <- do.call(rbind, mget(names_test[grep("[1-9]", names_test)]))
rownames(data_general_two) <- NULL
data_general_two <- data_general_two %>% filter(replications != 0)

rm(list = ls()[(ls() == 'data_general_one' | ls() == 'data_general_two') == FALSE])
colnames(data_general_two)[1] <- c('plasticity')

#########################################################
path_general3 <- here('output/files_to_results/stable_fast_envi_without_cost')
path_files <- list.files(path_general3, pattern = 'csv')
path_complete <- paste0(path_general3, "/", path_files)
data_general_three <- vector('list', length(path_complete))

for(i in seq_along(path_complete)){
  data_general_three[[i]] <- read.csv2(path_complete[i])
  assign(paste0('data', i), data_general_three[[i]])
}

all_data <- unlist(eapply(.GlobalEnv, is.data.frame))
names_test <- names(all_data)[all_data]
data_general_three <- do.call(rbind, mget(names_test[grep("[1-9]", names_test)]))
rownames(data_general_three) <- NULL
data_general_three <- data_general_three %>% filter(replications != 0)

rm(list = ls()[(ls() == 'data_general_one' | ls() == 'data_general_two' | ls() == 'data_general_three') == FALSE])
colnames(data_general_three)[1] <- c('plasticity')

#########################################################
path_general4 <- here('output/files_to_results/stable_fast_envi_with_cost')
path_files <- list.files(path_general4, pattern = 'csv')
path_complete <- paste0(path_general4, "/", path_files)
data_general_four<- vector('list', length(path_complete))

for(i in seq_along(path_complete)){
  data_general_four[[i]] <- read.csv2(path_complete[i])
  assign(paste0('data', i), data_general_four[[i]])
}

all_data <- unlist(eapply(.GlobalEnv, is.data.frame))
names_test <- names(all_data)[all_data]
data_general_four <- do.call(rbind, mget(names_test[grep("[1-9]", names_test)]))
rownames(data_general_four) <- NULL
data_general_four <- data_general_four %>% filter(replications != 0)

rm(list = ls()[(ls() == 'data_general_one' | ls() == 'data_general_two' | ls() == 'data_general_three' | ls() == 'data_general_four') == FALSE])
colnames(data_general_four)[1] <- c('plasticity')

data_general_one$cost <- "without"
data_general_two$cost <- "with"
data_general_three$cost <- "without"
data_general_four$cost <- "with"
new_data <- rbind(data_general_one, data_general_two,
                  data_general_three, data_general_four)
new_data$enviroment_type[new_data$enviroment_type == "stable_fast"] <- "fast"
new_data$enviroment_type[new_data$enviroment_type == "stable_low"] <- "low"

plasticity_which_species_die <- new_data %>%
  filter(alive_spec == 0) %>%
  select(plasticity, cost, enviroment_type) %>%
  group_by(plasticity, cost, enviroment_type) %>%
  summarise(count = n()) %>%
  filter(count == 10) %>% 
  select(plasticity, cost, enviroment_type)

plasticity_which_species_die <- as.numeric(unlist(plasticity_which_species_die))

data_with_alives <- new_data[!new_data$plasticity %in% plasticity_which_species_die, ]

new_data[is.na(new_data)] <- 0

result <- new_data %>%
  as_tibble() %>%
  filter(timesimulation > 40) %>%
  select(-c(timesimulation, abundance, alive_spec, occupancy)) %>%
  rename('Trait evolution' = traitevolution,
         Diversification = diversif,
         Speciation = speciation,
         Extinction = extinction) %>%
  group_by(plasticity, enviroment_type, cost, replications) %>%
  summarize_all(mean) %>%
  pivot_longer(col = -c(plasticity, enviroment_type, cost, replications)) %>%
  mutate(plasticity = as.factor(plasticity)) %>%
  ggplot(aes(x = plasticity, y = value, color=enviroment_type, size=cost)) +
  geom_boxplot() +
  scale_size_manual(values = c("with" = 1, "without" = 0)) +
  #geom_jitter() +
  facet_wrap(~name, scales = "free_y") +
  xlab("Plasticity") +
  #ggtitle('Effect of plasticity on adaptive evolution - climatic changes & cost context') +
  theme_bw() +
  scale_color_manual(values=c("#542788", "#b35806"),
                     name = "climatic change") +
  theme(plot.title =
          element_text(size = 14,
                       face = 2,
                       hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.position = "bottom") 

tiff(filename = file.path(here('output'), paste0("plot", "_", result[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 800, 
     height = 600,
     units = "px",
     res = 100)
print(result)
dev.off()

mechanism <- new_data %>%
  as_tibble() %>%
  filter(timesimulation > 40, abundance > 0) %>%
  select(-c(timesimulation, traitevolution, diversif, speciation, extinction, alive_spec)) %>%
  rename(Abundance = abundance,
         Occupancy = occupancy) %>%
  group_by(plasticity, enviroment_type, cost, replications) %>%
  summarize_all(mean) %>%
  pivot_longer(col = -c(plasticity, enviroment_type, cost, replications)) %>%
  mutate(plasticity = as.factor(plasticity)) %>%
  ggplot(aes(x = plasticity, y = value, color=enviroment_type, size=cost)) +
  geom_boxplot() +
  #geom_jitter() +
  scale_size_manual(values = c("with" = 1, "without" = 0)) +
  facet_wrap(~name, scales = "free_y") +
  xlab("Plasticity") +
  #ggtitle('Mechanisms: abundance and occupancy in climatic changes & cost context') +
  theme_bw() +
  scale_color_manual(values=c("#542788", "#b35806"),
                     name = "climatic change") +
  theme(plot.title =
          element_text(size = 14,
                       face = 2,
                       hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.position = "bottom") 

tiff(filename = file.path(here('output'), paste0("mechanism_plot", "_", result[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 800, 
     height = 400, 
     units = "px",
     res = 100)
print(mechanism)
dev.off()

#### testing cost differences ####
new_data %>%
  as_tibble() %>%
  filter(timesimulation > 40) %>%
  select(-c(timesimulation, abundance, alive_spec, occupancy)) %>%
  rename('Trait evolution' = traitevolution,
         Diversification = diversif,
         Speciation = speciation,
         Extinction = extinction) %>%
  group_by(plasticity, enviroment_type, cost, replications) %>%
  summarize_all(mean) %>%
  pivot_longer(col = -c(plasticity, enviroment_type, cost, replications)) %>%
  mutate(plasticity = as.factor(plasticity)) %>%
  ggplot(aes(x = cost, y = value)) + #cost in x as well
  geom_boxplot() +
  scale_size_manual(values = c("with" = 1, "without" = 0)) +
  #geom_jitter() +
  facet_wrap(~name, scales = "free_y") +
  xlab("Plasticity") +
  #ggtitle('Effect of plasticity on adaptive evolution - climatic changes & cost context') +
  theme_bw() +
  #scale_color_manual(values=c("#542788", "#b35806"),
                   #  name = "climatic change") +
  theme(plot.title =
          element_text(size = 14,
                       face = 2,
                       hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.position = "bottom") 

########################################################
#### ANOVA EXECUTION ####
#### data to analysis ####
path_general <- here('output/files_to_results/stable_slow_envi_without_cost')
#path_general <- here('output/files_to_results/stable_slow_envi_with_cost')
#path_general <- here('output/files_to_results/stable_fast_envi_without_cost')
#path_general <- here('output/files_to_results/stable_fast_envi_with_cost')
path_files <- list.files(path_general, pattern = 'csv')

if(length(path_files) >= 2){
  path_complete <- paste0(path_general, "/", path_files)
  data <- vector('list', length(path_complete))
  
  for(i in seq_along(path_complete)){
    data[[i]] <- read.csv2(path_complete[i])
    assign(paste0('data', i), data[[i]])
  }
  
  all_data <- unlist(eapply(.GlobalEnv, is.data.frame))
  data <- do.call(rbind, mget(names(all_data)[all_data]))
  rownames(data) <- NULL
  data <- data %>% filter(replications != 0)
} else {
  data <- read.csv2("output/rep_1_finalresult.csv")
  data <- data %>% filter(replications != 0)
  path <- here('output')
}
rm(list = ls()[(ls() == 'data') == FALSE])
colnames(data)[1] <- c('plasticity')


plasticity_which_species_die <- data %>%
  filter(alive_spec == 0) %>%
  select(plasticity) %>%
  group_by(plasticity) %>%
  summarise(count = n()) %>%
  filter(count == 10) %>% 
  select(plasticity)

plasticity_which_species_die <- as.numeric(unlist(plasticity_which_species_die))

data_with_alives <- data[!data$plasticity %in% plasticity_which_species_die, ]

data[is.na(data)] <- 0

result_aov <- data %>% 
  as_tibble() %>% 
  filter(timesimulation > 40) %>%
  select(-c(timesimulation, enviroment_type, abundance,
            occupancy, alive_spec)) %>% 
  group_by(plasticity, replications) %>% 
  summarize_all(mean) %>% 
  mutate(plasticity = as.factor(plasticity))

summary(aov(result_aov$extinction ~ result_aov$plasticity))
summary(aov(result_aov$traitevolution ~ result_aov$plasticity))
summary(aov(result_aov$speciation ~ result_aov$plasticity))
summary(aov(result_aov$diversif ~ result_aov$plasticity))

### tukey test ###
TukeyHSD(aov(result_aov$traitevolution ~ result_aov$plasticity))
TukeyHSD(aov(result_aov$speciation ~ result_aov$plasticity))
TukeyHSD(aov(result_aov$extinction ~ result_aov$plasticity))
TukeyHSD(aov(result_aov$diversif ~ result_aov$plasticity))

### cost anova ###
cost_aov <- new_data %>%
  as_tibble() %>%
  filter(timesimulation > 40) %>%
  select(-c(timesimulation, abundance, alive_spec, occupancy)) %>%
  rename('Trait evolution' = traitevolution,
         Diversification = diversif,
         Speciation = speciation,
         Extinction = extinction) %>%
  group_by(plasticity, enviroment_type, cost, replications) %>%
  summarize_all(mean) %>%
  pivot_longer(col = -c(plasticity, enviroment_type, cost, replications)) %>%
  mutate(plasticity = as.factor(plasticity))

spec <- cost_aov[cost_aov$name == "Speciation",]
Ext <- cost_aov[cost_aov$name == "Extinction",]
Div <- cost_aov[cost_aov$name == "Diversification",]
Tra <- cost_aov[cost_aov$name == "Trait evolution",]

summary(aov(spec$value ~ spec$cost))
summary(aov(Ext$value ~ spec$cost))
summary(aov(Div$value ~ spec$cost))
summary(aov(Tra$value ~ spec$cost))

#### BOX PLOT OF MEAN PER PLASTICITY WITH LIST ####
unique_plasticity <- unique(data$plasticity)
unique_replications <- unique(data$replications)

results_extinction <- data.frame(replications = rep(0, length(unique_replications) * length(unique_plasticity)),
                      plasticity = 0,
                      extinction = 0)

results_speciation <- data.frame(replications = rep(0, length(unique_replications) * length(unique_plasticity)),
                      plasticity = 0,
                      speciation = 0)

results_diversification <- data.frame(replications = rep(0, length(unique_replications) * length(unique_plasticity)),
                                     plasticity = 0,
                                     divers = 0)

results_trait <- data.frame(replications = rep(0, length(unique_replications) * length(unique_plasticity)),
                            plasticity = 0,
                            trait = 0)

count <- 0
for(i in seq_along(unique_plasticity)){
  for(j in seq_along(unique_replications)){
    count <- count + 1
    results_extinction$replications[count] <- unique_replications[j]
    results_extinction$plasticity[count] <- unique_plasticity[i]
    temp <- filter(data, plasticity == unique_plasticity[i], replications == unique_replications[j])
    results_extinction$extinction[count] <- temp$extinction %>% mean() 
  }
}

count <- 0
for(i in seq_along(unique_plasticity)){
  for(j in seq_along(unique_replications)){
    count <- count + 1
    results_speciation$replications[count] <- unique_replications[j]
    results_speciation$plasticity[count] <- unique_plasticity[i]
    temp <- filter(data, plasticity == unique_plasticity[i], replications == unique_replications[j])
    results_speciation$speciation[count] <- temp$speciation %>% mean()
  }
}

count <- 0
for(i in seq_along(unique_plasticity)){
  for(j in seq_along(unique_replications)){
    count <- count + 1
    results_diversification$replications[count] <- unique_replications[j]
    results_diversification$plasticity[count] <- unique_plasticity[i]
    temp <- filter(data, plasticity == unique_plasticity[i], replications == unique_replications[j])
    results_diversification$divers[count] <- temp$divers %>% mean()
  }
}

count <- 0
for(i in seq_along(unique_plasticity)){
  for(j in seq_along(unique_replications)){
    count <- count + 1
    results_trait$replications[count] <- unique_replications[j]
    results_trait$plasticity[count] <- unique_plasticity[i]
    temp <- filter(data, plasticity == unique_plasticity[i], replications == unique_replications[j])
    results_trait$trait[count] <- temp$trait %>% mean()
  }
}

ggplot(results_extinction, aes(plasticity, extinction, group = plasticity)) +
  geom_boxplot() +
  geom_jitter()

ggplot(results_speciation, aes(plasticity, speciation, group = plasticity)) +
  geom_boxplot() + 
  geom_jitter()

ggplot(results_diversification, aes(plasticity, divers, group = plasticity)) +
  geom_boxplot() + 
  geom_jitter()

ggplot(results_trait, aes(plasticity, trait, group = plasticity)) +
  geom_boxplot() + 
  geom_jitter()

#### EXPLORING RESULTS ####

# Min of max timesteps #
data %>%
  as.tibble() %>%
  group_by(plasticity, replications) %>%
  summarise(max = max(timesimulation)) %>%
  group_by(plasticity) %>%
  summarise(min = min(max)) %>%
  arrange(min)

# max timesteps #
data %>%
  as.tibble() %>%
  group_by(plasticity) %>%
  summarise(max = max(timesimulation)) %>%
  arrange(max)

# events of speciation, extinction, trait evolution, and diversification #
data %>%
  as.tibble() %>%
  filter(extinction > 0, timesimulation <= 172) %>%
  group_by(plasticity) %>%
  summarise(len = length(extinction) / 10) %>%
  arrange(len)

# mean, per plasticity, maximum of timesteps in all replications #
data %>%
  as.tibble() %>%
  group_by(plasticity, replications) %>%
  summarise(max = max(timesimulation)) %>%
  group_by(plasticity) %>%
  summarise(mean = mean(max)) %>%
  arrange(mean)

# max speciation, extinction, trait, and diversify #
data %>%
  as.tibble() %>%
  filter(timesimulation > 150) %>%
  group_by(plasticity) %>%
  summarise(max = max(speciation)) %>%
  arrange(max)
