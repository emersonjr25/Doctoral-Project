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

#### data ####
path_general <- here('output/files_to_results/')
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

# data %>%
#   filter(replications == 10) %>%
#   select(plasticity, timesimulation) %>%
#   group_by(plasticity) %>%
#   summarise(max = max(timesimulation)) %>%
#   arrange(max)

# data %>%
#   filter(extinction > 0, plasticity == 0.25, timesimulation <= 190) %>%
#   summarise(mean(extinction))

#### RESULT FINAL USING TIDYVERSE ####
result <- data %>% 
    as_tibble() %>% 
    filter(timesimulation > 100) %>%
    select(-timesimulation) %>% 
    rename('Trait evolution' = traitevolution,
           Diversification = diversif,
           Speciation = speciation,
           Extinction = extinction) %>%
    group_by(plasticity, replications) %>% 
    summarize_all(mean) %>% 
    pivot_longer(col = -c(plasticity, replications)) %>%
    mutate(plasticity = as.factor(plasticity)) %>% 
    ggplot(aes(x = plasticity, y = value)) + 
    geom_boxplot() + 
    geom_jitter() +
    facet_wrap(~name, scales = "free_y") + 
    xlab("Plasticity") + ylab('Value') +
    ggtitle('Effect of plasticity on adaptive evolution') + 
    theme_bw() +
    theme(plot.title = 
            element_text(size = 16, 
                        face = 2, 
                        hjust = 0.5), 
            axis.title.x = element_text(size = 14), 
            axis.title.y = element_text(size = 14))
    

tiff(filename = file.path(here('output'), paste0("plot", "_", result[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 800,
     height = 600,
     units = "px",
     res = 100)
print(result)
dev.off()

#### ANOVA EXECUTION ####
result_aov <- data %>% 
  as_tibble() %>% 
  filter(timesimulation > 100) %>%
  select(-timesimulation) %>% 
  group_by(plasticity, replications) %>% 
  summarize_all(mean) %>% 
  mutate(plasticity = as.factor(plasticity))

summary(aov(result_aov$speciation ~ result_aov$plasticity))
summary(aov(result_aov$extinction ~ result_aov$plasticity))
summary(aov(result_aov$diversif ~ result_aov$plasticity))
summary(aov(result_aov$traitevolution ~ result_aov$plasticity))

# manova_result <- manova(cbind(traitevolution, 
#             extinction,
#             diversif,
#             speciation) ~ plasticity,
#        data = result_aov)
# 
# summary(manova_result, tol = 0)

### tukey test ###
TukeyHSD(aov(result_aov$speciation ~ result_aov$plasticity))
TukeyHSD(aov(result_aov$extinction ~ result_aov$plasticity))
TukeyHSD(aov(result_aov$diversif ~ result_aov$plasticity))
TukeyHSD(aov(result_aov$traitevolution ~ result_aov$plasticity))

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
