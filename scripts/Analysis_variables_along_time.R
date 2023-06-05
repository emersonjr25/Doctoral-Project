#####################################################
########### PLASTICITY AND EVOLUTION STUDY ##########
#####################################################
### Main goal: Verify the effect of plasticity on adaptive evolution #####################
##### Diversification and Trait Evolution ~ Plasticity
##### Methods: Computational simulation #############
##### Script to analysis results from simulation ####

#### packages ####

library(dplyr)
library(ggplot2)
library(here)

#### data ####
path_general <- here('output/')
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
  data <- read.csv2("output/rep_4_finalresult.csv")
  data <- data %>% filter(replications != 0)
  path <- here('output')
}
colnames(data)[1] <- c('plasticity')

#### VISUALIZATION GRAPH BEFORE SAVE ####
colors <- c('#e7e1ef','#d4b9da','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')
ggplot(data = data, aes(timesimulation, speciation, colour = plasticity, group = plasticity)) + 
  labs(title = paste("Speciation along time with plasticity")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_color_gradientn(colours = colors, na.value = NA) +
  geom_smooth(se = FALSE)

ggplot(data = data, aes(timesimulation, extinction, colour = plasticity, group = plasticity)) + 
  labs(title = paste("Extinction along time with plasticity")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_color_gradientn(colours = colors, na.value = NA) +
  geom_smooth(se = FALSE)

ggplot(data = data, aes(timesimulation, diversif, colour = plasticity, group = plasticity)) + 
  labs(title = paste("Diversification along time with plasticity - all data")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_color_gradientn(colours = colors, na.value = NA) +
  geom_smooth(se = FALSE)

ggplot(data = data, aes(timesimulation, traitevolution, colour = plasticity, group = plasticity)) + 
  labs(title = paste("Trait evolution along time with plasticity - all data")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_color_gradientn(colours = colors, na.value = NA) +
  geom_smooth(se = FALSE)


#### plot speciation ####
colors <- c('#e7e1ef','#d4b9da','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')
general_plot_1 <- ggplot(data = data, aes(timesimulation, speciation, colour = plasticity, group = plasticity)) + 
  labs(title = paste("Speciation along time with plasticity")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_color_gradientn(colours = colors, na.value = NA) +
  geom_smooth(se = FALSE)


tiff(filename = file.path(path, paste0("plot", "_", general_plot_1[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 15,
     height = 15,
     units = "cm",
     res = 100)
print(general_plot_1)
dev.off()

#### plot extinction ####

general_plot_2 <- ggplot(data = data, aes(timesimulation, extinction, colour = plasticity, group = plasticity)) + 
  labs(title = paste("Extinction along time with plasticity")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_color_gradientn(colours = colors, na.value = NA) +
  geom_smooth(se = FALSE)

tiff(filename = file.path(path, paste0("plot", "_", general_plot_2[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 15,
     height = 15,
     units = "cm",
     res = 100)
print(general_plot_2)
dev.off()


#### plot diversification ####

general_plot_3 <- ggplot(data = data, aes(timesimulation, diversif, colour = plasticity, group = plasticity)) + 
  labs(title = paste("Diversification along time with plasticity - all data")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_color_gradientn(colours = colors, na.value = NA) +
  geom_smooth(se = FALSE)

tiff(filename = file.path(path, paste0("plot", "_", general_plot_3[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 15,
     height = 15,
     units = "cm",
     res = 100)
print(general_plot_3)
dev.off()

#### plot trait evolution ####

general_plot_4 <- ggplot(data = data, aes(timesimulation, traitevolution, colour = plasticity, group = plasticity)) + 
  labs(title = paste("Trait evolution along time with plasticity - all data")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  scale_color_gradientn(colours = colors, na.value = NA) +
  geom_smooth(se = FALSE)

tiff(filename = file.path(path, paste0("plot", "_", general_plot_4[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 15,
     height = 15,
     units = "cm",
     res = 100)
print(general_plot_4)
dev.off()
