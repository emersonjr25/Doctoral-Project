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
dados <- read.csv2("output/finalresult.csv")
dados <- dados %>% filter(replications != 0)
path <- here("output")

#### plot speciation ####
ggplot(data = dados, aes(timesimulation, speciation, colour = plasticidade, group = plasticidade)) + 
  labs(title = paste("Speciation along time with plasticity")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  geom_smooth(se = FALSE)


tiff(filename = file.path(path, paste0("plot", "_", general_plot_1[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 15,
     height = 15,
     units = "cm",
     res = 100)
print(general_plot_1)
dev.off()

for(i in unique(dados$plasticidade)){
  subdata <- dados %>% filter(plasticidade == i)
  tiff(filename = file.path(path, paste0("plot", "_", general_plot_1[["labels"]][["y"]], "_", "plas", i, ".tif")),
       width = 15,
       height = 15,
       units = "cm",
       res = 100)
  print(ggplot(subdata, aes(timesimulation, speciation)) +
          geom_point() + labs(title = paste("Speciation along time with plasticity- all data")) + theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(size = 16, hjust = 0.5),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)))
  dev.off()
}

#### plot extinction ####

ggplot(data = dados, aes(timesimulation, extinction, colour = plasticidade, group = plasticidade)) + 
  labs(title = paste("Extinction along time with plasticity")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  geom_smooth(se = FALSE)

tiff(filename = file.path(path, paste0("plot", "_", general_plot_2[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 15,
     height = 15,
     units = "cm",
     res = 100)
print(general_plot_2)
dev.off()

for(i in unique(dados$plasticidade)){
  subdata <- dados %>% filter(plasticidade == i)
  tiff(filename = file.path(path, paste0("plot", "_", general_plot_2[["labels"]][["y"]], "_", "plas", i, ".tif")),
       width = 15,
       height = 15,
       units = "cm",
       res = 100)
  print(ggplot(subdata, aes(timesimulation, extinction)) +
          geom_point() + labs(title = paste("Extinction along time with plasticity", i)) + theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(size = 16, hjust = 0.5),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)))
  dev.off()
}

#### plot diversification ####

ggplot(data = dados, aes(timesimulation, diversif, colour = plasticidade, group = plasticidade)) + 
  labs(title = paste("Diversification along time with plasticity - all data")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  geom_smooth(se = FALSE)

tiff(filename = file.path(path, paste0("plot", "_", general_plot_3[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 15,
     height = 15,
     units = "cm",
     res = 100)
print(general_plot_3)
dev.off()

for(i in unique(dados$plasticidade)){
  subdata <- dados %>% filter(plasticidade == i)
  tiff(filename = file.path(path, paste0("plot", "_", general_plot_3[["labels"]][["y"]], "_", "plas", i, ".tif")),
       width = 15,
       height = 15,
       units = "cm",
       res = 100)
  print(ggplot(subdata, aes(timesimulation, diversif)) +
          geom_point() + labs(title = paste("Diversification along time with plasticity", i)) + theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(size = 16, hjust = 0.5),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)))
  dev.off()
}

#### plot trait evolution ####
#subdata_final <- dados %>% filter(timesimulation >= 100)

ggplot(data = dados, aes(timesimulation, traitevolution, colour = plasticidade, group = plasticidade)) + 
  labs(title = paste("Trait evolution along time with plasticity - all data")) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  geom_smooth(se = FALSE)

tiff(filename = file.path(path, paste0("plot", "_", general_plot_4[["labels"]][["y"]], "_", "plas", "all", ".tif")),
     width = 15,
     height = 15,
     units = "cm",
     res = 100)
print(general_plot_4)
dev.off()

for(i in unique(dados$plasticidade)){
  subdata <- dados %>% filter(plasticidade == i)
  tiff(filename = file.path(path, paste0("plot", "_", general_plot_4[["labels"]][["y"]], "_", "plas", i, ".tif")),
       width = 15,
       height = 15,
       units = "cm",
       res = 100)
  print(ggplot(subdata, aes(timesimulation, traitevolution)) +
          geom_point() + labs(title = paste("Trait evolution along time with plasticity", i)) + theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                plot.title = element_text(size = 16, hjust = 0.5),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)))
  dev.off()
}
