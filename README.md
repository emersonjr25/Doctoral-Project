# Plasticity and Evolution Study Readme

## Article
This work has been published in the journal Ecological Modelling (https://doi.org/10.1016/j.ecolmodel.2024.110983). I invite you to read the conclusions of this article.

## Overview
This code is part of a computational simulation aimed at investigating the effect of plasticity on adaptive evolution. The study focuses on diversification and trait evolution in the context of plasticity. The simulation utilizes the `gen3sis` package and employs methods to input, modify, and run the model.

## Main Goal
The primary objective of this study is to verify the impact of plasticity on adaptive evolution. By exploring diversification and trait evolution through computational simulations, the code aims to provide insights into the dynamic relationship between plasticity and evolutionary processes.

## Dependencies
Ensure that the following R packages are installed before running the script:
- `gen3sis`
- `here`
- `dplyr`
- `raster`

You can install them using the following commands:
```R
install.packages("gen3sis")
install.packages("here")
install.packages("dplyr")
install.packages("raster")

# Simulation Configuration
The simulation is configured through the following parameters:
- `datapath`: Path to the data directory containing input files.
- `config`: Path to the configuration file (`config_worldcenter.R`) specifying simulation parameters.
- `landscape`: Path to the landscape data used in the simulation.
- `output_directory`: Output directory for simulation results.
- `timestep_restart`: Timestep at which to restart the simulation.
- `save_state`: Option to save the simulation state.
- `call_observer`: Type of observer to call during the simulation.
- `enable_gc`: Option to enable garbage collection.
- `verbose`: Verbosity level for simulation output.
Set these parameters according to your setup.

# PACKAGES
library(gen3sis)
library(here)
library(dplyr)
library(raster)

# SIMULATION
CARRYING CONFIGURATIONS AND PATHS
datapath <- here("data/raw/WorldCenter")
attach(loadNamespace('gen3sis'), name = 'gen3sis_all')
config = file.path(datapath, "config/config_worldcenter.R")
landscape = file.path(datapath, "landscape_new")
output_directory = NA
timestep_restart = NA
save_state = NA
call_observer = "all"
enable_gc = FALSE
verbose = 1

# Getting Started
1. Install the required R packages as mentioned in the "Dependencies" section.
2. Set the paths and configuration parameters in the script according to your setup.
3. Run all the script in file "Doctoral_project_chapter_1_definitive_simulations.R" in folder "scripts" to initiate the simulation.

# Notes
- Make sure that the required input files and data are available in the specified paths.

# Contact
For any inquiries or issues related to the code, please contact emersonjr25@hotmail.com.

# Acknowledgments
This code utilizes the `gen3sis` package and acknowledges its contribution to the simulation study.

Happy simulating! 🌐✨