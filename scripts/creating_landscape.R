# #####################################################
# ### Main goal: Verify the effect of plasticity on adaptive evolution #####################
# ##### Diversification and Trait Evolution ~ Plasticity
# ##### Methods: Computational simulation #############
# ##### Script to create a landscape with bigger timesteps #########
# 
# #### PACKAGES ####
# library(gen3sis)
# library(raster)
# library(here)
# 
# #### CREATING LANDSCAPES WITH SPECIFIC TIMESTEP TO USE ####
# datapath <- system.file(file.path("extdata", "WorldCenter"), package="gen3sis")
# 
# # create raster bricks
# temperature_brick <- brick(file.path(datapath, "input_rasters/temp_rasters.grd"))
# aridity_brick <- brick(file.path(datapath, "input_rasters/arid_rasters.grd"))
# area_brick <- brick(file.path(datapath, "input_rasters/area_rasters.grd"))
# 
# s_temp <- stack(temperature_brick, temperature_brick, temperature_brick, temperature_brick,
#                 temperature_brick, temperature_brick, temperature_brick, temperature_brick,
#                 temperature_brick, temperature_brick, temperature_brick, temperature_brick,
#                 temperature_brick, temperature_brick, temperature_brick, temperature_brick,
#                 temperature_brick, temperature_brick, temperature_brick, temperature_brick)
# 
# s_aridity <- stack(aridity_brick, aridity_brick, aridity_brick, aridity_brick,
#                    aridity_brick, aridity_brick, aridity_brick, aridity_brick,
#                    aridity_brick, aridity_brick, aridity_brick, aridity_brick,
#                    aridity_brick, aridity_brick, aridity_brick, aridity_brick,
#                    aridity_brick, aridity_brick, aridity_brick, aridity_brick)
# 
# s_area <- stack(area_brick, area_brick, area_brick, area_brick,
#                 area_brick, area_brick, area_brick, area_brick,
#                 area_brick, area_brick, area_brick, area_brick,
#                 area_brick, area_brick, area_brick, area_brick,
#                 area_brick, area_brick, area_brick, area_brick)
# 
# nlayers(s_temp)
# nlayers(s_aridity)
# nlayers(s_area)
# 
# landscapes_list <- list(temp=NULL, arid=NULL, area=NULL)
# 
# for(i in 1:nlayers(s_temp)){
#   landscapes_list$temp <- c(landscapes_list$temp, s_temp[[i]])
#   landscapes_list$arid <- c(landscapes_list$arid, s_aridity[[i]])
#   landscapes_list$area <- c(landscapes_list$area, s_area[[i]])
# }
# # define cost function, crossing water as double as land sites
# cost_function_water <- function(source, habitable_src, dest, habitable_dest) {
#   if(!all(habitable_src, habitable_dest)) {
#     return(2/1000)
#   } else {
#     return(1/1000)
#   }
# }
# 
# create_input_landscape(
#   landscapes = landscapes_list,
#   cost_function = cost_function_water,
#   output_directory = file.path(here("data", "raw", "landscape_new")),
#   directions = 8, # surrounding sites for each site
#   timesteps = paste0(round(1020:1,2), "Ma"),
#   crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
#   calculate_full_distance_matrices = FALSE) # full distance matrix
