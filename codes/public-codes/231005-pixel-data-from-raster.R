# Function to get raster pixel coordinates, values, and region names from a shapefile

get_raster_coordinates_and_regions <- function(raster_data, shapefile_data) {
  
  # Convert raster data to sf point dataframe with longitude and latitude
  raster_points <- rasterToPoints(raster_data) %>% data.frame() 
  
  raster_points_sp <- st_as_sf(raster_points, coords = c("x", "y"), crs = st_crs(shapefile_data)) %>% as("Spatial")

  # Find corresponding region for each point
  points_regions <- sp::over(raster_points_sp, shapefile_data)

  # Add coordinates, values, and NAMES to the data frame
  result_df <- points_regions %>% 
    select(NAME_1, NAME_2, NAME_3) %>%
    mutate(longitude = raster_points[,2], latitude = raster_points[,3], values = raster_points[,1])


  return(result_df)
}

#+++FUNCTION USAGE+++

# Load necessary libraries
library(raster)
library(sf)
library(geodata)
library(tidyverse)
library(sp)

# Set the directory for data
output_dir <- "data/"

# Load raster data
elevation_ETH <- elevation_30s(country = "ETH", path = output_dir)

elevation_ETH_raster <- raster(elevation_ETH)

# Load shapefile data
ethiopia_boundary <- gadm("ETH", level = 3, path = output_dir) %>% st_as_sf() %>% as("Spatial") 

# Example usage of the function with the provided inputs
raster_data <-  elevation_ETH_raster # Replace with your raster file
shapefile_data <- ethiopia_boundary      # Replace with your shapefile

output_data <- get_raster_coordinates_and_regions(raster_data, shapefile_data)
print(output_data)
