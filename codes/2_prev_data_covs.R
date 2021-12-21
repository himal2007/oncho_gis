setwd("C:/Users/user/OneDrive - LA TROBE UNIVERSITY/oncho_git/onchocerciasis/espen_ethiopia_data/analysis/210723_analysis_writing/")
# Import libraries --------------------------------------------------------


suppressMessages({
  library(plotly);
  library(rgeos);
  library(raster);
  library(ggplot2);
  library(rgdal);
  library(leaflet);
  library(sf);
  library(dplyr);
  library(climateStability);
                 })

# Import data -------------------------------------------------------------
eth_prev <- read.csv("data/eth_espen_nodule.csv")
eth_prev2 <- eth_prev %>% dplyr::select(CASES, N, LONG, LAT)
eth_prev2 <- eth_prev2 %>% filter(LONG !="") %>% na.omit()
# eth_st <- getData(name = "GADM", country = "ETH", level = 3)
# eth <- eth_st %>% st_as_sf()
covariates <- stack("../../data/covariates/covariates_5km_new.grd")
# covariates_normalised <- stack("../../data/covariates/covariates_1km_normalised.grd")

# Rescale covariates ------------------------------------------------------
# divide by 10
temp_covariates <- c("annual_mean_temp", "annual_diurnal_range", "maxtemp_warmest_mnth",
                     "mintemp_coldest_mnth", "temp_annual_range", "meantemp_wettestquart",
                     "mean_temp_driestquart", "meantemp_warmestquarter", "meantemp_coldestquart")
covariates_temp <- covariates[[temp_covariates]]/10

# divide by 1000
covs_1000 <- c("temp_seasonality", "dist_river_DIVA", "dist_river")
covstack_1000 <- covariates[[covs_1000]]/1000

# covariates as is
covs_as_is <- c("alt", "isothermality", "dist_river_WP", c(names(covariates)[13:20]), 
                "wetness_index", "Soil_moisture", "slope", "Population_density", "night_lights") 
covstack_as_is <- covariates[[covs_as_is]]

# rescaled to percentile
covs_percent <- c("NDVI_2003_11", "Flow_accumulation")
NDVI_2003_11 <- rescale0to1(covariates[["NDVI_2003_11"]]) * 100
Flow_accumulation <- rescale0to1(covariates[["Flow_accumulation"]]) * 100
mean_housing_2000_15 <- covariates[["mean_housing_2000_15"]] * 100

covariates_final <- stack(covariates_temp, covstack_1000, covstack_as_is, NDVI_2003_11,
                          Flow_accumulation, mean_housing_2000_15)

# Prepare selected covariates ---------------------------------------------

selected_covs <- c("alt", "annual_diurnal_range", "isothermality",
                   "precp_wettest_quart", "precp_seasonality", "precp_warmest_quart",
                   "precp_coldest_quarter", "NDVI_2003_11",
                   "dist_river_DIVA", "Flow_accumulation", "wetness_index",
                   "slope", "Soil_moisture", "Population_density", "mean_housing_2000_15", "night_lights")

selected_covariates <- covariates_final[[selected_covs]]  -> covariates

# Extract the data from the raster--------------------------------------------------------
pred_eth <- raster::extract(covariates, eth_prev2[c("LONG","LAT")], na.rm = TRUE, df = TRUE)
# pred_eth2 <- raster::extract(covariates, eth_prev[c("LONG","LAT")], buffer = 5000, fun = max, df = TRUE)
# pred_eth3 <- raster::extract(covariates, eth_prev[c("LONG","LAT")], buffer = 5000, fun = mode, df = TRUE)
eth_prev2 <- as.data.frame(cbind(eth_prev2, pred_eth))
# eth_prev2 <- na.omit(eth_prev2)
# write.csv(eth_prev2, "data/eth_espen_nodule_covs.csv")

# save(eth_prev2, selected_covs, selected_covariates, file = "data/selected_covs_files.RData")
