setwd("C:/Users/User/OneDrive - LA TROBE UNIVERSITY/oncho_git/onchocerciasis/espen_ethiopia_data/analysis/210723_analysis_writing/")

# Import libraries --------------------------------------------------------
suppressMessages({
  library(INLA);
  library(INLAutils);
  library(raster);
  library(sf);
  library(dplyr);
  library(tictoc);
  library(leaflet);
  library(ggregplot);
  library(tidyverse);
  library(tmap);
                 })

# Custom functions --------------------------------------------------------
source("../../data/custom_functions.R")

# Import data -------------------------------------------------------------
load("data/selected_covs_files_1km.RData")
data <- eth_prev2
data$CASES <- round(data$CASES)
data$N <- round(data$N)

m <- getData(name = "GADM", country = "ETH", level = 0)
covariates <- selected_covariates
varlist <- c("slope", "isothermality", "precp_seasonality", "NDVI_2003_11",
             "dist_river_DIVA", "Population_density")# , "Flow_accumulation", "Soil_moisture") # FAc and SM not included because of their WAIC and DIC scores
# pred_data <- covariates[[varlist]]
pred_data <- aggregate(covariates[[varlist]], fact = 5, fun = mean, na.rm = TRUE)
# save(pred_data, file = "data/pred_covs.RData")

ra <- pred_data
dp <- data.frame(rasterToPoints(ra))

# Mesh construction -------------------------------------------------------

coords <-  cbind(data$LONG, data$LAT)
bdry <- inla.sp2segment(m)
bdry$loc <- inla.mesh.map(bdry$loc)
mesh1 <- inla.mesh.2d(
  loc = coords, boundary = bdry, max.edge = c(0.5, 5),
  cutoff = 0.03 # mesh 8, initially 0.05 for mesh 5
)
mesh1$n

spde <- inla.spde2.matern(mesh1, alpha=2)
indexs <- inla.spde.make.index(name = "spatial.field", spde$n.spde)
length(indexs)

# Make projection matrix --------------------------------------------------

A <- inla.spde.make.A(mesh=mesh1,loc=as.matrix(coords));dim(A)
coop <- cbind(dp$x, dp$y)
Ap <- inla.spde.make.A(mesh = mesh1, loc = coop);dim(Ap)

# Model formula ---------------------------------------------------------
predterms <- as.formula(paste("y ~ 0 + b0 +", paste(varlist, collapse =  "+"), "+ f(spatial.field, model = spde)"))

datastack <- stack_data(data = data, dp = dp, cov_list = varlist)

# for fitting only
# datastack <- stack_data_fit(data= data, cov_list = varlist)

tic()
res <- inla(predterms,
            family = "zeroinflated.binomial.1", Ntrials = numtrials,
            data = inla.stack.data(datastack, spde = spde),
            control.family = list(link = "logit"),
            control.compute = list(dic = TRUE, waic = TRUE,
                                   cpo = TRUE, config = TRUE,
                                   openmp.strategy="huge"),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(datastack)
            )
)
toc()

# saveRDS(res, "docs/210902_inla_point_res_5km_mesh8_no_FAc_SM.rds")

# Model results -----------------------------------------------------------
res <- readRDS("docs/210902_inla_point_res_5km_mesh8_no_FAc_SM.rds")
summary(res)
reg_coff <- res$summary.fixed %>% arrange(desc(mean))
reg_coff %>% write.csv("docs/210902_reg_coff.csv")
res$summary.hyperpar %>% data.frame() # %>% write.csv("docs/hyper_par_theta.csv")
exp(res$summary.fixed)

p <- autoplot(res, plot.prior = TRUE)
# INLAutils::autoplot(res, priors = TRUE) # if INLAutils installed correctly
p[[1]] + theme_bw(base_family = "Verdana") + geom_vline(xintercept = 0, color = "blue", linetype = "dashed") + ylab("Probability Density") + xlab("Values")
p[[2]]

Efxplot(list(res))

# Mapping prediction ------------------------------------------------------
index <- inla.stack.index(stack = datastack, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]
prev_sd <-  res$summary.fitted.values[index, "sd"]
prev_med <- res$summary.fitted.values[index, "0.5quant"]

summary(res$summary.fitted.values[index,])

r <- getData(name = "alt", country = "ETH", mask = T)
raster_prev <- function(prev_mean){
  r_prev_mean <- rasterize(
    x = coop, y = ra, field = prev_mean,
    fun = mean)
  # r_prev_mean <- reproj(r_prev_mean,r)
}

r_prev_mean <- raster_prev(prev_mean = prev_mean)
r_prev_ul <- raster_prev(prev_ul)
r_prev_ll <- raster_prev(prev_ll)
r_pred_error <- (r_prev_ul - r_prev_ll)/2
r_prev_sd <- raster_prev(prev_sd)
r_prev_med <- raster_prev(prev_med)


# Export predicted raster -------------------------------------------------

# pred_raster <- stack(r_prev_mean, r_prev_ul, r_prev_ll, r_prev_sd)
# names(pred_raster) <- c("Mean_prevalence", "Upper_Limit", "Lower_Limit", "Prevalence_SD")
# writeRaster(pred_raster, "output/pred_raster.grd")

plotraster(r_prev_mean * 100, "Mean Prevalence")
p1 <- plotraster(r_prev_ul*100, "Upper Limit")
p2 <- plotraster(r_prev_ll*100, "Lower Limit")

spplot(r_prev_mean*100, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))
spplot(r_prev_med*100, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))

spplot(r_pred_error, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))
spplot(r_prev_sd, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))
spplot(r_prev_ul*100, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))
spplot(r_prev_ll*100, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))


# Plot predicted raster ---------------------------------------------------



# Mapping exceedance probabilities ----------------------------------------

index <- inla.stack.index(stack = datastack, tag = "pred")$data
marg <- res$marginals.fitted.values[index][[1]]
1 - inla.pmarginal(q = 0.20, marginal = marg)
excprob_10 <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.10, marginal = marg)})
raster_excprob_10 <- raster_prev(excprob_10)
spplot(raster_excprob_10, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))


excprob_35 <- sapply(res$marginals.fitted.values[index],
                     FUN = function(marg){1-inla.pmarginal(q = 0.35, marginal = marg)})
raster_excprob_35 <- raster_prev(excprob_35)
spplot(raster_excprob_35, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))

### Hyper endemic
excprob_60 <- sapply(res$marginals.fitted.values[index],
                     FUN = function(marg){1-inla.pmarginal(q = 0.60, marginal = marg)})
raster_excprob_60 <- raster_prev(excprob_60)
spplot(raster_excprob_60, col.regions = hcl.colors(1000, palette = "YlOrRd", rev = TRUE, alpha = .7))

# Export exceedance probability raster -------------------------------------------------

exc_raster <- stack(raster_excprob_10, raster_excprob_35, raster_excprob_60)
names(exc_raster) <- c("exc_prob_hypo", "exc_prob_meso", "exc_prob_hyper")
# writeRaster(exc_raster, "docs/exc_raster.grd", overwrite = TRUE)

# Plot the effect of the covariates ---------------------------------------

data_p <- dp[, varlist]
data_p$pred_prev <- prev_mean * 100
data_p$pred_prev_UL <- prev_ul
data_p$pred_prev_LL <- prev_ll
data_p$pred_prev_med <- prev_med

data_p_long <- data_p %>% pivot_longer(cols = names(data_p)[1:8],
                                       names_to = "covariates",
                                       values_to = "values") %>% 
  mutate(covariates = as_factor(covariates))    

levels(data_p_long$covariates) <- c("Slope", "Isothermality", "Precipitation seasonality",
                                    "Vegetation index", "Distance to the nearest river (km)", "Flow accumulation", "Soil moisture (mm)",
                                    "Population density")
p <- ggplot(data_p_long, aes(values, pred_prev)) +
  geom_smooth(method = "gam") + facet_wrap(~covariates, scales = "free", nrow = 2) + # default method is gam
  theme_bw(base_family = "Verdana", base_size = 12) +
  xlab("Values") + ylab("Onchocerciasis prevalence") + ylim(c(-20, 50))

# ggsave(plot = p, filename = "docs/var_relation_2.png", device = "png", dpi = 1000, width = 12, height = 6, units = "in" )

# Plot the spatial field - source: our coding club ------------------------
points.em <- mesh1$loc

stepsize <- 0.008                           # This is given in coordinates unit (in this case this is straightforward and correspond to 1km)
east.range <- diff(range(points.em[,1]))  # calculate the length of the Easting range
north.range <- diff(range(points.em[,2])) # calculate the length of the Northing range

nxy <- round(c(east.range, north.range)/stepsize)  # Calculate the number of cells in the x and y ranges

# Project the spatial field on the mesh vertices using the inla.mesh.projector() function
projgrid <- inla.mesh.projector(mesh1,
                                xlim = range(points.em[,1]),
                                ylim = range(points.em[,2]),
                                dims = nxy)
xmean <- inla.mesh.project(projgrid,
                           res$summary.random$spatial.field$mean)
xsd <- inla.mesh.project(projgrid,
                         res$summary.random$spatial.field$sd)

xmean2 <- t(xmean)
xmean3 <- xmean2[rev(1:length(xmean2[,1])),]
xmean_ras <- raster(xmean3,
                    xmn = range(projgrid$x)[1], xmx = range(projgrid$x)[2],
                    ymn = range(projgrid$y)[1], ymx = range(projgrid$y)[2],
                    crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))

xsd2 <- t(xsd)
xsd3 <- xsd2[rev(1:length(xsd2[,1])),]
xsd_ras <- raster(xsd3,
                  xmn = range(projgrid$x)[1], xmx =range(projgrid$x)[2],
                  ymn = range(projgrid$y)[1], ymx =range(projgrid$y)[2],
                  crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
xmean_ras2 <- mask(xmean_ras, m)
xsd_ras2 <- mask(xsd_ras, m)

library(RColorBrewer)
my.palette.post <- rev(brewer.pal(n = 9, name = "YlGnBu"))
my.palette.var <- brewer.pal(n = 9, name = "BuPu")
my.palette.heat <- rev(brewer.pal(n = 9, name = "RdYlBu"))
my.palette.blues <- (brewer.pal(n = 9, name = "Blues"))

par(mfrow = c(1,2), mar = c(2,2, 1,1))
plot(xmean_ras2, asp = 1, col = my.palette.heat, main = "Spatial effect (Mean)")
# points(coords, pch = 16, cex = 0.5)
plot(m, add = T)

plot(xsd_ras2, asp = 1, col = my.palette.var, main = "Spatial effect (SD)")
# points(coords, pch = 16, cex = 0.5)
plot(m, add = T)


## Plot spatial effects ----------------------------------------------------


p_mean <- tm_shape(xmean_ras2) + 
  tm_raster(title="Spatial field (Mean)", alpha = 1, palette = my.palette.heat, legend.is.portrait = T, style = "cont", midpoint = NA)+
  tm_shape(m)+
  tm_borders(lty = 1, lwd = 2, col = "black", alpha = .9) +
  tm_compass(type = "arrow", position = c(.1, .75)) +
  # tm_scale_bar(breaks = c(0, 100, 200, 400), text.size = .8, position = c(.6, .1)) +
  tm_layout(legend.outside = FALSE, legend.position = c(.61, .65), frame = FALSE, fontfamily = "Verdana")
p_mean

p_var <- tm_shape(xsd_ras2) + 
  tm_raster(title="Spatial field (SD)", alpha = 1, palette = my.palette.var, legend.is.portrait = T, style = "cont", midpoint = NA)+
  tm_shape(m)+
  tm_borders(lty = 1, lwd = 2, col = "black", alpha = .9) +
  # tm_compass(type = "arrow", position = c(.1, .75)) +
  tm_scale_bar(breaks = c(0, 100, 200, 400), text.size = .8, position = c(.6, .12)) +
  tm_layout(legend.outside = FALSE, legend.position = c(.61, .65), frame = FALSE, fontfamily = "Verdana")
p_var
spatial_field <- tmap_arrange(p_mean, p_var, ncol = 2, widths = c(1,1), heights = 1)
tmap_save(tm = spatial_field, filename = "docs/spatial_field.png", dpi = 1000, width = 10, height = 6, units = "in" )

# Calculate range for the spatial field -----------------------------------

spde.est <- inla.spde2.result(inla = res, name = "spatial.field", spde = spde,
                              do.transf = T)
#Kappa
inla.zmarginal(spde.est$marginals.kappa[[1]])
#variance
inla.zmarginal(spde.est$marginals.variance.nominal[[1]])
# inla.emarginal(function(x) x, spde.est$marginals.variance.nominal[[1]])
# inla.hpdmarginal(0.95, spde.est$marginals.variance.nominal[[1]])
#Range
inla.zmarginal(spde.est$marginals.range.nominal[[1]])
# inla.emarginal(function(x) x, spde.est$marginals.range.nominal[[1]])  
# inla.hpdmarginal(0.95, spde.est$marginals.range.nominal[[1]])

##Punama method
# spde.est$summary.log.range.nominal
# spde.est$summary.log.variance.nominal
# plot(spde.est$marginals.range.nominal[[1]], type = "l", xlab = "range nominal", ylab = "Density")

hyper_par <- data.frame(rbind(inla.zmarginal(spde.est$marginals.kappa[[1]])
, inla.zmarginal(spde.est$marginals.range.nominal[[1]])
, inla.zmarginal(spde.est$marginals.variance.nominal[[1]])))
hyper_par$parameters <- c("Kappa", "Range", "Variance")
write.csv(as.matrix(hyper_par), "docs/hyper_par.csv")
save(mesh1, r_prev_mean, r_prev_sd, data_p, file = "docs/INLA_model_run_5km.RData")


# Observed and predicted values -------------------------------------------
data$PREV <- data$CASES/data$N  *100
pred_stack <- stack(r_prev_mean, r_prev_med, r_prev_ll, r_prev_ul)
names(pred_stack) <- c("predicted_mean_prev", "pred_median_prev", "pred_LL_prev", "pred_UL_prev")

pred_eth <- raster::extract(pred_stack, eth_prev2[c("LONG","LAT")], na.rm = TRUE, df = TRUE)


data2 <- as.data.frame(cbind(data, pred_eth)) 
data2$predicted_mean_prev <- data2$predicted_mean_prev *100
cor(data$PREV, pred_eth$predicted_mean_prev, method = "spearman")
cor(data$PREV, pred_eth$pred_median_prev)

data2 <- data2 %>% dplyr::select(PREV, predicted_mean_prev,
                                                         pred_median_prev, pred_LL_prev,
                                                         pred_UL_prev) %>% mutate(abs_difference = (abs(PREV-predicted_mean_prev))) # %>% 
# %>% 
  pivot_longer(cols = predicted_mean_prev:pred_UL_prev,
               names_to = "variables",
               values_to = "values") %>% mutate(variables = as_factor(variables),
                                                abs_difference = GGally::rescale01(abs(PREV-values)))
library(RColorBrewer)
mypalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

p <-ggplot(data2) +
  geom_point(mapping = aes(PREV , predicted_mean_prev, col = abs_difference)) +
  geom_abline(slope=1, intercept=0, linetype = "dashed") +
  ggpubr::stat_cor(mapping = aes(PREV, predicted_mean_prev), data = data2, method = "spearman", label.x = 50, label.y = 5) +
  scale_color_gradientn(colours = mypalette(100)) + 
  theme_bw(base_family = "Verdana", base_size = 12) +
  xlab("Observed prevalence(%)") +
  ylab("Predicted prevalence (%)") + labs(color = "Absolute error") +
  theme(legend.position = "bottom", legend.key.size = unit(1, 'cm'),                                                                       #change legend key size
        legend.key.height = unit(.5, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm')) #+ facet_wrap(~variables)

ggsave(plot = p, filename = "docs/corr_obs_pred.png", device = "png", dpi = 1000, width = 5, height = 5, units = "in" )



# Variogram ---------------------------------------------------------------
library(gstat)
coordinates(data) <- ~LONG + LAT
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")

gs <- gstat(formula = PREV ~ 1, locations = data)
v <- variogram(gs, width = 20)
head(v)
plot(v)

fve <- fit.variogram(v, vgm(c("Mat")), fit.kappa = TRUE)
fve
plot(variogramLine(fve, 1000), type = "l", ylim = c(0, .03))
points(v[, 2:3], pch = 20, col = "red")

plot(v, fve)
