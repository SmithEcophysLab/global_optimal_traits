# global_optimal_traits.R
## script to predict global optimal leaf traits
## still a bit of work to do!

## load functions
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## load libraries
library(raster)
library(RColorBrewer)
library(maps)
library(tidyverse)
library(factoextra)

## source optimality model and related functions
source('../optimal_vcmax_r/calc_optimal_vcmax.R')
sourceDirectory('../optimal_vcmax_r/functions')

## test optimality model
calc_optimal_vcmax()

## read in global environmental data
cru_growingseason_data <- read.csv('data/cru_growingseason/cru_growingseason.csv')[,-1]
global_z_data <- read.csv('data/watch_elevation/z_globe.csv')[,-1]

## read in phenology data

## combine datasets by lat/lon
global_data_all <- left_join(cru_growingseason_data, global_z_data)

## subset out places with tmp < 10C (lowest leaf temp for Bernacchi 2003)
nrow(subset(global_data_all, tmp > 10))/nrow(global_data_all) # 67% of all data
global_data <- subset(global_data_all, tmp > 10)

## read in MODIS land cover data (to filter out non-veg sites)
modis_2001 = raster('data/landCoverMODIS/LC_hd_global_2001.tif')
modis_2002 = raster('data/landCoverMODIS/LC_hd_global_2002.tif')
modis_2003 = raster('data/landCoverMODIS/LC_hd_global_2003.tif')
modis_2004 = raster('data/landCoverMODIS/LC_hd_global_2004.tif')
modis_2005 = raster('data/landCoverMODIS/LC_hd_global_2005.tif')
modis_2006 = raster('data/landCoverMODIS/LC_hd_global_2006.tif')
modis_2007 = raster('data/landCoverMODIS/LC_hd_global_2007.tif')
modis_2008 = raster('data/landCoverMODIS/LC_hd_global_2008.tif')
modis_2009 = raster('data/landCoverMODIS/LC_hd_global_2009.tif')
modis_2010 = raster('data/landCoverMODIS/LC_hd_global_2010.tif')
modis_2011 = raster('data/landCoverMODIS/LC_hd_global_2011.tif')
modis_2012 = raster('data/landCoverMODIS/LC_hd_global_2012.tif')

modis = overlay(modis_2001, modis_2002, modis_2003, modis_2004, modis_2005, modis_2006, 
                modis_2007, modis_2008, modis_2009, modis_2010, modis_2011, modis_2012, fun = mean)
modis[modis == 16] <- 0 #barren
modis[modis > 0] <- 1 # vegetated

## create global data raster
global_data_raster <- rasterFromXYZ(cbind(global_data$lon, global_data$lat, global_data[,3:7]))
global_data_veg_raster <- global_data_raster * modis
global_data_veg <- as.data.frame(rasterToPoints(global_data_veg_raster))
colnames(global_data_veg) <- c('lon', 'lat', 'par', 'tmp', 'vpd', 'f', 'z')

## read isotope values, calculate beta distributions
### read isotope data
iso_data <- read.csv('data/iso_data/iso_data.csv') # from cheaib et al. (2025; file = data_clean_beta.csv)
head(iso_data)

### separate into c3 and c4
iso_data_c3 <- subset(iso_data, PS_pathway == 'C3')
iso_data_c4 <- subset(iso_data, PS_pathway == 'C4')
nrow(iso_data_c3)
nrow(iso_data_c4)

### calculate beta for c3
iso_data_c3$gammastar_pa <- calc_gammastar_pa(temp = iso_data_c3$tmp, z = iso_data_c3$z)
iso_data_c3$km_pa <- calc_km_pa(temp = iso_data_c3$tmp, z = iso_data_c3$z)
iso_data_c3$nstar <- calc_nstar(temp = iso_data_c3$tmp, z = iso_data_c3$z)
iso_data_c3$vpd_kpa <- calc_vpd(temp = iso_data_c3$tmp, z = iso_data_c3$z, vpdo = iso_data_c3$vpd)
iso_data_c3$ca <- iso_data_c3$CO2 * 1e-6 * calc_patm(iso_data_c3$z)
iso_data_c3$a_frac <- 4.4
iso_data_c3$b_frac <- 28
iso_data_c3$f_frac <- 12
iso_data_c3$chi <- (iso_data_c3$big_D13 - (iso_data_c3$a_frac + (iso_data_c3$f_frac * (iso_data_c3$gammastar_pa/iso_data_c3$ca))))/
  (iso_data_c3$b_frac - iso_data_c3$a_frac)
# hist(iso_data_c3$chi)
iso_data_c3$beta <- 1.6 * iso_data_c3$nstar * iso_data_c3$vpd_kpa * 1000 *
  (((iso_data_c3$chi - (iso_data_c3$gammastar_pa/iso_data_c3$ca))^2)/
                                                 (((1- iso_data_c3$chi)^2) * (iso_data_c3$km_pa + iso_data_c3$gammastar_pa)))
hist(subset(iso_data_c3, chi < 0.95 & chi > 0.2)$beta)
hist(log(subset(iso_data_c3, chi < 0.95 & chi > 0.2)$beta))

### calculate beta for c4
iso_data_c4$gammastar_pa <- calc_gammastar_pa(temp = iso_data_c4$tmp, z = iso_data_c4$z)
iso_data_c4$km_pa <- calc_km_pa(temp = iso_data_c4$tmp, z = iso_data_c4$z)
iso_data_c4$kp_pa <- calc_kp_temp_pa(temp = iso_data_c4$tmp, z = iso_data_c4$z)
iso_data_c4$nstar <- calc_nstar(temp = iso_data_c4$tmp, z = iso_data_c4$z)
iso_data_c4$vpd_kpa <- calc_vpd(temp = iso_data_c4$tmp, z = iso_data_c4$z, vpdo = iso_data_c4$vpd)
iso_data_c4$ca <- iso_data_c4$CO2 * 1e-6 * calc_patm(iso_data_c4$z)
iso_data_c4$a_frac <- 4.4
iso_data_c4$b_frac <- -5.7+0.2*30
iso_data_c4$f_frac <- 12
iso_data_c4$chi <- (iso_data_c4$big_D13 - (iso_data_c4$a_frac + (iso_data_c4$f_frac * (iso_data_c4$gammastar_pa/iso_data_c4$ca))))/
  (iso_data_c4$b_frac - iso_data_c4$a_frac)
# hist(iso_data_c4$chi)
hist(subset(iso_data_c4, chi > 0)$chi)
iso_data_c4$beta <- 1.6 * iso_data_c4$nstar * iso_data_c4$vpd_kpa * 1000 *
  (((iso_data_c4$chi)^2)/
     (((1- iso_data_c4$chi)^2) * (iso_data_c4$kp_pa)))
hist(subset(iso_data_c4, chi < 0.95 & chi > 0.1)$beta)
hist(log(subset(iso_data_c4, chi < 0.95 & chi > 0.1)$beta))

### calculate mean and stdev beta for c3 and c4
beta_c3_mean <- mean(log(subset(iso_data_c3, chi < 0.95 & chi > 0.2)$beta))
beta_c3_sd <- sd(log(subset(iso_data_c3, chi < 0.95 & chi > 0.2)$beta))
beta_c4_mean <- mean(log(subset(iso_data_c4, chi < 0.95 & chi > 0.1)$beta))
beta_c4_sd <- sd(log(subset(iso_data_c4, chi < 0.95 & chi > 0.1)$beta))

## read in and clean neon data
### read data
neon_data <- read.csv('data/neon/neon_core_terrestrial_metadata.csv')

### calculate nearest latitude
latitude_values <- global_data_veg$lat
neon_closest_lat <- c()
for(i in 1:length(neon_data$latitude)){
  
  temp_lat <- neon_data$latitude[i]
  closest_lat_position <- which(abs(latitude_values - temp_lat) == min(abs(latitude_values - temp_lat)))[1]
  closest_lat <- latitude_values[closest_lat_position]
  neon_closest_lat <- c(neon_closest_lat, closest_lat)
  
}
neon_data$closest_latitude <- neon_closest_lat

### calculone nearest longitude
longitude_values <- global_data_veg$lon
neon_closest_lon <- c()
for(i in 1:length(neon_data$longitude)){
  
  temp_lon <- neon_data$longitude[i]
  closest_lon_position <- which(abs(longitude_values - temp_lon) == min(abs(longitude_values - temp_lon)))[1]
  closest_lon <- longitude_values[closest_lon_position]
  neon_closest_lon <- c(neon_closest_lon, closest_lon)
  
}
neon_data$closest_longitude <- neon_closest_lon

### add climate to neon dataset
neon_data_clim <- left_join(neon_data, global_data_veg, by = c('closest_latitude' = 'lat', 'closest_longitude' = 'lon'))


#############################
### Run various models ##
#############################

## run model for c3 deciduous plants
global_optimal_traits_c3_deciduous <- calc_optimal_vcmax(pathway = 'C3',
                                                         deciduous = 'yes',
                                                         tg_c = global_data_veg$tmp, 
                                                         vpdo = global_data_veg$vpd,
                                                         paro = global_data_veg$par,
                                                         z = global_data_veg$z,
                                                         f = global_data_veg$f)

## add lat/lon
global_optimal_traits_c3_deciduous$lat <- global_data_veg$lat
global_optimal_traits_c3_deciduous$lon <- global_data_veg$lon

## run model for c3 evergreen plants
global_optimal_traits_c3_evergreen <- calc_optimal_vcmax(pathway = 'C3',
                                                         deciduous = 'no',
                                                         tg_c = global_data_veg$tmp, 
                                                         vpdo = global_data_veg$vpd,
                                                         paro = global_data_veg$par,
                                                         z = global_data_veg$z,
                                                         f = global_data_veg$f)

## add lat/lon
global_optimal_traits_c3_evergreen$lat <- global_data_veg$lat
global_optimal_traits_c3_evergreen$lon <- global_data_veg$lon

## run model for c4 deciduous plants
global_optimal_traits_c4_deciduous <- calc_optimal_vcmax(pathway = 'C4',
                                                        deciduous = 'yes',
                                                        tg_c = global_data_veg$tmp, 
                                                        vpdo = global_data_veg$vpd,
                                                        paro = global_data_veg$par,
                                                        z = global_data_veg$z,
                                                        f = global_data_veg$f)

## add lat/lon
global_optimal_traits_c4_deciduous$lat <- global_data_veg$lat
global_optimal_traits_c4_deciduous$lon <- global_data_veg$lon

## run model for neon sites with varying beta values
cbind(neon_data_clim$site_id, neon_data_clim$dominant_nlcd_classes)
### representative sites (2 of each type)
# c4 sites = cper[3], konz[6]
# c3 deciduous sites = scbi[12], unde[17]
# c3 evergreen sites = puum[11], wref[19]
# c3 mixed = harv[5], tall[15]

### cper
global_optimal_traits_cper <- data.frame()
beta_cper <- exp(rnorm(1000, beta_c4_mean, beta_c4_sd))
tmp_cper <- rnorm(1000, neon_data_clim$tmp[3], neon_data_clim$tmp[3] * 0.1)
vpd_cper <- rnorm(1000, neon_data_clim$vpd[3], neon_data_clim$vpd[3] * 0.1)
par_cper <- rnorm(1000, neon_data_clim$par[3], neon_data_clim$par[3] * 0.1)
for(i in 1:length(beta_cper)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C4',
                                       deciduous = 'yes',
                                       tg_c = tmp_cper[i], 
                                       vpdo = vpd_cper[i],
                                       paro = par_cper[i],
                                       z = neon_data_clim$z[3],
                                       f = neon_data_clim$f[3],
                                       beta = beta_cper[i])
  
  global_optimal_traits_cper <- rbind(global_optimal_traits_cper, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_cper$lat <- neon_data_clim$closest_latitude[3]
global_optimal_traits_cper$lon <- neon_data_clim$closest_longitude[3]

### konz
global_optimal_traits_konz <- data.frame()
beta_konz <- exp(rnorm(1000, beta_c4_mean, beta_c4_sd))
tmp_konz <- rnorm(1000, neon_data_clim$tmp[6], neon_data_clim$tmp[6] * 0.1)
vpd_konz <- rnorm(1000, neon_data_clim$vpd[6], neon_data_clim$vpd[6] * 0.1)
par_konz <- rnorm(1000, neon_data_clim$par[6], neon_data_clim$par[6] * 0.1)
for(i in 1:length(beta_konz)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C4',
                                       deciduous = 'yes',
                                       tg_c = tmp_konz[i], 
                                       vpdo = vpd_konz[i],
                                       paro = par_konz[i],
                                       z = neon_data_clim$z[6],
                                       f = neon_data_clim$f[6],
                                       beta = beta_konz[i])
  
  global_optimal_traits_konz <- rbind(global_optimal_traits_konz, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_konz$lat <- neon_data_clim$closest_latitude[6]
global_optimal_traits_konz$lon <- neon_data_clim$closest_longitude[6]

### scbi
global_optimal_traits_scbi <- data.frame()
beta_scbi <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
tmp_scbi <- rnorm(1000, neon_data_clim$tmp[12], neon_data_clim$tmp[12] * 0.1)
vpd_scbi <- rnorm(1000, neon_data_clim$vpd[12], neon_data_clim$vpd[12] * 0.1)
par_scbi <- rnorm(1000, neon_data_clim$par[12], neon_data_clim$par[12] * 0.1)
for(i in 1:length(beta_scbi)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = tmp_scbi[i], 
                                       vpdo = vpd_scbi[i],
                                       paro = par_scbi[i],
                                       z = neon_data_clim$z[12],
                                       f = neon_data_clim$f[12],
                                       beta = beta_scbi[i])
  
  global_optimal_traits_scbi <- rbind(global_optimal_traits_scbi, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_scbi$lat <- neon_data_clim$closest_latitude[12]
global_optimal_traits_scbi$lon <- neon_data_clim$closest_longitude[12]

### unde
global_optimal_traits_unde <- data.frame()
beta_unde <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
tmp_unde <- rnorm(1000, neon_data_clim$tmp[17], neon_data_clim$tmp[17] * 0.1)
vpd_unde <- rnorm(1000, neon_data_clim$vpd[17], neon_data_clim$vpd[17] * 0.1)
par_unde <- rnorm(1000, neon_data_clim$par[17], neon_data_clim$par[17] * 0.1)
for(i in 1:length(beta_unde)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = tmp_unde[i], 
                                       vpdo = vpd_unde[i],
                                       paro = par_unde[i],
                                       z = neon_data_clim$z[17],
                                       f = neon_data_clim$f[17],
                                       beta = beta_unde[i])
  
  global_optimal_traits_unde <- rbind(global_optimal_traits_unde, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_unde$lat <- neon_data_clim$closest_latitude[17]
global_optimal_traits_unde$lon <- neon_data_clim$closest_longitude[17]

### puum
global_optimal_traits_puum <- data.frame()
beta_puum <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
tmp_puum <- rnorm(1000, neon_data_clim$tmp[11], neon_data_clim$tmp[11] * 0.1)
vpd_puum <- rnorm(1000, neon_data_clim$vpd[11], neon_data_clim$vpd[11] * 0.1)
par_puum <- rnorm(1000, neon_data_clim$par[11], neon_data_clim$par[11] * 0.1)
for(i in 1:length(beta_puum)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'no',
                                       tg_c = tmp_puum[i], 
                                       vpdo = vpd_puum[i],
                                       paro = par_puum[i],
                                       z = neon_data_clim$z[11],
                                       f = neon_data_clim$f[11],
                                       beta = beta_puum[i])
  
  global_optimal_traits_puum <- rbind(global_optimal_traits_puum, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_puum$lat <- neon_data_clim$closest_latitude[11]
global_optimal_traits_puum$lon <- neon_data_clim$closest_longitude[11]

### sjer
global_optimal_traits_sjer <- data.frame()
beta_sjer <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
tmp_sjer <- rnorm(1000, neon_data_clim$tmp[13], neon_data_clim$tmp[13] * 0.1)
vpd_sjer <- rnorm(1000, neon_data_clim$vpd[13], neon_data_clim$vpd[13] * 0.1)
par_sjer <- rnorm(1000, neon_data_clim$par[13], neon_data_clim$par[13] * 0.1)
for(i in 1:length(beta_sjer)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'no',
                                       tg_c = tmp_sjer[i], 
                                       vpdo = vpd_sjer[i],
                                       paro = par_sjer[i],
                                       z = neon_data_clim$z[13],
                                       f = neon_data_clim$f[13],
                                       beta = beta_sjer[i])
  
  global_optimal_traits_sjer <- rbind(global_optimal_traits_sjer, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_sjer$lat <- neon_data_clim$closest_latitude[13]
global_optimal_traits_sjer$lon <- neon_data_clim$closest_longitude[13]

### harv_deciduous
global_optimal_traits_harv_deciduous <- data.frame()
beta_harv_deciduous <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
tmp_harv_deciduous <- rnorm(1000, neon_data_clim$tmp[5], neon_data_clim$tmp[5] * 0.1)
vpd_harv_deciduous <- rnorm(1000, neon_data_clim$vpd[5], neon_data_clim$vpd[5] * 0.1)
par_harv_deciduous <- rnorm(1000, neon_data_clim$par[5], neon_data_clim$par[5] * 0.1)
for(i in 1:length(beta_harv_deciduous)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = tmp_harv_deciduous[i], 
                                       vpdo = vpd_harv_deciduous[i],
                                       paro = par_harv_deciduous[i],
                                       z = neon_data_clim$z[5],
                                       f = neon_data_clim$f[5],
                                       beta = beta_harv_deciduous[i])
  
  global_optimal_traits_harv_deciduous <- rbind(global_optimal_traits_harv_deciduous, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_harv_deciduous$lat <- neon_data_clim$closest_latitude[5]
global_optimal_traits_harv_deciduous$lon <- neon_data_clim$closest_longitude[5]

### harv_evergreen
global_optimal_traits_harv_evergreen <- data.frame()
beta_harv_evergreen <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
tmp_harv_evergreen <- rnorm(1000, neon_data_clim$tmp[5], neon_data_clim$tmp[5] * 0.1)
vpd_harv_evergreen <- rnorm(1000, neon_data_clim$vpd[5], neon_data_clim$vpd[5] * 0.1)
par_harv_evergreen <- rnorm(1000, neon_data_clim$par[5], neon_data_clim$par[5] * 0.1)
for(i in 1:length(beta_harv_evergreen)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'no',
                                       tg_c = tmp_harv_evergreen[i], 
                                       vpdo = vpd_harv_evergreen[i],
                                       paro = par_harv_evergreen[i],
                                       z = neon_data_clim$z[5],
                                       f = neon_data_clim$f[5],
                                       beta = beta_harv_evergreen[i])
  
  global_optimal_traits_harv_evergreen <- rbind(global_optimal_traits_harv_evergreen, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_harv_evergreen$lat <- neon_data_clim$closest_latitude[5]
global_optimal_traits_harv_evergreen$lon <- neon_data_clim$closest_longitude[5]

### tall_deciduous
global_optimal_traits_tall_deciduous <- data.frame()
beta_tall_deciduous <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
tmp_tall_deciduous <- rnorm(1000, neon_data_clim$tmp[15], neon_data_clim$tmp[15] * 0.1)
vpd_tall_deciduous <- rnorm(1000, neon_data_clim$vpd[15], neon_data_clim$vpd[15] * 0.1)
par_tall_deciduous <- rnorm(1000, neon_data_clim$par[15], neon_data_clim$par[15] * 0.1)
for(i in 1:length(beta_tall_deciduous)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = tmp_tall_deciduous[i], 
                                       vpdo = vpd_tall_deciduous[i],
                                       paro = par_tall_deciduous[i],
                                       z = neon_data_clim$z[15],
                                       f = neon_data_clim$f[15],
                                       beta = beta_tall_deciduous[i])
  
  global_optimal_traits_tall_deciduous <- rbind(global_optimal_traits_tall_deciduous, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_tall_deciduous$lat <- neon_data_clim$closest_latitude[15]
global_optimal_traits_tall_deciduous$lon <- neon_data_clim$closest_longitude[15]

### tall_evergreen
global_optimal_traits_tall_evergreen <- data.frame()
beta_tall_evergreen <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
tmp_tall_evergreen <- rnorm(1000, neon_data_clim$tmp[15], neon_data_clim$tmp[15] * 0.1)
vpd_tall_evergreen <- rnorm(1000, neon_data_clim$vpd[15], neon_data_clim$vpd[15] * 0.1)
par_tall_evergreen <- rnorm(1000, neon_data_clim$par[15], neon_data_clim$par[15] * 0.1)
for(i in 1:length(beta_tall_evergreen)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'no',
                                       tg_c = tmp_tall_evergreen[i], 
                                       vpdo = vpd_tall_evergreen[i],
                                       paro = par_tall_evergreen[i],
                                       z = neon_data_clim$z[15],
                                       f = neon_data_clim$f[15],
                                       beta = beta_tall_evergreen[i])
  
  global_optimal_traits_tall_evergreen <- rbind(global_optimal_traits_tall_evergreen, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_tall_evergreen$lat <- neon_data_clim$closest_latitude[15]
global_optimal_traits_tall_evergreen$lon <- neon_data_clim$closest_longitude[15]

#############################
### global sim PCA ##
#############################

##global pft plots
### combine model outputs with new category
global_optimal_traits_c3_deciduous$pft <- 'c3_deciduous'
global_optimal_traits_c3_evergreen$pft <- 'c3_evergreen'
global_optimal_traits_c4_deciduous$pft <- 'c4_deciduous'

global_optimal_traits_all <- rbind(global_optimal_traits_c3_deciduous, 
                                   global_optimal_traits_c3_evergreen, 
                                   global_optimal_traits_c4_deciduous)

global_optimal_traits_all_select <- as.data.frame(select(subset(global_optimal_traits_all, par > 0 & vpd > 0 & tg_c > 0), 
                                                         lma, Anet, wue, gsw, chi, nue, nphoto, narea, nmass, vcmax25, jmax25, rd25,
                                                         vpd, tg_c, par, pft))
global_optimal_traits_all_select_nona <- na.omit(global_optimal_traits_all_select)

### fit pca
global_optimal_traits_all_pca <- prcomp(global_optimal_traits_all_select_nona[,c(1:7,12)], scale = T, center = T)
summary(global_optimal_traits_all_pca)
global_optimal_traits_all_pca$rotation[,1:3]

### plot results
arrow_scale_all <- 5 # Scale factor for arrows and labels to extend from origin
label_scale_all <- 5.5

loadings_pca_all <- as.data.frame(global_optimal_traits_all_pca$rotation[, 1:3]) 
loadings_pca_all$trait <- rownames(loadings_pca_all) ## ad trait column
loadings_pca_all$label_x <- with(loadings_pca_all, PC1 * label_scale_all)
loadings_pca_all$label_y <- with(loadings_pca_all, PC2 * label_scale_all)

pca_scores_all <- as.data.frame(global_optimal_traits_all_pca$x) # get scores
pca_scores_all$pft <- global_optimal_traits_all_select_nona$pft

global_optimal_traits_all_pca_plot_PC1PC2 <- ggplot(pca_scores_all, 
                                                    aes(x = PC1, y = PC2, group = pft, color = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_color_manual(values = c("lightgreen", "darkgreen", "brown")) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, alpha = 0.3) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_all,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_all, yend = PC2 * arrow_scale_all, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_all,
            aes(x = label_x, y = label_y, label = trait, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  labs(x = "PC1 (51%)", y = "PC2 (28%)") +
  guides(fill = guide_colorbar(title = "Density level"))

global_optimal_traits_all_pca_plot_PC1PC3 <- ggplot(pca_scores_all, aes(x = PC1, y = PC3, group = pft, color = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_color_manual(values = c("lightgreen", "darkgreen", "brown")) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.3) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_all,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_all, yend = PC2 * arrow_scale_all, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_all,
            aes(x = label_x, y = label_y, label = trait, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  labs(x = "PC1 (51%)", y = "PC3 (16%)") +
  guides(fill = guide_colorbar(title = "Density level"))

jpeg('results/plots/global_optimal_traits_all_pca_plot_PC1PC2.jpeg', width = 10, height = 10, units = 'in', res = 600)
plot(global_optimal_traits_all_pca_plot_PC1PC2)
dev.off()

jpeg('results/plots/global_optimal_traits_all_pca_plot_PC1PC3.jpeg', width = 10, height = 10, units = 'in', res = 600)
plot(global_optimal_traits_all_pca_plot_PC1PC3)
dev.off()

#############################
### global sim trait histograms ##
#############################

global_optimal_traits_all_select_nona$c3c4[global_optimal_traits_all_select_nona$pft=='c3_deciduous'] <- 'c3'
global_optimal_traits_all_select_nona$c3c4[global_optimal_traits_all_select_nona$pft=='c3_evergreen'] <- 'c3'
global_optimal_traits_all_select_nona$c3c4[global_optimal_traits_all_select_nona$pft=='c4_deciduous'] <- 'c4'

global_optimal_traits_all_hist_Anet <- ggplot(data = global_optimal_traits_all_select_nona, 
                                              aes(x = Anet, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('green', 'brown'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression(italic('A')[net] * ' (µmol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '')

global_optimal_traits_all_hist_gsw <- ggplot(data = global_optimal_traits_all_select_nona, 
                                              aes(x = gsw, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('green', 'brown'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression(italic('g')[sw] * ' (mol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '')

global_optimal_traits_all_hist_wue <- ggplot(data = global_optimal_traits_all_select_nona, 
                                             aes(x = wue, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('green', 'brown'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression('iWUE' * ' (µmol mol'^'-1'*')'), y = 'Density', fill = '')

global_optimal_traits_all_hist_rd25 <- ggplot(data = global_optimal_traits_all_select_nona, 
                                              aes(x = rd25, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('green', 'brown'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression(italic('R')[d25] * ' (µmol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '')

global_optimal_traits_all_hist_chi <- ggplot(data = global_optimal_traits_all_select_nona, 
                                              aes(x = chi, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('green', 'brown'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression(italic('χ') * ' (mol mol'^'-1'*')'), y = 'Density', fill = '')

global_optimal_traits_all_hist_nphoto <- ggplot(data = global_optimal_traits_all_select_nona, 
                                             aes(x = nphoto, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('green', 'brown'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression(italic('N')[photo] * ' (g m'^'-2'*')'), y = 'Density', fill = '')

global_optimal_traits_all_hist_lma <- ggplot(data = global_optimal_traits_all_select_nona, 
                                                aes(x = lma, fill = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("lightgreen", "darkgreen", "brown"), 
                    labels = c(expression('C'[3] * 'deciduous'), expression('C'[3] * 'evergreen'), expression('C'[4] * 'deciduous'))) +
  labs(x = expression('LMA' * ' (g m'^'-2'*')'), y = 'Density', fill = '')

global_optimal_traits_all_hist_nue <- ggplot(data = global_optimal_traits_all_select_nona, 
                                             aes(x = nue, fill = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("lightgreen", "darkgreen", "brown"), 
                    labels = c(expression('C'[3] * 'deciduous'), expression('C'[3] * 'evergreen'), expression('C'[4] * 'deciduous'))) +
  labs(x = expression('NUE' * ' (µmol gN'^'-1' * ' s'^'-1' * ')'), y = 'Density', fill = '')

jpeg(filename = "results/plots/global_optimal_traits_hist_all.jpeg", 
     width = 26, height = 14, units = 'in', res = 600)
multiplot(global_optimal_traits_all_hist_Anet,
          global_optimal_traits_all_hist_lma,
          global_optimal_traits_all_hist_gsw,
          global_optimal_traits_all_hist_wue,
          global_optimal_traits_all_hist_nphoto,
          global_optimal_traits_all_hist_nue,
          global_optimal_traits_all_hist_rd25,
          cols = 4)
dev.off()


##site plots
### combine model outputs with new category
global_optimal_traits_harv$site <- 'harv'
global_optimal_traits_ornl$site <- 'ornl'
global_optimal_traits_tall$site <- 'tall'
global_optimal_traits_unde$site <- 'unde'
global_optimal_traits_scbi$site <- 'scbi'
global_optimal_traits_cper$site <- 'cper'
global_optimal_traits_konz$site <- 'konz'
global_optimal_traits_onaq$site <- 'onaq'
global_optimal_traits_puum$site <- 'puum'
global_optimal_traits_sjer$site <- 'sjer'

global_optimal_traits_sites <- rbind(global_optimal_traits_harv, 
                                   global_optimal_traits_ornl,
                                   global_optimal_traits_scbi,
                                   global_optimal_traits_tall,
                                   global_optimal_traits_unde,
                                   global_optimal_traits_cper,
                                   global_optimal_traits_konz,
                                   global_optimal_traits_onaq,
                                   global_optimal_traits_puum,
                                   global_optimal_traits_sjer)

global_optimal_traits_sites_select <- as.data.frame(select(subset(global_optimal_traits_sites, par > 0 & vpd > 0 & tg_c > 0), 
                                                         lma, Anet, wue, gsw, chi, nue, nphoto, narea, nmass, vcmax25, jmax25, rd25,
                                                         vpd, tg_c, par, site))
global_optimal_traits_sites_select_nona <- na.omit(global_optimal_traits_sites_select)

### fit pca
global_optimal_traits_sites_pca <- prcomp(global_optimal_traits_sites_select_nona[,c(1:4,6:7)], scale = T, center = T)
summary(global_optimal_traits_sites_pca)
global_optimal_traits_sites_pca$rotation[,1:3]

### plot results
arrow_scale_sites <- 5 # Scale factor for arrows and labels to extend from origin
label_scale_sites <- 5.5

loadings_pca_sites <- as.data.frame(global_optimal_traits_sites_pca$rotation[, 1:3]) 
loadings_pca_sites$trait <- rownames(loadings_pca_sites) ## ad trait column
loadings_pca_sites$label_x <- with(loadings_pca_sites, PC1 * label_scale_sites)
loadings_pca_sites$label_y <- with(loadings_pca_sites, PC2 * label_scale_sites)

pca_scores_sites <- as.data.frame(global_optimal_traits_sites_pca$x) # get scores
pca_scores_sites$site <- global_optimal_traits_sites_select_nona$site

global_optimal_traits_sites_pca_plot_PC1PC2 <- ggplot(pca_scores_sites, 
                                                    aes(x = PC1, y = PC2, group = site, color = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.5, size = 0.5) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, alpha = 0.1) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_sites,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_sites, yend = PC2 * arrow_scale_sites, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_sites,
            aes(x = label_x, y = label_y, label = trait, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  labs(x = "PC1", y = "PC2") +
  guides(fill = guide_colorbar(title = "Density level"))

global_optimal_traits_sites_pca_plot_PC1PC2_nolines <- ggplot(pca_scores_sites, aes(x = PC1, y = PC2, group = site, color = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.5, size = 0.5) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
  scale_fill_viridis_c(option = "turbo") +
  labs(x = "PC1", y = "PC2") +
  guides(fill = guide_colorbar(title = "Density level"))

global_optimal_traits_sites_pca_plot_PC2PC3 <- ggplot(pca_scores_sites, aes(x = PC2, y = PC3, group = site, color = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.5, size = 0.5) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_sites,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_sites, yend = PC2 * arrow_scale_sites, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_sites,
            aes(x = label_x, y = label_y, label = trait, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE) +
  labs(x = "PC2", y = "PC3") +
  guides(fill = guide_colorbar(title = "Density level"))

# jpeg('results/plots/global_optimal_traits_sites_pca_plot_PC1PC2.jpeg', width = 10, height = 10, units = 'in', res = 600)
# plot(global_optimal_traits_sites_pca_plot_PC1PC2)
# dev.off()

# jpeg('results/plots/global_optimal_traits_sites_pca_plot_PC1PC2_nolines.jpeg', width = 10, height = 10, units = 'in', res = 600)
# plot(global_optimal_traits_sites_pca_plot_PC1PC2_nolines)
# dev.off()

################
## future sims##
################

## run model for c3 deciduous plants under future environments
global_optimal_traits_c3_deciduous_fut <- calc_optimal_vcmax(pathway = 'C3',
                                                         deciduous = 'yes',
                                                         tg_c = global_data_veg$tmp+3, 
                                                         vpdo = global_data_veg$vpd,
                                                         paro = global_data_veg$par,
                                                         z = global_data_veg$z,
                                                         f = global_data_veg$f,
                                                         cao = 1000)

## add lat/lon
global_optimal_traits_c3_deciduous_fut$lat <- global_data_veg$lat
global_optimal_traits_c3_deciduous_fut$lon <- global_data_veg$lon

## pca
### select and scale traits
# global_optimal_traits_c3_deciduous_fut_scale <- scale(select(global_optimal_traits_c3_deciduous_fut, 
#                                                          lma, chi, gsw, vcmax25, jmax25, Al, nphoto, rd25, narea, nmass, tg_c, vpd, par))
# 
# ### fit pca
# global_optimal_traits_c3_deciduous_fut_pca <- princomp(na.omit(global_optimal_traits_c3_deciduous_fut_scale))
# summary(global_optimal_traits_c3_deciduous_fut_pca)
# global_optimal_traits_c3_deciduous_fut_pca$loadings[, 1:2]
# 
# ### plot results
# global_optimal_traits_c3_deciduous_fut_pca_lineplot <- fviz_pca_biplot(global_optimal_traits_c3_deciduous_fut_pca, 
#                                                                    col.var = "red",
#                                                                    alpha.ind = 0.01,
#                                                                    geom = c("point"))
# 
# # jpeg('results/plots/global_optimal_traits_c3_deciduous_fut_pca_lineplot.jpeg', width = 10, height = 10, units = 'in', res = 600)
# # plot(global_optimal_traits_c3_deciduous_fut_pca_lineplot)
# # dev.off()

## run model for c3 evergreen plants under future environments
global_optimal_traits_c3_evergreen_fut <- calc_optimal_vcmax(pathway = 'C3',
                                                         deciduous = 'no',
                                                         tg_c = global_data_veg$tmp+3, 
                                                         vpdo = global_data_veg$vpd,
                                                         paro = global_data_veg$par,
                                                         z = global_data_veg$z,
                                                         f = global_data_veg$f,
                                                         cao = 1000)

## add lat/lon
global_optimal_traits_c3_evergreen_fut$lat <- global_data_veg$lat
global_optimal_traits_c3_evergreen_fut$lon <- global_data_veg$lon

# ## pca
# ### select and scale traits
# global_optimal_traits_c3_evergreen_fut_scale <- scale(select(global_optimal_traits_c3_evergreen_fut, 
#                                                          lma, chi, gsw, vcmax25, jmax25, Al, nphoto, rd25, narea, nmass, tg_c, vpd, par))
# 
# ### fit pca
# global_optimal_traits_c3_evergreen_fut_pca <- princomp(na.omit(global_optimal_traits_c3_evergreen_fut_scale))
# summary(global_optimal_traits_c3_evergreen_fut_pca)
# global_optimal_traits_c3_evergreen_fut_pca$loadings[, 1:2]
# 
# ### plot results
# global_optimal_traits_c3_evergreen_fut_pca_lineplot <- fviz_pca_biplot(global_optimal_traits_c3_evergreen_fut_pca, 
#                                                                    col.var = "red",
#                                                                    alpha.ind = 0.01,
#                                                                    geom = c("point"))
# 
# # jpeg('results/plots/global_optimal_traits_c3_evergreen_fut_pca_lineplot.jpeg', width = 10, height = 10, units = 'in', res = 600)
# # plot(global_optimal_traits_c3_evergreen_fut_pca_lineplot)
# # dev.off()

## run model for c4 deciduous plants
global_optimal_traits_c4_deciduous_fut <- calc_optimal_vcmax(pathway = 'C4',
                                                         deciduous = 'yes',
                                                         tg_c = global_data_veg$tmp+3, 
                                                         vpdo = global_data_veg$vpd,
                                                         paro = global_data_veg$par,
                                                         z = global_data_veg$z,
                                                         f = global_data_veg$f,
                                                         cao = 1000)

## add lat/lon
global_optimal_traits_c4_deciduous_fut$lat <- global_data_veg$lat
global_optimal_traits_c4_deciduous_fut$lon <- global_data_veg$lon

# ## pca
# ### select and scale traits
# global_optimal_traits_c4_deciduous_fut_scale <- scale(select(global_optimal_traits_c4_deciduous_fut, 
#                                                          lma, chi, gsw, vpmax25, vcmax25, jmax25, Al, nphoto, rd25, narea, nmass, tg_c, vpd, par))
# 
# ### fit pca
# global_optimal_traits_c4_deciduous_fut_pca <- princomp(na.omit(global_optimal_traits_c4_deciduous_fut_scale))
# summary(global_optimal_traits_c4_deciduous_fut_pca)
# global_optimal_traits_c4_deciduous_fut_pca$loadings[, 1:2]
# 
# ### plot results
# global_optimal_traits_c4_deciduous_fut_pca_lineplot <- fviz_pca_biplot(global_optimal_traits_c4_deciduous_fut_pca, 
#                                                                    col.var = "red",
#                                                                    alpha.ind = 0.01,
#                                                                    geom = c("point"))
# 
# # jpeg('results/plots/global_optimal_traits_c4_deciduous_fut_pca_lineplot.jpeg', width = 10, height = 10, units = 'in', res = 600)
# # plot(global_optimal_traits_c4_deciduous_fut_pca_lineplot)
# # dev.off()


#############################
### run model all together ##
#############################

## combine model outputs with new category
global_optimal_traits_c3_deciduous_fut$pft <- 'c3_deciduous'
global_optimal_traits_c3_evergreen_fut$pft <- 'c3_evergreen'
global_optimal_traits_c4_deciduous_fut$pft <- 'c4_deciduous'

global_optimal_traits_all_fut <- rbind(global_optimal_traits_c3_deciduous_fut, global_optimal_traits_c3_evergreen_fut, global_optimal_traits_c4_deciduous_fut)

global_optimal_traits_all_fut_select <- as.data.frame(select(global_optimal_traits_all_fut, 
                                                                  lma, Anet, wue, gsw, chi, nue, nphoto, narea, nmass, vcmax25, jmax25, rd25,
                                                                  vpd, tg_c, par, pft))
global_optimal_traits_all_fut_select_nona <- na.omit(global_optimal_traits_all_fut_scale)

### fit pca
global_optimal_traits_all_fut_pca <- prcomp(global_optimal_traits_all_fut_select_nona[,c(1:4,6:7)], scale = T, center = T)
summary(global_optimal_traits_all_fut_pca)
global_optimal_traits_all_fut_pca$rotation[,1:3]

### plot results
arrow_scale_all_fut <- 5 # Scale factor for arrows and labels to extend from origin
label_scale_all_fut <- 5.5

loadings_pca_all_fut <- as.data.frame(global_optimal_traits_all_fut_pca$rotation[, 1:3]) 
loadings_pca_all_fut$trait <- rownames(loadings_pca_all_fut) ## ad trait column
loadings_pca_all_fut$label_x <- with(loadings_pca_all_fut, PC1 * label_scale_all_fut)
loadings_pca_all_fut$label_y <- with(loadings_pca_all_fut, PC2 * label_scale_all_fut)

pca_scores_all_fut <- as.data.frame(global_optimal_traits_all_fut_pca$x) # get scores
pca_scores_all_fut$pft <- global_optimal_traits_all_fut_scale_nona$pft

global_optimal_traits_all_fut_pca_plot_PC1PC2 <- ggplot(pca_scores_all_fut, aes(x = PC1, y = PC2, group = pft, color = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.5, size = 0.5) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_all_fut,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_all_fut, yend = PC2 * arrow_scale_all_fut, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_all_fut,
            aes(x = label_x, y = label_y, label = trait, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE) +
  labs(x = "PC1", y = "PC2") +
  guides(fill = guide_colorbar(title = "Density level"))

global_optimal_traits_all_fut_pca_plot_PC2PC3 <- ggplot(pca_scores_all_fut, aes(x = PC2, y = PC3, group = pft, color = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.5, size = 0.5) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_all_fut,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_all_fut, yend = PC2 * arrow_scale_all_fut, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_all_fut,
            aes(x = label_x, y = label_y, label = trait, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE) +
  labs(x = "PC2", y = "PC3") +
  guides(fill = guide_colorbar(title = "Density level"))

# jpeg('results/plots/global_optimal_traits_all_fut_pca_plot_PC1PC2.jpeg', width = 10, height = 10, units = 'in', res = 600)
# plot(global_optimal_traits_all_fut_pca_plot_PC1PC2)
# dev.off()


## create map plots for each trait for each plant type
### create color palette for maps
pale = colorRampPalette(c('white', rev(brewer.pal(10,'Spectral'))))
cols = pale(28)

### c3 deciduous - Anet
#### Anet
hist(global_optimal_traits_c3_deciduous$Anet)
# arg_c3_deciduous_Anet = list(at = seq(0, 40, 4), labels = seq(0, 40, 4))
c3_deciduous_Anet_raster <- rasterFromXYZ(cbind(global_optimal_traits_c3_deciduous$lon,
                                                 global_optimal_traits_c3_deciduous$lat,
                                                 global_optimal_traits_c3_deciduous$Anet))
plot(c3_deciduous_Anet_raster)
# plot(c3_deciduous_Anet_raster, col=cols, breaks = seq(0, 210, 10), cex.axis=1.5, yaxt = 'n', xaxt = 'n', 
#      lab.breaks = seq(0, 210, 10), ylim = c(-90, 90), 
#      legend.args=list(text=expression(italic('V'*"'")[cmax25]*' (µmol m'^'-2'*' s'^'-1'*')'), 
#                       line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg_c3_deciduous_Anet)
maps::map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
# axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
# axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)

#### nue
hist(global_optimal_traits_c3_deciduous$nue)
# arg_c3_deciduous_nue = list(at = seq(15, 40, 5), labels = seq(15, 40, 5))
c3_deciduous_nue_raster <- rasterFromXYZ(cbind(global_optimal_traits_c3_deciduous$lon,
                                                 global_optimal_traits_c3_deciduous$lat,
                                                 global_optimal_traits_c3_deciduous$nue))
plot(c3_deciduous_nue_raster)
# plot(c3_deciduous_nue_raster, col=cols, breaks = seq(0, 200, 10), cex.axis=1.5, yaxt = 'n', xaxt = 'n', 
#      lab.breaks = seq(0, 200, 10), ylim = c(-90, 90), 
#      legend.args=list(text=expression(italic('M'*"'")[area]), 
#                       line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg_c3_deciduous_nue)
maps::map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
# axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
# axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)

#### gsw
hist(global_optimal_traits_c3_deciduous$gsw)
# arg_c3_deciduous_gsw = list(at = seq(0, 0.02, 0.001), labels = seq(0, 0.02, 0.001))
c3_deciduous_gsw_raster <- rasterFromXYZ(cbind(global_optimal_traits_c3_deciduous$lon,
                                               global_optimal_traits_c3_deciduous$lat,
                                               global_optimal_traits_c3_deciduous$gsw))
plot(c3_deciduous_gsw_raster)
# plot(c3_deciduous_gsw_raster, col=cols, breaks = seq(0, 0.02, 0.001), cex.axis=1.5, yaxt = 'n', xaxt = 'n', 
#      lab.breaks = seq(0, 0.02, 0.001), ylim = c(-90, 90), 
#      legend.args=list(text=expression(italic('M'*"'")[area]), 
#                       line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg_c3_deciduous_gsw)
maps::map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
# axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
# axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)













## create color palette for maps
pale = colorRampPalette(c('white', rev(brewer.pal(10,'Spectral'))))
cols = pale(28)
arg_vcmax = list(at = seq(0, 135, 15), labels = seq(0, 135, 15))
arg_lma = list(at = seq(0, 20000, 2000), labels = seq(0, 20000, 2000))

## create and plot rasters for different traits
### vcmax
vcmax_raster <- rasterFromXYZ(cbind(global_optimal_traits$lon, global_optimal_traits$lat, global_optimal_traits$vcmax))
vcmax_raster_mod <- vcmax_raster * modis
plot(vcmax_raster_mod, col=cols, breaks = seq(0, 140, 5), cex.axis=1.5, yaxt = 'n', xaxt = 'n', 
     lab.breaks = seq(0, 140, 5), ylim = c(-90, 90), 
     legend.args=list(text=expression(italic('V'*"'")[cmax]*' (µmol m'^'-2'*' s'^'-1'*')'), 
                      line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg_vcmax)
map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)

### lma
lma_raster <- rasterFromXYZ(cbind(global_optimal_traits$lon, global_optimal_traits$lat, global_optimal_traits$lma))
lma_raster_mod <- lma_raster * modis
plot(lma_raster_mod, col=cols, breaks = seq(0, 20000, 2000), cex.axis=1.5, yaxt = 'n', xaxt = 'n', 
     lab.breaks = seq(0, 140, 5), ylim = c(-90, 90), 
     legend.args=list(text=expression(italic('V'*"'")[cmax]*' (µmol m'^'-2'*' s'^'-1'*')'), 
                      line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg_lma)
map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)




