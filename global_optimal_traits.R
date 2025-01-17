# global_optimal_traits.R
## script to predict global optimal leaf traits

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
global_tmp_data <- read.csv('data/cru_growingseason/cru_tmp_climExtract_growingseason_globe.csv')[,-1]
global_vpd_data <- read.csv('data/cru_growingseason/cru_vpd_climExtract_growingseason_globe.csv')[,-1]
global_par_data <- read.csv('data/cru_growingseason/cru_par_climExtract_growingseason_globe.csv')[,-1]
global_z_data <- read.csv('data/watch_elevation/z_globe.csv')[,-1]

## read in phenology data

## combine datasets by lat/lon
global_data <- left_join(global_tmp_data, global_vpd_data)
global_data <- left_join(global_data, global_par_data)
global_data <- left_join(global_data, global_z_data)

## run model
global_optimal_traits <- calc_optimal_vcmax(tg_c = global_data$tmp, 
                                            vpdo = global_data$vpd,
                                            paro = global_data$par,
                                            z = global_data$z)

## add lat/lon
global_optimal_traits$lat <- global_data$lat
global_optimal_traits$lon <- global_data$lon

## pca
### select and scale traits
global_optimal_traits_scale <- scale(select(global_optimal_traits, lma, chi, vcmax, jmax, Al, nphoto))

### fit pca
global_optimal_traits_pca <- princomp(na.omit(global_optimal_traits_scale))
summary(global_optimal_traits_pca)
global_optimal_traits_pca$loadings[, 1:2]

### plot results
global_optimal_traits_pca_lineplot <- fviz_pca_var(global_optimal_traits_pca, col.var = "black")

jpeg('results/plots/global_optimal_traits_pca_lineplot.jpeg')
plot(global_optimal_traits_pca_lineplot)
dev.off()


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




