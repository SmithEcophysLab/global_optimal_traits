# global_optimal_traits.R
## script to predict global optimal leaf traits
## still a bit of work to do!

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
global_data <- left_join(cru_growingseason_data, global_z_data)

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

## pca
### select and scale traits
global_optimal_traits_c3_deciduous_scale <- scale(select(global_optimal_traits_c3_deciduous, 
                                                         lma, chi, gsw, vcmax25, jmax25, Al, nphoto, rd25, narea, nmass, tg_c, vpd, par))

### fit pca
global_optimal_traits_c3_deciduous_pca <- princomp(na.omit(global_optimal_traits_c3_deciduous_scale))
summary(global_optimal_traits_c3_deciduous_pca)
global_optimal_traits_c3_deciduous_pca$loadings[, 1:2]

### plot results
global_optimal_traits_c3_deciduous_pca_lineplot <- fviz_pca_biplot(global_optimal_traits_c3_deciduous_pca, 
                                                                   col.var = "red",
                                                                   alpha.ind = 0.01,
                                                                   geom = c("point"))

# jpeg('results/plots/global_optimal_traits_c3_deciduous_pca_lineplot.jpeg')
# plot(global_optimal_traits_c3_deciduous_pca_lineplot)
# dev.off()

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

## pca
### select and scale traits
global_optimal_traits_c3_evergreen_scale <- scale(select(global_optimal_traits_c3_evergreen, 
                                                         lma, chi, vcmax25, jmax25, Al, nphoto, rd25, narea, nmass, tg_c, vpd, par))

### fit pca
global_optimal_traits_c3_evergreen_pca <- princomp(na.omit(global_optimal_traits_c3_evergreen_scale))
summary(global_optimal_traits_c3_evergreen_pca)
global_optimal_traits_c3_evergreen_pca$loadings[, 1:2]

### plot results
global_optimal_traits_c3_evergreen_pca_lineplot <- fviz_pca_biplot(global_optimal_traits_c3_evergreen_pca, 
                                                                   col.var = "red",
                                                                   alpha.ind = 0.01,
                                                                   geom = c("point"))

# jpeg('results/plots/global_optimal_traits_c3_evergreen_pca_lineplot.jpeg')
# plot(global_optimal_traits_c3_evergreen_pca_lineplot)
# dev.off()

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

## pca
### select and scale traits
global_optimal_traits_c4_deciduous_scale <- scale(select(global_optimal_traits_c4_deciduous, 
                                                         lma, chi, vpmax25, vcmax25, jmax25, Al, nphoto, rd25, narea, nmass, tg_c, vpd, par))

### fit pca
global_optimal_traits_c4_deciduous_pca <- princomp(na.omit(global_optimal_traits_c4_deciduous_scale))
summary(global_optimal_traits_c4_deciduous_pca)
global_optimal_traits_c4_deciduous_pca$loadings[, 1:2]

### plot results
global_optimal_traits_c4_deciduous_pca_lineplot <- fviz_pca_biplot(global_optimal_traits_c4_deciduous_pca, 
                                                                   col.var = "red",
                                                                   alpha.ind = 0.01,
                                                                   geom = c("point"))

# jpeg('results/plots/global_optimal_traits_c4_deciduous_pca_lineplot.jpeg')
# plot(global_optimal_traits_c4_deciduous_pca_lineplot)
# dev.off()


## create map plots for each trait for each plant type
### create color palette for maps
pale = colorRampPalette(c('white', rev(brewer.pal(10,'Spectral'))))
cols = pale(28)

### c3 deciduous - vcmax25, nmass, chi, lma
#### vcmax25
hist(global_optimal_traits_c3_deciduous$vcmax25)
arg_c3_deciduous_vcmax25 = list(at = seq(0, 210, 10), labels = seq(0, 210, 10))
c3_deciduous_vcmax25_raster <- rasterFromXYZ(cbind(global_optimal_traits_c3_deciduous$lon,
                                                 global_optimal_traits_c3_deciduous$lat,
                                                 global_optimal_traits_c3_deciduous$vcmax2525))
plot(c3_deciduous_vcmax25_raster, col=cols, breaks = seq(0, 210, 10), cex.axis=1.5, yaxt = 'n', xaxt = 'n', 
     lab.breaks = seq(0, 210, 10), ylim = c(-90, 90), 
     legend.args=list(text=expression(italic('V'*"'")[cmax25]*' (µmol m'^'-2'*' s'^'-1'*')'), 
                      line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg_c3_deciduous_vcmax25)
maps::map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)

#### lma
hist(global_optimal_traits_c3_deciduous$lma)
arg_c3_deciduous_lma = list(at = seq(0, 200, 10), labels = seq(0, 200, 10))
c3_deciduous_lma_raster <- rasterFromXYZ(cbind(global_optimal_traits_c3_deciduous$lon,
                                                 global_optimal_traits_c3_deciduous$lat,
                                                 global_optimal_traits_c3_deciduous$lma))
plot(c3_deciduous_lma_raster, col=cols, breaks = seq(0, 200, 10), cex.axis=1.5, yaxt = 'n', xaxt = 'n', 
     lab.breaks = seq(0, 200, 10), ylim = c(-90, 90), 
     legend.args=list(text=expression(italic('M'*"'")[area]), 
                      line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg_c3_deciduous_lma)
maps::map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)

#### nmass
hist(global_optimal_traits_c3_deciduous$nmass)
arg_c3_deciduous_nmass = list(at = seq(0, 0.02, 0.001), labels = seq(0, 0.02, 0.001))
c3_deciduous_nmass_raster <- rasterFromXYZ(cbind(global_optimal_traits_c3_deciduous$lon,
                                               global_optimal_traits_c3_deciduous$lat,
                                               global_optimal_traits_c3_deciduous$nmass))
plot(c3_deciduous_nmass_raster, col=cols, breaks = seq(0, 0.02, 0.001), cex.axis=1.5, yaxt = 'n', xaxt = 'n', 
     lab.breaks = seq(0, 0.02, 0.001), ylim = c(-90, 90), 
     legend.args=list(text=expression(italic('M'*"'")[area]), 
                      line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg_c3_deciduous_nmass)
maps::map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)













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




