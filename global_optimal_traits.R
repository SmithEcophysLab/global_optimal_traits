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

# ## pca
# ### select and scale traits
# global_optimal_traits_c3_deciduous_scale <- scale(select(global_optimal_traits_c3_deciduous, 
#                                                          lma, Anet, wue, nue, nphoto, narea, tg_c, vpd, par))
# 
# ### fit pca
# global_optimal_traits_c3_deciduous_pca <- prcomp(na.omit(global_optimal_traits_c3_deciduous_scale))
# summary(global_optimal_traits_c3_deciduous_pca)
# global_optimal_traits_c3_deciduous_pca$rotation[,1:3]
# 
# ### plot results
# arrow_scale_c3_deciduous <- 5 # Scale factor for arrows and labels to extend from origin
# label_scale_c3_deciduous <- 5.5
# loadings_pca_c3_deciduous <- as.data.frame(global_optimal_traits_c3_deciduous_pca$rotation[, 1:2]) 
# loadings_pca_c3_deciduous$trait <- rownames(loadings_pca_c3_deciduous) ## ad trait column
# pca_scores_c3_deciduous <- as.data.frame(global_optimal_traits_c3_deciduous_pca$x) # get scores
# loadings_pca_c3_deciduous$label_x <- with(loadings_pca_c3_deciduous, PC1 * label_scale_c3_deciduous)
# loadings_pca_c3_deciduous$label_y <- with(loadings_pca_c3_deciduous, PC2 * label_scale_c3_deciduous)
# 
# global_optimal_traits_c3_deciduous_pca_plot <- ggplot(pca_scores_c3_deciduous, aes(x = PC1, y = PC2)) +
#   theme_minimal(base_size = 14) +
#   theme(axis.title = element_text(size = 18, face = "bold"),
#     axis.text = element_text(size = 14),
#     axis.line = element_line(color = "black", linewidth = 0.6),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()) +
#   coord_equal() +
#   geom_point(alpha = 0.03, size = 0.5) +
#   stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
#   scale_fill_viridis_c(option = "turbo") +
#   geom_segment(data = loadings_pca_c3_deciduous,
#                aes(x = 0, y = 0, xend = PC1 * arrow_scale_c3_deciduous, yend = PC2 * arrow_scale_c3_deciduous),
#                arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
#   geom_text(data = loadings_pca_c3_deciduous,
#             aes(x = label_x, y = label_y, label = trait),
#             size = 4, fontface = "bold", parse = TRUE) +
#   labs(x = "PC1", y = "PC2") +
#   guides(fill = guide_colorbar(title = "Density level"))
# 
# # jpeg('results/plots/global_optimal_traits_c3_deciduous_pca_plot.jpeg', width = 10, height = 10, units = 'in', res = 600)
# # plot(global_optimal_traits_c3_deciduous_pca_plot)
# # dev.off()

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

# ## pca
# ### select and scale traits
# global_optimal_traits_c3_evergreen_scale <- scale(select(global_optimal_traits_c3_evergreen, 
#                                                          lma, Anet, wue, nue, nphoto, narea, tg_c, vpd, par))
# 
# ### fit pca
# global_optimal_traits_c3_evergreen_pca <- prcomp(na.omit(global_optimal_traits_c3_evergreen_scale))
# summary(global_optimal_traits_c3_evergreen_pca)
# global_optimal_traits_c3_evergreen_pca$rotation[,1:3]
# 
# ### plot results
# arrow_scale_c3_evergreen <- 5 # Scale factor for arrows and labels to extend from origin
# label_scale_c3_evergreen <- 5.5
# loadings_pca_c3_evergreen <- as.data.frame(global_optimal_traits_c3_evergreen_pca$rotation[, 1:2]) 
# loadings_pca_c3_evergreen$trait <- rownames(loadings_pca_c3_evergreen) ## ad trait column
# pca_scores_c3_evergreen <- as.data.frame(global_optimal_traits_c3_evergreen_pca$x) # get scores
# loadings_pca_c3_evergreen$label_x <- with(loadings_pca_c3_evergreen, PC1 * label_scale_c3_evergreen)
# loadings_pca_c3_evergreen$label_y <- with(loadings_pca_c3_evergreen, PC2 * label_scale_c3_evergreen)
# 
# global_optimal_traits_c3_evergreen_pca_plot <- ggplot(pca_scores_c3_evergreen, aes(x = PC1, y = PC2)) +
#   theme_minimal(base_size = 14) +
#   theme(axis.title = element_text(size = 18, face = "bold"),
#         axis.text = element_text(size = 14),
#         axis.line = element_line(color = "black", linewidth = 0.6),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() +
#   geom_point(alpha = 0.03, size = 0.5) +
#   stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
#   scale_fill_viridis_c(option = "turbo") +
#   geom_segment(data = loadings_pca_c3_evergreen,
#                aes(x = 0, y = 0, xend = PC1 * arrow_scale_c3_evergreen, yend = PC2 * arrow_scale_c3_evergreen),
#                arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
#   geom_text(data = loadings_pca_c3_evergreen,
#             aes(x = label_x, y = label_y, label = trait),
#             size = 4, fontface = "bold", parse = TRUE) +
#   labs(x = "PC1", y = "PC2") +
#   guides(fill = guide_colorbar(title = "Density level"))
# 
# # jpeg('results/plots/global_optimal_traits_c3_evergreen_pca_plot.jpeg', width = 10, height = 10, units = 'in', res = 600)
# # plot(global_optimal_traits_c3_evergreen_pca_plot)
# # dev.off()
# 
# ## run model for c4 deciduous plants
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
# global_optimal_traits_c4_deciduous_scale <- scale(select(global_optimal_traits_c4_deciduous, 
#                                                         lma, chi, gsw, vpmax25, vcmax25, jmax25, Anet, wue, nue, nphoto, rd25, narea, nmass, tg_c, vpd, par))
# 
# ### fit pca
# global_optimal_traits_c4_deciduous_pca <- prcomp(na.omit(global_optimal_traits_c4_deciduous_scale))
# summary(global_optimal_traits_c4_deciduous_pca)
# global_optimal_traits_c4_deciduous_pca$rotation[,1:3]
# 
# ### plot results
# arrow_scale_c4_deciduous <- 5 # Scale factor for arrows and labels to extend from origin
# label_scale_c4_deciduous <- 5.5
# loadings_pca_c4_deciduous <- as.data.frame(global_optimal_traits_c4_deciduous_pca$rotation[, 1:2]) 
# loadings_pca_c4_deciduous$trait <- rownames(loadings_pca_c4_deciduous) ## ad trait column
# pca_scores_c4_deciduous <- as.data.frame(global_optimal_traits_c4_deciduous_pca$x) # get scores
# loadings_pca_c4_deciduous$label_x <- with(loadings_pca_c4_deciduous, PC1 * label_scale_c4_deciduous)
# loadings_pca_c4_deciduous$label_y <- with(loadings_pca_c4_deciduous, PC2 * label_scale_c4_deciduous)
# 
# global_optimal_traits_c4_deciduous_pca_plot <- ggplot(pca_scores_c4_deciduous, aes(x = PC1, y = PC2)) +
#   theme_minimal(base_size = 14) +
#   theme(axis.title = element_text(size = 18, face = "bold"),
#         axis.text = element_text(size = 14),
#         axis.line = element_line(color = "black", linewidth = 0.6),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_equal() +
#   geom_point(alpha = 0.03, size = 0.5) +
#   stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
#   scale_fill_viridis_c(option = "turbo") +
#   geom_segment(data = loadings_pca_c4_deciduous,
#                aes(x = 0, y = 0, xend = PC1 * arrow_scale_c4_deciduous, yend = PC2 * arrow_scale_c4_deciduous),
#                arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
#   geom_text(data = loadings_pca_c4_deciduous,
#             aes(x = label_x, y = label_y, label = trait),
#             size = 4, fontface = "bold", parse = TRUE) +
#   labs(x = "PC1", y = "PC2") +
#   guides(fill = guide_colorbar(title = "Density level"))
# 
# # jpeg('results/plots/global_optimal_traits_c4_deciduous_pca_plot.jpeg', width = 10, height = 10, units = 'in', res = 600)
# # plot(global_optimal_traits_c4_deciduous_pca_plot)
# # dev.off()

## run model for neon sites with varying beta values
cbind(neon_data_clim$site_id, neon_data_clim$dominant_nlcd_classes)
# c4 sites = CPER, KONZ

### bona
global_optimal_traits_bona <- data.frame()
beta_bona <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
for(i in 1:length(beta_bona)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = neon_data_clim$tmp[1], 
                                       vpdo = neon_data_clim$vpd[1],
                                       paro = neon_data_clim$par[1],
                                       z = neon_data_clim$z[1],
                                       f = neon_data_clim$f[1],
                                       beta = beta_bona[i])
  
  global_optimal_traits_bona <- rbind(global_optimal_traits_bona, optimal_traits)
  
}

### clbj
global_optimal_traits_clbj <- data.frame()
beta_clbj <- exp(rnorm(1000, beta_c3_mean, beta_c3_sd))
for(i in 1:length(beta_clbj)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = neon_data_clim$tmp[2], 
                                       vpdo = neon_data_clim$vpd[2],
                                       paro = neon_data_clim$par[2],
                                       z = neon_data_clim$z[2],
                                       f = neon_data_clim$f[2],
                                       beta = beta_clbj[i])
  
  global_optimal_traits_clbj <- rbind(global_optimal_traits_clbj, optimal_traits)
  
}

#######
## need to finish out rest of sites
#######


#############################
### PCA all together ##
#############################

## combine model outputs with new category
global_optimal_traits_c3_deciduous$pft <- 'c3_deciduous'
global_optimal_traits_c3_evergreen$pft <- 'c3_evergreen'
global_optimal_traits_c4_deciduous$pft <- 'c4_deciduous'

global_optimal_traits_all <- rbind(global_optimal_traits_c3_deciduous, global_optimal_traits_c3_evergreen, global_optimal_traits_c4_deciduous)

global_optimal_traits_all_select <- as.data.frame(select(subset(global_optimal_traits_all, par > 0 & vpd > 0 & tg_c > 0), 
                                                         lma, Anet, wue, gsw, chi, nue, nphoto, narea, nmass, vcmax25, jmax25, rd25,
                                                         vpd, tg_c, par, pft))
global_optimal_traits_all_select_nona <- na.omit(global_optimal_traits_all_select)

### fit pca
global_optimal_traits_all_pca <- prcomp(global_optimal_traits_all_select_nona[,c(1:4,6:7)], scale = T, center = T)
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

global_optimal_traits_all_pca_plot_PC1PC2 <- ggplot(pca_scores_all, aes(x = PC1, y = PC2, group = pft, color = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.5, size = 0.5) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_all,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_all, yend = PC2 * arrow_scale_all, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_all,
            aes(x = label_x, y = label_y, label = trait, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE) +
  labs(x = "PC1", y = "PC2") +
  guides(fill = guide_colorbar(title = "Density level"))

global_optimal_traits_all_pca_plot_PC1PC2_nolines <- ggplot(pca_scores_all, aes(x = PC1, y = PC2, group = pft, color = pft)) +
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

global_optimal_traits_all_pca_plot_PC2PC3 <- ggplot(pca_scores_all, aes(x = PC2, y = PC3, group = pft, color = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.5, size = 0.5) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.5) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_all,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_all, yend = PC2 * arrow_scale_all, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_all,
            aes(x = label_x, y = label_y, label = trait, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE) +
  labs(x = "PC2", y = "PC3") +
  guides(fill = guide_colorbar(title = "Density level"))

# jpeg('results/plots/global_optimal_traits_all_pca_plot_PC1PC2.jpeg', width = 10, height = 10, units = 'in', res = 600)
# plot(global_optimal_traits_all_pca_plot_PC1PC2)
# dev.off()

# jpeg('results/plots/global_optimal_traits_all_pca_plot_PC1PC2_nolines.jpeg', width = 10, height = 10, units = 'in', res = 600)
# plot(global_optimal_traits_all_pca_plot_PC1PC2_nolines)
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




