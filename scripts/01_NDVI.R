#### ELISE GALLOIS  - elise.gallois94@gmail.com 

#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(esquisse)
library(sf)
library(raster)
library(viridis)
library(readxl)
library(readr)
library(maps)
library(rgeos)
library(terra)
library(lsatTS)
library(purrr)
library(data.table)
library(stringr)
library(rgee)
library(strex)

# connect R to the GEE
ee_Initialize()
db <- 'CGIAR/SRTM90_V4'
image <- ee$Image(db)
image$bandNames()$getInfo()

#### 2 - MODIS: LOAD DATA ####
ndvi <- raster('data/MODIS_NDVI_slopes_Sval.tiff')
plot(ndvi/1000) # visualise 2000-2021 peak season ndvi change slopes
hist(ndvi/1000) # plot histogram of distribution

#### 3 - get coords of all plots included in polar alien hunters work
invasive_sites <- read_excel("data/invasive_plant_spp_svalbard_2016-7.xlsx",
                             sheet = "Spp per site")
str(invasive_sites)

(meanplot <- plot(ndvi/1000, main = "NDVI Change 2000-2021", col = viridis(100)))
#### 4 - extract NDVI change for each site
v <- vect(cbind(invasive_sites$X, invasive_sites$Y), crs="+proj=utm +zone=33X +datum=WGS84  +units=m")
v
y <- project(v, "+proj=longlat +datum=WGS84")
y

cord.dec = SpatialPoints(cbind(invasive_sites$X, invasive_sites$Y), proj4string=CRS("+proj=utm +zone=33X +datum=WGS84  +units=m"))
cord.UTM <- spTransform(cord.dec, crs(ndvi))

XY_sp <- SpatialPoints(cbind(invasive_sites$X, invasive_sites$Y), 
                       proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
cord.UTM <- spTransform(cord.dec , crs(ndvi))

crs(ndvi)
crs(cord.UTM)

## Now we just need to extract the raster values from these locations                
ndvi_point <- raster::extract(ndvi, cord.dec, method='simple',df=TRUE)
plot(ndvi)

#### 5 - LANDSAT: LOAD DATA ####
landsat <- read_csv("data/NDVI_yards.csv")
str(landsat)

# get DOY values
landsat$doy <- as.numeric(format(landsat$date_edit, "%j"))
# get year values
landsat$year <- as.numeric(format(landsat$date_edit, "%Y"))
# get month values
landsat$month <- as.numeric(format(landsat$date_edit, "%m"))

#remove all NDVI values below o.2 to exclude negatives and snow
landsat <- subset(landsat, mean > 0) 
hist(landsat$mean)
# visualise greening curves for all years
(year_plot <- landsat %>%
    filter(month >= 5L & month <= 9L) %>% # filter july and aug
    ggplot() +
    aes(x = doy, y = mean) +
    geom_smooth() +
    geom_point(aes(color=mean))+
    scale_color_distiller(palette = "YlGn", direction=-1) +
    labs(y = "NDVI", x = "DOY (2014 - 2021)", fill = "year") +
    facet_wrap(vars(year)) +
    theme_classic() +
    theme(legend.position = "none"))

# visualise greening curves for all plots
(site_plot <- landsat %>%
    filter(month >= 5L & month <= 9L) %>% # filter july and aug
    ggplot() +
    aes(x = doy, y = mean) +
    geom_point(aes(color=mean))+
    geom_smooth(aes(color=mean))+
    scale_color_distiller(palette = "YlGn", direction=-1) +
    labs(y = "NDVI", x = "DOY (2014 - 2021)", fill = "name") +
    facet_wrap(vars(name)) +
    theme_classic() +
    theme(legend.position = "none"))

# see change in peak season NDVI over time
peak <- landsat %>%  
  filter(month >= 7L & month <= 8L) %>% # filter july and aug
  group_by(name) %>% # group by site
  mutate(mean_summer = mean(mean)) %>% # calculate average summer NDVI
  mutate(max_summer = max(mean)) %>% # calculate max summer NDVI
  ungroup() # no longer grouped by year




#### 6 - SENTINEL: LOAD DATA ####
sentinel <- read_csv("data/NDVI_sentinel.csv")
str(sentinel)

# get DOY values
sentinel$doy <- as.numeric(format(sentinel$date_edit, "%j"))
# get year values
sentinel$year <- as.numeric(format(sentinel$date_edit, "%Y"))
# get month values
sentinel$month <- as.numeric(format(sentinel$date_edit, "%m"))

#remove all NDVI values below o.2 to exclude negatives and snow
sentinel <- subset(sentinel, mean > 0) 
hist(sentinel$mean)
# visualise greening curves for all years
(year_plot_sen <- sentinel %>%
    filter(month >= 5L & month <= 10L) %>% # filter july and aug
    ggplot() +
    aes(x = doy, y = mean) +
    geom_smooth() +
    geom_point(aes(color=mean))+
    scale_color_distiller(palette = "YlGn", direction=-1) +
    labs(y = "NDVI", x = "DOY (2014 - 2021)", fill = "year") +
    facet_wrap(vars(year)) +
    theme_classic() +
    theme(legend.position = "none"))

# visualise greening curves for all plots
(site_plot_sen <- sentinel %>%
    #filter(month >= 5L & month <= 10L) %>% # filter july and aug
    ggplot() +
    aes(x = doy, y = mean) +
    geom_point(aes(color=mean))+
    geom_smooth(aes(color=mean))+
    scale_color_distiller(palette = "YlGn", direction=-1) +
    labs(y = "NDVI", x = "DOY (2014 - 2021)", fill = "name") +
    facet_wrap(vars(name)) +
    theme_classic() +
    theme(legend.position = "none"))


#### 7 - Logan's Landsat Tool ####
#install.packages("/Users/elisegallois/Desktop/lsatTS", repos=NULL, type="source")
library('lsatTS')
library(geojsonio)

#### 7a - LSAT Trial Polygons ####
# trial with greendog yards
# Specify a region 
green_dog_a <- st_polygon(
  list(matrix(c(15.94988, 78.17300,
                15.94988, 78.17384,
                15.95391, 78.17384,
                15.95391, 78.17300,
                15.94988, 78.17300),
              ncol = 2, byrow = TRUE)))
greendog_a_sf <- st_sfc(green_dog_a, crs = 4326) %>% st_sf()



# Use lsat_get_pixel_centers to retrieve pixel centers and plot to a file that can be added to this documentation.
# We set plot_map to a file path (or just TRUE) to view 
greendog_a_poly <- lsat_get_pixel_centers(greendog_a_sf, plot_map = "figures/lsat_get_pixel_centers.png")

green_dog_b <- st_polygon(
  list(matrix(c(15.94731, 78.17141,
                15.94731, 78.17266,
                15.9784, 78.17266,
                15.9784, 78.17141,
                15.94731, 78.17141),
              ncol = 2, byrow = TRUE)))
greendog_b_sf <- st_sfc(green_dog_b, crs = 4326) %>% st_sf()

# Use lsat_get_pixel_centers to retrieve pixel centers and plot to a file that can be added to this documentation.
# We set plot_map to a file path (or just TRUE) to view 
greendog_b_poly <- lsat_get_pixel_centers(greendog_b_sf, plot_map = "figures/lsat_get_pixel_centers.png")


green_dog_c <- st_polygon(
  list(matrix(c(15.9396, 78.1753,
                15.9396, 78.1670,
                15.9784, 78.1670,
                15.9784, 78.1753,
                15.9396, 78.1753),
              ncol = 2, byrow = TRUE)))
greendog_c_sf <- st_sfc(green_dog_c, crs = 4326) %>% st_sf()

# Use lsat_get_pixel_centers to retrieve pixel centers and plot to a file that can be added to this documentation.
# We set plot_map to a file path (or just TRUE) to view 
greendog_c_poly <- lsat_get_pixel_centers(greendog_c_sf, plot_map = "figures/lsat_get_pixel_centers.png")

test_regions_sf <- st_sfc(green_dog_a, green_dog_b, green_dog_c, crs = 4326) %>% st_sf() %>%
  mutate(region = c("green_a", "green_b", "green_c"))

# Export time-series using lsat_export_ts()
task_list <- lsat_export_ts(pixel_coords_sf = greendog_c_poly, startJulian = 152, endJulian = 273,
                            file_prefix = 'greendog_c', drive_export_dir = 'earth_engine/lsat_greendog_c')

# Create a list of data files exported from GEE and then read them in to R as a data.table object 
data.files <- list.files('data/gee_output/earth_engine-lsat_greendog_c', full.names = T, pattern = 'greendog_c')
lsat.dt <- do.call("rbind", lapply(data.files,data.table::fread))

# Format the exported data
lsat.dt <- lsat_general_prep(lsat.dt)

# Clean the data by filtering out clouds, snow, and water, as well as radiometric and geometric errors
lsat.dt <- lsat_clean_data(lsat.dt)

# Summarize the availability of Landsat data for each pixel
lsat_summarize_data_avail(lsat.dt)
ggsave('figures/figure_greendoc_c_observation_density.jpg', width = 6, height = 4, units = 'in', dpi = 400)

# Compute the Normalized Difference Vegetation Index (NDVI)
lsat.dt <- lsat_calc_spec_index(lsat.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using random forest models and overwrite data in the NDVI column  
lsat.dt <- lsat_calibrate_rf(lsat.dt, band.or.si = 'ndvi', doy.rng = 152:273, 
                             train.with.highlat.data = T, outdir = 'output/ndvi_xcal_smry/', overwrite.col = T)

# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, si = 'ndvi')
# Fit phenological models (cubic splines) to each time series
ggsave('figures/figure_greendog_c_phenological_curves.jpg', width = 9, height = 7, units = 'in', dpi = 400)

# Summarize vegetation index for the "growing season", including estimating annual max vegetation index
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, si = 'ndvi', min.frac.of.max = 0.75)

# Evaluate estimates of annual maximum NDVI
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, si = 'ndvi', min.obs = 5, reps = 2, min.frac.of.max = 0.75)
ggsave('figures/figure_greendog_c_ndvi_max_evaluation.jpg', width = 6, height = 4, units = 'in', dpi = 400)


# compute trends in annual greening
lsat.trnds.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', 2006:2021, sig = 0.1)

# Write out data.table with growing season summaries
write(lsat.gs.dt, 'data/lsat_greendogc_annual_growing_season_summaries.csv')

# Convert to simple feature and write out shapefile
lsat.trend.sf <- lsat.trend.dt %>% st_as_sf(coords=c('longitude', 'latitude'), crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
lsat.trend.sf <- lsat.trend.sf %>% st_transform(crs = 3413)
st_write(lsat.trend.sf, dsn = 'data/lsat_ndvimax_trends.shp',append=FALSE)

# Create time series figure for each trend class
lsat.gs.dt <- lsat.gs.dt[lsat.trend.dt, on = 'sample.id']

lsat.trend.cls.yrly.dt <- lsat.gs.dt[, .(ndvi.max.avg = mean(ndvi.max,na.rm = T), ndvi.max.sd = sd(ndvi.max, na.rm = T),
                                         n = .N), by = c('trend.cat','year')]

lsat.trend.cls.yrly.dt[, ndvi.max.se := ndvi.max.sd/sqrt(n)]


lsat.trend.cls.yrly.dt[, trend.cat := factor(trend.cat, levels = c('browning','no_trend','greening'), 
                                             labels = c('browning','no trend','greening'))]
trend.cols <- c('darkgoldenrod4','ivory3','darkgreen')

ggplot(lsat.trend.cls.yrly.dt[year != 2017], aes(year, ndvi.max.avg, group = trend.cat, color = trend.cat)) + 
  ylim(0.40,0.65) + 
  labs(y='Landsat NDVImax', x='Year') + 
  geom_ribbon(aes(ymin=ndvi.max.avg-ndvi.max.se,ymax=ndvi.max.avg+ndvi.max.se, fill=trend.cat),alpha=0.3, linetype=0)+
  geom_line(aes(color = trend.cat), alpha = 1, size=1) + 
  scale_fill_manual(values = trend.cols, name = 'Trend class') + 
  scale_color_manual(values = trend.cols, name = 'Trend class')+
  theme_bw() +
  theme(legend.position=c(0.8, 0.2), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14),
        plot.title=element_text(hjust = 0.5))



ggplot(lsat.gs.dt) +
  aes(x = year, y = ndvi.max) +
  geom_point(size = 1L) +
  geom_smooth(method = "lm") +
  scale_fill_manual(values = "darkgreen") +               
  scale_colour_manual(values = "darkgreen") +   
  labs(x = "Year", y = "NDVI-max", title = "NDVImax at the Green Dog Yard \"C\"") +
  theme_classic()

#### 7b - LSAT Trial Points ####

# import all sites
animal_sites_svalbard <- read_csv("data/animal_sites_svalbard.csv")

# rename site to 'sample id'
animal_sites_svalbard$sample_id <- animal_sites_svalbard$site

# convert to spatialdataframe
test_points_sf <- st_as_sf(animal_sites_svalbard,coords=c("x","y"), crs = 4326)
crs(test_points_sf)
# make sure points are correct
plot(test_points_sf$geometry)

# Export time-series using lsat_export_ts()
task_list <- lsat_export_ts(test_points_sf, startJulian = 120, endJulian = 273,
                            file_prefix = 'animal_point', drive_export_dir = 'earth_engine/lsat_points')

# check progress of GEE extractions
ee_monitoring() 

data.files <- list.files('data/gee_output/earth_engine-lsat_points', full.names = T, pattern = 'animal')
# prepare data table
lsat.dt <- do.call("rbind", lapply(data.files, data.table::fread))

# setnames(lsat.dt, 'my_unique_location_column','sample.id') 
lsat.dt <- lsat_general_prep(lsat.dt)

# clean surface level reflectance data of snow, cloud and water using GEE masks
lsat.dt <- lsat_clean_data(lsat.dt)

subset(df, date>4 & date<6)


# calculate 'neighbourhood mean' of 3x3 pixels to get surroundings around yards
lsat.dt <- lsat_neighborhood_mean(lsat.dt)

# get summary information about total observations over time period
data.summary.dt <- lsat_summarize_data_avail(lsat.dt)
data.summary.dt

# Compute NDVI or other vegetation index
lsat.dt <- lsat_calc_spec_index(lsat.dt, si = 'ndvi')



# Cross-calibrate NDVI among sensors using an approach based on Random Forest machine learning
lsat.dt <- lsat_calibrate_rf(lsat.dt, band.or.si = 'ndvi', doy.rng = 121:273, train.with.highlat.data = T, outdir = 'gee_outputoutput/ndvi_xcal_smry/', overwrite.col = T)


# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, si = 'ndvi')
ggsave('figures/figure_points_phenological_curves.jpg', width = 9, height = 7, units = 'in', dpi = 400)

# Summarize vegetation index for the "growing season", including estimating annual max vegetation index
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, si = 'ndvi', min.frac.of.max = 0.75)

# Evaluate estimates of annual maximum NDVI
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, si = 'ndvi', min.obs = 5, reps = 2, min.frac.of.max = 0.75)
ggsave('figures/figure_greendog_c_ndvi_max_evaluation.jpg', width = 6, height = 4, units = 'in', dpi = 400)


# compute trends in annual greening
lsat.trnds.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', 2006:2021, sig = 0.1)


#### 7c - Doughnut Polygons ####
yard <- vect("data/yard_ring.shp")
plot(yard)

# create buffer
buffer <-  terra::buffer(yard, width = 50, quadsegs = 5)
plot(buffer)

# create doughnut shapes
doughnut <- erase(buffer, yard)
plot(doughnut)

# buffer check
outfile <- 'data/yard_ring.shp'
writeVector(doughnut, outfile, overwrite=TRUE) # yep - plotted it in qgis and we have beautiful beautiful doughnuts

# now re-read in file as a sf object and send to GEE
dog_poly <- read_sf('data/yard_ring.shp')

dog_poly <- rename(dog_poly, yard = sample.id)
plot(dog_poly)
dog_poly$region <- paste0("unique_id_", 1:nrow(dog_poly))

dog_poly$region <- as.factor(dog_poly$region)

#reorder columns
dog_poly <- dog_poly[, c(3,4,2)]

dog_poly <- dog_poly %>% st_transform(crs = 4326)
crs(dog_poly)

# split
DY_10 <- subset(dog_poly, yard == "DY_10")
DY_11 <- subset(dog_poly, yard == "DY_11")
DY_12 <- subset(dog_poly, yard == "DY_12")
DY_13 <- subset(dog_poly, yard == "DY_13")
DY_3 <- subset(dog_poly, yard == "DY_3")
DY_4 <- subset(dog_poly, yard == "DY_4")
DY_5 <- subset(dog_poly, yard == "DY_5")
DY_6 <- subset(dog_poly, yard == "DY_6")
DY_7 <- subset(dog_poly, yard == "DY_7")
DY_8 <- subset(dog_poly, yard == "DY_8")
DY_9 <- subset(dog_poly, yard == "DY_9")
DY_bar <- subset(dog_poly, yard == "DY_bar")
HH_barcows <- subset(dog_poly, yard == "HH_barcows")
HH_dogs <- subset(dog_poly, yard == "HH_dogs")
HH_lyr <- subset(dog_poly, yard == "HH_lyr")
ST_lyr <- subset(dog_poly, yard == "ST_lyr")

# Use lsat_get_pixel_centers to retrieve pixel centers and plot to a file that can be added to this documentation.
# We set plot_map to a file path (or just TRUE) to view 
DY_10_poly <- lsat_get_pixel_centers(DY_10, buffer = 0, plot_map = T)
DY_11_poly <- lsat_get_pixel_centers(DY_11, buffer = 0, plot_map = T)
DY_12_poly <- lsat_get_pixel_centers(DY_12, buffer = 0, plot_map = T)
DY_13_poly <- lsat_get_pixel_centers(DY_13, buffer = 0, plot_map = T)
DY_3_poly <- lsat_get_pixel_centers(DY_3, buffer = 0, plot_map = T)
DY_4_poly <- lsat_get_pixel_centers(DY_4, buffer = 0, plot_map = T)
DY_5_poly <- lsat_get_pixel_centers(DY_5, buffer = 0, plot_map = T)
DY_6_poly <- lsat_get_pixel_centers(DY_6, buffer = 0, plot_map = T)
DY_7_poly <- lsat_get_pixel_centers(DY_7, buffer = 0, plot_map = T)
DY_8_poly <- lsat_get_pixel_centers(DY_8, buffer = 0, plot_map = T)
DY_9_poly <- lsat_get_pixel_centers(DY_9, buffer = 0, plot_map = T)
DY_bar_poly <- lsat_get_pixel_centers(DY_bar, buffer = 0, plot_map = T)
HH_barcows_poly <- lsat_get_pixel_centers(HH_barcows, buffer = 0, plot_map = T)
HH_dogs_poly <- lsat_get_pixel_centers(HH_dogs, buffer = 0, plot_map = T)
HH_lyr_poly <- lsat_get_pixel_centers(HH_lyr, buffer = 0, plot_map = T)
ST_lyr_poly <- lsat_get_pixel_centers(ST_lyr, buffer = 0, plot_map = T)

# add "_key" to sample_ID
DY_10_poly$sample_id <- paste((DY_10_poly$sample_id), "_DY_10")
DY_11_poly$sample_id <- paste((DY_11_poly$sample_id), "_DY_11")
DY_12_poly$sample_id <- paste((DY_12_poly$sample_id), "_DY_12")
DY_13_poly$sample_id <- paste((DY_13_poly$sample_id), "_DY_13")
DY_3_poly$sample_id <- paste((DY_3_poly$sample_id), "_DY_3")
DY_4_poly$sample_id <- paste((DY_4_poly$sample_id), "_DY_4")
DY_5_poly$sample_id <- paste((DY_5_poly$sample_id), "_DY_5")
DY_6_poly$sample_id <- paste((DY_6_poly$sample_id), "_DY_6")
DY_7_poly$sample_id <- paste((DY_7_poly$sample_id), "_DY_7")
DY_8_poly$sample_id <- paste((DY_8_poly$sample_id), "_DY_8")
DY_9_poly$sample_id <- paste((DY_9_poly$sample_id), "_DY_9")
DY_bar_poly$sample_id <- paste((DY_bar_poly$sample_id), "_DY_bar")
HH_barcows_poly$sample_id <- paste((HH_barcows_poly$sample_id), "_HH_barcows")
HH_dogs_poly$sample_id <- paste((HH_dogs_poly$sample_id), "_HH_dogs")
HH_lyr_poly$sample_id <- paste((HH_lyr_poly$sample_id), "_HH_lyr")
ST_lyr_poly$sample_id <- paste((ST_lyr_poly$sample_id), "_ST_lyr")


# Bind Rows

yard_ring_pix <- bind_rows(DY_10_poly, DY_11_poly, DY_12_poly, DY_13_poly,
                           DY_3_poly, DY_4_poly, DY_5_poly, DY_6_poly,
                           DY_7_poly, DY_8_poly, DY_9_poly, DY_bar_poly,
                           HH_barcows_poly,HH_dogs_poly,HH_lyr_poly,ST_lyr_poly)



# Extract a time-series of Landsat surface reflectance measurements for each Landsat pixel
task_list <- lsat_export_ts(pixel_coords_sf = yard_ring_pix, startJulian = 120, endJulian = 273,
                            file_prefix = 'doughnut', drive_export_dir = 'earth_engine/lsat_yardring')

# check status
ee_monitoring()




# Create a list of data files exported from GEE and then read them in to R as a data.table object 
data.files <- list.files('data/earth_engine-lsat_yardring', full.names = T, pattern = 'doughnut')
lsat.dt <- do.call("rbind", lapply(data.files, fread))

# Format the exported data
lsat.dt <- lsat_general_prep(lsat.dt)

# Clean the data by filtering out clouds, snow, and water, as well as radiometric and geometric errors
lsat.dt <- lsat_clean_data(lsat.dt, geom.max = 15, cloud.max = 80, sza.max = 130, 
                           filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)

# Summarize the availability of Landsat data for each pixel
lsat_summarize_data_avail(lsat.dt)
ggsave('figures/figure_yardring_observation_density.jpg', width = 6, height = 4, units = 'in', dpi = 400)


# Compute the Normalized Difference Vegetation Index (NDVI)
lsat.dt <- lsat_calc_spec_index(lsat.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using random forest models and overwrite data in the NDVI column  
lsat.dt <- lsat_calibrate_rf(lsat.dt, band.or.si = 'ndvi', doy.rng = 120:270, 
                             train.with.highlat.data = T, outdir = 'output/ndvi_xcal_smry/', overwrite.col = T)

# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, si = 'ndvi', test.run = F)
ggsave('figures/figure_yardring_phenological_curves.jpg', width = 9, height = 7, units = 'in', dpi = 400)

# Summarize vegetation index for the "growing season", including estimating annual max vegetation index
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, si = 'ndvi', min.frac.of.max = 0.75)


# Evaluate estimates of annual maximum NDVI
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, si = 'ndvi', min.obs = 5, reps = 2, min.frac.of.max = 0.75)
ggsave('figures/figure_yardring_ndvi_max_evaluation.jpg', width = 6, height = 4, units = 'in', dpi = 400)

# Write out data.table with growing season summaries
fwrite(lsat.gs.dt, 'output/lsat_annual_growing_season_summaries.csv')

# Compute temporal trends in NDVImax
lsat.trend.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', yrs = 2004:2021, legend.position = c(0.66,0.93))
ggsave('figures/figure_yardring_ndvi_max_trend_distribution.jpg', width = 6, height = 8, units = 'in', dpi = 400)

# Convert trend data table to simple feature and write out shapefile
lsat.trend.sf <- lsat.trend.dt %>% st_as_sf(coords=c('longitude', 'latitude'), crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
lsat.trend.sf <- lsat.trend.sf %>% st_transform(crs = 3413)
st_write(lsat.trend.sf, dsn = 'data/lsat_ndvimax_trends_doughnut.shp')



lsat.gs.dt <- lsat.gs.dt %>%
  mutate(yard = str_after_nth(sample.id, "_", 2)) # new column with yard ID

# rename DY_10 _DY_10 to DY_10 
levels(lsat.gs.dt$yard)[levels(lsat.gs.dt$yard)=="DY_10 _DY_10"] <- "DY_10"

# get mean ndvi.max for each plot
lsat.gs.dt2 <- lsat.gs.dt %>%
  group_by(yard) %>%
  summarise(Mean.nd = mean(ndvi.max))

# NDVImax histogram for each site
(hist <- ggplot(lsat.gs.dt) +
  aes(x = ndvi.max) +
  geom_histogram(bins = 30L, fill = "#1f9e89") +
  geom_vline(data = lsat.gs.dt2, mapping = aes(xintercept = Mean.nd), colour = "blue", linetype = "dashed") +
  theme_classic() +
  facet_wrap(vars(yard)))

hist + geom_vline(data=lsat.gs.dt, aes(yintercept=ndvi.max, color=yard),
                  linetype="dashed")

# NDVI Changes over time
(ndvi_change <- ggplot(lsat.gs.dt) +
  aes(x = year, y = ndvi.max, colour = ndvi.max) +
  geom_point(size = 1L) +
  geom_smooth(method=lm, color="gold") +
  scale_color_distiller(palette = "Greens", direction = 1) +
  labs(x = "Year", y = "NDVI Max", color = "NDVI Max") +
  theme_classic() +
  facet_wrap(vars(yard)))


lsat.trend.dt <- lsat.trend.dt %>%
  mutate(yard = str_after_nth(sample.id, "_", 2)) # new column with yard ID

lsat.pheno.dt <- lsat.pheno.dt %>%
  mutate(yard = str_after_nth(sample.id, "_", 2)) # new column with yard ID

# greening curves by yard
(curve_yard <- ggplot(lsat.pheno.dt) +
  aes(x = doy, y = ndvi, colour = year) +
  geom_point(size = 1L) +
  scale_color_viridis_c(option = "viridis") +
    labs(y= 'NDVI ', x='Day of Year') + 
    theme_classic() +
  facet_wrap(vars(yard)))

ggplot(lsat.trend.dt) +
  aes(x = longitude, y = latitude, colour = total.change.pcnt) +
  geom_point(size = 3L, shape = 15) +
  scale_color_distiller(palette = "BrBG", direction = -2) +
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")

ggplot(lsat.gs.dt) +
  aes(x = longitude, y = latitude, colour = ndvi.max) +
  geom_point(size = 3L, shape = 15) +
  cale_color_gradientn(name = 'NDVI max',  colours = c('gold','grey','green')) + 
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")



