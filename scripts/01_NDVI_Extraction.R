#### ELISE GALLOIS  - elise.gallois94@gmail.com 

# This script uses Uses lsatTS package to extract landsat 30x30m NDVI values between 1985-2021
# Uses lsatTS package to calculate NDVI max and date of NDVI max across all sites
# We save an NDVImax and Phenology csv for all sites in ' Step 6'

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
library(geojsonio)

#### 2 - connect R to the GEE ####
ee_Initialize()
db <- 'CGIAR/SRTM90_V4'
image <- ee$Image(db)
image$bandNames()$getInfo()


#### 3 - LSAT Tool- Doughnut Polygons - ACTIVE YARDS ####
yard <- vect("data/yard_para.shp")
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

# Reproducibility note: this will generate a data folder in your google drive. 
# Zip this folder and save it in the "data/" folder of your SvalbardDogHub repo
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
lsat.dt <- lsat_calibrate_rf(lsat.dt, band.or.si = 'ndvi', doy.rng = 120:240, 
                             train.with.highlat.data = T, outdir = 'output/ndvi_xcal_smry/', overwrite.col = T)
fwrite(lsat.dt, 'output/Apheno_dog.csv') # "raw" data aka calibrated but no pheno cubic splines fitted

# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, si = 'ndvi', window.min.obs = 4, test.run = F)
ggsave('figures/figure_yardring_phenological_curves.jpg', width = 9, height = 7, units = 'in', dpi = 400)
fwrite(lsat.pheno.dt, 'output/pheno_dog.csv')

# Summarize vegetation index for the "growing season", including estimating annual max vegetation index
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, si = 'ndvi', min.frac.of.max = 0.75)


# Evaluate estimates of annual maximum NDVI
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, si = 'ndvi', min.obs = 5, reps = 2, min.frac.of.max = 0.75)
ggsave('figures/figure_yardring_ndvi_max_evaluation.jpg', width = 6, height = 4, units = 'in', dpi = 400)

# Write out data.table with growing season summaries
fwrite(lsat.gs.dt, 'output/lsat_annual_growing_season_summaries.csv')

# Compute temporal trends in NDVImax
lsat.trend.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', yr.tolerance = 1, 
                                 yrs = 2000:2021,  legend.position = c(0.66,0.93))
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
  aes(x = longitude, y = latitude, colour = slope) +
  geom_point(size = 3L, shape = 15) +
  scale_color_distiller(palette = "BrBG", direction = -2) +
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")

ggplot(lsat.gs.dt) +
  aes(x = longitude, y = latitude, colour = ndvi.max) +
  geom_point(size = 3L, shape = 15) +
  scale_color_gradientn(name = 'NDVI max',  colours = c('gold','grey','green')) + 
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")


#### 4 - LSAT Tool-  Polygons - REF SITES ####
ref <- vect("data/ref_sites.shp")
plot(ref)

# create buffer
buffer <-  terra::buffer(ref, width = 50, quadsegs = 5)
plot(buffer)


# buffer check
outfile <- 'data/ref_ring.shp'
writeVector(buffer, outfile, overwrite=TRUE) # yep - plotted it in qgis and looks good!

# now re-read in file as a sf object and send to GEE
ref_poly <- read_sf('data/ref_ring.shp')

plot(ref_poly$geometry)
ref_poly$region <- paste0("unique_id_", 1:nrow(ref_poly))

ref_poly$region <- as.factor(ref_poly$region)
colnames(ref_poly)

#reorder columns
ref_poly <- ref_poly[, c(1,4,5)]

# make sure all in lat/lon
ref_poly <- ref_poly %>% st_transform(crs = 4326)

# split
REF_fest <- subset(ref_poly, site == "REF_fest")
REF_tem <- subset(ref_poly, site == "REF_tem")
REF_bjor <- subset(ref_poly, site == "REF_bjor")
BC_Fjort <- subset(ref_poly, site == "BC_Fjort")
BC_Oss <- subset(ref_poly, site == "BC_Oss")
BC_stu <- subset(ref_poly, site == "BC_stu")
BC_skans <- subset(ref_poly, site == "BC_skans")
BC_alkhor <- subset(ref_poly, site == "BC_alkhor")

# Use lsat_get_pixel_centers to retrieve pixel centers and plot to a file that can be added to this documentation.
# We set plot_map to a file path (or just TRUE) to view 
REF_fest_poly <- lsat_get_pixel_centers(REF_fest, buffer = 0, plot_map = T)
REF_tem_poly <- lsat_get_pixel_centers(REF_tem, buffer = 0, plot_map = T)
REF_bjor_poly <- lsat_get_pixel_centers(REF_bjor, buffer = 0, plot_map = T)
BC_Fjort_poly <- lsat_get_pixel_centers(BC_Fjort, buffer = 0, plot_map = T)
BC_Oss_poly <- lsat_get_pixel_centers(BC_Oss, buffer = 0, plot_map = T)
BC_stu_poly <- lsat_get_pixel_centers(BC_stu, buffer = 0, plot_map = T)
BC_skans_poly <- lsat_get_pixel_centers(BC_skans, buffer = 0, plot_map = T)
BC_alkhor_poly <- lsat_get_pixel_centers(BC_alkhor, buffer = 0, plot_map = T)


# add "_key" to sample_ID
REF_fest_poly$sample_id <- paste((REF_fest_poly$sample_id), "_REF_fest")
REF_tem_poly$sample_id <- paste((REF_tem_poly$sample_id), "_REF_tem")
REF_bjor_poly$sample_id <- paste((REF_bjor_poly$sample_id), "_REF_bjor")
BC_Fjort_poly$sample_id <- paste((BC_Fjort_poly$sample_id), "_BC_Fjort")
BC_Oss_poly$sample_id <- paste((BC_Oss_poly$sample_id), "_BC_Oss")
BC_stu_poly$sample_id <- paste((BC_stu_poly$sample_id), "_BC_stu")
BC_skans_poly$sample_id <- paste((BC_skans_poly$sample_id), "_BC_skans")
BC_alkhor_poly$sample_id <- paste((BC_alkhor_poly$sample_id), "_BC_alkhor")

# Bind Rows
ref_ring_pix <- bind_rows(REF_fest_poly, REF_tem_poly, REF_bjor_poly, BC_Fjort_poly,
                          BC_Oss_poly, BC_stu_poly, BC_skans_poly, BC_alkhor_poly)


# Extract a time-series of Landsat surface reflectance measurements for each Landsat pixel
task_list <- lsat_export_ts(pixel_coords_sf = ref_ring_pix, startJulian = 120, endJulian = 273,
                            file_prefix = 'ref_bird', drive_export_dir = 'earth_engine/lsat_refring')

# check status
ee_monitoring()


# Create a list of data files exported from GEE and then read them in to R as a data.table object 
data.files <- list.files('data/earth_engine-lsat_yardring', full.names = T, pattern = 'ref_bird')
lsat.dt <- do.call("rbind", lapply(data.files, fread))


# Format the exported data
lsat.dt <- lsat_general_prep(lsat.dt)

# Clean the data by filtering out clouds, snow, and water, as well as radiometric and geometric errors
lsat.dt <- lsat_clean_data(lsat.dt, geom.max = 15, cloud.max = 80, sza.max = 130, 
                           filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)

# Summarize the availability of Landsat data for each pixel
lsat_summarize_data_avail(lsat.dt)
ggsave('figures/figure_refring_observation_density.jpg', width = 6, height = 4, units = 'in', dpi = 400)


# Compute the Normalized Difference Vegetation Index (NDVI)
lsat.dt <- lsat_calc_spec_index(lsat.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using random forest models and overwrite data in the NDVI column  
lsat.dt <- lsat_calibrate_rf(lsat.dt, band.or.si = 'ndvi', doy.rng = 120:270, 
                             train.with.highlat.data = T, outdir = 'output/ndvi_xcal_smry_ref/', overwrite.col = T)
fwrite(lsat.dt, 'output/Apheno_ref.csv') # "raw" data aka calibrated but no pheno cubic splines fitted

# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, si = 'ndvi', test.run = F)
ggsave('figures/figure_refring_phenological_curves.jpg', width = 9, height = 7, units = 'in', dpi = 400)
fwrite(lsat.pheno.dt, 'output/pheno_ref.csv')

# Summarize vegetation index for the "growing season", including estimating annual max vegetation index
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, si = 'ndvi', min.frac.of.max = 0.75)


# Evaluate estimates of annual maximum NDVI
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, si = 'ndvi', min.obs = 5, reps = 2, min.frac.of.max = 0.75)
ggsave('figures/figure_refring_ndvi_max_evaluation.jpg', width = 6, height = 4, units = 'in', dpi = 400)

# Write out data.table with growing season summaries
fwrite(lsat.gs.dt, 'output/lsat_annual_growing_season_summaries_REF.csv')

# Compute temporal trends in NDVImax
lsat.trend.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', yr.tolerance = 1, 
                                 yrs = 1990:2021, nyr.min.frac = 0.4,
                                 legend.position = c(0.66,0.93))
ggsave('figures/figure_refring_ndvi_max_trend_distribution.jpg', width = 6, height = 8, units = 'in', dpi = 400)


# Convert trend data table to simple feature and write out shapefile
lsat.trend.sf <- lsat.trend.dt %>% st_as_sf(coords=c('longitude', 'latitude'), crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
lsat.trend.sf <- lsat.trend.sf %>% st_transform(crs = 3413)
st_write(lsat.trend.sf, dsn = 'data/lsat_ndvimax_trends_doughnut.shp')

lsat.gs.dt <- lsat.gs.dt %>%
  mutate(yard = str_after_nth(sample.id, "_", 2)) # new column with yard ID


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
  aes(x = longitude, y = latitude, colour = slope) +
  geom_point(size = 3L, shape = 15) +
  scale_color_distiller(palette = "BrBG", direction = -2) +
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")

ggplot(lsat.gs.dt) +
  aes(x = longitude, y = latitude, colour = ndvi.max) +
  geom_point(size = 3L, shape = 15) +
  scale_color_gradientn(name = 'NDVI max',  colours = c('gold','grey','green')) + 
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")


#### 5 - LSAT Tool-  Polygons - EXTRA SITES ####
pig <- vect("data/pyr.shp") # added later - pig yard coordinates for pyramiden
plot(pig)

# create buffer
buffer <-  terra::buffer(pig, width = 50, quadsegs = 5)
plot(buffer)


# buffer check
outfile <- 'data/pig_ring.shp'
writeVector(buffer, outfile, overwrite=TRUE) # yep - plotted it in qgis and looks good!

# now re-read in file as a sf object and send to GEE
pig_poly <- read_sf('data/pig_ring.shp')

plot(pig_poly$geometry)
pig_poly$region <- paste0("unique_id_", 1:nrow(pig_poly))

pig_poly$region <- as.factor(pig_poly$region)


# make sure all in lat/lon
pig_poly <- pig_poly %>% st_transform(crs = 4326)

# get pixel centres
pig_poly <- lsat_get_pixel_centers(pig_poly, buffer = 0, plot_map = T)


# add "_key" to sample_ID
pig_poly$sample_id <- paste((pig_poly$sample_id), "_PYR_pig")

# Extract a time-series of Landsat surface reflectance measurements for each Landsat pixel
task_list <- lsat_export_ts(pixel_coords_sf = pig_poly, startJulian = 120, endJulian = 273,
                            file_prefix = 'ref_pig', drive_export_dir = 'earth_engine/lsat_pigring')

# check status
ee_monitoring()

# Create a list of data files exported from GEE and then read them in to R as a data.table object 
data.files <- list.files('data/earth_engine-lsat_yardring', full.names = T, pattern = 'ref_pig')
lsat.dt <- do.call("rbind", lapply(data.files, fread))


# Format the exported data
lsat.dt <- lsat_general_prep(lsat.dt)

# Clean the data by filtering out clouds, snow, and water, as well as radiometric and geometric errors
lsat.dt <- lsat_clean_data(lsat.dt, geom.max = 15, cloud.max = 80, sza.max = 130, 
                           filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)

# Summarize the availability of Landsat data for each pixel
lsat_summarize_data_avail(lsat.dt)
ggsave('figures/figure_refpig_observation_density.jpg', width = 6, height = 4, units = 'in', dpi = 400)


# Compute the Normalized Difference Vegetation Index (NDVI)
lsat.dt <- lsat_calc_spec_index(lsat.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using random forest models and overwrite data in the NDVI column  
lsat.dt <- lsat_calibrate_rf(lsat.dt, band.or.si = 'ndvi', doy.rng = 120:270, 
                             train.with.highlat.data = T, outdir = 'output/ndvi_xcal_smry_pig/', overwrite.col = T)
fwrite(lsat.dt, 'output/Apheno_pig.csv') # "raw" data aka calibrated but no pheno cubic splines fitted

# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, si = 'ndvi', test.run = F)
ggsave('figures/figure_refpig_phenological_curves.jpg', width = 9, height = 7, units = 'in', dpi = 400)
fwrite(lsat.pheno.dt, 'output/pheno_pig.csv')

# Summarize vegetation index for the "growing season", including estimating annual max vegetation index
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, si = 'ndvi', min.frac.of.max = 0.75)


# Evaluate estimates of annual maximum NDVI
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, si = 'ndvi', min.obs = 5, reps = 2, min.frac.of.max = 0.75)
ggsave('figures/figure_refpig_ndvi_max_evaluation.jpg', width = 6, height = 4, units = 'in', dpi = 400)

# Write out data.table with growing season summaries
fwrite(lsat.gs.dt, 'output/lsat_annual_growing_season_summaries_PIG.csv')

# Compute temporal trends in NDVImax
lsat.trend.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', yr.tolerance = 1, 
                                 yrs = 2004:2021, nyr.min.frac = 0.4,
                                 legend.position = c(0.66,0.93))
ggsave('figures/figure_refpig_ndvi_max_trend_distribution.jpg', width = 6, height = 8, units = 'in', dpi = 400)


lsat.gs.dt <- lsat.gs.dt %>%
  mutate(yard = str_after_nth(sample.id, "_", 2)) # new column with yard ID


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
  aes(x = longitude, y = latitude, colour = slope) +
  geom_point(size = 3L, shape = 15) +
  scale_color_distiller(palette = "BrBG", direction = -2) +
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")

ggplot(lsat.gs.dt) +
  aes(x = longitude, y = latitude, colour = ndvi.max) +
  geom_point(size = 3L, shape = 15) +
  scale_color_gradientn(name = 'NDVI max',  colours = c('gold','grey','green')) + 
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")


## Extra ref sites ##
ref <- vect("data/extra_ref.shp") # extra reference sites
plot(ref)

# create buffer
buffer <-  terra::buffer(ref, width = 50, quadsegs = 5)
plot(buffer)


# buffer check
outfile <- 'data/extra_ring.shp'
writeVector(buffer, outfile, overwrite=TRUE) # yep - plotted it in qgis and looks good!

# now re-read in file as a sf object and send to GEE
ref_poly <- read_sf('data/extra_ring.shp')

plot(ref_poly$geometry)
ref_poly$site <- paste("REF_",  (ref_poly$reference_))

ref_poly$region <- paste0("unique_id_", 1:nrow(ref_poly))

ref_poly$region <- as.factor(ref_poly$region)

colnames(ref_poly)
#reorder columns
ref_poly <- ref_poly[, c(8,7,9,2)]

# make sure all in lat/lon
ref_poly <- ref_poly %>% st_transform(crs = 4326)

# split
REF_col <- subset(ref_poly, site == "REF_ Colesbukta")
REF_odin <- subset(ref_poly, site == "REF_ Odindalen")
REF_gron <- subset(ref_poly, site == "REF_ Grøndalen")
REF_rein <- subset(ref_poly, site == "REF_ Reindalen")
REF_bolt <- subset(ref_poly, site == "REF_ Bolterdalen")
REF_lein <- subset(ref_poly, site == "REF_ Leinstrandodden")


# Use lsat_get_pixel_centers to retrieve pixel centers and plot to a file that can be added to this documentation.
# We set plot_map to a file path (or just TRUE) to view 
REF_col_poly <- lsat_get_pixel_centers(REF_col, buffer = 0, plot_map = T, lsat_WRS2_scene_bounds = "data/WRS-2_bound_world_0.kml")
REF_odin_poly <- lsat_get_pixel_centers(REF_odin, buffer = 0, plot_map = T, lsat_WRS2_scene_bounds = "data/WRS-2_bound_world_0.kml")
REF_gron_poly <- lsat_get_pixel_centers(REF_gron, buffer = 0, plot_map = T, lsat_WRS2_scene_bounds = "data/WRS-2_bound_world_0.kml")
REF_rein_poly <- lsat_get_pixel_centers(REF_rein, buffer = 0, plot_map = T, lsat_WRS2_scene_bounds = "data/WRS-2_bound_world_0.kml")
REF_bolt_poly <- lsat_get_pixel_centers(REF_bolt, buffer = 0, plot_map = T, lsat_WRS2_scene_bounds = "data/WRS-2_bound_world_0.kml")
REF_lein_poly <- lsat_get_pixel_centers(REF_lein, buffer = 0, plot_map = T, lsat_WRS2_scene_bounds = "data/WRS-2_bound_world_0.kml")



# add "_key" to sample_ID
REF_col_poly$sample_id <- paste((REF_col_poly$sample_id), "_REF_col")
REF_odin_poly$sample_id <- paste((REF_odin_poly$sample_id), "_REF_odin")
REF_gron_poly$sample_id <- paste((REF_gron_poly$sample_id), "_REF_gron")
REF_rein_poly$sample_id <- paste((REF_rein_poly$sample_id), "_REF_rein")
REF_bolt_poly$sample_id <- paste((REF_bolt_poly$sample_id), "_REF_bolt")
REF_lein_poly$sample_id <- paste((REF_lein_poly$sample_id), "_REF_lein")


# Bind Rows
ref_ring_pix <- bind_rows(REF_col_poly, REF_odin_poly, REF_gron_poly, REF_rein_poly,
                          REF_bolt_poly, REF_lein_poly)


# Extract a time-series of Landsat surface reflectance measurements for each Landsat pixel
task_list <- lsat_export_ts(pixel_coords_sf = ref_ring_pix, startJulian = 120, endJulian = 273,
                            file_prefix = 'ref_extra', drive_export_dir = 'earth_engine/lsat_refextraring')

# check status
ee_monitoring()


# Create a list of data files exported from GEE and then read them in to R as a data.table object 
data.files <- list.files('data/earth_engine-lsat_yardring', full.names = T, pattern = 'ref_extra')
lsat.dt <- do.call("rbind", lapply(data.files, fread))


# Format the exported data
lsat.dt <- lsat_general_prep(lsat.dt)

# Clean the data by filtering out clouds, snow, and water, as well as radiometric and geometric errors
lsat.dt <- lsat_clean_data(lsat.dt, geom.max = 15, cloud.max = 80, sza.max = 130, 
                           filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)

# Summarize the availability of Landsat data for each pixel
lsat_summarize_data_avail(lsat.dt)
ggsave('figures/figure_refring_observation_density.jpg', width = 6, height = 4, units = 'in', dpi = 400)


# Compute the Normalized Difference Vegetation Index (NDVI)
lsat.dt <- lsat_calc_spec_index(lsat.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using random forest models and overwrite data in the NDVI column  
lsat.dt <- lsat_calibrate_rf(lsat.dt, band.or.si = 'ndvi', doy.rng = 120:270, 
                             train.with.highlat.data = T, outdir = 'output/ndvi_xcal_smry_ref/', overwrite.col = T)
fwrite(lsat.dt, 'output/Apheno_refextra.csv') # "raw" data aka calibrated but no pheno cubic splines fitted

# Fit phenological models (cubic splines) to each time series
lsat.pheno.dt <- lsat_fit_phenological_curves(lsat.dt, si = 'ndvi', test.run = F)
ggsave('figures/figure_refring_phenological_curves.jpg', width = 9, height = 7, units = 'in', dpi = 400)
fwrite(lsat.pheno.dt, 'output/pheno_refextra.csv')

# Summarize vegetation index for the "growing season", including estimating annual max vegetation index
lsat.gs.dt <- lsat_summarize_growing_seasons(lsat.pheno.dt, si = 'ndvi', min.frac.of.max = 0.75)


# Evaluate estimates of annual maximum NDVI
lsat.gs.eval.dt <- lsat_evaluate_phenological_max(lsat.pheno.dt, si = 'ndvi', min.obs = 5, reps = 2, min.frac.of.max = 0.75)
ggsave('figures/figure_refring_ndvi_max_evaluation.jpg', width = 6, height = 4, units = 'in', dpi = 400)

# Write out data.table with growing season summaries
fwrite(lsat.gs.dt, 'output/lsat_annual_growing_season_summaries_REFextra.csv')

# Compute temporal trends in NDVImax
lsat.trend.dt <- lsat_calc_trend(lsat.gs.dt, si = 'ndvi.max', yr.tolerance = 1, 
                                 yrs = 2004:2021, nyr.min.frac = 0.4,
                                 legend.position = c(0.66,0.93))
ggsave('figures/figure_extraring_ndvi_max_trend_distribution.jpg', width = 6, height = 8, units = 'in', dpi = 400)


# Convert trend data table to simple feature and write out shapefile
lsat.trend.sf <- lsat.trend.dt %>% st_as_sf(coords=c('longitude', 'latitude'), crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
lsat.trend.sf <- lsat.trend.sf %>% st_transform(crs = 3413)
st_write(lsat.trend.sf, dsn = 'data/lsat_ndvimax_trends_extra.shp')

lsat.gs.dt <- lsat.gs.dt %>%
  mutate(yard = str_after_nth(sample.id, "_", 2)) # new column with yard ID


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
  aes(x = longitude, y = latitude, colour = slope) +
  geom_point(size = 3L, shape = 15) +
  scale_color_distiller(palette = "BrBG", direction = -2) +
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")

ggplot(lsat.gs.dt) +
  aes(x = longitude, y = latitude, colour = ndvi.max) +
  geom_point(size = 3L, shape = 15) +
  scale_color_gradientn(name = 'NDVI max',  colours = c('gold','grey','green')) + 
  theme_classic() +
  facet_wrap(vars(yard), scales = "free")

#### 6 - MERGE ALL SITES TOGETHER ####
# add in all the gs datasets
dog <- read_csv("output/lsat_annual_growing_season_summaries.csv")
ref <- read_csv("output/lsat_annual_growing_season_summaries_REF.csv")
ref_x <- read_csv("output/lsat_annual_growing_season_summaries_REFextra.csv")
pig <- read_csv("output/lsat_annual_growing_season_summaries_PIG.csv")

# merge into one big one with all sites
sval_green <- rbind(dog, ref, ref_x, pig)

sval_green <- sval_green %>%
  mutate(site = str_after_nth(sample.id, "_", 2)) # new column with site ID

(all_ndvimax <- ggplot(sval_green) +
  aes(x = year, y = ndvi.max, colour = ndvi.max) +
  geom_point(size = 1L) +
  geom_smooth(method=lm, color="gold") +
  scale_color_viridis_c(option = "viridis") +
  labs(x = "Year", y = "NDVImax", title = "NDVI max trends: 1985-2021", color = "NDVImax") +
  theme_classic() +
  facet_wrap(vars(site)))
ggsave('figures/all_sites_ndvimax_trend.jpg', width = 9, height = 9, units = 'in', dpi = 400)

# save csv
fwrite(sval_green, 'output/lsat_NDVImx_all.csv')

# add in all the pheno datasets WITH fitted pheno cubic splines
dog <- read_csv("output/pheno_dog.csv")
ref <- read_csv("output/pheno_ref.csv")
ref_x <- read_csv("output/pheno_refextra.csv")
pig <- read_csv("output/pheno_pig.csv")

# merge into one big one with all sites
sval_phen <- rbind(dog, ref,ref_x, pig)

sval_phen <- sval_phen %>%
  mutate(site = str_after_nth(sample.id, "_", 2)) # new column with site ID

(curve_yard <- ggplot(sval_phen) +
    aes(x = doy, y = ndvi, colour = year) +
    geom_point(size = 1L) +
    scale_color_viridis_c(option = "viridis") +
    labs(y= 'NDVI ', x='Day of Year') + 
    theme_classic() +
    facet_wrap(vars(site)))

ggsave('figures/all_sites_curve_yard.jpg', width = 9, height = 9, units = 'in', dpi = 400)

# save csv
fwrite(sval_phen, 'output/lsat_phen_all.csv') 



# add in all the pheno datasets WITHOUT cubic splines fitted
dog <- read_csv("output/Apheno_dog.csv")
ref <- read_csv("output/Apheno_ref.csv")
ref_x <- read_csv("output/Apheno_refextra.csv")
pig <- read_csv("output/Apheno_pig.csv")

# merge into one big one with all sites
sval_phen <- rbind(dog, ref,ref_x, pig)

sval_phen <- sval_phen %>%
  mutate(site = str_after_nth(sample.id, "_", 2)) # new column with site ID

(curve_yard <- ggplot(sval_phen) +
    aes(x = doy, y = ndvi, colour = year) +
    geom_point(size = 1L) +
    scale_color_viridis_c(option = "viridis") +
    labs(y= 'NDVI ', x='Day of Year') + 
    theme_classic() +
    facet_wrap(vars(site)))

ggsave('figures/all_sites_curve_yard.jpg', width = 9, height = 9, units = 'in', dpi = 400)

# save csv
fwrite(sval_phen, 'output/raw_phen_all.csv')

