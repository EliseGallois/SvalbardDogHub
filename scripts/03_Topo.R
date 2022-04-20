#### ELISE GALLOIS  - elise.gallois94@gmail.com 



#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(esquisse)
library(sf)
library(raster)
library(rgdal)
library(viridis)
library(readxl)
library(spatstat)  
library(maptools)  
library(dynatopmodel)
library(ggeffects)
library(insol) 


#### 2 - LOAD DATA ####
longyear_dem <- raster('data/longyearbyen_arcticdem_2m.tiff')
plot(longyear_dem)

barents_dem <- raster('data/barentsburg_arcticdem_2m.tiff')
plot(barents_dem)

#### 3 - Make contour map ####
contour <- rasterToContour(longyear_dem)
class(contour)
plot(contour)
plot(longyear_dem, add=TRUE)


#### Get slope, hillshade, aspect ####

# tifs of slope
long_slope <- raster::terrain(longyear_dem, opt = 'slope', unit = 'degrees')  #calculate slope
plot(long_slope)
# write new slope raster
raster::writeRaster(long_slope, 'data/longyear_slope.tif', format = 'GTiff', overwrite = TRUE)

barents_slope <- raster::terrain(barents_dem, opt = 'slope', unit = 'degrees')  #calculate slope
plot(barents_slope)
# write new slope raster
raster::writeRaster(barents_slope, 'data/barents_slope_slope.tif', format = 'GTiff', overwrite = TRUE)

# tifs of aspect
long_aspect <- raster::terrain(longyear_dem, opt = 'aspect', unit = 'degrees') #calculate aspect
plot(long_aspect)
# write new aspect raster
raster::writeRaster(long_aspect, 'data/longyear_aspect.tif', format = 'GTiff', overwrite = TRUE)

barents_aspect <- raster::terrain(barents_dem, opt = 'aspect', unit = 'degrees') #calculate aspect
plot(barents_aspect)
# write new aspect raster
raster::writeRaster(barents_aspect, 'data/barents_aspect.tif', format = 'GTiff', overwrite = TRUE)









#### Topographic Wetness Index ####
# resample raster to ensure x and y res has same res
r <- raster(ext=extent(longyear_dem), res = 2) 
dem <- raster::resample(longyear_dem, r, method = "bilinear")


newproj <- "+proj=utm +zone=33X +datum=WGS84  +units=m"
dem <- projectRaster(longyear_dem, crs=newproj, res = 2)

# calculate wetness index using dynatopmodel package
twi <- upslope.area(longyear_dem, log = TRUE, atb = TRUE, deg = 0.1, fill.sinks = FALSE)
sp::plot(twi$atb, main=c("Topographic Wetness Index log (m^2/m)"))
hist(twi$atb)

#save mean twi as rasters
writeRaster(twi$atb, 'data/spatial/twi.tif', format = 'GTiff')
# code to load data
twi <- raster('data/spatial/twi.tif')

