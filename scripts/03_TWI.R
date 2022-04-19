#### ELISE GALLOIS  - elise.gallois94@gmail.com 



#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(esquisse)
library(sf)
library(raster)
library(viridis)
library(readxl)
library(spatstat)  
library(maptools)  
library(dynatopmodel)
library(ggeffects)

#### 2 - LOAD DATA ####
longyear_dem <- raster('data/longyearbyen_arcticdem_2m.tiff')
plot(longyear_dem)

#### 3 - Make contour map ####
contour <- rasterToContour(longyear_dem)
class(contour)
plot(contour)
plot(longyear_dem, add=TRUE)

#### Topographic Wetness Index ####
# resample raster to ensure x and y res has same res
r <- raster(ext=extent(longyear_dem), res = 2) 
dem <- raster::resample(longyear_dem, r, method = "bilinear")
# calculate wetness index using dynatopmodel package
twi <- upslope.area(dem, log = TRUE, atb = TRUE, deg = 0.1, fill.sinks = FALSE)
sp::plot(twi$atb, main=c("Topographic Wetness Index log (m^2/m)"))
hist(twi$atb)

#save mean twi as rasters
writeRaster(twi$atb, 'data/spatial/twi.tif', format = 'GTiff')
# code to load data
twi <- raster('data/spatial/twi.tif')

