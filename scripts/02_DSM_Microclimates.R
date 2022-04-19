#### ELISE GALLOIS  - elise.gallois94@gmail.com 



#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(esquisse)
library(sf)
library(raster)
library(viridis)
library(readxl)

#### 2 - LOAD DATA ####
longyear_dem <- raster('data/longyearbyen_arcticdem_2m.tiff')
plot(dem)
