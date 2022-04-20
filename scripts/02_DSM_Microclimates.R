#### ELISE GALLOIS  - elise.gallois94@gmail.com 



#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(esquisse)
library(sf)
library(raster)
library(viridis)
library(readxl)
library(microclima)
library(NicheMapR)

#### 2 - LOAD DATA ####
# Get DEM for Longyearbyen Dog Yards
r <- get_dem(lat = 78.17162812, long = 15.95069815, resolution = 30,xdims = 200, ydims = 200)
plot(r)

# Takes ~ c. 5 minutes to run
temps <- runauto(r, "30/07/2019", "15/08/2019", hgt = 0.1, l = NA, x = NA,
                 habitat = "Barren or sparsely vegetated")

mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "orange",
                            "red"))(255)
meantemp <- temps$tmean
extent(meantemp) <- temps$e
par(mfrow = c(1, 1))
plot(meantemp, main = "Mean temperature", col = mypal)
# Interactive 3D plot
require(plotly)
zrange<-list(range = c(0, 3000))
plot_ly(z = ~is_raster(r)) %>%
  add_surface(surfacecolor = ~is_raster(meantemp)) %>%
  layout(scene = list(zaxis = zrange))