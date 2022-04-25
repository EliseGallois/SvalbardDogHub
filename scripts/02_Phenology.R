#### ELISE GALLOIS  - elise.gallois94@gmail.com 


#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(esquisse)
library(sf)
library(raster)
library(viridis)
library(readxl)
library(readr)

#### 2 - LOAD DATA ####
pheno <- read_csv("output/lsat_phen_all.csv") # phenology (large file)
ndvi_max <- read_csv("output/lsat_NDVImx_all.csv") # ndvimax only (streamlined file)

# plot all 'date of 'NDVImax'
ggplot(ndvi_max) +
  aes(x = year, y = ndvi.max.doy, colour = ndvi.max) +
  geom_point(size = 1L) +
  geom_smooth(method=lm, color="gold") +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  facet_wrap(vars(site))

# plot greening curves for all sites
(curve_yard <- ggplot(pheno) +
    aes(x = doy, y = ndvi, colour = year) +
    geom_point(size = 1L) +
    geom_smooth(color="red") +
    scale_color_viridis_c(option = "viridis") +
    labs(y= 'NDVI ', x='Day of Year') + 
    theme_classic() +
    facet_wrap(vars(site)))





