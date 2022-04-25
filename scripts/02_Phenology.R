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
  labs(y= 'DOY of NDVImax ', x='Year') + 
  theme_classic() +
  facet_wrap(vars(site))

# plot greening curves for all sites
pheno$year <- as.factor(pheno$year)
                        
(curve_yard <- ggplot(pheno) +
    geom_point(aes(x=doy, y=ndvi), 
                alpha = 0.2, size = 0.1,colour = "grey") +
    geom_smooth(aes(x=doy, y=ndvi, colour = year, group = sample.id), 
                span = 0.2, alpha = 0.5, size = 0.5, se = F, show.legend=T) +
    scale_color_viridis_d(option = "magma") +
    labs(y= 'NDVI ', x='Day of Year') + 
    theme_classic()) +
    facet_wrap(vars(site))
  
  
ggsave('figures/all_sites_pheno_trend.jpg', width = 9, height = 9, units = 'in', dpi = 400)





