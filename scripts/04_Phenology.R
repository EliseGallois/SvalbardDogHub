#### ELISE GALLOIS  - elise.gallois94@gmail.com 


#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(ggplot2)
library(viridis)
library(gridExtra)
library(tidyverse)
library(plyr)


#### 2 - LOAD DATA ####
pheno <- read.csv("output/lsat_phen_all.csv") # phenology (large file)
as_tibble(head(pheno))

#### 3 - NDVI THRESHOLD ####


# calulate NDVI RATIO
thresh <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>% 
  group_by(site, year) %>% 
  mutate(threshold = as.numeric(((ndvi - (min(ndvi)))/(max(ndvi)-min(ndvi)))))
  
thresh <- thresh %>% 
  group_by(site, year) %>% 
  filter(threshold >= 0.24 & threshold <= 0.26) %>% 
  slice(1) 



# yearly trends per site
(annual.plot.greenup <- ggplot(thresh) +
    geom_point(aes(x = year, y = doy, colour = site), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year, y = doy)) + 
    ylab("Green-Up DOY (25% of Maximum Curviture) \n") +
    xlab("Year") +
    scale_color_viridis_d(option = "viridis") +
    scale_fill_viridis_d(option = "viridis") +
    theme_classic())

