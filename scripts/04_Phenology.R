#### ELISE GALLOIS  - elise.gallois94@gmail.com 


#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(ggplot2)
library(viridis)
library(gridExtra)
library(tidyverse)
library(plyr)
library(lubridate)
library(gridExtra)


#### 2 - LOAD DATA ####
pheno <- read.csv("output/lsat_phen_all.csv") # phenology (large file)
as_tibble(head(pheno))

#### 3 - NDVI THRESHOLD ####

# rename incorrect names
pheno$site <- recode(pheno$site, PYR_pig = 'HH_pyrpig')
pheno$site <- recode(pheno$site, REF_tem = 'BC_tem')

# calulate NDVI RATIO
thresh25 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(site,year)) %>% 
  dplyr::group_by(Year, site) %>% 
  dplyr::mutate(threshold = as.numeric(((ndvi - 0)/(max(ndvi) - 0)))) %>%  # calculate NDVI ratio (NDVImin set as 0 == snow cover)
  mutate(ratio = round(threshold, digits = 2)) %>% 
  filter(ratio == 0.25) %>% #get only the first row
  slice(1) %>%
  distinct(unique_id, .keep_all= TRUE) %>% # remove duplicate rows
  ungroup() 

thresh50 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(site,year)) %>% 
  dplyr::group_by(Year, site) %>% 
  dplyr::mutate(threshold = as.numeric(((ndvi - 0)/(max(ndvi) - 0)))) %>%  # calculate NDVI ratio (NDVImin set as 0 == snow cover)
  mutate(ratio = round(threshold, digits = 2)) %>% 
  filter(ratio == 0.5) %>% #get only the first row
  slice(1) %>%
  distinct(unique_id, .keep_all= TRUE) %>% # remove duplicate rows
  ungroup() 

thresh75 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(site,year)) %>% 
  dplyr::group_by(Year, site) %>% 
  dplyr::mutate(threshold = as.numeric(((ndvi - 0)/(max(ndvi) - 0)))) %>%  # calculate NDVI ratio (NDVImin set as 0 == snow cover)
  mutate(ratio = round(threshold, digits = 2)) %>% 
  filter(ratio == 0.75) %>% #get only the first row
  slice(1) %>%
  distinct(unique_id, .keep_all= TRUE) %>% # remove duplicate rows
  ungroup() 


# MEAN DOY GREEN-UP STATS
mean(thresh25$doy) # 25% = 181 aka 1st July
mean(thresh50$doy) # 50% = 196 aka 15th July
mean(thresh75$doy) # 75% = 211 aka 30th July


# yearly trends per site 
(greenup25 <- ggplot(thresh25) +
    geom_point(aes(x = year, y = doy, colour = type), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year, y = doy)) + 
    ylab("Green-Up DOY (25% of Maximum Curviture) \n") +
    xlab("Year") +
    scale_color_viridis_d(option = "viridis") +
    scale_fill_viridis_d(option = "viridis") +
    theme_classic())

(greenup50 <- ggplot(thresh50) +
    geom_point(aes(x = year, y = doy, colour = type), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year, y = doy)) + 
    ylab("Green-Up DOY (50% of Maximum Curviture) \n") +
    xlab("Year") +
    scale_color_viridis_d(option = "viridis") +
    scale_fill_viridis_d(option = "viridis") +
    theme_classic())

(greenup75 <- ggplot(thresh75) +
    geom_point(aes(x = year, y = doy, colour = type), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year, y = doy)) + 
    ylab("Green-Up DOY (75% of Maximum Curviture) \n") +
    xlab("Year") +
    scale_color_viridis_d(option = "viridis") +
    scale_fill_viridis_d(option = "viridis") +
    theme_classic())

(greenup <- grid.arrange(greenup25, greenup50, greenup75, ncol = 1))


