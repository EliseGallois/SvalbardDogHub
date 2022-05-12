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


#### 50% ####
splines50 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(sample.id,year)) %>% 
  dplyr::group_by(Year, sample.id) %>% 
  mutate(ratio = round(spl.fit, digits = 2)) %>% 
  filter(ratio == 0.38)


splines50a <- splines50 %>% 
  mutate(site_doy_id  = paste(unique_id,doy)) %>% 
  distinct(site_doy_id, .keep_all= TRUE) %>% # remove duplicate rows
  subset(ave(unique_id, unique_id, FUN = length) > 1) %>% 
  dplyr::group_by(unique_id) %>% 
  arrange(doy) %>% 
  dplyr::mutate(
    greenup = dplyr::first(doy),
    senescence = dplyr::last(doy)) 
  
  
# wide data  
thresh50 <-  pivot_longer(splines50a,cols = c("greenup", "senescence"),
               names_to = "phase",
               values_to = "phase_doy")  

# identify duplicate rows
thresh50$site_doy_id <- paste(thresh50$site_doy_id,thresh50$phase)

# remove duplicate rows
thresh50 = thresh50[!duplicated(thresh50$site_doy_id),]
   
#### 25% ####
splines25 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(sample.id,year)) %>% 
  dplyr::group_by(Year, sample.id) %>% 
  mutate(ratio = round(spl.fit, digits = 2)) %>% 
  filter(ratio == 0.19)


splines25a <- splines25 %>% 
  mutate(site_doy_id  = paste(unique_id,doy)) %>% 
  distinct(site_doy_id, .keep_all= TRUE) %>% # remove duplicate rows
  subset(ave(unique_id, unique_id, FUN = length) > 1) %>% 
  dplyr::group_by(unique_id) %>% 
  arrange(doy) %>% 
  dplyr::mutate(
    greenup = dplyr::first(doy),
    senescence = dplyr::last(doy)) 


# long data  
thresh25 <-  pivot_longer(splines25a,cols = c("greenup", "senescence"),
                          names_to = "phase",
                          values_to = "phase_doy")  

# identify duplicate rows
thresh25$site_doy_id <- paste(thresh25$site_doy_id,thresh25$phase)

# remove duplicate rows
thresh25 = thresh25[!duplicated(thresh25$site_doy_id),]

# yearly trends per site 

(greenup50 <- ggplot(thresh50) +
    geom_point(aes(x = year, y = phase_doy, colour = phase), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = phase_doy, colour = phase))) + 
    ylab("Phenophase DOY (50% of Maximum Curviture) \n") +
    xlab("Year") +
    scale_color_viridis_d(option = "viridis") +
    scale_fill_viridis_d(option = "viridis") +
    theme_classic() +
  facet_wrap(vars(type))
ggsave('figures/50_sen_ts.jpg', width = 9, height = 5, units = 'in', dpi = 400)


(greenup25 <- ggplot(thresh25) +
    geom_point(aes(x = year, y = phase_doy, colour = phase), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = phase_doy, colour = phase))) + 
  ylab("Phenophase DOY (25% of Maximum Curviture) \n") +
  xlab("Year") +
  scale_color_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis") +
  theme_classic() +
  facet_wrap(vars(type))
ggsave('figures/25_sen_ts.jpg', width = 9, height = 5, units = 'in', dpi = 400)




