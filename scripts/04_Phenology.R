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

# remove sites not suitable for pheno analysis (no greening curves)
pheno <- pheno %>% 
  filter(!(site %in% c("REF_fest", "REF_lein", "REF_odin", "REF_rein")))


#### 50% ####
splines50 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(sample.id,year)) %>% 
  dplyr::group_by(unique_id) %>% 
  dplyr::mutate(ratio = max(ndvi)*0.5) %>% 
  mutate(ratio = round(ratio, digits = 1)) %>% 
  mutate(ndvi = round(ndvi, digits = 1)) %>%
  filter(ndvi == ratio) 


splines50a <- splines50 %>% 
  mutate(site_doy_id  = paste(unique_id,doy)) %>% 
  distinct(site_doy_id, .keep_all= TRUE) %>% # remove duplicate rows
  dplyr::group_by(unique_id) %>% 
  arrange(doy) %>% 
  dplyr::mutate(
    greenup = dplyr::first(doy),
    senescence = dplyr::last(doy)) 
  
# range between greenup and senescence
splines50a$range<-(splines50a$senescence-splines50a$greenup)

# remove any range <5
splines50a <- splines50a[which(splines50a$range >= 6),]

  
# long data  
thresh50 <-  pivot_longer(splines50a,cols = c("greenup", "senescence"),
               names_to = "phase",
               values_to = "phase_doy")  

# identify duplicate rows
thresh50$site_doy_id <- paste(thresh50$site_doy_id,thresh50$phase)

# remove duplicate rows
thresh50 = thresh50[!duplicated(thresh50$site_doy_id),]
   
#### 30% ####
splines30 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(sample.id,year)) %>% 
  dplyr::group_by(unique_id) %>% 
  dplyr::mutate(ratio = max(ndvi)*0.3) %>% 
  mutate(ratio = round(ratio, digits = 1)) %>% 
  mutate(ndvi = round(ndvi, digits = 1)) %>%
  filter(ndvi == ratio) 

splines30a <- splines30 %>% 
  mutate(site_doy_id  = paste(unique_id,doy)) %>% 
  distinct(site_doy_id, .keep_all= TRUE) %>% # remove duplicate rows
  dplyr::group_by(unique_id) %>% 
  arrange(doy) %>% 
  dplyr::mutate(
    greenup = dplyr::first(doy),
    senescence = dplyr::last(doy)) 

# range between greenup and senescence
splines30a$range<-(splines30a$senescence-splines30a$greenup)

# remove any range <5
splines30a <- splines30a[which(splines30a$range >= 6),]


# long data  
thresh30 <-  pivot_longer(splines30a,cols = c("greenup", "senescence"),
                          names_to = "phase",
                          values_to = "phase_doy")  

# identify duplicate rows
thresh30$site_doy_id <- paste(thresh30$site_doy_id,thresh30$phase)

# remove duplicate rows
thresh30 = thresh30[!duplicated(thresh30$site_doy_id),]


##### 5 - PLOT TRENDS BY TYPE ####

# yearly trends per site 
(greenup50 <- ggplot(thresh50) +
    geom_point(aes(x = year, y = phase_doy, colour = phase), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = phase_doy, fill = phase, colour = phase))) + 
    ylab("Phenophase DOY (50% of Maximum Curviture) \n") +
    xlab("Year") +
    scale_color_manual(values = c("#09E88F", "#E0C707")) +
   scale_fill_manual(values = c("#09E88F", "#E0C707")) +
  theme_classic() +
  facet_wrap(vars(type))
ggsave('figures/50_sen_ts.jpg', width = 9, height = 5, units = 'in', dpi = 400)



(greenup30 <- ggplot(thresh30) +
    geom_point(aes(x = year, y = phase_doy, colour = phase), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = phase_doy, colour = phase))) + 
  ylab("Phenophase DOY (25% of Maximum Curviture) \n") +
  xlab("Year") +
  scale_color_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis") +
  theme_classic() +
  facet_wrap(vars(type))
ggsave('figures/25_sen_ts.jpg', width = 9, height = 5, units = 'in', dpi = 400)

(phen <- grid.arrange(greenup25, greenup50, ncol = 1))


