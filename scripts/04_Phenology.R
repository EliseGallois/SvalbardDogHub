#### ELISE GALLOIS  - elise.gallois94@gmail.com 


#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(ggplot2)
library(viridis)
library(gridExtra)
library(tidyverse)
library(plyr)
library(lubridate)
library(lme4)
library(ggeffects)
library(sjPlot)
library(scales)
        
#### 2 - LOAD DATA ####
pheno <- read.csv("output/raw_phen_all.csv") # phenology (large file)
as_tibble(head(pheno))

#### 3 - NDVI THRESHOLD ####

# rename incorrect names
pheno$site <- recode(pheno$site, PYR_pig = 'HH_pyrpig')
pheno$site <- recode(pheno$site, REF_tem = 'BC_tem')

# remove sites not suitable for pheno analysis (no greening curves)
#pheno <- pheno %>% 
 # filter(!(site %in% c("REF_fest", "REF_lein", "REF_odin", "REF_rein")))


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

phen50 <- thresh50 %>% 
  mutate(site_doy_un  = paste(unique_id,phase_doy)) %>% 
  distinct(site_doy_un, .keep_all= TRUE) 


   
#### 25% ####
splines25 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(sample.id,year)) %>% 
  dplyr::group_by(unique_id) %>% 
  dplyr::mutate(ratio = max(ndvi)*0.25) %>% 
  mutate(ratio = round(ratio, digits = 2)) %>% 
  mutate(ndvi = round(ndvi, digits = 2)) %>%
  filter(ndvi == ratio) 


splines25a <- splines25 %>% 
  mutate(site_doy_id  = paste(unique_id,doy)) %>% 
  distinct(site_doy_id, .keep_all= TRUE) %>% # remove duplicate rows
  dplyr::group_by(unique_id) %>% 
  arrange(doy) %>% 
  dplyr::mutate(
    greenup = dplyr::first(doy),
    senescence = dplyr::last(doy)) 

# range between greenup and senescence
splines25a$range<-(splines25a$senescence-splines25a$greenup)

# remove any range <5
splines25a <- splines25a[which(splines25a$range >= 6),]


# long data  
thresh25 <-  pivot_longer(splines25a,cols = c("greenup", "senescence"),
                          names_to = "phase",
                          values_to = "phase_doy")  

# identify duplicate rows
thresh25$site_doy_id <- paste(thresh25$site_doy_id,thresh25$phase)

# remove duplicate rows
thresh25 = thresh25[!duplicated(thresh25$site_doy_id),]

##### 5 - PLOT TRENDS BY TYPE ####

# yearly trends per site 
(greenup50 <- ggplot(phen50) +
    geom_point(aes(x = year, y = phase_doy, colour = phase), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = phase_doy, fill = phase, colour = phase))) + 
    ylab("Phenophase DOY (50% of Maximum Curviture) \n") +
    xlab("Year") +
  # xlim(2004,2021) +
    scale_color_manual(values = c("#09E88F", "#E0C707")) +
   scale_fill_manual(values = c("#09E88F", "#E0C707")) +
  theme_classic() +
  facet_wrap(vars(type))
ggsave('figures/50_sen_ts.jpg', width = 9, height = 5, units = 'in', dpi = 400)


##### 6 - MODEL PHENOLOGY CHANGE ####
thresh50$year <- as.numeric(thresh50$year)
thresh50$year.scaled <- scale(I(thresh50$year - 1984), center = T)


# ndvi max doy change over time
# year as random effect??????????????
ndvi_max_m <- lmer(phase_doy ~ year.scaled*type + phase + (1|site),
                   data = thresh50)
summary(ndvi_max_m)
anova(ndvi_max_m)
tab_model(ndvi_max_m)

# Visualises random effects
(re.effects <- plot_model(ndvi_max_m, type = "re", show.values = TRUE))


# Using model predict for the tern_m model outputs and plotting the predictions
ggpredict(ndvi_max_m, terms = c("year.scaled","phase", "type"),type = "re") %>% plot()

predictions <- ggpredict(ndvi_max_m, terms = c("year.scaled", "phase","type"), type = "re") 

# rescale year
predictions$year <- rescale(predictions$x, to = c(1985, 2021))

(ggplot(predictions) + 
    geom_line(aes(x = year, y = predicted, colour = group), size = 1) + 
    geom_ribbon(data = predictions, aes(ymin = conf.low, ymax = conf.high, x = year, fill = group), 
    alpha = 0.2) + 
    geom_point(data = phen50, 
               aes(x = year, y = phase_doy, colour = phase),  alpha = 0.1) +
    scale_colour_manual(values = c("#09E88F", "#E0C707"))+
    scale_fill_manual(values = c("#09E88F", "#E0C707"))+
    labs(x = "\nYear", y = "Phenophase DOY (50% of Maximum Curviture) \n", title = "Phenology Change over Time (1985-2021)\n", colour = "Phenophase") +
    theme_classic() +
    facet_wrap(vars(type)) +
    theme(plot.title = element_text(size = 16, hjust = 0.5), # formats the title, font size 15, centres it, and moves it upwards from the graph
          plot.margin = unit(c(1,1,1,1), units = , "cm"),  # adds 1 cm margin around the figure) 
          legend.position = "right")) # positions the legend below the graph


#### NDVImax from Raw ####

#### 100% ####
splinesmax <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(sample.id,year)) %>% 
  dplyr::group_by(unique_id) %>% 
  dplyr::mutate(ratio = max(ndvi)) %>% 
  mutate(ratio = round(ratio, digits = 1)) %>% 
  mutate(ndvi = round(ndvi, digits = 1)) %>%
  filter(ndvi == ratio) 


splinesmax <- splinesmax %>% 
  mutate(site_doy_id  = paste(unique_id,doy)) %>% 
  distinct(site_doy_id, .keep_all= TRUE) %>% # remove duplicate rows
  dplyr::group_by(unique_id) %>% 
  arrange(doy) %>% 
  dplyr::mutate(
    ndvimax = dplyr::first(doy)) 



# identify duplicate rows
splinesmax$site_doy_id <- paste(splinesmax$site_doy_id)

maxi <- splinesmax %>% 
  mutate(site_doy_un  = paste(unique_id,doy)) %>% 
  distinct(site_doy_un, .keep_all= TRUE) 


(max_ndvi <- ggplot(maxi) +
    geom_point(aes(x = year, y = ndvi, colour = ndvi), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = ndvi), colour = "gold")) + 
  ylab("NDVImax \n") +
  xlab("Year") +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  facet_wrap(vars(site))
ggsave('figures/50_sen_ts.jpg', width = 9, height = 5, units = 'in', dpi = 400)





