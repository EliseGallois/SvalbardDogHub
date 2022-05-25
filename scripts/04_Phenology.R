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
library('ggthemes')

        
#### 2 - LOAD DATA ####
pheno <- read.csv("output/raw_phen_all.csv") # phenology (large file)
as_tibble(head(pheno))

#### 3 - NDVI THRESHOLD ####

# rename incorrect names
pheno$site <- recode(pheno$site, PYR_pig = 'HH_pyrpig')
pheno$site <- recode(pheno$site, REF_tem = 'BC_tem')


# plot raw pheno data
ggplot(pheno) +
  aes(x = doy, y = ndvi, colour = ndvi) +
  geom_point(size = 0.2, alpha = 0.4) +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  facet_wrap(vars(site))

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
  facet_wrap(vars(site))
ggsave('figures/50_sen_ts_site.jpg', width = 9, height = 5, units = 'in', dpi = 400)


##### 6 - MODEL PHENOLOGY CHANGE ####
thresh50$year <- as.numeric(thresh50$year)
thresh50$year.scaled <- scale(I(thresh50$year - 1984), center = 0)
thresh50$yearfactor <- as.factor(thresh50$year)

# filter for greenup
green <- thresh50 %>% 
  filter(phase %in% "greenup")

# filter for senescence
senesc <- thresh50 %>% 
  filter(phase %in% "senescence")


stab <- maxi %>% 
  filter(type %in% "Active stable")


# get finite years
green$yearno <- I(green$year - 1984)
senesc$yearno <- I(senesc$year - 1984)

# ndvi max doy change over time
greenup_m <- lmer(phase_doy ~ -1 + yearno*type +  (1|site) + (1|yearfactor),
                   data = green)
summary(greenup_m)
tab_model(greenup_m, digits = 2)

# Visualises random effects
(re.effects <- plot_model(greenup_m,  show.values = TRUE))

# ndvi max doy change over time
senesc_m <- lmer(phase_doy ~ -1 + yearno*type +  (1|site) + (1|yearfactor),
                  data = senesc)
summary(senesc_m)
tab_model(senesc_m, digits = 2)

# Visualises random effects
(re.effects <- plot_model(senesc_m,  show.values = TRUE))


#### Greenup plot ####
# Using model predict for the tern_m model outputs and plotting the predictions
ggpredict(greenup_m, terms = c("yearno", "type")) %>% plot()

# make predictions a data frame
predictions <- ggpredict(greenup_m, terms = c("yearno","type"))
predictions2 <- ggpredict(senesc_m, terms = c("yearno","type"))

# rename group to match type
predictions$type <- predictions$group
predictions2$type <- predictions$group


# rename factors
green$type <- recode_factor(green$type, BC = "Bird cliffs", 
                           DY = "Active dogyards",
                           HH = "Historic animal husbandry",
                           RE = "Reference sites",
                           ST = "Active stable")

senesc$type <- recode_factor(senesc$type, BC = "Bird cliffs", 
                            DY = "Active dogyards",
                            HH = "Historic animal husbandry",
                            RE = "Reference sites",
                            ST = "Active stable")

predictions$type <- recode_factor(predictions$type, BC = "Bird cliffs", 
                                  DY = "Active dogyards",
                                  HH = "Historic animal husbandry",
                                  RE = "Reference sites",
                                  ST = "Active stable")

predictions2$type <- recode_factor(predictions2$type, BC = "Bird cliffs", 
                                  DY = "Active dogyards",
                                  HH = "Historic animal husbandry",
                                  RE = "Reference sites",
                                  ST = "Active stable")

# Using model predict for the model outputs and plotting the predictions
ggplot(green) +
  geom_point(data = green, 
             aes(x = year, y = phase_doy,  group = type),  alpha = 0.4, colour = "#09E88F") +
  geom_line(data = predictions,aes(x = (x+1984), y = predicted, group = group), size = 1, colour = "#09E88F") +
  
  geom_ribbon(data = predictions, aes(ymin = conf.low, ymax = conf.high, x = (x+1984), group = group), 
              alpha = 0.3, fill = "#09E88F") +  
  geom_point(data = senesc, 
             aes(x = year, y = phase_doy,  group = type),  alpha = 0.4, colour = "#E0C707") +
  geom_line(data = predictions2,aes(x = (x+1984), y = predicted, group = group), size = 1, colour = "#E0C707") +
  
  geom_ribbon(data = predictions2, aes(ymin = conf.low, ymax = conf.high, x = (x+1984), group = group), 
              alpha = 0.3, fill = "#E0C707") +  
  theme_classic() +
  facet_wrap(vars(type)) +
  labs(x = "\nYear", y = "Phenophase Day of Year \n", colour = "Phenophase") +
  theme_classic() +
  facet_wrap(vars(type)) +
  theme(plot.title = element_text(size = 16, hjust = 0.5), 
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

ggsave('figures/phenology_model.jpg', width = 9, height = 5, units = 'in', dpi = 400)


#### NDVImax from Raw ####

#### 100% ####
splinesmax <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(sample.id,year)) %>% 
  dplyr::group_by(unique_id) %>% 
  dplyr::mutate(ratio = max(ndvi)) %>% 
  #mutate(ratio = round(ratio, digits = 1)) %>% 
  #mutate(ndvi = round(ndvi, digits = 1)) %>%
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
  facet_wrap(vars(site), scales='free') +
  theme_classic() +
  theme(axis.line=element_line()) +
  scale_x_continuous(limits=c(1985,2021)) + scale_y_continuous(limits=c(0,0.85))




(max_ndvi <- ggplot(maxi) +
    geom_point(aes(x = year, y = ndvi, colour = ndvi), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = ndvi), colour = "gold")) + 
  ylab("NDVImax \n") +
  xlab("Year") +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  scale_x_continuous(limits=c(1985,2021)) + scale_y_continuous(limits=c(0,0.85)) +
  facet_wrap(vars(type))

ggsave('figures/100_max_ts.jpg', width = 9, height = 5, units = 'in', dpi = 400)


# scale year variable (can be rescaled later)
maxi$year <- as.numeric(maxi$year)
maxi$year.scaled <- scale(I(maxi$year - 1984), center = T)
maxi$yearfactor <- as.factor(maxi$year)
scaling_attributes <- data.frame(predictor = c("year"),
                                 scale = c(attributes(maxi$year.scaled)$"scaled:scale"))

# scaling attribute for year:
10.81378

# ndvi max doy change over time
ndvi_max_m <- lmer(ndvi ~ -1 + year.scaled*type + (1|site) + (1|yearfactor),
                   data = maxi)

summary(ndvi_max_m)
tab_model(ndvi_max_m, digits = 6)
plot_model(ndvi_max_m)

# SLOPES BY SITE
# (baseline slope + site slope) x (number of days / scaling attribute for year)
(0.064891)*(36/10.81378) # BC:0.2160277
(0.064891+0.070355)*(36/10.81378) #DY:0.4502455 ***
(0.064891+0.051280)*(36/10.81378) #HH: 0.3867432 ***
(0.064891+0.013073)*(36/10.81378) #RE: 0.2595488 ***
(0.064891+0.034627)*(36/10.81378) #ST: 0.3313039 ***

# ERROR BY SITE
(0.014900)*(36/10.81378) # BC:0.04960338
(0.002991)*(36/10.81378) #DY:0.009957295 ***
(0.003711)*(36/10.81378) #HH: 0.01235424 ***
(0.003582)*(36/10.81378) #RE:  0.01192478 ***
(0.005769)*(36/10.81378) #ST: 0.0192055 ***

# Visualises random effects
(re.effects <- plot_model(ndvi_max_m,  show.values = TRUE))
ggpredict(ndvi_max_m, se=TRUE,interactive=TRUE)
str(maxi)


# Using model predict for the tern_m model outputs and plotting the predictions
ggpredict(ndvi_max_m, terms = c("year.scaled", "type")) %>% plot()

# make predictions a data frame
predictions <- ggpredict(ndvi_max_m, terms = c("year.scaled","type"))

# rescale back to year range
predictions$x <- as.integer(rescale(predictions$x, to = c(1985, 2021), 
                         from = range(predictions$x, na.rm = TRUE, finite = FALSE)))
# rename group to match type
predictions$type <- predictions$group


               
# rename factors
maxi$type <- recode_factor(maxi$type, BC = "Bird cliffs", 
                             DY = "Active dogyards",
                           HH = "Historic animal husbandry",
                           RE = "Reference sites",
                           ST = "Active stable")

predictions$type <- recode_factor(predictions$type, BC = "Bird cliffs", 
                           DY = "Active dogyards",
                           HH = "Historic animal husbandry",
                           RE = "Reference sites",
                           ST = "Active stable")


ggplot(maxi) +
  geom_point(data = maxi, 
             aes(x = year, y = ndvi, colour = ndvi, group = type),  alpha = 0.4) +
  scale_color_viridis_c(option = "viridis") +
  geom_line(data = predictions,aes(x = x, y = predicted,group = group), size = 1, colour = "gold") +

  geom_ribbon(data = predictions, aes(ymin = conf.low, ymax = conf.high, x = x, group = group), 
              alpha = 0.6, fill = "gold") +  
  theme_classic() +
  facet_wrap(vars(type)) +
  labs(x = "\nYear", y = "NDVImax (100% of Maximum Annual Curviture) \n", colour = "NDVImax per pixel") +
  theme_classic() +
  ylim(0,1) +
  facet_wrap(vars(type)) +
  theme(plot.title = element_text(size = 16, hjust = 0.5), 
        plot.margin = unit(c(1,1,1,1), units = , "cm"))
ggsave('figures/ndvi_max_model.jpg', width = 9, height = 5, units = 'in', dpi = 400)

fggs