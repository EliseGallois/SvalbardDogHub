#### ELISE GALLOIS  - elise.gallois94@gmail.com 
# This script extracts greenup, senescence and peak NDVI dates and values
# This script contains statistical models to analyse the trends over time of these extracted values

#### 1 - LOAD PACKAGES ####
library(tidyverse)
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


# plot raw pheno data (get sense of 'greening curves')
ggplot(pheno) +
  aes(x = doy, y = ndvi, colour = ndvi) +
  geom_point(size = 0.2, alpha = 0.4) +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  facet_wrap(vars(site))


#### 4 - EXTRACT GREENUP AND SENSESCENCE DATES ####
splines50 <- pheno %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>%    # get land use type classification 
  mutate(Year = lubridate::ymd(year, truncated = 2L)) %>%  # set year to date type
  mutate(unique_id  = paste(sample.id,year)) %>% 
  dplyr::group_by(unique_id) %>% 
  dplyr::mutate(ratio = max(ndvi)*0.5) %>% 
  mutate(ratio = round(ratio, digits = 1)) %>% 
  mutate(ndvi = round(ndvi, digits = 1)) %>%
  filter(ndvi == ratio) 

# extract senescence and greenup dates
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

# keep only distinct values
maxi <- splinesmax %>% 
  mutate(site_doy_un  = paste(unique_id,doy)) %>% 
  distinct(site_doy_un, .keep_all= TRUE) 

# get NDVImax trends over time for all sites
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
ggsave('figures/all_ndvi_max.jpg', width = 9, height = 9, units = 'in', dpi = 400)


# filter above for case study plots
case_study <- maxi %>% 
  filter(site %in% c("DY_12", "HH_lyr","HH_dogs", "BC_alkhor","REF_fest")) 

case_study$site <- recode(case_study$site, DY_12 = 'Active Dog yard')
case_study$site <- recode(case_study$site, HH_dogs = 'Abandoned Dog yard')
case_study$site <- recode(case_study$site, BC_alkhor = 'Bird Cliff')
case_study$site <- recode(case_study$site, REF_fest = 'Natural tundra')
case_study$site <- recode(case_study$site, HH_lyr = 'Abandoned Livestock')


(case_study_ndvi <- ggplot(case_study) +
    geom_point(aes(x = year, y = ndvi, colour = ndvi), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = ndvi), colour = "gold")) + 
  ylab("NDVImax \n") +
  xlab("Year") +
  scale_color_viridis_c(option = "viridis") +
  facet_wrap(vars(site), scales='free') +
  theme_classic() +
  theme(axis.line=element_line()) +
  scale_x_continuous(limits=c(1985,2021)) + scale_y_continuous(limits=c(0,0.85))
ggsave('figures/case_study_ndvi_max.jpg', width = 10, height = 5, units = 'in', dpi = 1000)

# individual plots
abandoned_dy <- case_study %>% 
  filter(site == "Abandoned Dog yard") 

(ab_dy_ndvi <- ggplot(abandoned_dy) +
    geom_point(aes(x = year, y = ndvi, colour = ndvi), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = ndvi), colour = "gold")) + 
  ylab("NDVImax \n") +
  xlab("Year") +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  theme(axis.line=element_line()) +
  scale_x_continuous(limits=c(1985,2021)) + scale_y_continuous(limits=c(0,0.85))
ggsave('figures/abandoned_dy_ndvi_max.jpg', width = 6, height = 3, units = 'in', dpi = 1000)

active_dy <- case_study %>% 
  filter(site == "Active Dog yard") 

(ac_dy_ndvi <- ggplot(active_dy) +
    geom_point(aes(x = year, y = ndvi, colour = ndvi), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = ndvi), colour = "gold")) + 
  ylab("NDVImax \n") +
  xlab("Year") +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  theme(axis.line=element_line()) +
  scale_x_continuous(limits=c(1985,2021)) + scale_y_continuous(limits=c(0,0.85))
ggsave('figures/active_dy_ndvi_max.jpg', width = 6, height = 3, units = 'in', dpi = 1000)

bc <- case_study %>% 
  filter(site == "Bird Cliff") 

(ac_dy_ndvi <- ggplot(bc) +
    geom_point(aes(x = year, y = ndvi, colour = ndvi), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = ndvi), colour = "gold")) + 
  ylab("NDVImax \n") +
  xlab("Year") +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  theme(axis.line=element_line()) +
  scale_x_continuous(limits=c(1985,2021)) + scale_y_continuous(limits=c(0,0.85))
ggsave('figures/bird_ndvi_max.jpg', width = 6, height = 3, units = 'in', dpi = 1000)


natural <- case_study %>% 
  filter(site == "Natural tundra") 

(nat_ndvi <- ggplot(natural) +
    geom_point(aes(x = year, y = ndvi, colour = ndvi), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = ndvi), colour = "gold")) + 
  ylab("NDVImax \n") +
  xlab("Year") +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  theme(axis.line=element_line()) +
  scale_x_continuous(limits=c(1985,2021)) + scale_y_continuous(limits=c(0,0.85))
ggsave('figures/natural_tundra_ndvi_max.jpg', width = 6, height = 3, units = 'in', dpi = 1000)

ab_livestock <- case_study %>% 
  filter(site == "Abandoned Livestock") 

(ab_livestock <- ggplot(ab_livestock) +
    geom_point(aes(x = year, y = ndvi, colour = ndvi), alpha = 0.5, size = 2) +
    geom_smooth(method=lm, aes(x = year,  y = ndvi), colour = "gold")) + 
  ylab("NDVImax \n") +
  xlab("Year") +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  theme(axis.line=element_line()) +
  scale_x_continuous(limits=c(1985,2021)) + scale_y_continuous(limits=c(0,0.85))
ggsave('figures/abandoned_livestock_ndvi_max.jpg', width = 6, height = 3, units = 'in', dpi = 1000)

# get NDVImax trends over time grouped by land use types
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
maxi$year.scaled <- scale(I(maxi$year - 1984), center = T) # scale so model converges
maxi$yearfactor <- as.factor(maxi$year) # make a factor version of year for use as random effect

# ndvi max doy change over time
ndvi_max_m <- lmer(ndvi ~ -1 + year.scaled*type + (1|site) + (1|yearfactor),
                   data = maxi)

summary(ndvi_max_m)
tab_model(ndvi_max_m, digits = 6)
plot_model(ndvi_max_m)

# CALCULATE SLOPES BY SITE
# (baseline slope + site slope) x (number of days / scaling attribute for year)
(0.064891)*(36/10.81378) # BC:0.2160277
(0.064891+0.070355)*(36/10.81378) #DY:0.4502455 ***
(0.064891+0.051280)*(36/10.81378) #HH: 0.3867432 ***
(0.064891+0.013073)*(36/10.81378) #RE: 0.2595488 ***
(0.064891+0.034627)*(36/10.81378) #ST: 0.3313039 ***

# ERROR BY SITE
(0.014900)*(36/10.81378) # BC:0.04960338
(0.002991)*(36/10.81378) #DY: 0.009957295 ***
(0.003711)*(36/10.81378) #HH: 0.01235424 ***
(0.003582)*(36/10.81378) #RE:  0.01192478 ***
(0.005769)*(36/10.81378) #ST: 0.0192055 ***

# Visualises random effects
(re.effects <- plot_model(ndvi_max_m,  show.values = TRUE))

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

