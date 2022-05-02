#### ELISE GALLOIS  - elise.gallois94@gmail.com 


#### 1 - LOAD PACKAGES ####
library(tidyverse)
library(esquisse)
library(sf)
library(raster)
library(viridis)
library(readxl)
library(readr)
library(ggrepel)
library(plyr)


#### 2 - LOAD DATA ####
pheno <- read_csv("output/lsat_phen_all.csv") # phenology (large file)
ndvi_max <- read_csv("output/lsat_NDVImx_all.csv") # ndvimax only (streamlined file)


### 3 - EXPLORATORY PLOTS ####
# rename pyr_pig to HH_pyrpig
ndvi_max$site <- recode(ndvi_max$site, PYR_pig = 'HH_pyrpig')
ndvi_max$site <- recode(ndvi_max$site, PYR_pig = 'HH_pyrpig')

# subset into types
ndvi_group <- ndvi_max %>% 
  mutate(type = substr(site, start = 1, stop = 2)) %>% 
  mutate(quantile = ntile(ndvi.max, 3)) %>% 
  mutate(quantile = as.factor(quantile)) %>% 
  filter(ndvi.max.doy >= 121 & ndvi.max.doy <= 244) # between START OF MAY and END OF AUGUST


# calculate quantiles of NDVImax
ndvi_group$quantile <- recode(ndvi_group$quantile, '1' = 'Low NDVImax')
ndvi_group$quantile <- recode(ndvi_group$quantile, '2' = 'Mid NDVImax')
ndvi_group$quantile <- recode(ndvi_group$quantile, '3' = 'High NDVImax')


# plot ndvi max change by quantile
ggplot(ndvi_group) +
  aes(x = year, y = ndvi.max, colour = quantile) +
  geom_point(size = 1L) +
  geom_smooth(method=lm,  aes(colour = quantile)) +
  scale_color_viridis_d(option = "viridis") +
  labs(y= 'Annual NDVImax ', x='Year', title = "Annual NDVImax Change by Quantile") + 
  theme_classic() 


# plot ndvi max change by land use type
ggplot(ndvi_group) +
  aes(x = year, y = ndvi.max, colour = type) +
  geom_point(size = 1L) +
  geom_smooth(method=lm,  aes(colour = type)) +
  scale_colour_manual(values = c("#65CF7E", "#EB9C1D", "#1BC4DE", "#D997F7", "#BD403E"))+
  labs(y= 'Annual NDVImax ', x='Year', title = "Annual NDVImax Change by Land Use Type") + 
  theme_classic() 



# plot all 'date of 'NDVImax'
ggplot(ndvi_group) +
  aes(x = year, y = ndvi.max.doy, colour = ndvi.max) +
  geom_point(size = 1L) +
  geom_smooth(method=lm, color="gold") +
  scale_color_viridis_c(option = "viridis") +
  labs(y= 'DOY of NDVImax ', x='Year') + 
  theme_classic() 


ggplot(ndvi_group) +
  aes(x = year, y = ndvi.max, colour = ndvi.max) +
  geom_point(size = 1L) +
  geom_smooth(method=lm, color="gold") +
  scale_color_viridis_c(option = "viridis") +
  labs(y= 'Annual NDVImax ', x='Year') + 
  theme_classic()


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


# get mean ndvi.max doy for each plot
lsat.gs.dt2 <- ndvi_group %>%
  group_by(site) %>%
  summarise(Mean.nd = median(ndvi.max.doy))

(hist <- ggplot(ndvi_group) +
    aes(x = ndvi.max.doy) +
    geom_histogram(bins = 30L, fill = "#1f9e89") +
    geom_vline(data = lsat.gs.dt2, mapping = aes(xintercept = Mean.nd), colour = "blue", linetype = "dashed") +
    theme_classic() +
    facet_wrap(vars(site)))

#### 4 - PHENOLOGY TRENDS PLOTS ####
# get mean ndvi.max  for each site and year
pheno_av <- ndvi_group %>%
  group_by(site, year) %>%
  summarise_at(vars(ndvi.max), list(name = mean_ndvixmax))
yearly.means <- ddply(ndvi_group,c("site","year"),summarise,mean=mean(ndvi.max))

# yearly trends per site
(annual.plot <- ggplot(yearly.means) +
  geom_point(aes(x = year, y = mean, colour = site), alpha = 0.5, size = 2) +
  geom_line(aes(x = year, y = mean, colour = site), alpha = 0.5, size = 0.5) +
  geom_smooth(method=lm, aes(x = year, y = mean, colour = site, fill = site), alpha = 0.2, show.legend=F) + 
  ylab("Annual Max. NDVI\n") +
  xlab("Year") +
    scale_color_viridis_d(option = "viridis") +
    scale_fill_viridis_d(option = "viridis") +
  theme_classic() +
    facet_wrap(vars(site)))
ggsave('figures/all_sites_ndvimax_annual.jpg', width = 9, height = 9, units = 'in', dpi = 400)

# get mean ndvi.max DOY  for each site and year

pheno.yearly.means <- ddply(ndvi_group,c("site","year"),summarise,mean=mean(ndvi.max.doy))

# yearly trends per site
(annual.plot.pheno <- ggplot(pheno.yearly.means) +
    geom_point(aes(x = year, y = mean, colour = site), alpha = 0.5, size = 2) +
    geom_line(aes(x = year, y = mean, colour = site), alpha = 0.5, size = 0.5) +
    geom_smooth(method=lm, aes(x = year, y = mean, colour = site, fill = site), alpha = 0.2, show.legend=F) + 
    ylab("Annual Max. NDVI DOY\n") +
    xlab("Year") +
    scale_color_viridis_d(option = "viridis") +
    scale_fill_viridis_d(option = "viridis") +
    theme_classic() +
    facet_wrap(vars(site)))
ggsave('figures/all_sites_ndvidoy_annual.jpg', width = 9, height = 9, units = 'in', dpi = 400)

ndvi_group$year <- as.factor(ndvi_group$year)

ndvi_group %>%
  filter(site %in% c("DY_4")) %>%
  ggplot() +
  aes(x = year, y = ndvi.max.doy, fill = year, alpha = 0.6) +
  geom_boxplot() +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2, aes(color = year)) +
  scale_fill_hue() +
  scale_color_hue() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(vars(site)) 


#### 5 - PHENOLOGY TRENDS MODELS ####


