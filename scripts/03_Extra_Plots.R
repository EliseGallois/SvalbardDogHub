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
library(lme4)
library(sjPlot)
library(stargazer)
library(report)
library(scales)
library(phenex)


#### 2 - LOAD DATA ####
pheno <- read_csv("output/lsat_phen_all.csv") # phenology (large file)
ndvi_max <- read_csv("output/lsat_NDVImx_all.csv") # ndvimax only (streamlined file)


### 3 - NDVI MAX EXPLORATORY PLOTS ####
# rename pyr_pig to HH_pyrpig, and 
ndvi_max$site <- recode(ndvi_max$site, PYR_pig = 'HH_pyrpig')
ndvi_max$site <- recode(ndvi_max$site, REF_tem = 'BC_tem')
pheno$site <- recode(pheno$site, REF_tem = 'BC_tem')


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
  
ggsave('figures/extra/all_sites_pheno_trend_splines.jpg', width = 9, height = 9, units = 'in', dpi = 400)


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
ggsave('figures/extra/all_sites_ndvimax_annual_splines.jpg', width = 9, height = 9, units = 'in', dpi = 400)

#### 5 - heterogeneity ####
ggplot(ndvi_group) +
  aes(x = site, y = ndvi.max, fill = type) +
  geom_boxplot(adjust = 1L, scale = "area") +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.9), alpha = 0.1) +
  scale_fill_hue() +
  theme_minimal()
 


