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

#### 3 - EXTRACT PHENOPHASES#### 
# DOY Function ----
DOY.data <- function(year, DOY, NDVI) {
  y <- cbind.data.frame(DOY, NDVI)
  names(y) <- c("DOY", "NDVI")
  if(leapYears(year[1])==TRUE){
    x <- cbind.data.frame(seq(1, 366))
    names(x) <- c("NDVI")
    z <- merge(x, y, all.x = TRUE, all.y = TRUE)
  } else {
    x <- cbind.data.frame(seq(1, 365))
    names(x) <- c("DOY")
    z <- merge(x, y, all.x = TRUE, all.y = TRUE)
  }
  return(z)
}

# pheno data cleaning
pheno_data <- pheno %>% 
  filter(!is.na(ndvi)) %>% 
  mutate(NDVI=ndvi/1000) %>% 
  filter(year < 2022)

# Reformat MODIS data
pheno.data <- pheno_data %>% 
  dplyr::select(site = sample.id, NDVI, DOY = doy, year)

# Add zeros
num.years <- max(pheno.data$year) - min(pheno.data$year) + 1
num.sites <- length(unique(pheno.data$site))

DOY <- cbind.data.frame(DOY = rep(c(rep(1, num.years), 
                                    rep(32, num.years), 
                                    rep(60, num.years), 
                                    rep(305, num.years), 
                                    rep(335, num.years)), 
                                  num.sites))

NDVI <- cbind.data.frame(NDVI = rep(0, length(DOY$DOY)))
year <- cbind.data.frame(year = rep(rep(seq(min(pheno.data$year), max(pheno.data$year)), 5), num.sites))
site <- cbind.data.frame(site = rep(rep(as.character(unique(pheno.data$site)), 5), num.years))
site <- arrange(site, site)
zero.data <- cbind.data.frame(site, year, DOY, NDVI)

pheno.data <- rbind(pheno.data, zero.data)
pheno.data$DOY <- as.numeric(pheno.data$DOY)
pheno.data$NDVI_percent <- as.numeric(pheno.data$NDVI)

# run DOY function to create vectors for modelNDVI
greenup <- pheno.data %>% group_by(site, year) %>% 
  do(., DOY.data(.$year, .$DOY, .$NDVI_percent))

# run modelNDVI function - takes a while
greenup.all <- greenup %>% 
  # group by site and year
  group_by(site, year) %>% 
  # apply modelNDVI function
  do(., ndvi.values = unlist(modelNDVI(ndvi.values=.$NDVI, year.int=.$year, correction='none', method="DLogistic", MARGIN=2, multipleSeasons=FALSE, doParallel=FALSE, silent=TRUE))) %>%
  # Save just the DOY and modelledValues in individual columns of the tibble for plotting
  mutate(modVal = ifelse(leapYears(year[1])==TRUE, 
                         list(cbind(DOY = seq(1:366), modVal = modelledValues(ndvi.values[[1]]))), 
                         list(cbind(DOY = seq(1:365), modVal = modelledValues(ndvi.values[[1]])))
  )) %>% 
  # create phenoPhase dates and max NDVI
  mutate(
    greenup.date.05 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.05)$mean, 
    greenup.date.50 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.50)$mean, 
    greenup.date.95 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.95)$mean, 
    ndvi.date.max = phenoPhase(ndvi.values[[1]], phase="maxval", method="local")$mean, 
    senescence.date.05 = phenoPhase(ndvi.values[[1]], phase="senescence", method="local", threshold=0.05)$mean, 
    senescence.date.50 = phenoPhase(ndvi.values[[1]], phase="senescence", method="local", threshold=0.50)$mean, 
    senescence.date.95 = phenoPhase(ndvi.values[[1]], phase="senescence", method="local", threshold=0.95)$mean,
    gs.length.05 = senescence.date.05 - greenup.date.05,
    gs.length.50 = senescence.date.50 - greenup.date.50,
    gs.length.95 = senescence.date.95 - greenup.date.95,
    test = list(ndvi.values[[1]])
    #, integrateTimeserie = integrateTimeserie(ndvi.values[[1]][[1]], start=greenup.date.50, end=senescence.date.50, n=1000)
  )

save(greenup.all, file = "data/greenup.all.RData")

greenup.greenup <- greenup.all %>% dplyr::select(site, year, greenup.date.05, greenup.date.50, greenup.date.95)
greenup.max <- greenup.all %>% dplyr::select(site, year, ndvi.date.max)
greenup.scen <- greenup.all %>% dplyr::select(site, year, senescence.date.05, senescence.date.50, senescence.date.95)
greenup.GSL <- greenup.all %>% dplyr::select(site, year, greenup.date.50, senescence.date.50, gs.length.50)
greenup.modVal <- greenup.all %>% dplyr::select(site, year, modVal) %>% group_by(site, year) %>% mutate(DOY = list(modVal[[1]][,1]), modVal = list(modVal[[1]][,2])) %>% unnest()
#greenup.modVal <-  merge(greenup, greenup.modVal)

#  greenup trends - by site ----
ggplot(greenup.greenup) +
  geom_point(aes(x = year, y = greenup.date.05), colour = "#c1e3a3", alpha = 0.5, size = 2) +
  geom_point(aes(x = year, y = greenup.date.50), colour = "#64B91B", alpha = 0.5, size = 2) +
  geom_point(aes(x = year, y = greenup.date.95), colour = "#468112", alpha = 0.5, size = 2) +
  geom_line(aes(x = year, y = greenup.date.05), colour = "#c1e3a3", alpha = 0.5, size = 0.5) +
  geom_line(aes(x = year, y = greenup.date.50), colour = "#64B91B", alpha = 0.5, size = 0.5) +
  geom_line(aes(x = year, y = greenup.date.95), colour = "#468112", alpha = 0.5, size = 0.5) +
  facet_wrap( ~ site, ncol=5) +
  geom_smooth(method=lm, aes(x = year, y = greenup.date.05), colour = "#c1e3a3", fill = "#c1e3a3", alpha = 0.2, show.legend=F) +
  geom_smooth(method=lm, aes(x = year, y = greenup.date.50), colour = "#64B91B", fill = "#64B91B", alpha = 0.2, show.legend=F) +  
  geom_smooth(method=lm, aes(x = year, y = greenup.date.95), colour = "#468112", fill = "#468112", alpha = 0.2, show.legend=F) +
  ylab("Landsat NDVI green up dates\n") +
  xlab("Time (years)") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"))
ggsave("figures/landsat_greenup.png", width = 7, height = 7, units = "in", dpi = 300)





### 4 - EXPLORATORY PLOTS ####
# rename pyr_pig to HH_pyrpig, and 
ndvi_max$site <- recode(ndvi_max$site, PYR_pig = 'HH_pyrpig')
ndvi_max$site <- recode(ndvi_max$site, REF_tem = 'BC_tem')


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

#### 5 - PHENOLOGY TRENDS PLOTS ####
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

ndvi_group$year_fac <- as.factor(ndvi_group$year)

ndvi_group %>%
  filter(site %in% c("DY_4")) %>%
  ggplot() +
  aes(x = year_fac, y = ndvi.max.doy, fill = year, alpha = 0.6) +
  geom_boxplot() +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2, aes(color = year)) +
  scale_fill_hue() +
  scale_color_hue() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(vars(site)) 


#### 6 - PHENOLOGY TRENDS MODELS ####
ndvi_group$year <- as.numeric(ndvi_group$year)
ndvi_group$year.scaled <- scale(I(ndvi_group$year - 1984), center = T)


# ndvi max doy change over time
ndvi_max_m <- lmer(ndvi.max.doy ~ year.scaled*type + (1|site),
                   data = ndvi_group)
summary(ndvi_max_m)
anova(ndvi_max_m)


stargazer(ndvi_max_m, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
tab_model(ndvi_max_m)


# Visualises random effects
(re.effects <- plot_model(ndvi_max_m, type = "re", show.values = TRUE))


# Using model predict for the tern_m model outputs and plotting the predictions
ggpredict(ndvi_max_m, terms = c("year.scaled", "type"),type = "re") %>% plot()

predictions <- ggpredict(ndvi_max_m, terms = c("year.scaled", "type"), type = "re") 

# rescale year
predictions$year <- rescale(predictions$x, to = c(1985, 2021))

(ggplot(predictions) + 
   geom_line(aes(x = year, y = predicted, colour = group), size = 1) + 
   #geom_ribbon(data = predictions, aes(ymin = conf.low, ymax = conf.high, x = year, fill = group), 
               #alpha = 0.2) + 
   geom_point(data = ndvi_group, 
              aes(x = year, y = ndvi.max.doy, colour = type), alpha = 0.5) +
    scale_colour_manual(values = c("#65CF7E", "#EB9C1D", "#1BC4DE", "#D997F7", "#BD403E"))+
    scale_fill_manual(values = c("#65CF7E", "#EB9C1D", "#1BC4DE", "#D997F7", "#BD403E"))+
   labs(x = "\nYear", y = "NDVImax DOY (Landsat Surveys 30m Resolution) \n", title = "NDVImax DOY Change over Time (1985-2021)\n", colour = "Land Use Type") +
   theme_classic() +
   theme(plot.title = element_text(size = 16, hjust = 0.5), # formats the title, font size 15, centres it, and moves it upwards from the graph
         plot.margin = unit(c(1,1,1,1), units = , "cm"),  # adds 1 cm margin around the figure) 
         legend.position = "right")) # positions the legend below the graph



# ndvi max doy change over time
greening_m <- lmer(ndvi.max ~ year.scaled*type + (1|site),
                   data = ndvi_group)
summary(greening_m)
anova(greening_m)

tab_model(greening_m)


# Visualises random effects
(re.effects <- plot_model(greening_m, type = "re", show.values = TRUE))


# Using model predict for the tern_m model outputs and plotting the predictions
ggpredict(greening_m, terms = c("year.scaled", "type"),type = "re") %>% plot()

predictions <- ggpredict(greening_m, terms = c("year.scaled", "type"), type = "re") 

# rescale year
predictions$year <- rescale(predictions$x, to = c(1985, 2021))

(ggplot(predictions) + 
    geom_line(aes(x = year, y = predicted, colour = group), size = 1) + 
    #geom_ribbon(data = predictions, aes(ymin = conf.low, ymax = conf.high, x = year, fill = group), 
    #alpha = 0.2) + 
    geom_point(data = ndvi_group, 
               aes(x = year, y = ndvi.max, colour = type), alpha = 0.5) +
    scale_colour_manual(values = c("#65CF7E", "#EB9C1D", "#1BC4DE", "#D997F7", "#BD403E"))+
    scale_fill_manual(values = c("#65CF7E", "#EB9C1D", "#1BC4DE", "#D997F7", "#BD403E"))+
    labs(x = "\nYear", y = "NDVImax  (Landsat Surveys 30m Resolution) \n", title = "NDVImax Change over Time (1985-2021)\n", colour = "Land Use Type") +
    theme_classic() +
    theme(plot.title = element_text(size = 16, hjust = 0.5), # formats the title, font size 15, centres it, and moves it upwards from the graph
          plot.margin = unit(c(1,1,1,1), units = , "cm"),  # adds 1 cm margin around the figure) 
          legend.position = "right")) # positions the legend below the graph




