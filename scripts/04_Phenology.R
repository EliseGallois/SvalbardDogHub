#### ELISE GALLOIS  - elise.gallois94@gmail.com 


#### 1 - LOAD PACKAGES ####
library(phenex)
library(tidyverse)
library(ggplot2)
library(viridis)
library(gridExtra)
library(tidyverse)
library(plyr)


#### 2 - LOAD DATA ####
pheno <- read.csv("output/lsat_phen_all.csv") # phenology (large file)
as_tibble(head(pheno))


pheno$ndvi <- as.numeric(as.character(pheno$ndvi))
pheno$year <- as.integer(as.character(pheno$year))

#### 3 - EXTRACT PHENOPHASES#### 
# DOY Function ----
DOY.data <- function(year, DOY, NDVI) {
  y <- cbind.data.frame(DOY, NDVI)
  names(y) <- c("DOY", "NDVI")
  if(leapYears(year[1]) == TRUE){
    x <- cbind.data.frame(seq(1, 366))
    names(x) <- c("DOY")
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
  mutate(type = substr(site, start = 1, stop = 2)) %>% 
  mutate(NDVI=ndvi) %>% 
  filter(site %in% c("DY_12", "DY_3"))


# Reformat MODIS data
pheno.data <- pheno_data %>% 
  dplyr::select(site = site, NDVI, DOY = doy, year)

pheno.data <- ddply(pheno.data,c("site","year","DOY"),
                    summarise,NDVI=mean(NDVI))

pheno.data$NDVI <- as.numeric(as.character(pheno.data$NDVI))
pheno.data$year <- as.numeric(as.character(pheno.data$year))


# Reformat MODIS data
pheno.data <- cbind.data.frame(site = as.character(pheno.data$site), 
                               year = pheno.data$year,
                               DOY = pheno.data$DOY,
                               NDVI = as.numeric(pheno.data$NDVI)) 

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

# run DOY function to create vectors for modelNDVI
greenup <- pheno.data %>% group_by(site, year) %>% 
  do(., DOY.data(.$year, .$DOY, .$NDVI))

greenup$year <- as.integer(greenup$year)
greenup$NDVI <- as.numeric(as.character(greenup$NDVI))

glimpse(greenup)

# run modelNDVI function - takes a while
greenup.all <- greenup %>% 
  # group by site and year
  group_by(site, year) %>% 
  # apply modelNDVI function
  do(., ndvi.values = unlist(modelNDVI(ndvi.values=.$NDVI, year.int=.$year, correction='none', 
                                       method="DLogistic", MARGIN=2, multipleSeasons=FALSE, doParallel=FALSE, silent=TRUE))) %>%
  # create phenoPhase dates and max NDVI
  mutate(
    greenup.date.40 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.40)$mean, 
    greenup.date.95 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.95)$mean
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



## PHENO ##
# Libraries ----
library(phenex)
library(tidyverse)
library(ggplot2)
library(viridis)
library(gridExtra)
library(tidyverse)

# DOY Function ----

DOY.data <- function(year, DOY, NDVI) {
  y <- cbind.data.frame(DOY, NDVI)
  names(y) <- c("DOY", "NDVI")
  if(leapYears(year[1]) == TRUE){
    x <- cbind.data.frame(seq(1, 366))
    names(x) <- c("DOY")
    z <- merge(x, y, all.x = TRUE, all.y = TRUE)
  } else {
    x <- cbind.data.frame(seq(1, 365))
    names(x) <- c("DOY")
    z <- merge(x, y, all.x = TRUE, all.y = TRUE)
  }
  return(z)
}

# Read data ----
MODIS <- read.csv("output/lsat_phen_all.csv")
as_tibble(head(MODIS))

MODIS <- MODIS %>% mutate(doy = as.integer(MODIS$doy))
MODIS$year <- as.numeric(as.character(MODIS$year))
MODIS$ndvi <- as.numeric(as.character(MODIS$ndvi))

MODIS <- MODIS %>% filter(!is.na(ndvi)) 

# Reformat MODIS data
MODIS.data <- cbind.data.frame(site = as.character(MODIS$site), 
                               NDVI = as.numeric(MODIS$ndvi), 
                               DOY = MODIS$doy, year = MODIS$year) 

MODIS.data <- ddply(MODIS.data,c("site","year","DOY"),
                    summarise,NDVI=mean(NDVI))

# Remove negative values
# MODIS.data <- MODIS.data %>% filter(MODIS.data$NDVI > 0)

# Add zeros
num.years <- max(MODIS.data$year) - min(MODIS.data$year) + 1
num.sites <- length(unique(MODIS.data$site))

# Phenex calculations ----

DOY <- cbind.data.frame(DOY = rep(c(rep(1, num.years), rep(32, num.years), rep(60, num.years), rep(305, num.years), rep(335, num.years)), num.sites))
NDVI <- cbind.data.frame(NDVI = rep(0, length(DOY$DOY)))
year <- cbind.data.frame(year = rep(rep(seq(min(MODIS.data$year),max(MODIS.data$year)), 5), num.sites))
site <- cbind.data.frame(site = rep(rep(as.character(unique(MODIS.data$site)), 5), num.years))
site <- arrange(site, site)
zero.data <- cbind.data.frame(site, year, DOY, NDVI)

MODIS.data <- rbind(MODIS.data, zero.data)


# run DOY function to create vectors for modelNDVI
greenup <- MODIS.data %>% group_by(site, year) %>%
  do(., DOY.data(.$year, .$DOY, .$NDVI))

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
    greenup.date.40 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.40)$mean, 
    greenup.date.60 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.)$mean, 
    greenup.date.95 = phenoPhase(ndvi.values[[1]], phase="greenup", method="local", threshold=0.95)$mean, 
    ndvi.date.max = phenoPhase(ndvi.values[[1]], phase="maxval", method="local")$mean
  )

# Save the greenup df for faster loading
save(greenup.all, file="data/greenup.qhi.kanger.RData")

# Make data.frames for plotting
greenup.greenup <- greenup.all %>% dplyr::select(site, year, greenup.date.05, greenup.date.50, greenup.date.95)
greenup.max <- greenup.all %>% dplyr::select(site, year, ndvi.date.max)
greenup.scen <- greenup.all %>% dplyr::select(site, year, senescence.date.05, senescence.date.50, senescence.date.95)
greenup.GSL <- greenup.all %>% dplyr::select(site, year, greenup.date.50, senescence.date.50, gs.length.50)
greenup.modVal <- greenup.all %>% dplyr::select(site, year, modVal) %>% group_by(site, year) %>% mutate(DOY = list(modVal[[1]][,1]), modVal = list(modVal[[1]][,2])) %>% unnest()
#greenup.modVal <-  merge(greenup, greenup.modVal)

# MODIS greenup trends - by site ----
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
  ylab("MODIS6 NDVI green up dates\n") +
  xlab("Time (years)") +
  #annotate("text", label = ("A. MODIS6 Greenup Dates"), x = 2004, y = 200, size = 5) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"))
ggsave("figures/MODIS_NDVI_phenology_plots/MODIS_greenup.png", width = 7, height = 7, units = "in", dpi = 300)


