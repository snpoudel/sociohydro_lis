setwd("C:/SANDEEP/SocioHydro_Model/Data/Surge.threshold_vs_mean.elev")
library(tidyverse)

#input files
f1 <- read.csv("mean_elev.csv")
f2 <- read.csv("SH_Parameters.csv")
f2 <- f2 %>% select(census_tract, Surge_threshold) %>% arrange(census_tract)

#make both files have same census tracts by filtering
f1 <- f1 %>% filter(FIPS %in% f2$census_tract) %>% arrange(FIPS) 

#merge to have surge threshold and mean elevation in same file
f <- data.frame(census_tract = f2$census_tract,surge_threshold = f2$Surge_threshold,mean_elevation = f1$Mean.Elevation)

#plot
theme_set(theme_bw())
f %>% ggplot(aes(x= mean_elevation, y = surge_threshold))+
    geom_smooth(method = lm, se =T)+
    geom_point(size = 0.5)
