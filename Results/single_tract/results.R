library(Hmisc)
library(corrplot)
library(tidyverse)


setwd("C:/SANDEEP/SocioHydro_Model/Results/single_tract")

para <- read.csv("SH_Parameters.csv")
para <- para %>% dplyr::select(-census_tract)
SH <- read.csv("SH_Results.csv")
SH <- SH %>% dplyr::select(-census_tract)
rmse <- read.csv("SH_RMSE.csv")
#Barplot of all parameters
hist.data.frame(para, nclass = 10)

#Correlation matrix of all parameters
corrplot(cor(para), method = "number")

#Plot of SH Results



