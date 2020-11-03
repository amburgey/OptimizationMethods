#### Nov 2, 2020; S. Amburgey ####
## Read in data (survey info and captures) from numerous trapping projects
## Select what projects to use
## Format for SCR analysis

rm(list=ls())

library(lubridate); library(reshape2); library(dplyr); library(sp); library(ggplot2); library(rgdal); library(raster); library(sf); library(tidyverse); library(plyr)


###### SURVEY DATA ######
allsurv <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/TRAPinfo.csv")

## Study areas deemed suitable
subsurv <- subset(allsurv, SITE %in% c("ASOI","NWFN","NWFO","HMUI","HMUR","HMU","RTSI"))

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subsurv <- subset(subsurv, PROJECTCODE %in% c("HMU TOX DROP 2 TRAP", "TOX DROP TRAP","KB TRAP","POST KB TRAP 1","POST KB TRAP 2", "NWFN TRAP 1","NWFN TRAP 2 LINVIS","NWFN TRAP 3","NWFN TRAP 4 LCM","NWFN VISTRAP TRAP","PRE BT1 TRAP","POST BT2 TRAP","NWFO TRAP 1","NWFO TRAP 2","NWFO TRAP 3","NWFO TRAP 4","GNWR SSP MOUSE/BIRD"))

## 