##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area

## CP was created in 2004 and has been used in several projects, resulting in a rich time series with surveys occurring at various densities of snakes

rm(list=ls())

source("Select&PrepTrapData.R")   ## Creation of subcap and subsurv (cleaned up)
source("Trapping/DataPrep/DataPrepCP.R")              ## Functions to reshape survey and capture data

library(secr); library(reshape2); library(jagsUI)

## Subset capture data (subcap) and effort/survey data (subsurv)
CPcaps <- subset(subcap, SITE == "NWFN")
CPsurv <- subset(subsurv, SITE == "NWFN")

## Subset to specific NWFN project
CPcaps <- subset(CPcaps, PROJECTCODE == "NWFN TRAP 1")
CPsurv <- subset(CPsurv, PROJECTCODE == "NWFN TRAP 1")

## SECIFY TIME FRAME
time <- c("06","07")
time2 <- c("2004-05-01","2004-08-31")


##### SPECIFY DIMENSIONS OF CP #####
## Make study area grid to ensure correct size
## 27 transects (VIS and TRAP) with 13 points each, 8m from VIS to TRAP and 16m between points
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
J <- nrow(locs)

## Define state-space of point process. (i.e., where animals live).
## Don't need to estimate state-space since we know it (5 ha/50000 m2 enclosed pop) but do need this to help make integration grid below
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Area of CP
A <- (Xu-Xl)*(Yu-Yl)


##### USE CATEGORICAL GRID CELL LOCATIONS #####
## Surveys locations
X <- as.matrix(locs)

#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected (V = visual, T = trap)
capPROJ <- subSnk(SITEcaps=CPcaps, type=c("TRAPTYPE"), info=c("M"))
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=time)  ## this is using 2 months (Feb - Mar)
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=CPsurv, time=time)
## Check data to make sure no missing effort or captured snakes were on survey dates (throws error if dim mismatch)
checkDims(SCReff, SCRcaps)