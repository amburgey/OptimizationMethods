##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area

## CP was created in 2004 and has been used in several projects, resulting in a rich time series with surveys occurring at various densities of snakes

rm(list=ls())

source("Select&PrepTrapData.R")   ## Creation of subcap and subsurv (cleaned up)
source("Trapping/DataPrep/DataPrepCP.R")    ## Functions to reshape survey and capture data
source("Trapping/DataPrep/OverlayCPGrid.R")

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
cellsize <- c(5,5)  ## dimensions of integration grid cell
CPspecs <- overlayCP(CPcaps, cellsize)  ## ignore warnings, all about projections
## Area (5 ha/50,000 m2): 
A <- 50000


##### USE CATEGORICAL GRID CELL LOCATIONS #####
## Surveys locations
X <- as.matrix(CPspecs$tran[,-1])[,2:3]
J <- nrow(X)

#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected (V = visual, T = trap)
capPROJ <- subSnk(SITEcaps=CPcaps, type=c("TRAPTYPE"), info=c("M"))
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=time)  ## this is using 2 months (Feb - Mar)
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=CPsurv, time=time)
## Check data to make sure no missing effort or captured snakes were on survey dates (throws error if dim mismatch)
checkDims(SCReff, SCRcaps)


## Inits for activity centers, can't take mean grid cell location where each snake was found as snakes found all over CP
## Instead take first cell location where captured for each individual
locs <- CPspecs$tran
colnames(locs)[2] <- c("CellID")
vsst <- list()
for(i in 1:nrow(dat$y)){
  vsst[i] <- apply(dat$y,1,function(x) which(x==1))[[i]][1]
  vsst <- unlist(vsst)
}