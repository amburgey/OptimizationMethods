##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area
##### VIS 2 - 2006 project

source("Real Data Analysis/Visual surveys/DataPrep/DataPrepCP_VIS2.R")   ## Functions to reshape survey and capture data

  
## Subset capture data (subcap) and effort/survey data (subsurv)
CPcaps <- subset(subcap, SITE == "NWFN")
CPsurv <- subset(subsurv, SITE == "NWFN")

## Subset to specific NWFN project
CPcaps <- subset(CPcaps, PROJECTCODE == "NWFN VIS 2")
CPsurv <- subset(CPsurv, PROJECTCODE == "NWFN VIS 2")

## SECIFY TIME FRAME
time <- c("02","03")
time2 <- c("2006-01-01","2006-05-31")


##### SPECIFY DIMENSIONS OF CP - DO HERE AND USE FOR ALL STUDIES IN CP #####
cellsize <- c(10,10)  ## dimensions of integration grid cell
CPspecs <- overlayCP(CPcaps, cellsize)  ## ignore warning about Z-dimension
## Area (5 ha/50,000 m2): 
A <- sum(CPspecs$area)


##### USE CATEGORICAL GRID CELL LOCATIONS #####
## Surveys locations
fullX <- CPspecs$tran
X <- as.matrix(CPspecs$tran[,-1])[,2:3]
J <- nrow(X)

#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected (V = visual, T = trap)
capPROJ <- subSnk(SITEcaps=CPcaps, type=c("TRAPTYPE"), info=c("V"))
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=time)  ## this is using 2 months
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=CPsurv, time=time)
## Check no duplicates surveys being retained
SCRcaps <- checkSnks(SCRcaps=SCRcaps)
## Check data to make sure no missing effort or captured snakes were on survey dates (throws error if dim mismatch)
checkDims(SCReff, SCRcaps)

#### FORMAT DATA FOR TRADITIONAL SCR ANALYSIS ####
## Add GridID to captures so sorting using that instead of location
colnames(fullX)[1] <- c("Point")
SCRcaps <- merge(SCRcaps, fullX[,1:2], by = c("Point"))
dat <- prepSCR(SCRcaps, SCReff, grid = fullX)

## Observations, already in order of CellID locations
y <- dat$y

## Uniquely marked individuals
nind <- nrow(y)

## Get sizes of individuals
snsz <- getSizeman(capPROJ, SCRcaps, subcap, time=time2)[,2] ## if some snake sizes are missing than expand window of time
## Categorize by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
snsz <- ifelse(snsz < 850, 1,
               ifelse(snsz >= 850 & snsz < 950, 2,
                      ifelse(snsz >= 950 & snsz < 1150, 3,
                             ifelse(snsz >= 1150, 4, -9999))))
if(max(snsz) == -9999) stop('snake size incorrect')
L <- length(unique(snsz))
ngroup <- as.vector(table(snsz))

## Active/not active for when transects run, already in order of 1-351 CellID locations
act <- as.matrix(dat$act[,-1])
colnames(act) <- NULL
K <- rowSums(act)

## Number of survey occasions
nocc <- ncol(act)

## Inits for activity centers, can't take mean grid cell location where each snake was found as snakes found all over CP
## Instead take first cell location where captured for each individual
locs <- CPspecs$tran
colnames(locs)[2] <- c("CellID")
vsst <- list()
for(i in 1:nrow(dat$y)){
  vsst[i] <- apply(dat$y,1,function(x) which(x>=1))[[i]][1]
  vsst <- unlist(vsst)
}
