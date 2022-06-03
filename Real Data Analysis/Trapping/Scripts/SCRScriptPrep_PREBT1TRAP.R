##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area
##### PREBT1 - 2008 project
##### The purpose of this code is to prepare this dataset for spatial capture recapture analysis

source("Real Data Analysis/Trapping/DataPrep/DataPrepCPPREBT1TRAP.R")    ## Functions to reshape survey and capture data

## Subset capture data (subcap) and effort/survey data (subsurv)
CPcaps <- subset(subcap, SITE == "NWFN")
CPsurv <- subset(subsurv, SITE == "NWFN")

## Subset to specific NWFN project
CPcaps <- subset(CPcaps, PROJECTCODE == "PRE BT1 TRAP")
CPsurv <- subset(CPsurv, PROJECTCODE == "PRE BT1 TRAP")

## SECIFY TIME FRAME
time <- c("07","08")
time2 <- c("2008-06-01","2008-09-30")


##### SPECIFY DIMENSIONS OF CP #####
## Done for POSTBT2 so use for all projects

#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected (V = visual, M = trap)
capPROJ <- subSnk(SITEcaps=CPcaps, type=c("TRAPTYPE"), info=c("M"))
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=time)  ## this is using 2 months (Feb - Mar)
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=CPsurv, time=time)
## Check no duplicates surveys being retained, no erroneous snake captures
SCRcaps <- checkSnks(SCRcaps=SCRcaps)
## Check data to make sure no missing effort or captured snakes were on survey dates (throws error if dim mismatch)
checkDims(SCReff, SCRcaps)

##### USE CATEGORICAL GRID CELL LOCATIONS #####
## Surveys locations
slocs <- paste(rep(unique(SCReff$TRANSECT),each=13),1:13,sep="")
fullX <- subset(CPspecs$tran, TranID %in% slocs)
X <- as.matrix(fullX[,-1])[,2:3]
J <- nrow(X)

#### FORMAT DATA FOR TRADITIONAL SCR ANALYSIS ####
## Add GridID to captures so sorting using that instead of location
colnames(fullX)[1] <- c("Point")
SCRcaps <- merge(SCRcaps, fullX[,1:2], by = c("Point"))
## If error and need to do manual
dat <- prepSCRman(SCRcaps, SCReff, grid = fullX)

## Observations, already in order of CellID locations
y <- dat$y

## Uniquely marked individuals
nind <- nrow(y)

## Get sizes of individuals
# snsz <- getSize(capPROJ, SCRcaps, subcap)[,2]  ## if all snakes have a measurement during that project
snsz <- getSizeman(capPROJ, SCRcaps, subcap, time=time2)[,2] ## if some snake sizes are missing than expand window of time
## Categorize by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
snsz <- ifelse(snsz < 850, 1,
               ifelse(snsz >= 850 & snsz < 950, 2,
                      ifelse(snsz >= 950 & snsz < 1150, 3,
                             ifelse(snsz >= 1150, 4, -9999))))
if(min(snsz) == -9999) stop('snake size incorrect')
L <- length(unique(snsz))
ngroup <- as.vector(table(snsz))

## Active/not active for when transects run, already in order of 1-J CellID locations
act <- as.matrix(dat$act[,-1])
colnames(act) <- NULL
K <- rowSums(act)

## Number of survey occasions
nocc <- ncol(act)

## Inits for activity centers, can't take mean grid cell location where each snake was found as snakes found all over CP
## Instead take first cell location where captured for each individual
locs <- subset(CPspecs$tran, TranID %in% slocs)
vsst <- list()
for(i in 1:nrow(dat$y)){
  vsst[i] <- apply(dat$y,1,function(x) which(x==1))[[i]][1]
  vsst <- unlist(vsst)
}
