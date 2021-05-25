##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area

## CP was created in 2004 and has been used in several projects, resulting in a rich time series with surveys occurring at various densities of snakes

rm(list=ls())

source("Select&PrepTrapData.R")   ## Creation of subcap and subsurv (cleaned up)
source("Trapping/DataPrep/DataPrepCP.R")    ## Functions to reshape survey and capture data
source("Trapping/DataPrep/OverlayCPGrid.R")

library(nimble); library(reshape2); library(jagsUI)

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
cellsize <- c(10,10)  ## dimensions of integration grid cell
CPspecs <- overlayCP(CPcaps, cellsize)  ## ignore warnings, all about projections
## Area (5 ha/50,000 m2): 
A <- 50000


#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected (V = visual, T = trap)
capPROJ <- subSnk(SITEcaps=CPcaps, type=c("TRAPTYPE"), info=c("M"))
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=time)  ## this is using 2 months (Feb - Mar)
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=CPsurv, time=time)
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
# dat <- prepSCR(SCRcaps, SCReff, grid = fullX)
## If error and need to do manual
dat <- prepSCRman(SCRcaps, SCReff, grid = fullX)

## Observations, already in order of CellID locations
y <- dat$y

## Uniquely marked individuals
nind <- nrow(y)

## Get sizes of individuals
snsz <- getSize(capPROJ, SCRcaps, subcap)[,2]  ## if all snakes have a measurement during that project
# snsz <- getSizeman(capPROJ, SCRcaps, subcap, time=time2)[,2] ## if some snake sizes are missing than expand window of time
## Categorize by size (1 = <850, 2 = 850-<950, 3 = 950-<1150, 1150 and >)
snsz <- ifelse(snsz < 850, 1,
               ifelse(snsz >= 850 & snsz < 950, 2,
                      ifelse(snsz >= 950 & snsz < 1150, 3,
                             ifelse(snsz >= 1150, 4, -9999))))
if(max(snsz) == -9999) stop('snake size incorrect')
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

#### FORMAT DATA FOR SEMI-COMPLETE LIKELIHOOD SCR ANALYSIS ####

e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

## Integration grid
Ggrid <- cellsize                                #spacing (check sensitivity to spacing)
G <- CPspecs$intgrd[,2:3]
Gpts <- dim(G)[1]                         #number of integration points
a <- Ggrid[1]*Ggrid[2] #CPspecs$area                              #area of each integration grid
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")
points(X, pch=16, col="red")


########################################################
##NIMBLE code for a King et al 2016 semicomplete likelihood

NimModel <- nimbleCode({

  sigma ~ dunif(0,100)
  alpha1 <- 1/(2*sigma*sigma)

  for(l in 1:L){   # size categories
    #prior for intercept
    p0[l] ~ dunif(0,1)
    alpha0[l] <- logit(p0[l])
    
    # Posterior conditional distribution for N-n (and hence N):
    n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures by size category
    Ngroup[l] <- ngroup[l] + n0[l]
  }
  
  N <- sum(Ngroup[1:L])
  
  for(l in 1:L){
    #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
    one_minus_detprob[l,1:Gpts,1:J] <- 1 - p0[l]*exp(-alpha1*Gdist[1:Gpts,1:J]*Gdist[1:Gpts,1:J])*K[1:J]
    #Prob of failure to detect each size category across entire study area and time period
    pdot.temp[l,1:Gpts] <- 1 - prod(one_minus_detprob[l,1:Gpts,1:J])
    #pdot.temp is very close to zero and will lock model up with out this
    pdot[l,1:Gpts] <- max(pdot.temp[l,1:Gpts], 1.0E-10)  
    #prob of detecting a size category at least once in S (a=area of each integration grid)
    # pstar.temp[l] <- (sum(pdot[l,1:Gpts])*a)/A
    # pstar[l] <- CorrectPstar(ps = pstar.temp[l])
  
    # Zero trick for initial 1/pstar^n
    loglikterm[l] <- -ngroup[l] * log(pstar[l])
    lambda[l] <- -loglikterm[l] + 10000
    dummy[l] ~ dpois(lambda[l])
  } #L

  # prior prob for each grid cell   
  pi[1:Gpts] ~ ddirch(b[1:Gpts])
 
  for(i in 1:n){  ## observed n
    ## Activity center grid cell
    s[i] ~ dcat(pi[1:Gpts])
    
    # Model for capture histories of observed individuals:
    for(j in 1:J){  ## number of traps
      y[i,j] ~ dbin(p[i,j],K[j])
      p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
    }#J
  }#I
  
  #derived proportion in each size class
    piGroup[1:L] <- Ngroup[1:L]/N
})

#######################################################

## NIMBLE functions
# CorrectPstar <- nimbleFunction(
#   run = function(ps = double(0)){
#     returnType(double(1))
#     if(ps>=1) return(1)
#     if(ps<1) return(ps)
#   }
# )

## MCMC settings
nchains <- 3; nAdapt=200; nburnin <- 100; niter <- 100+nburnin; nthin <- 1

## Data and constants
constants <- list(J=J, A=A, Gpts=Gpts, nocc=nocc, a=a, n=nind, L=L, size=snsz)

data <- list(y=y, Gdist=Gdist, ngroup=ngroup, dummy=rep(0,L), b=rep(1,Gpts), K=K, pstar=c(0.7,0.5,0.2,0.2))

inits <- list (sigma=runif(1,40,50), n0=(nind+30), s=vsst, p0=runif(L,.002,.003))

parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0","Ngroup","piGroup","pdot","one_minus_detprob")

## Compile and run in NIMBLE
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=data, inits=inits, check=FALSE)
conf <- configureMCMC(Rmodel, monitors=parameters, control = list(adaptInterval = nAdapt), thin=nthin) 

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)

out <- runMCMC(Cmcmc, niter = niter , nburnin = nburnin , nchains = nchains, inits=inits,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
end.time<-Sys.time()
end.time-start.time


