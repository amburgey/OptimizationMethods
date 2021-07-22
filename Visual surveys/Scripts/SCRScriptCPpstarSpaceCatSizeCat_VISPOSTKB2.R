##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area

## CP was created in 2004 and has been used in several projects, resulting in a rich time series with surveys occurring at various densities of snakes

rm(list=ls())

source("Select&PrepVisualData.R")   ## Creation of subcap and subsurv (cleaned up)
source("Visual surveys/DataPrep/DataPrepCP_VISPOSTKB2.R")              ## Functions to reshape survey and capture data
source("Visual surveys/DataPrep/OverlayCPGrid.R")

library(secr); library(reshape2); library(jagsUI)

## Subset capture data (subcap) and effort/survey data (subsurv)
CPcaps <- subset(subcap, SITE == "NWFN")
CPsurv <- subset(subsurv, SITE == "NWFN")

## Subset to specific NWFN project
CPcaps <- subset(CPcaps, PROJECTCODE == "POST KB VIS 2")
CPsurv <- subset(CPsurv, PROJECTCODE == "POST KB VIS 2")

## SECIFY TIME FRAME
time <- c("01","02")
time2 <- c("2011-11-01","2012-03-30")


##### SPECIFY DIMENSIONS OF CP #####
cellsize <- c(10,10)  ## dimensions of integration grid cell
CPspecs <- overlayCP(CPcaps, cellsize)  ## ignore warnings, all about projections
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
SCRcaps <- subYr(SITEcaps=capPROJ, time=time)  ## this is using 2 months (Feb - Mar)
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
## If error and need to do manual
# dat <- prepSCRman(SCRcaps, SCReff, grid = fullX)

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
Ggrid <- cellsize                         #spacing
G <- CPspecs$intgrd[,2:3]
Gpts <- dim(G)[1]                         #number of integration points
a <- CPspecs$area                         #area of each integration grid
Gdist <- e2dist(G, X)                     #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")
points(X, pch=16, col="red")


########################################################
##Jags model for a King et al 2016 semicomplete likelihood

cat("
model {
  
  sigma ~ dunif(0,100)
  alpha1 <- 1/(2*sigma*sigma)
  
  for(l in 1:L){   # 4 size categories
    #prior for intercept
    p0[l] ~ dunif(0,1)
    alpha0[l] <- logit(p0[l])
    
    # Posterior conditional distribution for N-n (and hence N):
    n0[l] ~ dnegbin(pstar[l],ngroup[l])  # number of failures by size category
    Ngroup[l] <- ngroup[l] + n0[l]
  }
  
  N <- sum(Ngroup[1:L])  # successful observations plus failures to observe of each size = total N
  
  #Probability of capture for integration grid points
  #pdot = probability of being detected at least once (given location)
  
  for(l in 1:L){  # size category
    for(g in 1:Gpts){ # Gpts = number of points on integration grid
      for(j in 1:J){  # J = number of traps
        #Probability of an individual of size i being missed at grid cell g and trap j multiplied by total effort (K) at that trap
        miss_allK[l,g,j] <- pow((1 - p0[l]*exp(-alpha1*Gdist[g,j]*Gdist[g,j])),K[j])
      } #J
      pdot.temp[l,g] <- 1 - prod(miss_allK[l,g,]) #Prob of detect each size category across entire study area and time period
      pdot[l,g] <- max(pdot.temp[l,g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
    } #G
    pstar[l] <- (sum(pdot[l,1:Gpts]*a[1:Gpts]))/A #prob of detecting a size category at least once in S (a=area of each integration grid, given as data)
    
    # Zero trick for initial 1/pstar^n
    loglikterm[l] <- -ngroup[l] * log(pstar[l])
    lambda[l] <- -loglikterm[l] + 10000
    dummy[l] ~ dpois(lambda[l]) # dummy = 0; entered as data
  } #L
  
  # prior prob for each grid cell (setting b[1:Gpts] = rep(1,Gpts) is a uniform prior across all cells)   
  pi[1:Gpts] ~ ddirch(b[1:Gpts])
  
  for(i in 1:n){  ## n = number of observed individuals
    ## For use when defining traps on a grid cell
    s[i] ~ dcat(pi[1:Gpts])
    
    # Model for capture histories of observed individuals:
    for(j in 1:J){  ## J = number of traps
      y[i,j] ~ dpois(p[i,j]*K[j])
      p[i,j] <- p0[size[i]]*exp(-alpha1*Gdist[s[i],j]*Gdist[s[i],j])
    }#J
  }#I
  
  #derived proportion in each size class
  for(l in 1:L){
    piGroup[l] <- Ngroup[l]/N
  }
}
",file = "Visual surveys/Models/SCRpstarCATsizeCAT_CP.txt")

#######################################################

## MCMC settings
nc <- 3; nAdapt=200; nb <- 100; ni <- 2500+nb; nt <- 1

## Data and constants
jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=rep(0,L), b=rep(1,Gpts), size=snsz, L=L, ngroup=ngroup) # ## semicomplete likelihood
# jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=0, b=rep(1,Gpts)) 

inits <- function(){
  list (sigma=runif(1,30,40), n0=c(ngroup+10), s=vsst, p0=runif(L,.002,.003))
}
# inits <- function(){
#   list (sigma=runif(1,30,40), n0=(nind+30), s=vsst, p0=runif(1,.002,.003))
# }

parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0","Ngroup","piGroup")
# parameters <- c("p0","sigma","pstar","alpha0","alpha1","N","n0")

out <- jags("Visual surveys/Models/SCRpstarCATsizeCAT_CP.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## might have to use "factories" to keep JAGS from locking up with large categorical distribution, will speed things up a little
# out <- jags("Archive/SCRpstarCAT_CP.txt", data=jags.data, inits=inits, parallel=TRUE,
#             n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters, factories = "base::Finite sampler FALSE") ## might have to use "factories" to keep JAGS from locking up with large categorical distribution, will speed things up a little


save(out, file="Visual surveys/Results/NWFNVISPOSTKB2_SCRpstarvisCATsizeCATdpois.Rdata")
# save(out, file="Visual surveys/Results/NWFNVISPOSTKB2_SCRpstarvisCATNOSIZEgrid10.Rdata")

