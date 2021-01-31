##### Begin analysis small - start with one site and build from there #####
## Create dataframes of captures (PITTAG by Date) and activity (TRANSECT & LOCATION by Date)

rm(list=ls())

source("Select&PrepVisualData.R")  ## Creation of subcap and subsurv
source("SpecifyStateSpace.R")   ## Creation of hmuSpace
source("OverlayGridObservations.R")  ## PrepDat lives here

library(tibble); library(jagsUI); library(abind) #library(nimble)

# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

## Start with HMU as that state space is a little tricky
## HMU EDGE EFFECT
caps <- tibble(subcap)
caps$seen <- c(1)
survs <- tibble(subsurv)

## Prepare data
data <- PrepDat(hmuEdge=hmuEdge, subsurv=subsurv, hmuSpace=hmuSpace)

## When these locations were "active" (i.e., surveys done)
act <- data$act[,-(1:2)] ## remove LOCATION info
## Number of nights surveys were conducted at each transect grid
K <- apply(act, 1, sum)  ## number of occasions a camera was active
nocc <- ncol(act)

## Location of captured snakes on grid
obs <- as.matrix(data$caps[,-1])  ## remove PITTAG ID
colnames(obs) <- NULL
rownames(obs) <- NULL
# Uniquely marked individuals
nind <- nrow(obs)

## Set up data augmentation for the encounter histories, add 0s to end of data
M <- 500  ## 250 hits the limit
y <- rbind(obs,array(0,dim=c((M-nrow(obs)),ncol(obs))))

## Location of visual searches on grid
locs <- as.matrix(data$locs[,2:3])
row.names(locs) <- data$locs[,1]
J <- nrow(locs)

## Habitat matrix needed to define state space due to orientation of study area PLUS one-way fence (snakes can leave but not enter)
habmat <- HMUhabmat   ## From SpecifyStateSpace.R
xu = nrow(habmat)
yu = ncol(habmat)
A = 550000  ## known area of HMU (55 ha -> 550,000m2)

## Initial values for activity centers
set.seed(182021)
SX <- rep(mean(locs[,1]),M) ## CHANGE IF HARD TIME CONVERGING
SY <- rep(mean(locs[,2]),M) ## CHANGE IF HARD TIME CONVERGING


###### JAGS Model ###### > phase to NIMBLE once JAGS works
cat(file="trickpatch_model.txt","
model {
  lam0 ~ dunif(0,5)
  sigma ~ dunif(0,100)
  psi ~ dunif(0,1)
     
  for (i in 1:M){           # loop through the augmented population
    z[i] ~ dbern(psi)       # state of individual i (real or imaginary)
    SX[i] ~ dunif(0,xu)     # priors for the activity centers for each individual
    SY[i] ~ dunif(0,yu)     # lower x coordinate = 0, xu is the upper x value
    pOK[i] <- habmat[trunc(SX[i]+1), trunc(SY[i]+1)]      # habitat check
    OK[i] ~ dbern(pOK[i])   # OK[i] = 1, the ones trick
    
    for(j in 1:J) {         # loop through the J camera trap locations
      d[i,j] <- pow(SX[i]-locs[j,1], 2) + pow(SY[i]-locs[j,2],2)
      p[i,j] <- z[i]*lam0*exp(-(d[i,j])/(2*sigma*sigma))
      y[i,j] ~ dpois(p[i,j]*K[j])
    } ## J
    
  } ## M
  
  N <- sum(z[1:M])   # derive number (check that N << M)
  D <- N / A         # derive density

}")


# MCMC settings
nc <- 3; nAdapt=1000; nb <- 1; ni <- 2000+nb; nt <- 1

# data and constants
jags.data <- list (y=y, locs=locs, M=M, J=J, xu=xu, yu=yu, A=A, K=K, habmat=habmat)

inits <- function(){
  list (sigma=runif(1,50,70), z=c(rep(1,nind),rep(0,M-nind)), SX=SX, SY=SY, psi=runif(1), lam0=runif(1,0.05,0.1))
}

parameters <- c("sigma","psi","N","D","lam0")

out <- jags("trickpatch_model.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters)

save(out, file="Results/hmuEdge_SCRvisM500.Rdata")  ## M = 250 (~4-5hrs)






## Next NCRU as that state space is entirely unfenced
## NCRI/NCRR EDGE EFFECT

ncrEE

## Prepare data
data <- PrepDat(hmuEdge=hmuEdge, subsurv=subsurv, hmuSpace=hmuSpace)

## When these locations were "active" (i.e., surveys done)
act <- data$act[,-(1:2)] ## remove LOCATION info
## Number of nights surveys were conducted at each transect grid
K <- apply(act, 1, sum)  ## number of occasions a camera was active
nocc <- ncol(act)

## Location of captured snakes on grid
obs <- as.matrix(data$caps[,-1])  ## remove PITTAG ID
colnames(obs) <- NULL
rownames(obs) <- NULL
# Uniquely marked individuals
nind <- nrow(obs)

## Set up data augmentation for the encounter histories, add 0s to end of data
M <- 500  ## 250 hits the limit
y <- rbind(obs,array(0,dim=c((M-nrow(obs)),ncol(obs))))

## Location of visual searches on grid
locs <- as.matrix(data$locs[,2:3])
row.names(locs) <- data$locs[,1]
J <- nrow(locs)

## Habitat matrix needed to define state space due to orientation of study area PLUS one-way fence (snakes can leave but not enter)
habmat <- HMUhabmat   ## From SpecifyStateSpace.R
xu = nrow(habmat)
yu = ncol(habmat)
A = 550000  ## known area of HMU (55 ha -> 550,000m2)

## Initial values for activity centers
set.seed(182021)
SX <- rep(mean(locs[,1]),M) ## CHANGE IF HARD TIME CONVERGING
SY <- rep(mean(locs[,2]),M) ## CHANGE IF HARD TIME CONVERGING


###### JAGS Model ###### > phase to NIMBLE once JAGS works
cat(file="trickpatch_model.txt","
model {
  lam0 ~ dunif(0,5)
  sigma ~ dunif(0,100)
  psi ~ dunif(0,1)
     
  for (i in 1:M){           # loop through the augmented population
    z[i] ~ dbern(psi)       # state of individual i (real or imaginary)
    SX[i] ~ dunif(0,xu)     # priors for the activity centers for each individual
    SY[i] ~ dunif(0,yu)     # lower x coordinate = 0, xu is the upper x value
    pOK[i] <- habmat[trunc(SX[i]+1), trunc(SY[i]+1)]      # habitat check
    OK[i] ~ dbern(pOK[i])   # OK[i] = 1, the ones trick
    
    for(j in 1:J) {         # loop through the J camera trap locations
      d[i,j] <- pow(SX[i]-locs[j,1], 2) + pow(SY[i]-locs[j,2],2)
      p[i,j] <- z[i]*lam0*exp(-(d[i,j])/(2*sigma*sigma))
      y[i,j] ~ dpois(p[i,j]*K[j])
    } ## J
    
  } ## M
  
  N <- sum(z[1:M])   # derive number (check that N << M)
  D <- N / A         # derive density

}")


# MCMC settings
nc <- 3; nAdapt=1000; nb <- 1; ni <- 2000+nb; nt <- 1

# data and constants
jags.data <- list (y=y, locs=locs, M=M, J=J, xu=xu, yu=yu, A=A, K=K, habmat=habmat)

inits <- function(){
  list (sigma=runif(1,50,70), z=c(rep(1,nind),rep(0,M-nind)), SX=SX, SY=SY, psi=runif(1), lam0=runif(1,0.05,0.1))
}

parameters <- c("sigma","psi","N","D","lam0")

out <- jags("trickpatch_model.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters)

save(out, file="Results/hmuEdge_SCRvisM500.Rdata")  ## M = 250 (~4-5hrs)

