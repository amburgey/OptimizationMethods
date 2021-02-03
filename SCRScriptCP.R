##### CP (Closed Pop, aka NWFN) is a 5-ha closed (fenced to entry and exit of snakes) study area

## CP was created in 2004 and has been used in several projects, resulting in a rich time series with surveys occurring at various densities of snakes

rm(list=ls())

source("Select&PrepVisualData.R")  ## Creation of subcap and subsurv
source("DataPrepCP.R")
# source("SpecifyStateSpace.R")   ## Creation of NWFN grid (not using right now due to issue with perimeter surveys and knowing effort)

library(secr); library(reshape2); library(jagsUI)

## Capture data (subcap) and effort/survey data (subsurv)
CPcaps <- subset(subcap, SITE == "NWFN")
CPsurv <- subset(subsurv, SITE == "NWFN")

## Subset to specific NWFN project (options = MWFM VIS 2, NWFN VIS HL 1, NWFN VIS HL 2, PRE NT2 VIS, POST BT2 VIS, POST KB VIS 1, POST KB VIS 2, POST KB VIS 3 EXTRA, POST KB VIS 3, NWFN VISPACE, NWFN SCENT VIS TRAIL)
CPcaps <- subset(CPcaps, PROJECTCODE == "NWFN VIS 2")
CPsurv <- subset(CPsurv, PROJECTCODE == "NWFN VIS 2")

##### SPECIFY DIMENSIONS OF CP #####
## Make study area grid to ensure correct size
## 27 transects (VIS and TRAP) with 13 points each, 8m from VIS to TRAP and 16m between points
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)
J <- nrow(locs)
X <- as.matrix(locs)
# X <- as.matrix(sapply(CPlocs[,-3], as.numeric))   ## functionalize to not call full SpecifyStateSpace.R
# scrtraps <- nrow(CPlocs)

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
## Don't need to estimate state-space since we know it (5 ha enclosed pop)
delta<- 11.874929
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
# delta<- 0
# Xl<-min(X[,1]) - delta
# Xu<-max(X[,1]) + delta
# Yl<-min(X[,2]) - delta
# Yu<-max(X[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)


#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected or the size of snakes involved
capPROJ <- subSnk(SITEcaps=CPcaps, type=c("TRAPTYPE"), info=c("V"))
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=c("02","04"))  ## this is using 4 months (Feb - April)
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=CPsurv, time=c("02","04"))
## Check data to make sure no missing effort or captured snakes were on survey dates (throws error if dim mismatch)
checkDims(SCReff, SCRcaps)

#### FORMAT DATA FOR ANALYSIS ####
dat <- prepSCR(SCRcaps, SCReff)

# Observations
y <- dat$y
colnames(y) <- NULL

# Uniquely marked individuals
nind <- nrow(y)

# Active/not active for when transects run (one less day surveyed for TL)
act <- as.matrix(dat$act[,-1])
colnames(act) <- NULL
K <- rowSums(act)

# Number of survey occasions
nocc <- ncol(act)


## Data augmentation
# M <- 150
# y <- rbind(y,array(0,dim=c((M-nrow(y)),ncol(y))))

## Initial values for activity centers
# set.seed(02022021)
# sst <- cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
# for(i in 1:nind){
#   sst[i,1] <- mean( X[y[i,]>0,1] )
#   sst[i,2] <- mean( X[y[i,]>0,2] )
# }


########################################################
##Jags model for a simple SCR model

# cat("
# model {
# 
#   lam0 ~ dunif(0,5)
#   sigma ~ dunif(0,100)
#   psi ~ dunif(0,1)
# 
#   for(i in 1:M){
#     z[i] ~ dbern(psi)
#     s[i,1] ~ dunif(Xl,Xu)
#     s[i,2] ~ dunif(Yl,Yu)
#   
#     for(j in 1:J) {
#       d[i,j] <- pow(s[i,1]-locs[j,1], 2) + pow(s[i,2]-locs[j,2],2)
#       p[i,j] <- z[i]*lam0*exp(-(d[i,j])/(2*sigma*sigma))
#       y[i,j] ~ dpois(p[i,j]*K[j])
#       } # locations
#     } # individuals
#   
#   N <- sum(z[])
#   D <- N/A
#   
# }
# ",file = "SCR_CP.txt")

########################################################
##Jags model for a King et al 2016 semicomplete likelihood

cat("
model {

  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)

  for(i in 1:n){
    s[i,1] ~ dunif(Xl,Xu)
    s[i,2] ~ dunif(Yl,Yu)
  }
  
  # Posterior conditional distribution for N-n (and hence N):
  n0 ~ dnegbin(pstar,n)
  N <- n + n0
  
  # pdot = probability of being detected at least once (given location)
  # calculate esa numerically using the integration grid
  
  for(i in 1:G){ # G = number of points on integration grid
    for(s in 1:S){  # WHAT IS THIS???
      for(k in 1:K){  # Assuming this is effort of some kind?
        one_minus_detprob[i,s,k] <- 1 - exp(-dist2[i,k] * tau/2)   ## WHAT IS DIST2???
      }
    }
    pdot.temp[i] <- 1 - prod(one_minus_detprob[i,,])
    pdot[i] <- max(pdot.temp[i], 1.0E-10)
  }
  
  esa <- sum(pdot[])*a # a = size of grid square in numerical integration
  pstar <- esa / A
  
  # Zero trick for initial 1/pstar^n
  loglikterm <- -n * log(pstar)
  lambda <- -loglikterm + 300  ## can play around with 1000 (setting upper limit of pop size)? Changed to 300 for CP
  dummy ~ dpois(lambda) # dummy = 0; entered as data, WHAAAAAAAAAAAAAAAAAAT?
  
  # Model for capture histories of observed individuals:
  for(i in 1:n){
    for(k in 1:K){
      for(s in 1:S){
        capthist[i,s,k] ~ dbern(detprob[i,s,k])
        detprob[i,s,k] <- exp(-r2[i,k] * tau/2 )
      }
      r2[i,k] <- pow(X[i] - traps[k,1], 2) + pow(Y[i] - traps[k,2], 2)
    }
  }
}
",file = "SCRpstar_CP.txt")

#######################################################

## MCMC settings
nc <- 3; nAdapt=10; nb <- 1; ni <- 100+nb; nt <- 1

## Data and constants
# jags.data <- list (y=y, locs=X, M=M, J=J, Xu=Xu, Xl=Xl, Yu=Yu, Yl=Yl, A=A, K=K)  ## data aug
jags.data <- list (y=y, G=nrow(X), J=J, Xu=Xu, Xl=Xl, Yu=Yu, Yl=Yl, A=A, K=K, a=(8*16), n=nind, dummy=0)  ## semicomplete likelihood

inits <- function(){
  list (sigma=runif(1,50,70), z=c(rep(1,nind),rep(0,M-nind)), s=sst, psi=runif(1), lam0=runif(1,0.05,0.1))
}

parameters <- c("sigma","psi","N","D","lam0")

out <- jags("SCR_CP.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters)

save(out, file="Results/NWFNVIS2_SCRpstarvis4mo.Rdata")  ## M = 150 (XXXXhrs)


