#### SIMULATION USING SCR FRAMEWORK ####
## Detect snakes based on different detection process for each method
## Use sigma info from "normal" snake density to start (23 snakes/ha; Rodda et al. 1999, Tyrrell et al. 2009, Christy et al. 2010)
## sigma ~ dgamma(274.69, 7.27) based on 2015 SCR analysis
## sigma mean (CI; SD) = 37.8 (33.6, 42.4; 2.3)

rm(list=ls())
library(secr)
library(nimble)
library(coda)
library(devtools)
set.seed(2018)

Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# MCMC settings
nc <- 3; nAdapt=500; nb <- 1000; ni <- 5000+nb; nt <- 1

## Define study area grid (random example currently)
locs <- as.matrix(secr::make.grid(nx = 20, ny = 20, spacex = 8, spacey = 8))

## Which parts of grid have traps
set.seed(922020)
a=sample(200, 15)
Yscr=locs[a,]
ntraps <- nrow(Yscr)

## Define state-space of point process. (i.e., where animals live).
## "delta" just adds a fixed buffer to the outer extent of the traps.
delta <- 20  ## will need to play with this to figure out delta where N doesn't change substantially
Xl<-min(locs[,1]) - delta
Xu<-max(locs[,1]) + delta
Yl<-min(locs[,2]) - delta
Yu<-max(locs[,2]) + delta
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)

## Establish parameters for detection
lam0<- 0.07  # Christy et al. (2010) mean detection prob (but could reach 0.18 under ideal situations)
sigma<- 37.8  # Based on Amburgey et al. (in review) and Gardner SSP analysis

## Number of snakes based on probable density on Guam
## 23 snakes/ha -> 0.0023 snakes/m2 (see notes at beginning)
N <- round(A*0.0023)

## Number of nights trapping
K <- 20

## Simulate activity centers
sx<-runif(N,Xl,Xu)
sy<-runif(N,Yl,Yu)
S<-cbind(sx,sy) 

## Distance between each individual AC and each each trap
## Function to calculate distance between two sets of (x,y) locations
e2dist <- function(A, B)  {
  xdif <- outer(A[, 1], B[, 1], "-")
  ydif <- outer(A[, 2], B[, 2], "-")
  sqrt(xdif^2 + ydif^2)
}

D <- e2dist(S,Yscr)
muy <- lam0*exp(-(D*D)/(2*sigma^2))






## Simulate observations of individuals by trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:nrow(Y)){
  Y[i,]<-rpois(ntraps,K*muy[i,])
}
datnam<-paste('Dat_', iter, '.R', sep='')
dput(Y, datnam)

## NOTE NOTE NOTE NOTE
## Y is a matrix of encounter frequencies of EACH individual in EACH trap
## As simulated here it includes the "all 0" observations.  We want
## to delete those to mimic real data.

Yscr=Y[,a]
totalcaps<-apply(Yscr,1,sum)
Yscr<-Yscr[totalcaps>0,]

## Data augmentation
M <- 200
y <- rbind(y,matrix(0,nrow=M-nind,ncol=ncol(y)))
z <- c(rep(1,nind),rep(0,M-nind))

## Starting values for s
sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
for(i in 1:nind){
  sst[i,1] <- mean( X[y[i,]>0,1] )
  sst[i,2] <- mean( X[y[i,]>0,2] )
}


## NIMBLE model is nearly identical to BUGS
code <- nimbleCode({
  alpha0 ~ dnorm(0,.1)
  logit(p0) <- alpha0
  alpha1 ~ dnorm(0,.1)
  sigma <- sqrt(1/(2*alpha1))
  psi ~ dunif(0,1)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d[i,1:J] <- pow(pow(s[i,1]-X[1:J,1],2) + pow(s[i,2]-X[1:J,2],2),0.5) #accepts vector assignments like R
    p[i,1:J] <- z[i]*p0*exp(- alpha1*d[i,1:J]*d[i,1:J])
    for(j in 1:J){
      y[i,j] ~ dbin(p[i,j],K)
    }#J
  }#i
  N <- sum(z[1:M]) #Must specify dimensions in NIMBLE
  D <- N/64
})#code


# Separate data and constants (constants appear only on right-hand side of formulas)
nim.data <- list (y=y)
constants <- list (X=X, K=K, M=M, J=J, xlim=xlim, ylim=ylim)

# Initial values (same as BUGS)
inits <- function(){
  list (alpha0=rnorm(1,-4,.4), alpha1=runif(1,1,2), s=sst, z=z)
}

# Parameters (same as BUGS)
parameters <- c("alpha0","alpha1","sigma","N","D")


## Nimble steps
start.time <- Sys.time()  
Rmodel <- nimbleModel(code=code, constants=constants, data=nim.data)
conf <- configureMCMC(Rmodel,monitors=parameters,control = list(adaptInterval = nAdapt))
Rmcmc <- buildMCMC(conf)  #produces an uncompiled R mcmc function
Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
samplesList <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits,
                       setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time <- Sys.time()
SCR0time<-end.time - start.time

#summarize
summaryList<-summary(samplesList)
outSummary<-cbind(summaryList$statistics[,c("Mean","SD")],summaryList$quantiles[,c("2.5%","50%", "97.5%")],gelman.diag(samplesList,multivariate=FALSE)$psrf[,1],effectiveSize(samplesList))
colnames(outSummary)[6:7]<-c("Rhat","n.eff")
round(outSummary,2)

#plot results
plot(samplesList[,"sigma"])