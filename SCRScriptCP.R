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
## Check area: 
A <- (Xu-Xl)*(Yu-Yl)


#### PREP DATA FOR SCR ANALYSIS ####
## Subset data based on how it was collected or the size of snakes involved
capPROJ <- subSnk(SITEcaps=CPcaps, type=c("TRAPTYPE"), info=c("V"))
## Subset data based on sampling time of interest and order by dates and sites
SCRcaps <- subYr(SITEcaps=capPROJ, time=c("03","04"))  ## specify month range
## Find effort for this set of snakes and time
SCReff <- effSnk(eff=CPsurv, time=c("03","04"))
## Check data to make sure no missing effort or captured snakes were on survey dates (throws error if dim mismatch)
checkDims(SCReff, SCRcaps)

#### FORMAT DATA FOR TRADITIONAL SCR ANALYSIS ####
dat <- prepSCR(SCRcaps, SCReff)

## Observations
y <- dat$y
colnames(y) <- NULL

## Uniquely marked individuals
nind <- nrow(y)

## Active/not active for when transects run (one less day surveyed for TL)
act <- as.matrix(dat$act[,-1])
colnames(act) <- NULL
K <- rowSums(act)

## Number of survey occasions
nocc <- ncol(act)


## Data augmentation
# M <- 150
# y <- rbind(y,array(0,dim=c((M-nrow(y)),ncol(y))))

## Initial values for activity centers
set.seed(02022021)
# sst <- cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
sst <- cbind(runif(nind,Xl,Xu),runif(nind,Yl,Yu))
for(i in 1:nind){
  sst[i,1] <- mean( X[y[i,]>0,1] )
  sst[i,2] <- mean( X[y[i,]>0,2] )
}


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

#######################################################

# ## MCMC settings
# nc <- 3; nAdapt=10; nb <- 1; ni <- 100+nb; nt <- 1
# 
# ## Data and constants
# # jags.data <- list (y=y, locs=X, M=M, J=J, Xu=Xu, Xl=Xl, Yu=Yu, Yl=Yl, A=A, K=K)  ## data aug
# jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, Xu=Xu, Xl=Xl, Yu=Yu, Yl=Yl, A=A, K=K, a=a, n=nind, dummy=0)  ## semicomplete likelihood
# 
# inits <- function(){
#   list (sigma=runif(1,50,70), z=c(rep(1,nind),rep(0,M-nind)), s=sst, psi=runif(1), lam0=runif(1,0.05,0.1))
# }
# 
# parameters <- c("sigma","psi","N","D","lam0")
# 
# out <- jags("SCR_CP.txt", data=jags.data, inits=inits, parallel=TRUE,
#             n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters)
# 
# save(out, file="Results/NWFNVIS2_SCRvisM1504mo.Rdata")


#### FORMAT DATA FOR SEMI-COMPLETE LIKELIHOOD SCR ANALYSIS ####

e2dist <- function (x, y) {
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#Integration grid
Ggrid <- 5                               #spacing (verify sensitivity to spacing)
Xlocs <- rep(seq(Xl, Xu, Ggrid), times = 47)
Ylocs <- rep(seq(Yl, Yu, Ggrid), each = 44)
G <- cbind(Xlocs, Ylocs)
# Xlocs <- seq(Yl,Yu,Ggrid)          
# G <- cbind(sort(rep(Xlocs,length(Xlocs))),rep(Xlocs,length(Xlocs))) #integration grid locations
Gpts <- dim(G)[1]                            #number of integration points
a <- Ggrid^2                                 #area of each integration grid
Gdist <- e2dist(G, X)                      #distance between integration grid locations and traps
plot(G, pch=16, cex=.5, col="grey")
points(X, pch=16, col="red")
points(Xl,Yu, pch=21, col="blue")  #check that state space area match for CP
points(Xl,Yl, pch=21, col="blue")
points(Xu,Yu, pch=21, col="blue")
points(Xu,Yl, pch=21, col="blue")

## compute the squared distance matrix to be used as DATA if using integratin points as trap locations  
# dsq <- matrix(NA, nPix, J)
# for(i in 1:nPix) {
#   dsq[i,] <-  ((G$x[i] - trapmat[,"x"])^2 + (G$y[i] - trapmat[,"y"])^2)
# }

########################################################
##Jags model for a King et al 2016 semicomplete likelihood

cat("
model {

  p0 ~ dunif(0,1)
  alpha0 <- logit(p0)
  sigma ~ dunif(0,100)
  alpha1 <- 1/(2*sigma*sigma)
  
  # Posterior conditional distribution for N-n (and hence N):
  n0 ~ dnegbin(pstar,n)  # number of failures
  N <- n + n0  # successful observations plus failures to observe = total N
  
  #Probability of capture for integration grid points
  #pdot = probability of being detected at least once (given location)
  
  ## K loop if any time-varying covariates need to be added in
  # for(g in 1:Gpts){ # Gpts = number of points on integration grid
  #   for(k in 1:nocc){  # K = number of occasions
  #     for(j in 1:J){  # J = number of traps
  #     #Probability of being missed at grid cell g at survey k and trap j
  #       one_minus_detprob[g,k,j] <- 1 - p0*exp(-alpha1*Gdist[g,j]*Gdist[g,j])*act[k,j] #Gdist given as data
  #     } #J
  #   } #K
  #   pdot.temp[g] <- 1 - prod(one_minus_detprob[g,,]) #Prob of failure to detect across entire study area and time period
  #   pdot[g] <- max(pdot.temp[g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
  # } #G

## Remove k loop (for decreased run time?)
  for(g in 1:Gpts){ # Gpts = number of points on integration grid
    for(j in 1:J){  # J = number of traps
      #Probability of being missed at grid cell g and trap j multiplied by total effort (K) at that trap
      one_minus_detprob[g,j] <- 1 - p0*exp(-alpha1*Gdist[g,j]*Gdist[g,j])*K[j] #Gdist given as data
    } #J
  pdot.temp[g] <- 1 - prod(one_minus_detprob[g,]) #Prob of failure to detect across entire study area and time period
  pdot[g] <- max(pdot.temp[g], 1.0E-10)  #pdot.temp is very close to zero and will lock model up with out this
} #G
  
  pstar <- (sum(pdot[1:Gpts])*a)/A   #prob of detecting an individual at least once in S (a=area of each integration grid, a given as data)
  
  ### NO CHANGING, TO MAKE JAGS/BUGS LIKELIHOOD FUNCTION PROPERLY ### 
  # Zero trick for initial 1/pstar^n
  loglikterm <- -n * log(pstar)
  lambda <- -loglikterm + 1000
  dummy ~ dpois(lambda) # dummy = 0; entered as data
  ### NO CHANGING, TO MAKE JAGS/BUGS LIKELIHOOD FUNCTION PROPERLY ###
  
  
  for(i in 1:n){  ## n = number of observed individuals
  ## For use when defining state space and traps traditionally
    # s[i,1] ~ dunif(Xl,Xu)
    # s[i,2] ~ dunif(Yl,Yu)
  ## For use when defining traps on a grid cell
    pi[1:Gpts] ~ ddirch(b[1:Gpts])
    s[i] ~ dcat(pi[1:Gpts])
    
    # Model for capture histories of observed individuals:
    for(j in 1:J){  ## J = number of traps
      y[i,j] ~ dbin(p[i,j],K[j])
      d[i,j] <- pow(pow(s[i,1]-locs[j,1],2) + pow(s[i,2]-locs[j,2],2),0.5)
      p[i,j] <- p0*exp(-alpha1*d[i,j]*d[i,j])
      # trapdist calculated outside model for dcat option, could use some integration grid and denote which are traps
      # g[i,j] <- exp(-dsq[s[i],j]/sigma2)  ## check if sigma or sigma squared
    }#J
  }#n
}
",file = "SCRpstar_CP.txt")

#######################################################

## MCMC settings
nc <- 3; nAdapt=1000; nb <- 1; ni <- 10000+nb; nt <- 1

## Data and constants
jags.data <- list (y=y, Gpts=Gpts, Gdist=Gdist, J=J, locs=X, Xu=Xu, Xl=Xl, Yu=Yu, Yl=Yl, A=A, K=K, nocc=nocc, a=a, n=nind, dummy=0, b=rep(1,Gpts), act=t(act)) # ## semicomplete likelihood

inits <- function(){
  list (sigma=runif(1,45,50), n0=(M-nind), s=sst[1:nind,], p0=runif(1,.002,.003)) #ran at 0.002 and 0.003 before
}

parameters <- c("p0","sigma","pstar","alpha0","alpha1","N")

out <- jags("SCRpstar_CP.txt", data=jags.data, inits=inits, parallel=TRUE,
            n.chains=nc, n.burnin=nb,n.adapt=nAdapt, n.iter=ni, parameters.to.save=parameters)#, factories = "base::Finite sampler FALSE") ## might have to use to keep JAGS from locking up with large categorical distribution, will speed things up a little

save(out, file="Results/NWFNVIS2_SCRpstarvistest2noK2months.Rdata")  ## M = 150 (XXXXhrs)


